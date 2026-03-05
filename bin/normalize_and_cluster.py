#!/usr/bin/env python3
"""
normalize_and_cluster.py – Merge per-sample k-mer counts, CPM-normalise,
and k-means cluster k-mers whose expression profiles suggest they derive
from the same TE transcript sequence (not locus).

Sub-commands
    merge   – combine per-sample k-mer TSVs into a matrix, apply CPM
    cluster – k-means on the CPM matrix; output cluster assignments +
              per-cluster aggregated expression + diagnostic plots
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# ────────────────────────────────────────────────────────────────────
#  MERGE
# ────────────────────────────────────────────────────────────────────
def cmd_merge(args):
    kmer_files = args.kmer_files.split()
    depths = pd.read_csv(args.depth_file, sep="\t")
    depths = depths.set_index("sample_id")["total_reads"].to_dict()

    frames = {}
    for fpath in kmer_files:
        # Derive sample id from filename: {sample}_kmer_counts.tsv.gz
        sid = os.path.basename(fpath).replace("_kmer_counts.tsv.gz", "")
        df = pd.read_csv(fpath, sep="\t")
        # Collapse across TE families for the matrix; keep family annotation separate
        df["kmer_te"] = df["kmer"] + "|" + df["te_family"]
        frames[sid] = df.set_index("kmer_te")["count"]

    # Build sample × kmer matrix
    mat = pd.DataFrame(frames).fillna(0).astype(int)
    mat.to_csv(args.out_raw, sep="\t", compression="gzip")

    # CPM normalisation: counts / total_reads * 1e6
    cpm = mat.copy().astype(float)
    for sid in cpm.columns:
        total = float(depths.get(sid, cpm[sid].sum()))
        if total == 0:
            total = 1.0
        cpm[sid] = cpm[sid] / total * 1e6
    cpm.to_csv(args.out_cpm, sep="\t", compression="gzip")

    print(f"[merge] {len(kmer_files)} samples, {mat.shape[0]} k-mers",
          file=sys.stderr)


# ────────────────────────────────────────────────────────────────────
#  CLUSTER
# ────────────────────────────────────────────────────────────────────
def cmd_cluster(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cpm = pd.read_csv(args.cpm_matrix, sep="\t", index_col=0)

    # Drop k-mers with zero variance (uninformative)
    variance = cpm.var(axis=1)
    cpm = cpm.loc[variance > 0]

    if cpm.shape[0] == 0:
        sys.exit("[cluster] No variable k-mers found. Aborting.")

    # Scale (z-score across samples per k-mer)
    scaler = StandardScaler()
    X = scaler.fit_transform(cpm.values)

    # ── Elbow plot to help choose k ─────────────────────────────────
    max_k = min(args.n_clusters * 2, X.shape[0])
    ks = list(range(2, max_k + 1, max(1, max_k // 20)))
    inertias = []
    for kk in ks:
        km = KMeans(n_clusters=kk, n_init=10, random_state=42, max_iter=300)
        km.fit(X)
        inertias.append(km.inertia_)

    plt.figure(figsize=(6, 4))
    plt.plot(ks, inertias, "o-")
    plt.xlabel("Number of clusters")
    plt.ylabel("Inertia")
    plt.title("Elbow plot – k-mer clustering")
    plt.tight_layout()
    plt.savefig(args.out_elbow, dpi=150)
    plt.close()

    # ── Final k-means ───────────────────────────────────────────────
    n_clust = min(args.n_clusters, X.shape[0])
    km = KMeans(n_clusters=n_clust, n_init=20, random_state=42, max_iter=500)
    labels = km.fit_predict(X)

    # Cluster assignments
    assign = pd.DataFrame({
        "kmer_te": cpm.index,
        "cluster": labels,
    })
    # Split kmer_te back to kmer + te_family
    split = assign["kmer_te"].str.rsplit("|", n=1, expand=True)
    assign["kmer"] = split[0]
    assign["te_family"] = split[1] if split.shape[1] > 1 else ""
    assign.to_csv(args.out_clusters, sep="\t", index=False)

    # Per-cluster mean expression (CPM) across samples
    cpm_with_label = cpm.copy()
    cpm_with_label["cluster"] = labels
    cluster_expr = cpm_with_label.groupby("cluster").mean()
    cluster_expr.to_csv(args.out_expression, sep="\t")

    # ── PCA visualisation ───────────────────────────────────────────
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X)

    plt.figure(figsize=(7, 6))
    scatter = plt.scatter(coords[:, 0], coords[:, 1],
                          c=labels, cmap="tab20", s=4, alpha=0.6)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title("PCA of k-mer expression – coloured by cluster")
    plt.colorbar(scatter, label="cluster")
    plt.tight_layout()
    plt.savefig(args.out_pca, dpi=150)
    plt.close()

    print(f"[cluster] {cpm.shape[0]} k-mers → {n_clust} clusters",
          file=sys.stderr)


# ────────────────────────────────────────────────────────────────────
#  CLI
# ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd")

    # merge
    m = sub.add_parser("merge")
    m.add_argument("--kmer-files", required=True,
                   help="Space-separated list of per-sample k-mer TSV.gz files")
    m.add_argument("--depth-file", required=True,
                   help="TSV with sample_id and total_reads columns")
    m.add_argument("--out-raw", required=True)
    m.add_argument("--out-cpm", required=True)

    # cluster
    c = sub.add_parser("cluster")
    c.add_argument("--cpm-matrix", required=True)
    c.add_argument("--n-clusters", type=int, default=50)
    c.add_argument("--out-clusters", required=True)
    c.add_argument("--out-expression", required=True)
    c.add_argument("--out-elbow", required=True)
    c.add_argument("--out-pca", required=True)

    args = parser.parse_args()
    if args.cmd == "merge":
        cmd_merge(args)
    elif args.cmd == "cluster":
        cmd_cluster(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
