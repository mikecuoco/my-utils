import pandas as pd
import pyranges as pr
from anndata import AnnData

from .rmsk import read_rmsk


def load_tetranscripts(
    te_counts: pd.DataFrame,
    coldata: pd.DataFrame,
    gene_gtf: str,
    rmsk_gtf: str,
    rmsk_out: str,
):
    assert (
        "condition" in coldata.columns
    ), "coldata must have a column named 'condition'"

    # gene info
    print(f"Reading gene GTF from {gene_gtf}...")
    gtf = pr.read_gtf(gene_gtf).df
    gtf = gtf[gtf.Feature == "gene"].set_index("gene_id")  # get genes only

    # te info
    # calculate age
    print(f"Reading RepeatMasker output from {rmsk_out}...")
    rmsk = read_rmsk(rmsk_out)
    rmsk["gene_id"] = rmsk["repName"] + ":" + rmsk["repFamily"] + ":" + rmsk["repClass"]
    avg_age = rmsk[["age", "gene_id"]].groupby("gene_id").mean()

    # read in rmsk gtf
    print(f"Reading TE GTF from {rmsk_gtf}...")
    rmsk = pr.read_gtf(rmsk_gtf).df
    rmsk["subfamily_id"] = rmsk["gene_id"]
    rmsk["gene_id"] = rmsk["gene_id"] + ":" + rmsk["family_id"] + ":" + rmsk["class_id"]
    n_copies = rmsk[["transcript_id", "gene_id"]].groupby("gene_id").size()
    rmsk = (
        rmsk[["gene_id", "subfamily_id", "family_id", "class_id"]]
        .drop_duplicates()
        .set_index("gene_id")
    )
    rmsk["n_copies"] = n_copies
    rmsk["avg_age"] = avg_age

    # combine gene and te info
    var_data = pd.concat([gtf, rmsk])

    # put everything together
    adata = AnnData(
        X=te_counts,
        dtype=te_counts.values.dtype,
        obs=coldata.loc[te_counts.index.values, :],
        var=var_data.loc[te_counts.columns.values, :],
    )
    adata.obs["total_reads"] = adata.X.sum(axis=1)

    return adata
