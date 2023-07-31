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
    """
    Load TEtranscripts results into an AnnData object
    :param te_counts: TEtranscripts counts table (columns are samples, rows are genes/TEs)
    :param coldata: TEtranscripts coldata table
    :param gene_gtf: GTF file with gene annotations
    :param rmsk_gtf: GTF file with TE annotations
    :param rmsk_out: RepeatMasker output file
    """

    assert (
        "condition" in coldata.columns
    ), "coldata must have a column named 'condition'"

    # gene info
    print(f"Reading gene GTF from {gene_gtf}...")
    gtf = pr.read_gtf(gene_gtf)

    # check if GENCODE
    if "gene" in gtf.Feature.unique():
        gtf = gtf.df
    else:
        genes = gtf.boundaries("gene_id").df
        genes["Feature"] = "gene"
        genes["Source"] = "refGene"
        gtf = pd.concat([gtf.df, genes]).sort_values(["Chromosome", "Start"])
        gtf = gtf.fillna(".")
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
    te_counts = te_counts.T
    adata = AnnData(
        X=te_counts,
        dtype=te_counts.values.dtype,
        obs=coldata.loc[te_counts.index.values, :],
        var=var_data.loc[te_counts.columns.values, :],
    )
    adata.obs["total_reads"] = adata.X.sum(axis=1)

    return adata
