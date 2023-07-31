import shutil
import tarfile
from gzip import GzipFile
from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd
import pyranges as pr
import requests


def read_remote_gzip_gtf(url: str, chrom: str = "chr2L"):
    """
    Download a gzipped gtf file, extract chromosome of interest, and read into a pyranges object
    :param url: url to download from
    :param chrom: chromosome to extract
    """

    response = requests.get(url, stream=True)
    with NamedTemporaryFile(suffix="gtf") as f:
        with GzipFile(fileobj=response.raw) as gz:
            shutil.copyfileobj(gz, f)
            gtf = pr.read_gtf(f.name, as_df=True)

    return gtf.loc[gtf.Chromosome == chrom]


if __name__ == "__main__":
    # setup data dir
    DATADIR = Path(__file__).parent / "data"
    DATADIR.mkdir(exist_ok=True)

    # RMSK GTF
    print(f"Saving RMSK GTF to {str(DATADIR / 'rmsk.gtf')}...")
    RMSK_GTF = "https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/test_data/testdata_small/dm3_rmsk_TE_trimmed.gtf.gz"
    rmsk_gtf = read_remote_gzip_gtf(RMSK_GTF)
    pr.PyRanges(rmsk_gtf).to_gtf(str(DATADIR / "rmsk.gtf"))
    rmsk_gtf["te_id"] = (
        rmsk_gtf["gene_id"] + ":" + rmsk_gtf["family_id"] + ":" + rmsk_gtf["class_id"]
    )

    # refgene GTF
    print(f"Saving refgene GTF to {str(DATADIR / 'refgene.gtf')}...")
    REFGENE_GTF = "https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/test_data/testdata_small/dm3_refGene_trimmed.gtf.gz"
    refgene_gtf = read_remote_gzip_gtf(REFGENE_GTF)
    pr.PyRanges(refgene_gtf).to_gtf(str(DATADIR / "refgene.gtf"))

    # TEtranscripts results
    print(f"Saving tetranscripts results to {str(DATADIR / 'res.cntTable')}...")
    RES = "https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/test_data/testdata_small/results/TEtranscripts_Galaxy_test.cntTable"
    res = pd.read_csv(RES, sep="\t", index_col=0)
    keep_te = res.index.isin(rmsk_gtf.te_id)
    keep_gene = res.index.isin(refgene_gtf.gene_id)
    res[keep_te | keep_gene].to_csv(DATADIR / "res.cntTable", sep="\t")

    # RMSK OUT
    print(f"Saving rmsk.out results to {str(DATADIR / 'rmsk.out')}...")
    RMSK = "https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromOut.tar.gz"
    response = requests.get(RMSK, stream=True)
    with NamedTemporaryFile(suffix="tar.gz", mode="wb") as f:
        shutil.copyfileobj(response.raw, f)

        with tarfile.open(f.name, mode="r:gz") as tar:
            member = tar.getmember("2L/chr2L.fa.out")
            tar.extract(member, path=DATADIR)

    shutil.move(DATADIR / "2L/chr2L.fa.out", DATADIR / "rmsk.out")
    shutil.rmtree(DATADIR / "2L")

    # Coldata
    print(f"Saving coldata to {str(DATADIR / 'coldata.csv')}...")
    df = pd.DataFrame(
        {
            "filename": [
                "testData_treatment_rep1_SE.bam.T",
                "testData_treatment_rep2_SE.bam.T",
                "testData_control_rep1_SE.bam.C",
                "testData_control_rep2_SE.bam.C",
            ],
            "condition": ["treatment", "treatment", "control", "control"],
            "replicate": ["rep1", "rep2", "rep1", "rep2"],
        }
    )
    df.to_csv(DATADIR / "coldata.csv", index=False)
