import pandas as pd

from myutils.tetranscripts import load_tetranscripts


def test_load_tetranscripts():
    df = pd.read_csv("tests/data/res.cntTable", sep="\t", index_col=0)
    coldata = pd.read_csv("tests/data/coldata.csv", index_col=0)

    load_tetranscripts(
        te_counts=df,
        coldata=coldata,
        gene_gtf="tests/data/refgene.gtf",
        rmsk_gtf="tests/data/rmsk.gtf",
        rmsk_out="tests/data/rmsk.out",
    )
