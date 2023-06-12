import pandas as pd
import pytest

from myutils.wrappers import bedtools_intersect


def make_test_bed():
    bed = pd.DataFrame(
        {
            "Chromosome": ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2"],
            "Start": [0, 10, 20, 0, 10, 20],
            "End": [10, 20, 30, 10, 20, 30],
        }
    )
    return bed


def run_bedtools_intersect(extra=""):
    return bedtools_intersect(make_test_bed(), make_test_bed(), extra)


def test_bedtools_intersect():
    inter = run_bedtools_intersect()
    assert inter.shape[1] == 6, "regular intersection must have 6 rows"


def test_bedtools_intersect_wb():
    inter = run_bedtools_intersect("-wb")
    assert inter.shape[1] == 12, "intersection with -wb must have 11 rows"


def test_bedtools_intersect_wao():
    inter = run_bedtools_intersect("-wao")
    assert inter.shape[1] == 13, "intersection with -wao must have 9 rows"
