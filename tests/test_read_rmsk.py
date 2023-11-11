import pytest

from myutils.rmsk import read_rmsk


def read(file):
    df = read_rmsk(file)
    assert "length" in df.columns, "length column not found"
    assert "age" in df.columns, "age column not found"
    assert "is_full_length" in df.columns, "is_full_length column not found"
    assert "has_promoter" in df.columns, "has_promoter column not found"


def test_ucsc_rmsk():
    with pytest.raises(AssertionError) as e:
        df = read_rmsk("tests/data/rmsk.ucsc.txt")


def test_read_rmsk():
    read("tests/data/rmsk.out")


def test_read_rmsk_gz():
    read("tests/data/rmsk.out.gz")


def test_read_rmsk_astrk():
    read("tests/data/rmsk_astrk.out")


def test_read_rmsk_url():
    # use fruit fly as a test case
    read("https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.out.gz")
