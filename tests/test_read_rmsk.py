import pytest

from myutils.rmsk import read_rmsk


def read(file):
    df = read_rmsk(file)
    assert df.shape[1] == 18, "Wrong number of columns"
    assert "repClass" in df.columns, "repClass column not found"
    assert "repFamily" in df.columns, "repFamily column not found"
    assert "length" in df.columns, "length column not found"
    assert "age" in df.columns, "age column not found"


def test_ucsc_rmsk():
    with pytest.raises(AssertionError) as e:
        df = read_rmsk("tests/rmsk.ucsc.txt")


def test_read_rmsk():
    read("tests/rmsk.out")
    read("tests/rmsk.out.gz")
    read("tests/rmsk_astrk.out")
