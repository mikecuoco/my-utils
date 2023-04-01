import pytest
from my_utils.rmsk import read_rmsk

def test_read_rmsk():
	df = read_rmsk("tests/rmsk.out")
	assert df.shape[1] == 18, "Wrong number of columns"
	assert "repClass" in df.columns, "repClass column not found"
	assert "repFamily" in df.columns, "repFamily column not found"
	assert "length" in df.columns, "length column not found"
	assert "age" in df.columns, "age column not found"

def test_ucsc_rmsk():
	with pytest.raises(AssertionError) as e:
		df = read_rmsk("tests/rmsk.ucsc.txt")
