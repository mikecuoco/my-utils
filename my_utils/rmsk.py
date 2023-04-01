import pandas as pd
from math import log

def calculate_age(milli_div, subsitution_rate = 2.2):
    p = milli_div / 1000 # The milliDiv column in the `rmsk.txt` file.
    p_part = (4 / 3) * p
    jc_dist = -0.75 * (log(1 - p_part))
    (jc_dist * 100) / (subsitution_rate * 2 * 100) * 1000
    return jc_dist

def read_rmsk(filename: str):

	# read first line to check if it is a valid rmsk file
	with open(filename, "rt") as f:
		line = f.readline()
		assert line == '   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat\n', "Not a valid rmsk file"

	# setup converter functions
	strand_conv = lambda x: "-" if x == "C" else "+"
	coord_conv = lambda x: int(x.rstrip(")").lstrip("("))
	perc_conv = lambda x: float(x)*10

	convs = {
		"milliDiv": perc_conv,
		"milliDel": perc_conv,
		"milliIns": perc_conv,
		"genoLeft": coord_conv,
		"strand": strand_conv,
		"repStart": coord_conv,
		"repLeft": coord_conv,
	}

	# read the rmsk file
	df = pd.read_csv(
		filename,
		skiprows=3,
		delim_whitespace=True,
		names=["swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClassFamily", "repStart", "repEnd", "repLeft", "id"],
		converters=convs
	)

	# split repClassFamily into repClass and repFamily on /
	df[["repClass", "repFamily"]] = df["repClassFamily"].str.split("/", expand=True)
	df.drop("repClassFamily", axis=1, inplace=True)

	# calculate length of each repeat
	df["length"] = df.apply(lambda x: x["repEnd"]-x["repLeft"] if x["strand"] == "-" else x["repEnd"]-x["repStart"], axis = 1)

	# calculate age of each repeat
	df["age"] = df["milliDiv"].apply(calculate_age)

	return df 
