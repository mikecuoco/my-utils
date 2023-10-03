import logging
from math import log

import pandas as pd


# From https://github.com/MarioniLab/CELLOseq/blob/93ce4f3b86014df9348abab3fc2bcaa9620c3ff4/manuscript/python_files/calculate_TE_JCage.py
# To calculate the TE age in million years of age (mya), we multiply the JC distance using the following formula:
# (JC_distance * 100) / (subsitution_rate * 2 * 100) * 1000
# For subsitution rate, we used 2.2 and 4.5 for human and mouse, according to Lander et al., 2001 and Waterston et al., 2002, respectively.
# Jukes-Cantor evolutionary distance https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)
def calculate_age(milli_div, subsitution_rate=2.2):
    p = (
        milli_div / 1000
    )  # The milliDiv column in the `rmsk.txt` file. Used to calculate proportion of sites that differ between this seq and the consensus seq
    p_part = (4 / 3) * p
    jc_dist = -0.75 * (log(1 - p_part))
    mya = (jc_dist * 100) / (subsitution_rate * 2 * 100) * 1000
    return mya


def has_promoter(x):
    "test if rmsk element has a promoter"
    assert x.repFamily in ["L1", "Alu"], "repFamily must be L1 or Alu"

    if x.repFamily == "L1":
        thresh = 125
    elif x.repFamily == "Alu":
        thresh = 70

    if x.strand == "+":
        return x.repStart < thresh
    elif x.strand == "-":
        return x.repLeft < thresh


def read_rmsk(filename: str):
    # read first line to check if it is a valid rmsk file
    line = pd.read_csv(filename, nrows=1, header=None).values[0][0]

    for w in [
        "perc",
        "query",
        "position in",
        "matching",
        "repeat",
    ]:
        assert w in line, f"Not a valid rmsk file: {w} not found in first line"

    # setup converter functions
    strand_conv = lambda x: "-" if x == "C" else "+"
    coord_conv = lambda x: int(x.rstrip(")").lstrip("("))
    perc_conv = lambda x: float(x) * 10

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
    logging.info(f"Reading RepeatMasker file: {filename}")
    df = pd.read_csv(
        filename,
        skiprows=3,
        delim_whitespace=True,
        names=[
            "swScore",
            "milliDiv",
            "milliDel",
            "milliIns",
            "genoName",
            "genoStart",
            "genoEnd",
            "genoLeft",
            "strand",
            "repName",
            "repClassFamily",
            "repStart",
            "repEnd",
            "repLeft",
            "id",
        ],
        converters=convs,
        engine="python",
        on_bad_lines=lambda x: x[:-1],
    )

    # split repClassFamily into repClass and repFamily on /
    if any([True for x in df["repClassFamily"].values if "/" in x]):
        df[["repClass", "repFamily"]] = df["repClassFamily"].str.split("/", expand=True)
        df.drop("repClassFamily", axis=1, inplace=True)

    # calculate length of each repeat
    df["length"] = df.apply(
        lambda x: x["repEnd"] - x["repLeft"]
        if x["strand"] == "-"
        else x["repEnd"] - x["repStart"],
        axis=1,
    )

    # calculate age of each repeat
    logging.info("Calculating evolutionary age of each repeat")
    df["age"] = df["milliDiv"].apply(calculate_age)

    # calculate promoter status
    logging.info("Calculating promoter status of L1 and Alu repeats")
    if "repFamily" in df.columns:
        df.loc[df["repFamily"].isin(["Alu", "L1"]), "has_promoter"] = df[
            df["repFamily"].isin(["Alu", "L1"])
        ].apply(has_promoter, axis=1)

    return df
