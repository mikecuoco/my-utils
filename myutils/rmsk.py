import logging
from math import log

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# From https://github.com/MarioniLab/CELLOseq/blob/93ce4f3b86014df9348abab3fc2bcaa9620c3ff4/manuscript/python_files/calculate_TE_JCage.py
# To calculate the TE age in million years of age (mya), we multiply the JC distance using the following formula:
# (JC_distance * 100) / (subsitution_rate * 2 * 100) * 1000
# For subsitution rate, we used 2.2 and 4.5 for human and mouse, according to Lander et al., 2001 and Waterston et al., 2002, respectively.
# Jukes-Cantor evolutionary distance https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)
# made it vectorized
def calculate_age(milli_div, subsitution_rate=2.2):
    """
    Vectorized calculation of evolutionary age (million years).
    """
    p = milli_div / 1000  # Proportion of sites differing
    p_part = (4 / 3) * p
    # Handle potential issues with log(0) by masking invalid values
    jc_dist = np.where(p_part < 1, -0.75 * np.log(1 - p_part), np.nan)
    mya = (jc_dist * 100) / (subsitution_rate * 2 * 100) * 1000
    return mya


def has_promoter(df, L1_cutoff=125, Alu_cutoff=4, SVA_cutoff=50):
    """
    Vectorized promoter status calculation.
    only works for human
    """
    conditions = [
        (df["repFamily"] == "L1")
        & (df["strand"] == "+")
        & (df["repStart"] <= L1_cutoff),
        (df["repFamily"] == "L1")
        & (df["strand"] == "-")
        & (df["repLeft"] <= L1_cutoff),
        (df["repFamily"] == "Alu")
        & (df["strand"] == "+")
        & (df["repStart"] <= Alu_cutoff),
        (df["repFamily"] == "Alu")
        & (df["strand"] == "-")
        & (df["repLeft"] <= Alu_cutoff),
        (df["repFamily"] == "SVA")
        & (df["strand"] == "+")
        & (df["repStart"] <= SVA_cutoff),
        (df["repFamily"] == "SVA")
        & (df["strand"] == "-")
        & (df["repLeft"] <= SVA_cutoff),
    ]
    return np.select(conditions, [True] * len(conditions), default=False)


def is_full_length(df, L1_cutoff=6000, Alu_cutoff=267, SVA_cutoff=1336):
    """
    Vectorized full-length status calculation.
    only works for human
    """
    conditions = [
        (df["repFamily"] == "L1") & (df["repEnd"] >= L1_cutoff),
        (df["repFamily"] == "Alu") & (df["repEnd"] >= Alu_cutoff),
        (df["repFamily"] == "SVA") & (df["repEnd"] >= SVA_cutoff),
    ]
    return np.select(conditions, [True] * len(conditions), default=False)


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

    logger.info(f"Reading RepeatMasker file: {filename}")

    # setup converter functions
    strand_conv = lambda x: "-" if x == "C" else "+"
    coord_conv = lambda x: int(x.rstrip(")").lstrip("("))
    perc_conv = lambda x: float(x) * 10

    # read the rmsk file
    df = pd.read_csv(
        filename,
        skiprows=3,
        sep="\s+",
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
        converters={
            "milliDiv": perc_conv,
            "milliDel": perc_conv,
            "milliIns": perc_conv,
            "genoLeft": coord_conv,
            "strand": strand_conv,
            "repStart": coord_conv,
            "repLeft": coord_conv,
        },
        engine="python",
        on_bad_lines=lambda x: x[:-1],
    )

    # split repClassFamily into repClass and repFamily on /
    if any([True for x in df["repClassFamily"].values if "/" in x]):
        df[["repClass", "repFamily"]] = df["repClassFamily"].str.split("/", expand=True)
    else:
        df["repClass"] = df["repClassFamily"]
        df["repFamily"] = None
    df.drop("repClassFamily", axis=1, inplace=True)

    # calculate length of each repeat
    logger.info("Calculating length of each repeat")
    df["length"] = np.where(
        df["strand"] == "+", df["repEnd"] - df["repStart"], df["repEnd"] - df["repLeft"]
    )

    # calculate age of each repeat
    logger.info("Calculating evolutionary age of each repeat")
    df["age"] = calculate_age(df["milliDiv"])

    # calculate promoter status
    logger.info(
        "Calculating full-length and promoter status of L1, Alu, and SVA repeats"
    )
    df["has_promoter"] = has_promoter(df)
    df["has_promoter"] = df["has_promoter"].astype(bool)
    df["is_full_length"] = is_full_length(df)
    df["is_full_length"] = df["is_full_length"].astype(bool)

    return df
