import os
from collections import defaultdict
from io import StringIO
from subprocess import DEVNULL, PIPE, CalledProcessError, Popen
from tempfile import NamedTemporaryFile

import pandas as pd
import pyranges as pr


def bedtools_intersect(a, b, extra=""):
    # check if bedtools is installed
    try:
        p = Popen("bedtools --version", shell=True, stdout=DEVNULL)
        p.wait()
    except CalledProcessError:
        raise ValueError("bedtools not installed")

    # recursively convert input to dictionary of dataframes/tempfiles
    def make_input(arg):
        if type(arg) == pd.DataFrame:
            f = NamedTemporaryFile()
            pr.PyRanges(arg).to_bed(f.name)
            return {"name": f.name, "file": f}
        elif type(arg) == str:
            assert os.path.exists(arg), f"{arg} does not exist"
            return {"name": arg}
        elif type(arg) == list:
            out = defaultdict(list)
            for a in arg:
                o = make_input(a)
                for i in o.keys():
                    out[i].append(o[i])
            return out
        else:
            raise ValueError(f"Invalid input type {type(arg)}")

    # if input is a dataframe, convert to bed and save to temporary file
    files = defaultdict(dict)
    for id, arg in zip(["a", "b"], [a, b]):
        files[id] = make_input(arg)

    # run bedtools intersect, collect output into dataframe
    cmd = f"bedtools intersect -a {files['a']['name']} "
    if type(files["b"]["name"]) == list:
        for f in files["b"]["name"]:
            cmd += f"-b {f} "
    else:
        cmd += f"-b {files['b']['name']} "

    cmd += extra
    print(f"Running command: {cmd}")
    p = Popen(cmd, shell=True, stdout=PIPE)
    with StringIO(p.stdout.read().decode()) as bed:
        df = pd.read_csv(bed, sep="\t", header=None)

    # close temporary files if they were created
    for id in files:
        if "file" in files[id].keys():
            if type(files[id]["file"]) == list:
                for f in files[id]["file"]:
                    f.close()
            else:
                files[id]["file"].close()

    return df
