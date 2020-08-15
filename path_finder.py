import argparse
import sys

import numpy as np

import path_apex as apex
import path_baseline as baseline
import path_curve as curve

# avoid "RecursionError: maximum recursion depth exceeded in comparison"
sys.setrecursionlimit(10000)

parser = argparse.ArgumentParser(description="Command line for Path Finder.")

parser.add_argument(
    "mode",
    type=str,
    help="Mode of path finder, baseline, apex and curve (dev), default is baseline",
)

parser.add_argument(
    "input_filename",
    type=str,
    help="Feature table generated from .mzTab, before feature filtering",
)

parser.add_argument(
    "outfile_name", type=str, help="Output file storing path information, could be .txt"
)

parser.add_argument(
    "intensity", type=float, help="intensity cut-off for feature filtering"
)

parser.add_argument(
    "intensity_ratio", type=float, help="intensity_ratio cut-off for filtering"
)

parser.add_argument("num_path", type=int, help="number of paths required")

parser.add_argument(
    "-infile_raw", type=str, help="Raw .mzTab file with only samples (curve mode only)",
)

parser.add_argument(
    "-intensity_accu",
    type=float,
    help="minimum requirement for feature intensity accumulation in a time range (apex and curve mode only)",
)

parser.add_argument(
    "-win_len", type=float, help="window length of baseline method (baseline mode only)"
)

parser.add_argument(
    "-restriction",
    type=float,
    nargs=2,
    help="restriction grid for clustering (curve mode only)",
)

parser.add_argument(
    "-isolation",
    type=float,
    help="isolation window in mz (baseline and apex mode only)",
)

parser.add_argument(
    "-delta",
    type=float,
    help="delay switching from feature to next feature (apex and curve mode only)",
)

args = parser.parse_args()

mode = args.mode
infile = args.input_filename
outfile = args.outfile_name
intensity = args.intensity
intensity_ratio = args.intensity_ratio
num_path = args.num_path

if mode == "apex":
    isolation = args.isolation
    intensity_accu = args.intensity_accu
    intensity_accu = np.exp(np.log(intensity_accu) + 2.5)
    delta = args.delta
    data = apex.ReadFile(infile)
    print("=============")
    print("Apex mode begin")
    print("=============")
    print("File Read")
    print("=============")

    data = apex.DataFilter(data, intensity, intensity_ratio)
    print("Begin Finding Path")
    print("=============")

    paths_rt, paths_mz, paths_charge, edge_intensity_dic = apex.PathGen(
        data, intensity_accu, num_path, delta
    )
    print("Paths Generated")
    print("=============")

    apex.WriteFile(outfile, paths_rt, paths_mz, paths_charge, edge_intensity_dic, isolation, delta)
    print("File Written")
    print("=============")

if mode == "baseline":
    isolation = args.isolation
    window_len = args.win_len
    data = baseline.ReadFile(infile)
    print("=============")
    print("Baseline mode begin")
    print("=============")
    print("File Read")
    print("=============")

    data = baseline.DataFilter(data, intensity, intensity_ratio)
    print("Begin Finding Path")
    print("=============")

    path = baseline.PathGen(data, window_len, num_path, isolation)
    print("Paths Generated")
    print("=============")

    baseline.WriteFile(outfile, path)
    print("File Written")
    print("=============")

if mode == "curve":
    intensity_accu = args.intensity_accu
    intensity_accu = np.exp(np.log(intensity_accu) + 2.5)
    infile_raw = args.infile_raw
    restriction = args.restriction
    delta = args.delta
    print("=============")
    print("Curve mode begin")
    indice_his = curve.PathGen(
        infile_raw,
        infile,
        intensity,
        intensity_ratio,
        intensity_accu,
        restriction,
        num_path,
        delta,
    )
    curve.WriteFile(outfile, indice_his, restriction, delta)
    print("File Written")
    print("=============")
