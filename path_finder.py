import argparse
import logging
import sys

import numpy as np

import path_apex as apex
import path_baseline as baseline
import path_curve as curve

# avoid "RecursionError: maximum recursion depth exceeded in comparison"
sys.setrecursionlimit(10000)

# set up log
logger = logging.getLogger('path_finder')
logger.setLevel(level=logging.INFO)
# Handler
handler = logging.FileHandler('path_finder.log')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

parser = argparse.ArgumentParser(description="Arguments for Path Finder.")

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
    "-delay",
    type=float,
    help="delay switching from feature to next feature (apex and curve mode only)",
)

parser.add_argument(
    "-min_scan",
    type=float,
    help="minimum scan time required (apex and curve mode)",
)

parser.add_argument(
    "-max_scan",
    type=float,
    help="maximum scan time required (apex and curve mode)",
)

args = parser.parse_args()

try:
    mode = args.mode
    infile = args.input_filename
    outfile = args.outfile_name
    intensity = args.intensity
    intensity_ratio = args.intensity_ratio
    num_path = args.num_path
    isolation = args.isolation # all
    intensity_accu = args.intensity_accu # curve and apex mode
    delay = args.delay # all
    window_len = args.win_len # baseline mode
    infile_raw = args.infile_raw # curve mode
    restriction = args.restriction # curve mode
    min_scan = args.min_scan # curve and apex mode
    max_scan = args.max_scan # curve and apex mode
except:
    logger.error("error in parsing args", exc_info=sys.exc_info())
    sys.exit()

if mode == "apex":
    try:
        intensity_accu = np.exp(np.log(intensity_accu) + 2.5)
    except:
        logger.error("intensity_accu argument is not valid", exc_info=sys.exc_info)
        sys.exit()
    if args.win_len is not None:
        logger.warning("win_len should not be input for apex mode")
    if args.infile_raw is not None:
        logger.warning("infile_raw should not be input for apex mode")
    if args.restriction is not None:
        logger.warning("restriction should not be input for apex mode")

    try:
        data = apex.ReadFile(infile)
    except:
        logger.error("error in reading data", exc_info=sys.exc_info())
        sys.exit()
    logger.info("=============")
    logger.info("Apex mode begin")
    logger.info("=============")
    logger.info("File Read")
    logger.info("=============")

    try:
        data = apex.DataFilter(data, intensity, intensity_ratio)
    except:
        logger.error("error in filtering data", exc_info=sys.exc_info())
        sys.exit()
    logger.info("Begin Finding Path")
    logger.info("=============")

    try:
        paths_rt, paths_mz, paths_charge, edge_intensity_dic = apex.PathGen(
            data, intensity_accu, num_path, delay, min_scan, max_scan
        )
    except:
        logger.error("error in generating path", exc_info=sys.exc_info())
        sys.exit()
    logger.info("Paths Generated")
    logger.info("=============")

    try:
        apex.WriteFile(outfile, paths_rt, paths_mz, paths_charge, edge_intensity_dic, isolation, delay, min_scan, max_scan)
    except:
        logger.error("error in generating path", exc_info=sys.exc_info())
        sys.exit()
    logger.info("File Written")
    logger.info("=============")

if mode == "baseline":
    if args.intensity_accu is not None:
        logger.warning("intensity_accu should not be input for baseline mode")
    if args.infile_raw is not None:
        logger.warning("infile_raw should not be input for baseline mode")
    if args.restriction is not None:
        logger.warning("restriction should not be input for baseline mode")
    if args.min_scan is not None:
        logger.warning("min_scan should not be input for baseline mode")
    if args.max_scan is not None:
        logger.warning("max_scan should not be input for baseline mode")

    try:
        data = baseline.ReadFile(infile)
    except:
        logger.error("error in reading data", exc_info=sys.exc_info())
        sys.exit()
    logger.info("=============")
    logger.info("Baseline mode begin")
    logger.info("=============")
    logger.info("File Read")
    logger.info("=============")

    try:
        data = baseline.DataFilter(data, intensity, intensity_ratio)
    except:
        logger.error("error in filtering data", exc_info=sys.exc_info())
        sys.exit()
    logger.info("Begin Finding Path")
    logger.info("=============")

    try:
        path = baseline.PathGen(data, window_len, num_path, isolation, delay)
    except:
        logger.error("error in generating path", exc_info=sys.exc_info())
        sys.exit()
    logger.info("Paths Generated")
    logger.info("=============")

    try:
        baseline.WriteFile(outfile, path)
    except:
        logger.error("error in generating path", exc_info=sys.exc_info())
        sys.exit()
    logger.info("File Written")
    logger.info("=============")

if mode == "curve":
    restriction[1] = max(restriction[1], isolation)
    try:
        intensity_accu = np.exp(np.log(intensity_accu) + 2.5)
    except:
        logger.error("intensity_accu argument is not valid", exc_info=sys.exc_info)
        sys.exit()
    if args.win_len is not None:
        logger.warning("win_len should not be input for apex mode")
    logger.info("=============")
    logger.info("Curve mode begin")
    logger.info("restriction: (%.4f, %.4f)", restriction[0], restriction[1])
    indice_his = curve.PathGen(
        infile_raw,
        infile,
        intensity,
        intensity_ratio,
        intensity_accu,
        restriction,
        num_path,
        delay,
        min_scan,
        max_scan,
    )
    try:
        curve.WriteFile(outfile, indice_his, restriction, delay, isolation, min_scan, max_scan)
    except:
        logger.error("error in writing to output", exc_info=sys.exc_info())
        exit()
    logger.info("File Written")
    logger.info("=============")
