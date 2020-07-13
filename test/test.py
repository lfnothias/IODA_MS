import pandas as pd
import sys

sys.path.insert(0, "..")

def test_mztab_import():
    from IODA_exclusion_workflow import *

    input_filename = 'https://drive.google.com/file/d/1LYk-PKsBWl4Pv7c1TlhQwaqwkF2T6sux/view?usp=sharing'
    min_intensity = 100
    rt_exclusion_margin = 5

    make_exclusion_from_mzTab(input_filename, min_intensity, rt_exclusion_margin)
