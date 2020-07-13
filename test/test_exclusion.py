import pandas as pd
import sys
from logzero import logger, logfile
import datetime
from subprocess import call

sys.path.insert(0, "..")

from IODA_run_OpenMS_exclusion import *
from IODA_exclusion_workflow import *

### IODA-exclusion-from-mzTab
def test_exclusion_from_mztab_googledrive():
    input_filename = 'https://drive.google.com/file/d/1LYk-PKsBWl4Pv7c1TlhQwaqwkF2T6sux/view?usp=sharing'
    min_intensity = 100
    rt_exclusion_margin = 5

    make_exclusion_from_mzTab(input_filename, min_intensity, rt_exclusion_margin)

def test_exclusion_from_mztab():

    input_filename = '../tests/Euphorbia/exclusion/ioda_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzTab'
    min_intensity = 100
    rt_exclusion_margin = 5

    make_exclusion_from_mzTab(input_filename, min_intensity, rt_exclusion_margin)
