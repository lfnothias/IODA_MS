import pandas as pd
import sys
from logzero import logger, logfile
import datetime
from subprocess import call

sys.path.insert(0, "..")

from IODA_run_OpenMS_targeted import *
from IODA_targeted_workflow import *

### IODA-targeted-from-mzTabs
def test_targeted_from_mztab_googledrive():
    input_filename = 'https://drive.google.com/file/d/1NGVzhrw-xZ4nMJserIQ7v4tYgcmraZ6g/view?usp=sharing'
    min_ratio_value = 5
    min_intensity_value = 1E5
    pretarget_rt_exclusion_time = 0.5
    posttarget_rt_exclusion_time = 3
    experiment_number = 3
    make_targeted_list_from_mzTab(input_filename, experiment_number, min_ratio_value, min_intensity_value, pretarget_rt_exclusion_time,posttarget_rt_exclusion_time,30)

def test_targeted_from_mztab_local():
    input_filename = '../tests/Euphorbia/Targeted/ioda_input/Euphorbia_rogers_latex_Blank_MS1_2uL_to_Euphorbia_rogers_latex_latex_MS1_2uL_mrgd.mzTab'
    min_ratio_value = 5
    min_intensity_value = 1E5
    pretarget_rt_exclusion_time = 0.5
    posttarget_rt_exclusion_time = 3
    experiment_number = 3
    make_targeted_list_from_mzTab(input_filename, experiment_number, min_ratio_value, min_intensity_value, pretarget_rt_exclusion_time,posttarget_rt_exclusion_time,30)
