import pandas as pd
import sys

sys.path.insert(0, "..")

from IODA_run_OpenMS_targeted import *
from IODA_targeted_workflow import *

### IODA-targeted-from-mzTabs
def test_targeted_from_mztab_googledrive():
    input_filename = 'https://drive.google.com/file/d/1NGVzhrw-xZ4nMJserIQ7v4tYgcmraZ6g/view?usp=sharing'
    ratio_value = 5
    min_intensity_value = 1E5
    experiment_number = 3

    make_targeted_list_from_mzTab(input_filename, experiment_number, ratio_value, min_intensity_value)

def test_targeted_from_mztab_local():
    input_filename = 'tests/Euphorbia/Targeted/ioda_input/Euphorbia_rogers_latex_Blank_MS1_2uL_to_Euphorbia_rogers_latex_latex_MS1_2uL_mrgd.mzTab'
    ratio_value = 5
    min_intensity_value = 1E5
    experiment_number = 3

    make_targeted_list_from_mzTab(input_filename, experiment_number, ratio_value, min_intensity_value)

### IODA-targeted-from-mzML
def test_targeted_from_mzml_local():
    input_BLANK = "tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
    input_SAMPLE = "tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
    ppm_error = 10
    feature_noise = 1E6

    IODA_targeted_workflow(input_BLANK,input_SAMPLE,ppm_error,feature_noise)

    ratio_value = 5
    min_intensity_value = 1E5
    experiment_number = 3

    make_targeted_list_from_mzTab('OpenMS_generated', experiment_number, ratio_value, min_intensity_value)

def test_targeted_from_mzml_google_drive():
    input_BLANK = "https://drive.google.com/file/d/11p2Jau2T-gCQb9KZExWdC7dy8AQWV__l/view?usp=sharing"
    input_SAMPLE = "https://drive.google.com/file/d/1_lOYEtsmEPAlfGVYbzJpLePPSitUp1yh/view?usp=sharing"
    ppm_error = 10
    feature_noise = 1E6

    IODA_targeted_workflow(input_BLANK,input_SAMPLE,ppm_error,feature_noise)

    ratio_value = 5
    min_intensity_value = 1E5
    experiment_number = 3

    make_targeted_list_from_mzTab('OpenMS_generated', experiment_number, ratio_value, min_intensity_value)

def test_targeted_from_mzml_massive():
    input_BLANK = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_blank_SPE_20181227092326.mzML"
    input_SAMPLE = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_F1-1_F2-1_NIST-1_To-1_20181227135238.mzML"
    ppm_error = 10
    feature_noise = 1E6

    IODA_targeted_workflow(input_BLANK,input_SAMPLE,ppm_error,feature_noise)

    ratio_value = 5
    min_intensity_value = 1E5
    experiment_number = 3

    make_targeted_list_from_mzTab('OpenMS_generated', experiment_number, ratio_value, min_intensity_value)
