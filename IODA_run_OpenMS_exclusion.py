# Author: Louis Felix Nothias, louisfelix.nothias@gmail.com, June 2020
import os
import subprocess
import sys
from xvfbwrapper import Xvfb
from logzero import logger, logfile
import datetime
import zipfile
from datetime import date
from IODA_exclusion_workflow import get_all_file_paths
from subprocess import call

def IODA_exclusion_workflow(input_mzML,ppm_error,narrow_noise_threshold,large_noise_threshold):
    #source_mzML = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/toppas_input/Blank.mzML"
    logfile('TOPPAS_Workflow/logfile_IODA_OpenMS_from_mzML.txt')

    TOPPAS_Pipeline = "toppas_Exclusion_workflow.toppas"
    TOPPAS_output_folder = "toppas_output"
    TOPPAS_folder = "TOPPAS_Workflow"
    os.system('rm download_results/IODA_OpenMS_results.zip')
    os.system('rm -r TOPPAS_Workflow/toppas_output/TOPPAS_out/')
    os.system('mkdir download_results')
    #large_noise = 5E5
    #narrow_noise = 1E5
    #ppm_error = 10

    #SOURCE_MZML_URL = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/toppas_input/Blank.mzML"
    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logger.info('STARTING the IODA-exclusion WORKFLOW with OpenMS')
    logger.info('======')
    logger.info('Getting the mzML, please wait ...')

    logger.info('This is the input: '+input_mzML)
    if input_mzML.startswith('http'):
        if 'google' in input_mzML:
            logger.info('This is the Google Drive download link:'+str(input_mzML))
            url_id = input_mzML.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_mzML = prefixe_google_download+url_id
            bashCommand1 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+"/toppas_input/Blank.mzML"
            cp1 = subprocess.run(bashCommand1,shell=True)
            cp1
        else:
            logger.info('This is the input file path: '+str(input_mzML))
            bashCommand2 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+"/toppas_input/Blank.mzML"
            cp2 = subprocess.run(bashCommand2,shell=True)
            cp2
    else:
        logger.info('This is the input file path: '+str(input_mzML))
        bashCommand3 = "cp "+input_mzML+" "+TOPPAS_folder+"/toppas_input/Blank.mzML"
        cp3 = subprocess.run(bashCommand3,shell=True)
        cp3

    logger.info('Copying the mzML to the OpenMS input folder')


    logger.info('======')
    logger.info('Changing variables of the OpenMS workflow')
    logger.info('   ppm error = '+str(ppm_error))
    logger.info('   narrow peak/feature noise threshold = '+str(narrow_noise_threshold))
    logger.info('   large peak/feature noise_threshold = '+str(large_noise_threshold))

    try:
        bashCommand0 = "wget https://github.com/lfnothias/IODA_MS/raw/test2/TOPPAS_Workflow/toppas_Exclusion_workflow.toppas -O TOPPAS_Workflow/toppas_Exclusion_workflow.toppas"
        cp0 = subprocess.run(bashCommand0,shell=True)
        cp0
        a_file = open("TOPPAS_Workflow/toppas_Exclusion_workflow.toppas", "r")
        list_of_lines = a_file.readlines()
    except:
        raise

    # Check format for variable
    try:
        float(large_noise_threshold)
    except ValueError:
        logger.info("== The noise level must be a float or an integer, such as 6.0e05 ==")

    try:
        float(narrow_noise_threshold)
    except ValueError:
        logger.info("== The noise level must be a float or an integer, such as 6.0e05 ==")

    try:
        float(ppm_error)
    except ValueError:
        logger.info("== The ppm error must be a float or an integer, such as 10 ppm ==")

    # Make string object for the noise level line FFM
    noise_line = '''            <ITEM name="noise_threshold_int" value="NOISE" type="double" description="Intensity threshold below which peaks are regarded as noise." required="false" advanced="false" />'''
    # Replace noise level for large features FFM
    list_of_lines[37] = noise_line.replace('NOISE',str(large_noise_threshold))
    # Replace noise level for narrow features FFM
    list_of_lines[95] = noise_line.replace('NOISE',str(narrow_noise_threshold))

    # Make string object for ppm error FFM
    ppm_line = '''            <ITEM name="mass_error_ppm" value="PPM_ERROR" type="double" description="Allowed mass deviation (in ppm)." required="false" advanced="false" />'''
    # Replace ppm error for large features FFM
    list_of_lines[42] = ppm_line.replace('PPM_ERROR',str(ppm_error))
    # Replace ppm error for narrow features FFM
    list_of_lines[100] = ppm_line.replace('PPM_ERROR',str(ppm_error))

    # Write out the file
    a_file = open("TOPPAS_Workflow/toppas_Exclusion_workflow.toppas", "w")
    a_file.writelines(list_of_lines)
    a_file.close()

    logger.info('======')
    logger.info('Initializing the OpenMS workflow')

    try:
        vdisplay = Xvfb()
        vdisplay.start()
    except:
        raise

    logger.info('======')
    logger.info('Running the OpenMS workflow, this usually takes less than a minute, please wait ...')

    bashCommand4 = "cd "+TOPPAS_folder+" && /openms-build/bin/ExecutePipeline -in "+TOPPAS_Pipeline+" -out_dir "+TOPPAS_output_folder
    try:
        cp4 = subprocess.run(bashCommand4,shell=True)
        cp4
    except CalledProcessError as e:
        logger.info("!!! There was an error with OpenMS workflow, please check your input files and parameters !!!")
        logger.info(e.output)
        raise

    vdisplay.stop()

    logger.info('======')
    logger.info('Completed the OpenMS workflow')
    logger.info('======')
    logger.info('Zipping up the OpenMS workflow results ...')
    get_all_file_paths('TOPPAS_Workflow/','download_results/IODA_OpenMS_results.zip')

    logger.info('======')
    logger.info('NOW CONTINUE WITH THE REST OF THE IODA-exclusion WORKFLOW')

if __name__ == "__main__":
    IODA_run_TOPPAS_exclusion(str(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
