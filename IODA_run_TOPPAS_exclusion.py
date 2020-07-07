# coding: utf-8
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

def IODA_exclusion_workflow(input_mzML,ppm_error,narrow_noise_threshold,large_noise_threshold):
    #source_mzML = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/toppas_input/Blank.mzML"
    TOPPAS_Pipeline = "toppas_Exclusion_workflow.toppas"
    TOPPAS_output_folder = "toppas_output"
    TOPPAS_folder = "TOPPAS_Workflow"
    os.system('mkdir results')
    os.system('mkdir download_results')
    #large_noise = 5E5
    #narrow_noise = 1E5
    #ppm_error = 10

    #SOURCE_MZML_URL = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/toppas_input/Blank.mzML"
    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logfile('results/logfile_IODA_from_mzML_'+str(today)+'.txt')
    print('======')
    print('Starting the IODA-Exclusion workflow from a mzML file')
    print('======')
    print('Getting the mzML, please wait ...')

    if input_mzML.startswith('http'):
        if 'google' in input_mzML:
            logger.info('This is the Google Drive download link:'+str(input_mzML))
            url_id = input_mzML.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_mzML = prefixe_google_download+url_id
            bashCommand1 = "wget -r "+input_mzML+" -O input.mzML"
            cp1 = subprocess.run(bashCommand1,shell=True)
            cp1
        else:
            logger.info('This is the input file path: '+str(input_mzML))
            bashCommand2 = "wget -r "+input_mzML+" -O input.mzML"
            cp2 = subprocess.run(bashCommand2,shell=True)
            cp2
    else:
        logger.info('This is the input file path: '+str(input_mzML))
        bashCommand2 = "cp "+input_mzML+" -O input.mzML"
        cp2 = subprocess.run(bashCommand2,shell=True)
        cp2

    bashCommand3 = "cp input.mzML "+TOPPAS_folder+"/"+TOPPAS_output_folder+'/Blank.mzML'
    cp3 = subprocess.run(bashCommand3,shell=True)
    cp3

    print('======')
    print('Copying the mzML to the TOPPAS/OpenMS input folder')


    print('======')
    print('Changing variables of the TOPPAS/OpenMS workflow')

    a_file = open("TOPPAS_Workflow/toppas_Exclusion_workflow.toppas", "r")
    list_of_lines = a_file.readlines()

    # Check format for variable
    try:
        float(large_noise_threshold)
    except ValueError:
        print("== The noise level must be a float or an integer, such as 6.0e05 ==")

    try:
        float(narrow_noise_threshold)
    except ValueError:
        print("== The noise level must be a float or an integer, such as 6.0e05 ==")

    try:
        float(ppm_error)
    except ValueError:
        print("== The ppm error must be a float or an integer, such as 10 ppm ==")

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

    print('======')
    print('Initializing the TOPPAS/OpenMS workflow')

    try:
        vdisplay = Xvfb()
        vdisplay.start()
    except:
        raise

    print('======')
    print('Running the TOPPAS/OpenMS workflow, this could take several minutes, please wait ...')

    bashCommand4 = "cd "+TOPPAS_folder+" && /openms-build/bin/ExecutePipeline -in "+TOPPAS_Pipeline+" -out_dir "+TOPPAS_output_folder
    print(bash_Command4)
    try:
        cp4 = subprocess.run(bashCommand4,shell=True)
        cp4
    except:
        raise

    vdisplay.stop()

    print('======')
    print('Completed the TOPPAS/OpenMS workflow')
    print('======')
    print('Zipping up the TOPPAS/OpenMS workflow files')
    get_all_file_paths('TOPPAS_Workflow/','download_results/IODA_exclusion_OpenMS_results.zip')

    print('======')
    print('Completed zipping up the TOPPAS/OpenMS workflow output files')

    print('======')
    print('You can continue the rest of the IODA workflow')

if __name__ == "__main__":
    IODA_run_TOPPAS_exclusion(str(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))# coding: utf-8
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

def IODA_exclusion_workflow(input_mzML,ppm_error,narrow_noise_threshold,large_noise_threshold):
    #source_mzML = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/toppas_input/Blank.mzML"
    TOPPAS_Pipeline = "toppas_Exclusion_workflow.toppas"
    TOPPAS_output_folder = "toppas_output"
    TOPPAS_folder = "TOPPAS_Workflow"
    os.system('mkdir results')
    os.system('mkdir download_results')
    #large_noise = 5E5
    #narrow_noise = 1E5
    #ppm_error = 10

    #SOURCE_MZML_URL = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/toppas_input/Blank.mzML"
    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logfile('results/logfile_IODA_from_mzML_'+str(today)+'.txt')
    print('======')
    print('Starting the IODA-Exclusion workflow from a mzML file')
    print('======')
    print('Getting the mzML, please wait ...')

    if input_mzML.startswith('http'):
        if 'google' in input_mzML:
            logger.info('This is the Google Drive download link:'+str(input_mzML))
            url_id = input_mzML.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_mzML = prefixe_google_download+url_id
            bashCommand1 = "wget -r "+input_mzML+" -O input.mzML"
            cp1 = subprocess.run(bashCommand1,shell=True)
            cp1
        else:
            logger.info('This is the input file path: '+str(input_mzML))
            bashCommand2 = "wget -r "+input_mzML+" -O input.mzML"
            cp2 = subprocess.run(bashCommand2,shell=True)
            cp2
    else:
        logger.info('This is the input file path: '+str(input_mzML))
        bashCommand2 = "cp "+input_mzML+" -O input.mzML"
        cp2 = subprocess.run(bashCommand2,shell=True)
        cp2

    bashCommand3 = "cp input.mzML "+TOPPAS_folder+"/"+TOPPAS_output_folder+'/Blank.mzML'
    cp3 = subprocess.run(bashCommand3,shell=True)
    cp3

    print('======')
    print('Copying the mzML to the TOPPAS/OpenMS input folder')


    print('======')
    print('Changing variables of the TOPPAS/OpenMS workflow')

    a_file = open("TOPPAS_Workflow/toppas_Exclusion_workflow.toppas", "r")
    list_of_lines = a_file.readlines()

    # Check format for variable
    try:
        float(large_noise_threshold)
    except ValueError:
        print("== The noise level must be a float or an integer, such as 6.0e05 ==")

    try:
        float(narrow_noise_threshold)
    except ValueError:
        print("== The noise level must be a float or an integer, such as 6.0e05 ==")

    try:
        float(ppm_error)
    except ValueError:
        print("== The ppm error must be a float or an integer, such as 10 ppm ==")

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

    print('======')
    print('Initializing the TOPPAS/OpenMS workflow')

    try:
        vdisplay = Xvfb()
        vdisplay.start()
    except:
        raise

    print('======')
    print('Running the TOPPAS/OpenMS workflow, this could take several minutes, please wait ...')

    bashCommand4 = "cd "+TOPPAS_folder+" && /openms-build/bin/ExecutePipeline -in "+TOPPAS_Pipeline+" -out_dir "+TOPPAS_output_folder
    print(bash_Command4)
    try:
        cp4 = subprocess.run(bashCommand4,shell=True)
        cp4
    except:
        raise

    vdisplay.stop()

    print('======')
    print('Completed the TOPPAS/OpenMS workflow')
    print('======')
    print('Zipping up the TOPPAS/OpenMS workflow files')
    get_all_file_paths('TOPPAS_Workflow/','download_results/IODA_exclusion_OpenMS_results.zip')

    print('======')
    print('Completed zipping up the TOPPAS/OpenMS workflow output files')

    print('======')
    print('You can continue the rest of the IODA workflow')

if __name__ == "__main__":
    IODA_run_TOPPAS_exclusion(str(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
