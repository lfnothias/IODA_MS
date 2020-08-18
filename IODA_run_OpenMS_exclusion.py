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
    # Test samples
        #source_mzML = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/toppas_input/Blank.mzML"
        #input_mzML = "tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #input_mzML = "https://drive.google.com/file/d/11p2Jau2T-gCQb9KZExWdC7dy8AQWV__l/view?usp=sharing"
        #input_mzML = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_blank_SPE_20181227092326.mzML"

    os.system('rm TOPPAS_Workflow/logfile_IODA_OpenMS_from_mzML.txt')
    logfile('TOPPAS_Workflow/logfile_IODA_OpenMS_from_mzML.txt')
    TOPPAS_Pipeline = "toppas_Exclusion_workflow.toppas"
    TOPPAS_output_folder = "toppas_output"
    TOPPAS_folder = "TOPPAS_Workflow"
    os.system('rm download_results/IODA_OpenMS_results.zip')
    os.system('rm -r TOPPAS_Workflow/toppas_input/*')
    os.system('rm -r TOPPAS_Workflow/toppas_output/TOPPAS_out/')
    os.system('mkdir download_results')

    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logger.info('STARTING the IODA-exclusion WORKFLOW with OpenMS')
    logger.info('======')
    logger.info('Getting the mzML, please wait ...')

    if input_mzML.startswith(('http','ftp')):
        if 'google' in input_mzML:
            logger.info('This is the Google Drive download link:'+str(input_mzML))
            logger.info('Downloading the mzML, please wait ...')
            url_id = input_mzML.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_mzML = prefixe_google_download+url_id
            bashCommand1 = "wget --no-check-certificate '"+input_mzML+"' -O "+TOPPAS_folder+"/toppas_input/Blank.mzML || rm -f "+TOPPAS_folder+"/toppas_input/Blank.mzML"
            cp1 = subprocess.run(bashCommand1,shell=True)
            try:
                cp1
            except subprocess.CalledProcessError:
                raise
        if 'massive.ucsd.edu' in input_mzML:
            logger.info('This is the MassIVE repository link: '+str(input_mzML))
            logger.info('Downloading the mzML, please wait ... ')
            bashCommand4 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+"/toppas_input/Blank.mzML || rm -f "+TOPPAS_folder+"/toppas_input/Blank.mzML"
            cp4 = subprocess.run(bashCommand4,shell=True)
            try:
                cp4
            except subprocess.CalledProcessError:
                raise
    else:
        #Check the file path is correct for local upload
        logger.info('This is the input file path: '+str(input_mzML))
        bashCommand3 = "cp "+input_mzML+" "+TOPPAS_folder+"/toppas_input/Blank.mzML"
        cp3 = subprocess.run(bashCommand3,shell=True)
        try:
            cp3
        except subprocess.CalledProcessError:
            raise
    # Error getting the file ! PLEASE VERY THE PATH TO THE FILE OR DOWNLOAD LINK ...
    try:
        f = open(TOPPAS_folder+'/toppas_input/Blank.mzML')
        f.close()
    except subprocess.CalledProcessError:
        logger.info('There was an error getting the file !')
    logger.info('The mzML file was found')

    logger.info('Copying the mzML to the OpenMS input folder. File will be renamed internally "Blank.mzML"')

    logger.info('======')
    logger.info('Changing variables of the OpenMS workflow')
    logger.info('   ppm error = '+str(ppm_error))
    logger.info('   narrow peak/feature noise threshold = '+str(narrow_noise_threshold))
    logger.info('   large peak/feature noise_threshold = '+str(large_noise_threshold))

    try:
        bashCommand0 = "wget https://github.com/lfnothias/IODA_MS/raw/master/TOPPAS_Workflow/toppas_Exclusion_workflow.toppas -O TOPPAS_Workflow/toppas_Exclusion_workflow.toppas"
        cp0 = subprocess.run(bashCommand0,shell=True)
        cp0
        a_file = open("TOPPAS_Workflow/toppas_Exclusion_workflow.toppas", "r")
        list_of_lines = a_file.readlines()
    except:
        raise

    # Check format for variable

    #"== The noise level must be a float or an integer, such as 6.0e05 =="
    try:
        float(large_noise_threshold)
    except subprocess.CalledProcessError:
        logger.info("== The noise level must be a float or an integer, such as 6.0e05 ==")

    #"== The noise level must be a float or an integer, such as 6.0e05 =="
    try:
        float(narrow_noise_threshold)
    except subprocess.CalledProcessError:
        logger.info("== The noise level must be a float or an integer, such as 6.0e05 ==")
    #"== The ppm error must be a float or an integer, such as 10 ppm =="
    try:
        float(ppm_error)
    except subprocess.CalledProcessError:
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
    except subprocess.CalledProcessError:
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

    # Error with the OpenMS workflow. No output files.
    try:
        f = open(TOPPAS_folder+'/toppas_output/TOPPAS_out/mzTab_Narrow/Blank.mzTab')
        f.close()
    except:
        logger.info('There was an issue with the OpenMS workflow ! See the log below.')
        f = open(TOPPAS_folder+'/toppas_output/TOPPAS.log', 'r')
        file_contents = f.read()
        logger.info(file_contents)
        raise

    logger.info('======')
    logger.info('Completed the OpenMS workflow')
    logger.info('======')
    logger.info('Zipping up the OpenMS workflow results ...')
    get_all_file_paths('TOPPAS_Workflow/','download_results/IODA_OpenMS_results.zip')

    logger.info('======')
    logger.info('NOW YOU CAN CONTINUE WITH THE REST OF THE WORKFLOW')

if __name__ == "__main__":
    IODA_run_TOPPAS_exclusion(str(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
