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

def IODA_targeted_workflow(blank_mzML,sample_mzML,ppm_tolerance,noise_level):
    #source_mzML1 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
    #source_mzML2 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
    logfile('TOPPAS_Workflow/logfile_IODA_OpenMS_from_mzML.txt')
    TOPPAS_Pipeline = "toppas_targeted_workflow_qOrbitrap_positive.toppas"
    TOPPAS_output_folder = "toppas_output"
    TOPPAS_input_folder = "toppas_input"
    TOPPAS_folder = "TOPPAS_Workflow"
    os.system('rm download_results/IODA_OpenMS_results.zip')
    os.system('rm -r TOPPAS_Workflow/'+TOPPAS_output_folder+'/TOPPAS_out/')
    os.system('mkdir download_results')
    #large_noise = 5E5
    #narrow_noise = 1E5
    #ppm_error = 10

    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logger.info('STARTING the IODA-targeted WORKFLOW with OpenMS')
    logger.info('======')
    logger.info('Path to the input files: ')
    logger.info('Blank: '+blank_mzML)
    logger.info('Sample: '+sample_mzML)

    # Collect the mzML and copy
    def download_copy_mzML(input_mzML, name_mzML):
        if input_mzML.startswith('http'):
            logger.info('Downloading the mzML files, please wait ...')
            if 'google' in input_mzML:
                logger.info('This is the Google Drive download link:'+str(input_mzML))
                url_id = input_mzML.split('/', 10)[5]
                prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
                input_mzML = prefixe_google_download+url_id
                bashCommand1 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML
                cp1 = subprocess.run(bashCommand1,shell=True)
                cp1
            else:
                #logger.info('The Google Drive file path is invalid: '+str(input_mzML))
                logger.info('This is the input file path: '+str(input_mzML))
                bashCommand2 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML
                cp2 = subprocess.run(bashCommand2,shell=True)
                cp2
        else:
            logger.info('Path to uploaded file: '+str(input_mzML))
            bashCommand2 = "cp "+input_mzML+" "+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+input_mzML.split('/', 10)[-1]
            cp2 = subprocess.run(bashCommand2,shell=True)
            cp2

    # Run the function for the two input smamples
    logger.info('Copying the mzML files ...')
    download_copy_mzML(blank_mzML, 'Blank.mzML' )
    download_copy_mzML(sample_mzML, 'Sample.mzML')

    logger.info('======')
    logger.info('Changing the variables of the OpenMS workflow ...')
    logger.info('   ppm error = '+str(ppm_tolerance))
    logger.info('   noise threshold = '+str(noise_level))

    try:
        bashCommand0 = "wget https://github.com/lfnothias/IODA_MS/raw/targeted_draft/"+TOPPAS_folder+'/'+TOPPAS_Pipeline+" -O "+TOPPAS_folder+'/'+TOPPAS_Pipeline
        cp0 = subprocess.run(bashCommand0,shell=True)
        cp0
        a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "r")
        list_of_lines = a_file.readlines()
    except:
        raise

    # Check format for variable
    try:
        float(noise_level)
    except ValueError:
        logger.info("== The noise level must be a float or an integer, such as 6.0e05 =")

    try:
        float(ppm_tolerance)
    except ValueError:
        logger.info("== The ppm error must be a float or an integer, such as 10 ppm =")

    # Preserve the original mzML file names in the OpenMS workflow for local files
    if sample_mzML.startswith('http'):
        if 'google' in sample_mzML:
            pass
    else:
        blank_filename = str(blank_mzML.split('/', 10)[-1])
        sample_filename = str(sample_mzML.split('/', 10)[-1])
        list_of_lines = [sub.replace('LISTITEM value="toppas_input/Blank.mzML', 'LISTITEM value="'+TOPPAS_input_folder+'/'+blank_filename) for sub in list_of_lines]
        list_of_lines = [sub.replace('LISTITEM value="toppas_input/Sample.mzML', 'LISTITEM value="'+TOPPAS_input_folder+'/'+sample_filename) for sub in list_of_lines]

    # Replace OpenMS workflow parameters
    list_of_lines = [sub.replace('NOISE', str(noise_level)) for sub in list_of_lines]
    list_of_lines = [sub.replace('PPM_ERROR', str(ppm_tolerance)) for sub in list_of_lines]
    # Write out the file
    a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "w")
    a_file.writelines(list_of_lines)
    a_file.close()

    logger.info('======')
    logger.info('Initializing the TOPPAS/OpenMS workflow')

    try:
        vdisplay = Xvfb()
        vdisplay.start()
    except:
        raise

    logger.info('======')
    logger.info('Running the TOPPAS/OpenMS workflow, this usually takes less than a minute, please wait ...')
    logger.info('If this takes longer, increase the noise_level value ...')

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
    logger.info('Zipping up the OpenMS workflow results ..')
    get_all_file_paths('TOPPAS_Workflow/','download_results/IODA_OpenMS_results.zip')

    logger.info('======')
    logger.info('Completed zipping up the OpenMS workflow result files')

    logger.info('======')
    logger.info('NOW CONTINUE WITH THE REST OF THE IODA-targeted WORKFLOW')

if __name__ == "__main__":
    IODA_targeted_workflow(str(sys.argv[1]),str(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
