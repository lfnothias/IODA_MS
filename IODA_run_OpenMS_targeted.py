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

def IODA_targeted_workflow(blank_mzML,sample_mzML,ppm_error,noise_threshold):
    #source_mzML1 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
    #source_mzML2 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
    TOPPAS_Pipeline = "toppas_targeted_workflow_qOrbitrap_positive.toppas"
    TOPPAS_output_folder = "toppas_output"
    TOPPAS_input_folder = "toppas_input"
    TOPPAS_folder = "TOPPAS_Workflow"
    os.system('rm download_results/IODA_OpenMS_results.zip')
    os.system('rm -r TOPPAS_Workflow/toppas_output/TOPPAS_out/')
    os.system('mkdir download_results')
    #large_noise = 5E5
    #narrow_noise = 1E5
    #ppm_error = 10

    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logfile('TOPPAS_Workflow/logfile_IODA_from_mzML_'+str(today)+'.txt')
    print('======')
    print('Starting the IODA-targeted workflow from a mzML file')
    print('======')
    print('Getting the mzML, please wait ...')

    logger.info('These are the path to the input files: ')
    logger.info('Blank: '+blank_mzML)
    logger.info('Blank: '+sample_mzML)

    # Collect the mzML and copy
    def download_copy_mzML(input_mzML, name_mzML):
        if input_mzML.startswith('http'):
            if 'google' in input_mzML:
                logger.info('This is the Google Drive download link:'+str(input_mzML))
                url_id = input_mzML.split('/', 10)[5]
                prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
                input_mzML = prefixe_google_download+url_id
                bashCommand1 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML
                print(bashCommand1)
                cp1 = subprocess.run(bashCommand1,shell=True)
                cp1
            else:
                #logger.info('The Google Drive file path is invalid: '+str(input_mzML))
                logger.info('This is the input file path: '+str(input_mzML))
                bashCommand2 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML
                cp2 = subprocess.run(bashCommand2,shell=True)
                cp2
        else:
            logger.info('This is the input file path: '+str(input_mzML))
            bashCommand2 = "cp "+input_mzML+" -O "+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML
            print(bashCommand2)
            cp2 = subprocess.run(bashCommand2,shell=True)
            cp2

    # Run the function for the two input smamples
    download_copy_mzML(blank_mzML, 'Blank.mzML' )
    download_copy_mzML(sample_mzML, 'Sample.mzML')

    print('======')
    print('Changing variables of the OpenMS workflow')
    logger.info('   ppm error = '+str(ppm_error))
    logger.info('   noise threshold = '+str(noise_threshold))

    try:
        bashCommand0 = "wget https://github.com/lfnothias/IODA_MS/raw/targeted_draft/"+TOPPAS_folder+'/'+TOPPAS_Pipeline+" -O "+TOPPAS_folder+'/'+TOPPAS_Pipeline
        print(bashCommand0)
        cp0 = subprocess.run(bashCommand0,shell=True)
        cp0
        a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "r")
        list_of_lines = a_file.readlines()
    except:
        raise

    # Check format for variable
    try:
        float(noise_threshold)
    except ValueError:
        print("== The noise level must be a float or an integer, such as 6.0e05 =")

    try:
        float(ppm_error)
    except ValueError:
        print("== The ppm error must be a float or an integer, such as 10 ppm =")

    print(list_of_lines[38])
    print(list_of_lines[43])

    list_of_lines = [sub.replace('NOISE', str(noise_threshold)) for sub in list_of_lines]
    list_of_lines = [sub.replace('PPM_ERROR', str(ppm_error)) for sub in list_of_lines]

    print(list_of_lines[38])
    print(list_of_lines[43])

    # Write out the file
    a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "w")
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
    get_all_file_paths('TOPPAS_Workflow/','download_results/IODA_OpenMS_results.zip')

    print('======')
    print('Completed zipping up the TOPPAS/OpenMS workflow output files')

    print('======')
    print('You can continue the rest of the IODA workflow')

if __name__ == "__main__":
    IODA_targeted_workflow(str(sys.argv[1]),str(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
