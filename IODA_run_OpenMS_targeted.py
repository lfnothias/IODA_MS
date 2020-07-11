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
    #source_mzML1 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
    #source_mzML2 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
    TOPPAS_Pipeline = "TOPPAS_Create_mzTab_for_Targeted_list_QE_positive.toppas"
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
    print('Starting the IODA-Exclusion workflow from a mzML file')
    print('======')
    print('Getting the mzML, please wait ...')

    logger.info('These are the path to the input files:)
    logger.info('Blank: '+blank_mzML)
    logger.info('Blank: '+sample_mzML)

    # Collect the mzML and copy
    def download_copy_mzML(input_mzML_path, name_mzML):
        if input_mzML.startswith('http'):
            if 'google' in input_mzML:
                logger.info('This is the Google Drive download link:'+str(input_mzML))
                url_id = input_mzML.split('/', 10)[5]
                prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
                input_mzML = prefixe_google_download+url_id
                bashCommand1 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+'/'+TOPPAS_output_folder+'/'+name_mzML
                cp1 = subprocess.run(bashCommand1,shell=True)
                cp1
            else:
                logger.info('The Google Drive file path is invalid: '+str(input_mzML))
                raise
        else:
            logger.info('This is the input file path: '+str(input_mzML))
            bashCommand2 = "cp "+input_mzML+" -O "+TOPPAS_folder+'/'+TOPPAS_output_folder+'/'+name_mzML
            cp2 = subprocess.run(bashCommand2,shell=True)
            cp2

    # Run the function for the two input smamples
    download_copy_mzML(blank_mzML, 'Blank.mzML' )
    download_copy_mzML(sample_mzML, 'Sample.mzML')

    # ==> Check the path are correct in the TOPPAS_Workflow
    # ==> Check the workflow location and name_mzML
    # ===> Check folder consistency OpenMS TOPPAS_out
    # Hardcore sample name TOPPAS_Workflow
    
    print('======')
    print('Changing variables of the OpenMS workflow')
    logger.info('   ppm error = '+str(ppm_error))
    logger.info('   narrow peak/feature noise threshold = '+str(narrow_noise_threshold))
    logger.info('   large peak/feature noise_threshold = '+str(large_noise_threshold))

    try:
        bashCommand0 = "wget https://github.com/lfnothias/IODA_MS/raw/master/TOPPAS_Workflow/TOPPAS_Create_mzTab_for_Targeted_list_QE_positive.toppas -O TOPPAS_Workflow/toppas_Exclusion_workflow.toppas"
        cp0 = subprocess.run(bashCommand0,shell=True)
        cp0
        a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "r")
        list_of_lines = a_file.readlines()
    except:
        raise

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
    IODA_run_TOPPAS_exclusion(str(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
