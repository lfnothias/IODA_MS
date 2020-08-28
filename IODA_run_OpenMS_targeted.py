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

def IODA_targeted_workflow(blank_mzML:str,sample_mzML:str,ppm_tolerance:float,noise_level:float):
    # Test samples
        #source_mzML1 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #source_mzML2 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
        #input_BLANK = "tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #input_SAMPLE = "tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
        #input_BLANK = "https://drive.google.com/file/d/11p2Jau2T-gCQb9KZExWdC7dy8AQWV__l/view?usp=sharing"
        #input_SAMPLE = "https://drive.google.com/file/d/1_lOYEtsmEPAlfGVYbzJpLePPSitUp1yh/view?usp=sharing"
        #input_BLANK = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_blank_SPE_20181227092326.mzML"
        #input_SAMPLE = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_F1-1_F2-1_NIST-1_To-1_20181227135238.mzML"

    os.system('rm TOPPAS_Workflow/logfile_IODA_OpenMS_from_mzML.txt')
    logfile('TOPPAS_Workflow/logfile_IODA_OpenMS_from_mzML.txt')
    TOPPAS_Pipeline = "toppas_targeted_workflow_qOrbitrap_positive.toppas"
    TOPPAS_output_folder = "toppas_output"
    TOPPAS_input_folder = "toppas_input"
    TOPPAS_folder = "TOPPAS_Workflow"
    os.system('rm download_results/IODA_OpenMS_results.zip')
    os.system('rm -r TOPPAS_Workflow/toppas_input/*')
    os.system('rm -r TOPPAS_Workflow/'+TOPPAS_output_folder+'/TOPPAS_out/')
    os.system('mkdir download_results')
    #large_noise = 5E5
    #narrow_noise = 1E5
    #ppm_error = 10

    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logger.info('STARTING the IODA-targeted WORKFLOW')
    logger.info('======')
    logger.info('Path to the input files: ')
    logger.info('Blank: '+blank_mzML)
    logger.info('Sample: '+sample_mzML)

    # Collect the mzML and copy
    def download_copy_mzML(input_mzML, name_mzML):
        if input_mzML.startswith(('http','ftp')):
            logger.info('Downloading the mzML files, please wait ...')
            if 'google' in input_mzML:
                logger.info('This is the Google Drive download link:'+str(input_mzML))
                logger.info('Downloading ...')
                url_id = input_mzML.split('/', 10)[5]
                prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
                input_mzML = prefixe_google_download+url_id
                bashCommand1 = "wget --no-check-certificate '"+input_mzML+"' -O "+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML+' || rm -f '+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML
                cp1 = subprocess.run(bashCommand1,shell=True)
                cp1

                # Error getting the file ! PLEASE, VERY THE DOWNLOAD LINK.
                try:
                    f = open(TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML)
                    f.close()
                except subprocess.CalledProcessError:
                    logger.info('There was an error downloading the mzML file(s) ! PLEASE, VERY THE PATH OR DOWNLOAD LINK')
                logger.info('The mzML file was download succesfully')

            if 'massive.ucsd.edu' in input_mzML:
                logger.info('This is the MassIVE repository link: '+str(input_mzML))
                logger.info('Downloading ... ')
                bashCommand4 = "wget -r "+input_mzML+" -O "+TOPPAS_folder+"/toppas_input/Blank.mzML || rm -f "+TOPPAS_folder+'/toppas_input/Blank.mzML'
                cp4 = subprocess.run(bashCommand4,shell=True)
                cp4
                # Error getting the file ! PLEASE, VERY THE DOWNLOAD LINK.
                try:
                    f = open(TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+name_mzML)
                    f.close()
                except subprocess.CalledProcessError:
                    logger.info('There was an error downloading the mzML file(s) ! PLEASE, VERY THE PATH OR DOWNLOAD LINK')
                logger.info('The mzML file was download succesfully')
        else:
            logger.info('Path to uploaded file: '+str(input_mzML))
            bashCommand2 = "cp "+input_mzML+" "+TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+input_mzML.split('/', 10)[-1]
            cp2 = subprocess.run(bashCommand2,shell=True)
            cp2
            # Error getting the file ! PLEASE, VERY THE FILE PATH OR DOWNLOAD LINK')
            try:
                f = open(TOPPAS_folder+'/'+TOPPAS_input_folder+'/'+input_mzML.split('/', 10)[-1])
                f.close()
            except subprocess.CalledProcessError:
                logger.info('There was an error getting the mzML file(s) ! PLEASE, VERY THE PATH OR DOWNLOAD LINK')
            logger.info('The mzML file was found')

    # Downloading the two input mzML
    logger.info('Copying the mzML files ...')
    download_copy_mzML(blank_mzML, 'Blank.mzML' )
    download_copy_mzML(sample_mzML, 'Sample.mzML')

    logger.info('======')
    logger.info('Changing the variables of the OpenMS workflow ...')
    logger.info('   ppm error = '+str(ppm_tolerance))
    logger.info('   noise threshold = '+str("{:.2e}".format(noise_level)))

    try:
        bashCommand0 = "wget https://github.com/lfnothias/IODA_MS/raw/MS2Planner_merge_w_master/"+TOPPAS_folder+'/'+TOPPAS_Pipeline+" -O "+TOPPAS_folder+'/'+TOPPAS_Pipeline
        cp0 = subprocess.run(bashCommand0,shell=True)
        cp0
        a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "r")
        list_of_lines = a_file.readlines()
    except subprocess.CalledProcessError:
        logger.info('ERROR getting the reference workflow')

    # Check format for variable
    #"== The noise level must be a float or an integer, such as 6.0e05 ="
    try:
        float(noise_level)
    except subprocess.CalledProcessError:
        logger.info("== The noise level must be a float or an integer, such as 6.0e05 =")

    #"== The ppm error must be a float or an integer, such as 10 ppm ="
    try:
        float(ppm_tolerance)
    except subprocess.CalledProcessError:
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
    list_of_lines = [sub.replace('"1E5"', str('"'+str(noise_level)+'"')) for sub in list_of_lines]
    list_of_lines = [sub.replace('"11"', str('"'+str(ppm_tolerance)+'"')) for sub in list_of_lines]
    # Write out the file
    a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "w")
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
    logger.info('If this takes longer, increase the noise treshold value ...')
    bashCommand4 = "cd "+TOPPAS_folder+" && /openms-build/bin/ExecutePipeline -in "+TOPPAS_Pipeline+" -out_dir "+TOPPAS_output_folder
    try:
        cp4 = subprocess.run(bashCommand4,shell=True)
        cp4
    except:
        logger.info('There was an issue with the OpenMS workflow ! See the log below.')
        f = open(TOPPAS_folder+'/toppas_output/TOPPAS.log', 'r')
        file_contents = f.read()
        logger.info(file_contents)
        raise

    vdisplay.stop()

    # Error with the OpenMS workflow. No output files.
    try:
        mzTab_file = os.listdir(TOPPAS_folder+"/toppas_output/TOPPAS_out/Targeted_MzTab/")[0]
        f = open(TOPPAS_folder+'/toppas_output/TOPPAS_out/Targeted_MzTab/'+mzTab_file)
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
    logger.info('Zipping up the OpenMS workflow results ..')
    get_all_file_paths('TOPPAS_Workflow/','download_results/IODA_OpenMS_results.zip')

    logger.info('======')
    logger.info('Completed zipping up the OpenMS workflow result files')

    logger.info('======')
    logger.info('NOW CONTINUE WITH THE REST OF THE IODA-targeted WORKFLOW')


def MS2Planner_Curve_OpenMS(sample_mzML:str):
    # Test samples
        #source_mzML1 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #source_mzML2 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
        #input_BLANK = "tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #input_SAMPLE = "tests/Euphorbia/Targeted/toppas_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
        #input_BLANK = "https://drive.google.com/file/d/11p2Jau2T-gCQb9KZExWdC7dy8AQWV__l/view?usp=sharing"
        #input_SAMPLE = "https://drive.google.com/file/d/1_lOYEtsmEPAlfGVYbzJpLePPSitUp1yh/view?usp=sharing"
        #input_BLANK = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_blank_SPE_20181227092326.mzML"
        #input_SAMPLE = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_F1-1_F2-1_NIST-1_To-1_20181227135238.mzML"

    TOPPAS_Pipeline = "MS1_MS2Planner_Curve_mzTab.toppas"
    TOPPAS_output_folder = "toppas_output"
    TOPPAS_input_folder = "toppas_input"
    TOPPAS_folder = "TOPPAS_Workflow"

    os.system('rm '+TOPPAS_folder+'/logfile_MS2Planner_OpenMS.txt')
    logfile(TOPPAS_folder+'/logfile_MS2Planner_OpenMS.txt')
    os.system('rm download_results/MS2Planner_OpenMS_results.zip')
    os.system('rm -r '+TOPPAS_folder+'/'+TOPPAS_output_folder+'/TOPPAS_out/MS2Planner*')
    os.system('mkdir download_results')
    os.system('mv '+TOPPAS_folder+'/'+TOPPAS_output_folder+'/TOPPAS.log '+TOPPAS_folder+'/'+TOPPAS_output_folder+'/TOPPAS_FFM.log')
    os.system('mkdir '+TOPPAS_folder+'/'+TOPPAS_output_folder+'/TOPPAS_out/MS2Planner_mzTab')


    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logger.info('STARTING the MS2Planner Curve processing')
    logger.info('======')
    logger.info('Path to the input files: ')
    logger.info('    Sample: '+sample_mzML)

    logger.info('======')
    logger.info('Download the latest version of the workflow from the repository ...')

    try:
        bashCommand0 = "wget https://github.com/lfnothias/IODA_MS/raw/MS2Planner_merge_w_master/"+TOPPAS_folder+'/'+TOPPAS_Pipeline+" -O "+TOPPAS_folder+'/'+TOPPAS_Pipeline
        logger.info(bashCommand0)
        cp0 = subprocess.run(bashCommand0,shell=True)
        cp0
        a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "r")
        list_of_lines = a_file.readlines()
    except subprocess.CalledProcessError:
        logger.info('ERROR getting the reference workflow')

    # Preserve the original mzML file names in the OpenMS workflow for local files
    if sample_mzML.startswith('http'):
        if 'google' in sample_mzML:
            pass
    else:
        sample_filename = str(sample_mzML.split('/', 10)[-1])
        list_of_lines = [sub.replace('LISTITEM value="toppas_input/Sample.mzML', 'LISTITEM value="'+TOPPAS_input_folder+'/'+sample_filename) for sub in list_of_lines]

    # Replace OpenMS workflow parameters
    # Write out the file
    a_file = open(TOPPAS_folder+'/'+TOPPAS_Pipeline, "w")
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
    logger.info('Running the OpenMS workflow, this usually takes couple minutes, please wait ...')

    bashCommand4 = "cd "+TOPPAS_folder+" && /openms-build/bin/ExecutePipeline -in "+TOPPAS_Pipeline+" -out_dir "+TOPPAS_output_folder
    try:
        cp4 = subprocess.run(bashCommand4,shell=True)
        cp4
    except subprocess.CalledProcessError as e:
        logger.info("!!! There was an error with OpenMS workflow, See the log !")
        logger.info(e.output)
        f = open(TOPPAS_folder+'/toppas_output/TOPPAS.log', 'r')
        file_contents = f.read()
        logger(file_contents)
        raise

    vdisplay.stop()

    # Error with the OpenMS workflow. No output files.
    try:
        mzTab_file = os.listdir(TOPPAS_folder+"/toppas_output/TOPPAS_out/MS2Planner_mzTab/")[0]
        f = open(TOPPAS_folder+'/toppas_output/TOPPAS_out/MS2Planner_mzTab/'+mzTab_file)
        f.close()
    except:
        logger.info('There was an issue with the OpenMS workflow ! See the log below. Not enough ressources might be available (memory, ...).')
        f = open(TOPPAS_folder+'/toppas_output/TOPPAS.log', 'r')
        file_contents = f.read()
        logger.info(file_contents)
        raise

    logger.info('======')
    logger.info('Completed the OpenMS workflow')
    logger.info('======')
    logger.info('Zipping up the OpenMS workflow results ..')

    try:
        cmdA = 'cp '+TOPPAS_folder+'/'+TOPPAS_output_folder+'/TOPPAS.log '+TOPPAS_folder+'/'+TOPPAS_output_folder+'/MS2Planner_mzTab/TOPPAS_PathFinder.log'
        os.system(cmdA)
        cmdB ='cp -r '+TOPPAS_folder+'/*MS2Planner*'+' '+TOPPAS_folder+'/'+TOPPAS_output_folder+'/MS2Planner_mzTab/'
        os.system(cmdB)
        get_all_file_paths(TOPPAS_folder+'/'+TOPPAS_output_folder+'/TOPPAS_out/MS2Planner_mzTab','download_results/MS2Planner_OpenMS_results.zip')
    except:
        raise

    logger.info('======')
    logger.info('Completed zipping up the OpenMS workflow result files')

    logger.info('======')
    logger.info('NOW CONTINUE WITH THE REST OF THE WORKFLOW')



if __name__ == "__main__":
    IODA_targeted_workflow(str(sys.argv[1]),str(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
