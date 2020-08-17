
# coding: utf-8
# Author: Louis Felix Nothias, louisfelix.nothias@gmail.com, June 2020
import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from io import StringIO
import warnings
from pandas.core.common import SettingWithCopyWarning
from format_to_qexactive_list import *
from zipfile import ZipFile
from logzero import logger, logfile
import datetime
from IODA_split_features import *
from format_to_qexactive_list import *
import subprocess
from subprocess import call
import pathlib

def convert_mzTab_to_table(input_filename: str,output_filename: str):
    """Take an mzTab containing two samples, output a table with mz, charge, rt, intensities."""
    df = pd.read_csv(input_filename, sep='\t', error_bad_lines=False, warn_bad_lines=False)

    # Get the metadata
    metadata = []
    start_row_consensus = 1
    for row in df['1.0.0']:
        metadata.append(row)
        start_row_consensus += 1

    # Change the type of the list
    metadata = [str(item) for item in metadata]
    [type(item) for item in metadata]

    # Get the filenames
    Filenames = []
    for x in metadata:
        if x.startswith('file:/'):
            x = x.split('/')[-1]
            #Remove duplicates
            if x not in Filenames:
                Filenames.append(x[:-5])

    #print('Filenames in the mzTab')
    #print(Filenames)

    Filename1 = Filenames[0]
    Filename2 = Filenames[1]

    # Display error message for additional samples
    for x in Filenames:
        if x == "ms_run[2]-location":
            Filename3 = Filenames[2]
            logger.info('Warning: There is more than two samples in that mzTab file. We support only two samples currently')
        if x == "ms_run[4]-location":
            Filename4 = Filenames[3]
            logger.info('Warning: There is more than three samples in that mzTab file. We support only two samples currently.')

    # Read and edit the table
    main_df = pd.read_csv(input_filename,sep='\t',index_col=0, skiprows=range(0,start_row_consensus))

    # Get the columns of interest
    feat_int = main_df.loc[:, main_df.columns.str.startswith('peptide_abundance_study_variable')]
    feat_mz = main_df.loc[:, main_df.columns.str.startswith('mass_to_charge')]
    feat_charge = main_df.loc[:, main_df.columns.str.startswith('charge')]
    feat_ret = main_df[['retention_time']]

    # Concat into a master table
    df_master = pd.concat([feat_mz,feat_ret,feat_charge,feat_int], axis=1)

    #Detection of blank
    #print('#Deducing the blank sample by comparing the sum of feature intensity between samples')
    column1_sum = df_master['peptide_abundance_study_variable[1]'].sum()
    logger.info('- For sample '+Filename1+' the sum of feature intensities is = '+str("{:.2e}".format(column1_sum)))
    column2_sum = df_master['peptide_abundance_study_variable[2]'].sum()
    logger.info('- For sample '+Filename2+' the sum of feature intensities = '+str("{:.2e}".format(column2_sum)))
    if column1_sum > column2_sum:
    #    logger.info('- The blank sample is assumed to be '+str(Filename2)+' in the mzTab-M')
    #    logger.info('- The samples is assumed to be '+str(Filename1)+' in the mzTab-M')
        df_master.rename(columns={'peptide_abundance_study_variable[1]':Filename2}, inplace=True)
        df_master.rename(columns={'peptide_abundance_study_variable[2]':Filename1}, inplace=True)
    if column1_sum < column2_sum:
    #    logger.info('- The blank sample is assumed to be '+str(Filename1)+' in the mzTab-M')
    #    logger.info('- The samples is assumed to be '+str(Filename2)+' in the mzTab-M')
        df_master.rename(columns={'peptide_abundance_study_variable[1]':Filename1}, inplace=True)
        df_master.rename(columns={'peptide_abundance_study_variable[2]':Filename2}, inplace=True)

    #Replace the sample headers for mandatory samples
    df_master.rename(columns={'mass_to_charge':"Mass [m/z]"}, inplace=True)

    # Replace the sample header for additional samples
    if 'peptide_abundance_study_variable[3]' in df_master.columns:
        df_master.rename(columns={'peptide_abundance_study_variable[3]':Filename3}, inplace=True)
    if 'peptide_abundance_study_variable[4]' in df_master.columns:
        df_master.rename(columns={'peptide_abundance_study_variable[4]':Filename4}, inplace=True)
    df_master = df_master.fillna(0).sort_values('retention_time')
    df_master.to_csv(output_filename, sep=',', index=False)
    return output_filename

def get_all_file_paths(directory,output_zip_path):
    # initializing empty file paths list
    file_paths = []

    # crawling through directory and subdirectories
    for root, directories, files in os.walk(directory):
        for filename in files:
            # join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)

            # writing files to a zipfile
    with ZipFile(output_zip_path,'w') as zip:
        # writing each file one by one
        for file in file_paths:
            zip.write(file)

    logger.info('All files zipped successfully!')

# Run the Path Finder workflow with baseline method
def run_path_finder_baseline_from_mzTab(input_filename:int, num_path:int, intensity_ratio:float, intensity_threshold:float, win_len:float, isolation:float, rt_margin:float):

    output_dir = 'results_targeted_pathfinder_baseline'
    os.system('rm -r '+output_dir)
    os.system('rm -r download_'+output_dir)
    os.system('mkdir '+output_dir)
    os.system('mkdir download_'+output_dir)
    logfile(output_dir+'/logfile.txt')

    logger.info('STARTING THE IODA targeted-from-mzTab WORKFLOW')
    if input_filename.startswith('http'):
        logger.info('File path was specified by the user')
        pass
    elif input_filename == 'OpenMS_generated':
        logger.info('The mzTab was generated with the IODA-OpenMS workflow')
        path_input_folder = "TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/"
        mzTab_file = os.listdir("TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/")[0]
        input_filename = path_input_folder+mzTab_file
    else:
        logger.info("the input_filename variable should be a valid path/download link or must be: 'OpenMS_generated', when using the OpenMS workflow online")

    now = datetime.datetime.now()
    logger.info(now)

    logger.info('======')
    logger.info('Getting the mzTab')
    if input_filename.startswith('http'):
        if 'google' in input_filename:
            logger.info('This is the Google Drive download link: '+str(input_filename))
            url_id = input_filename.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_filename = prefixe_google_download+url_id
            output_filename = output_dir+'/Converted_mzTab.csv'

        else:
            output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
            logger.info('This is the input file path: '+str(input_filename))
            logger.info('This is the output file path: '+str(output_filename))

    else:
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))


    # Convert the mzTab into a Table
    logger.info('======')
    logger.info('Converting mzTab to intermediate table format ...')
    convert_mzTab_to_table(input_filename,output_filename)
    logger.info('======')

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_filename)
    samplename = feature_table.columns[-1]
    logger.info('Assumed sample filename: '+samplename)
    blank_samplename = feature_table.columns[-2]
    logger.info('Assumed blank filename: ' +blank_samplename)
    logger.info('======')

    # User-defined parameters
    logger.info('User-defined parameters')
    ratio = intensity_ratio
    logger.info('Ratio between sample/blank for ion filtering = ' + str(ratio))
    min_intensity = intensity_threshold
    logger.info('Minimum intensity for ion filtering in sample = '+ str("{:.2e}".format(min_intensity)))
    logger.info('Retention time window (min.) for binning target ions = ' +str(win_len))
    logger.info('Isolation window (m/z) = ' +str(isolation))
    logger.info('Retention time margin (sec.) = ' +str(rt_margin))
    experiements = num_path
    logger.info('Number of iterative experiment(s) = ' + str(experiements))
    logger.info('======')

    # Running the table processing
    logger.info('Running Path Finder in Baseline mode ...')
    run_pathfinder_baseline(output_filename, output_filename[:-4]+'_PathFinder.csv', intensity_threshold, intensity_ratio, num_path, win_len, isolation)
    logger.info('======')

    Test_BestPath_Output = pathlib.Path(output_filename[:-4]+'_PathFinder.csv')
    try:
        if Test_BestPath_Output.exists ():
            logger.info("Best Path output found")
        else:
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            print("Problem when running Best Path !!!")
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            logger.info("Problem when running Best Path !!!")
    except:
        raise

    logger.info('Preparing results ...')
    transient_time = 0 #Hardcoded to keep the same def function with Curve mode. Parameter only used in Curve mode.
    logger.info('Preparing results ...')
    make_bestpath_targeted_lists_from_table(output_filename[:-4]+'_PathFinder.csv', rt_margin, transient_time)
    logger.info('======')

    logger.info('Cleaning and zipping workflow results files ...')

    # Cleaning files first

    #mkdir XCalibur
    os.system('mkdir '+output_dir+'/XCalibur')
    # mv files XCalibur
    os.system('mv '+output_dir+'/*formatted_QE* '+output_dir+'/XCalibur')

    #mkdir MQL
    os.system('mkdir '+output_dir+'/MaxQuantLive')

    # mv files MQL
    os.system('mv '+output_dir+'/*formatted_MQL* '+output_dir+'/MaxQuantLive')

    # mkdir intermediate files
    os.system('mkdir '+output_dir+'/intermediate_files')
    os.system('mkdir '+output_dir+'/plots')
    os.system('mkdir '+output_dir+'/log')

    # mv
    os.system('mv '+output_dir+'/*scatter_plot* '+output_dir+'/plots')
    os.system('mv '+output_dir+'/logfile.txt '+output_dir+'/log')
    os.system('mv '+output_dir+'/*.csv '+output_dir+'/intermediate_files')
    os.system('mv '+output_dir+'/*.txt '+output_dir+'/intermediate_files')

    get_all_file_paths(output_dir,'download_'+output_dir+'/IODA_Path_Finder_baseline_results.zip')

    logger.info('======')
    logger.info('END OF THE IODA-Path-Finder-Baseline-from-mzTab WORKFLOW')
    logger.info('======')
    print(' ')



# Run the Path Finder workflow with apex method
def run_path_finder_apex_from_mzTab(input_filename:int, num_path:int, intensity_ratio:float, intensity_threshold:float, intensity_accu:float, isolation:float, delta:float, rt_margin:float):

    output_dir = 'results_targeted_pathfinder_apex'
    os.system('rm -r '+output_dir)
    os.system('rm -r download_'+output_dir)
    os.system('mkdir '+output_dir)
    os.system('mkdir download_'+output_dir)
    logfile(output_dir+'/logfile.txt')

    logger.info('STARTING THE IODA targeted-from-mzTab WORKFLOW')
    if input_filename.startswith('http'):
        logger.info('File path was specified by the user')
        pass
    elif input_filename == 'OpenMS_generated':
        logger.info('The mzTab was generated with the IODA-OpenMS workflow')
        path_input_folder = "TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/"
        mzTab_file = os.listdir("TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/")[0]
        input_filename = path_input_folder+mzTab_file
    else:
        logger.info("the input_filename variable should be a valid path/download link or must be: 'OpenMS_generated', when using the OpenMS workflow online")

    now = datetime.datetime.now()
    logger.info(now)

    logger.info('======')
    logger.info('Getting the mzTab')
    if input_filename.startswith('http'):
        if 'google' in input_filename:
            logger.info('This is the Google Drive download link: '+str(input_filename))
            url_id = input_filename.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_filename = prefixe_google_download+url_id
            output_filename = output_dir+'/Converted_mzTab.csv'

        else:
            output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
            logger.info('This is the input file path: '+str(input_filename))
            logger.info('This is the output file path: '+str(output_filename))

    else:
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))


    # Convert the mzTab into a Table
    logger.info('======')
    logger.info('Converting mzTab to intermediate table format ...')
    convert_mzTab_to_table(input_filename,output_filename)
    logger.info('======')

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_filename)
    samplename = feature_table.columns[-1]
    logger.info('Assumed sample filename: '+samplename)
    blank_samplename = feature_table.columns[-2]
    logger.info('Assumed blank filename: ' +blank_samplename)
    logger.info('======')

    # User-defined parameters
    logger.info('User-defined parameters')
    ratio = intensity_ratio
    logger.info('Ratio between sample/blank for ion filtering = ' + str(ratio))
    min_intensity = intensity_threshold
    logger.info('Minimum intensity for ion filtering in sample = '+ str("{:.2e}".format(min_intensity)))
    logger.info('Precursor ion intensity to accumulate in the MS2 scan = ' +str("{:.2e}".format(intensity_accu)))
    logger.info('Isolation window (m/z) = ' +str(isolation))
    logger.info('Retention time margin (sec.) = ' +str(rt_margin))


    logger.info('Delay between targeted MS2 scans (sec)= ' +str(delta))
    experiements = num_path
    logger.info('Number of iterative experiment(s) = ' + str(experiements))
    logger.info('======')

    # Running the table processing
    logger.info('Running Path Finder in Apex mode ...')
    run_pathfinder_apex(output_filename, output_filename[:-4]+'_PathFinder.csv', intensity_threshold, intensity_ratio, num_path, intensity_accu, isolation, delta)
    logger.info('======')

    Test_BestPath_Output = pathlib.Path(output_filename[:-4]+'_PathFinder.csv')
    try:
        if Test_BestPath_Output.exists ():
            logger.info("Best Path output found")
        else:
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            print("Problem when running Best Path !!!")
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            logger.info("Problem when running Best Path !!!")
    except:
        raise

    transient_time = 0 #Hardcoded to keep the same def function with Curve mode. Parameter only used in Curve mode.
    logger.info('Preparing results ...')
    make_bestpath_targeted_lists_from_table(output_filename[:-4]+'_PathFinder.csv',rt_margin, transient_time)
    logger.info('======')

    logger.info('Cleaning and zipping workflow results files ...')

    # Cleaning files first

    #mkdir XCalibur
    os.system('mkdir '+output_dir+'/XCalibur')
    # mv files XCalibur
    os.system('mv '+output_dir+'/*formatted_QE* '+output_dir+'/XCalibur')

    #mkdir MQL
    os.system('mkdir '+output_dir+'/MaxQuantLive')

    # mv files MQL
    os.system('mv '+output_dir+'/*formatted_MQL* '+output_dir+'/MaxQuantLive')

    # mkdir intermediate files
    os.system('mkdir '+output_dir+'/intermediate_files')
    os.system('mkdir '+output_dir+'/plots')
    os.system('mkdir '+output_dir+'/log')

    # mv
    os.system('mv '+output_dir+'/*scatter_plot* '+output_dir+'/plots')
    os.system('mv '+output_dir+'/logfile.txt '+output_dir+'/log')
    os.system('mv '+output_dir+'/*.csv '+output_dir+'/intermediate_files')
    os.system('mv '+output_dir+'/*.txt '+output_dir+'/intermediate_files')

    get_all_file_paths(output_dir,'download_'+output_dir+'/IODA_Path_Finder_apex_results.zip')

    logger.info('======')
    logger.info('END OF THE IODA-Path-Finder-Apex-from-mzTab WORKFLOW')
    logger.info('======')
    print(' ')



# Run the Path Finder workflow with apex method
def run_path_finder_curve_from_mzTab(input_filename:int, num_path:int, intensity_ratio:float, intensity_threshold:float, input_filename_curve:int, intensity_accu:float, restriction:float, mz_accuracy:float, delta:float, rt_margin:float, transient_time:float):

    output_dir = 'results_targeted_pathfinder_curve'
    os.system('rm -r '+output_dir)
    os.system('rm -r download_'+output_dir)
    os.system('mkdir '+output_dir)
    os.system('mkdir download_'+output_dir)
    logfile(output_dir+'/logfile.txt')

    logger.info('STARTING THE IODA targeted-from-mzTab WORKFLOW')
    if input_filename.startswith('http'):
        logger.info('File path was specified by the user')
        pass
    elif input_filename == 'OpenMS_generated':
        logger.info('The mzTab was generated with the IODA-OpenMS workflow')
        path_input_folder = "TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/"
        mzTab_file = os.listdir("TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/")[0]
        input_filename = path_input_folder+mzTab_file
    else:
        logger.info("the input_filename variable should be a valid path/download link or must be: 'OpenMS_generated', when using the OpenMS workflow online")

    now = datetime.datetime.now()
    logger.info(now)

    logger.info('======')
    logger.info('Getting the mzTab')
    if input_filename.startswith('http'):
        if 'google' in input_filename:
            logger.info('This is the Google Drive download link: '+str(input_filename))
            url_id = input_filename.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_filename = prefixe_google_download+url_id
            output_filename = output_dir+'/Converted_mzTab.csv'

        else:
            output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
            logger.info('This is the input file path: '+str(input_filename))
            logger.info('This is the output file path: '+str(output_filename))

    else:
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))

    # Convert the mzTab into a Table
    logger.info('======')
    logger.info('Converting mzTab to intermediate table format ...')
    convert_mzTab_to_table(input_filename,output_filename)
    logger.info('======')

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_filename)
    samplename = feature_table.columns[-1]
    logger.info('Assumed sample filename: '+samplename)
    blank_samplename = feature_table.columns[-2]
    logger.info('Assumed blank filename: ' +blank_samplename)
    logger.info('======')

    # User-defined parameters
    logger.info('User-defined parameters')
    ratio = intensity_ratio
    logger.info('Ratio between sample/blank for ion filtering = ' + str(ratio))
    min_intensity = intensity_threshold
    logger.info('Minimum intensity for ion filtering in sample = '+ str("{:.2e}".format(min_intensity)))
    logger.info('Precursor ion intensity to accumulate in the MS2 scan = ' +str("{:.2e}".format(intensity_accu)))
    logger.info('Transient time and overhead (ms) = '+str(transient_time))
    logger.info('Input file for curve data : ' +str(input_filename_curve))
    logger.info('Retention time margin (sec.) = ' +str(rt_margin))
    logger.info('Restriction parameter : ' +str(restriction))
    logger.info('Mass accuracy (m/z): ' +str(mz_accuracy))
    logger.info('Delay between targeted MS2 scans (sec)= ' +str(delta))
    experiements = num_path
    logger.info('Number of iterative experiment(s) = ' + str(experiements))
    logger.info('======')

    #Checking user provided file for Path Finder
    mzTab_curve = pathlib.Path(input_filename_curve)
    try:
        if mzTab_curve.exists ():
            logger.info("mzTab for the Curve mode found")
        else:
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            print("Problem with the mzTab file or file path ! Please verify")
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            logger.info("Problem with the mzTab file or file path ! Please verify")
    except:
        raise

    #Running Path Finder
    logger.info('Running Path Finder in Curve mode ...')
    try:
        run_pathfinder_curve(output_filename, output_filename[:-4]+'_PathFinder.csv', intensity_threshold, intensity_ratio, num_path, input_filename_curve, intensity_accu, restriction, mz_accuracy, delta)
    except:
        raise

    Test_BestPath_Output = pathlib.Path(output_filename[:-4]+'_PathFinder.csv')
    try:
        if Test_BestPath_Output.exists ():
            logger.info("Best Path output found")
        else:
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            print("Problem when running Best Path !!!")
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            logger.info("Problem when running Best Path !!!")
    except:
        raise

    logger.info('======')
    logger.info('Preparing results ...')
    make_bestpath_targeted_lists_from_table(output_filename[:-4]+'_PathFinder.csv', rt_margin, transient_time)
    logger.info('======')

    logger.info('Cleaning and zipping workflow results files ...')

    # Cleaning files first

    #mkdir XCalibur
    os.system('mkdir '+output_dir+'/XCalibur')
    # mv files XCalibur
    os.system('mv '+output_dir+'/*formatted_QE* '+output_dir+'/XCalibur')

    #mkdir MQL
    os.system('mkdir '+output_dir+'/MaxQuantLive')

    # mv files MQL
    os.system('mv '+output_dir+'/*formatted_MQL* '+output_dir+'/MaxQuantLive')

    # mkdir intermediate files
    os.system('mkdir '+output_dir+'/intermediate_files')
    os.system('mkdir '+output_dir+'/plots')
    os.system('mkdir '+output_dir+'/log')

    # mv
    os.system('mv '+output_dir+'/*scatter_plot* '+output_dir+'/plots')
    os.system('mv '+output_dir+'/logfile.txt '+output_dir+'/log')
    os.system('mv '+output_dir+'/*.csv '+output_dir+'/intermediate_files')
    os.system('mv '+output_dir+'/*.txt '+output_dir+'/intermediate_files')

    get_all_file_paths(output_dir,'download_'+output_dir+'/IODA_Path_Finder_curve_results.zip')

    logger.info('======')
    logger.info('END OF THE IODA-Path-Finder-Apex-from-mzTab WORKFLOW')
    logger.info('======')
    print(' ')


### PathFinder
# This parse one line from the BestPath output and create a output table per path. The rows to skip define which line/path is parsed.
def bestpath_format(input_filename: str, output_filename: str, rows_to_skip:int):
    df_path = pd.read_csv(input_filename, sep=' ', header=None, skiprows=rows_to_skip, error_bad_lines=False, warn_bad_lines=False)

    #Make a list for the first row
    df_path_list = df_path.iloc[0].values.tolist()
    df_path_list.pop(0)
    nfeatures = int(len(df_path_list)/8)

    # Convert the list into a nested list
    target_list = []
    for entries in range(nfeatures):
            try:
                while len(df_path_list) > 8:
                    target_list.append(df_path_list[:8])
                    indexes = [0,1,2,3,4,5,6,7]
                    for index in sorted(indexes, reverse=True):
                        del df_path_list[index]
            except:
                continue
    #Make a dataframe
    target_table = pd.DataFrame(target_list)
    target_table = target_table.rename(columns={0: 'Mass [m/z]',1: 'mz_isolation',2: 'duration',3: 'rt_start',4: 'rt_end', 5: 'intensity', 6: 'rt_apex',7: 'charge'})
    #############################target_table = target_table[target_table['intensity'] > 0]

    logger.info('Valid target ions in path'+str(rows_to_skip+1)+' = '+str(target_table.shape[0]))

    target_table.to_csv(output_filename, sep=',', index=False)

# This parse BestPath output file and create output tables formatted for XCalibur and MaxQuant.live
def make_bestpath_targeted_lists_from_table(input_filename:str,rt_margin:float, transient_time:float):
    os.system("sed -i 's/\t/ /g' "+input_filename)
    logger.info('File processed: '+input_filename)
    logger.info('======')
    with open(input_filename) as file:
        counter = -1
        for line in file:
            try:
                counter += 1
                output_filename = input_filename[:-4]+"_"+str(counter+1)+'_formatted.txt'
                #Take the list and make a table
                bestpath_format(input_filename,output_filename, counter)
                logger.info('Formatting to XCalibur format ...')
                generate_QE_list_from_BestPath(output_filename, output_filename[:-4]+'_QE.csv',rt_margin)
                #Format for MaxQuant.Live targeted experiment
                logger.info('Formatting for MaxQuant.Live ...')
                generate_MQL_list_from_BestPath(output_filename, output_filename[:-4]+'_MQL.txt',rt_margin)
                generate_MQL_list_from_BestPath_MaxIT(output_filename, output_filename[:-4]+'_MQL_variableMaxIT.txt',transient_time)
                logger.info('=======')
            except:
                raise

    table_list_bestpath = []
    for x in range(0,counter+1):
        table_list_bestpath.append(input_filename[:-4]+"_"+str(x+1)+'_formatted.txt')

    logger.info('Plotting results ...')
    try:
        make_plot_bestpath1(table_list_bestpath,output_filename)
        make_plot_bestpath2(table_list_bestpath,output_filename)
    except:
        raise
    try:
        make_plot_bestpath3(table_list_bestpath,output_filename)
    except:
        pass


def run_pathfinder_baseline(input_filename:str, output_filename:str, intensity_threshold:float, intensity_ratio:float, num_path:int, win_len:float, isolation:float):
    cmd_baseline = ('python3 path_finder.py baseline '+input_filename+' '+output_filename+' '+str(intensity_threshold)+' '+str(intensity_ratio)+' '+str(num_path)+' -win_len '+str(win_len)+' -isolation '+str(isolation))
    logger.info('Command: '+cmd_baseline)
    os.system(cmd_baseline)

def run_pathfinder_apex(input_filename:str, output_filename:str, intensity_threshold:float, intensity_ratio:float, num_path:int, intensity_accu:float, isolation:float, delta:float):
    cmd_apex = ('python3 path_finder.py apex '+input_filename+' '+output_filename+' '+str(intensity_threshold)+' '+str(intensity_ratio)+' '+str(num_path)+' -intensity_accu '+str(intensity_accu)+' -isolation '+str(isolation)+' -delta '+str(delta))
    logger.info('Command: '+cmd_apex)
    os.system(cmd_apex)

def run_pathfinder_curve(input_filename:str, output_filename:str, intensity_threshold:float, intensity_ratio:float, num_path:int, input_filename_curve:str, intensity_accu:float, restriction:float, mz_accuracy:float, delta:float):
    cmd_curve = ('python3 path_finder.py curve '+input_filename+' '+output_filename+' '+str(intensity_threshold)+' '+str(intensity_ratio)+' '+str(num_path)+' -infile_raw '+str(input_filename_curve)+' -intensity_accu '+str(intensity_accu)+' -restriction '+str(restriction)+' '+str(mz_accuracy)+' -delta '+str(delta))
    logger.info('Command: '+cmd_curve)
    logger.info('Path Finder in Curve mode can take up to 10 minutes to complete ... please wait')
    try:
        cp0 = subprocess.run(cmd_curve,shell=True)
        cp0
    except subprocess.CalledProcessError:
        logger.info('ERROR running Path Finder ...')

#Best path generate mz / rt figures
def make_plot_bestpath1(table_list_bestpath, output_filename):
    Labels = []
    if len(table_list_bestpath) >= 1:
        table0 = pd.read_csv(table_list_bestpath[0], sep=',', header=0)
        plt.scatter('rt_apex','Mass [m/z]', data=table0, marker='o', color='blue',s=2, alpha=0.6)
        Label1 = ['Inj. 1, n = '+ str(table0.shape[0])+ ', median = '+ "{0:.2e}".format(table0['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table0['rt_apex'].mean())]
        Labels.append(Label1)

    if len(table_list_bestpath) >= 2:
        table1 = pd.read_csv(table_list_bestpath[1], sep=',', header=0)
        plt.scatter('rt_apex','Mass [m/z]',  data=table1, marker='o', color='violet',s=2, alpha=0.6)
        Label2 = ['Inj. 2, n = '+ str(table1.shape[0])+ ', median = '+ "{0:.2e}".format(table1['rt_apex'].median())  + ', mean = '+ "{0:.2e}".format(table1['rt_apex'].mean())]
        Labels.append(Label2)

    if len(table_list_bestpath) >= 3:
        table2 = pd.read_csv(table_list_bestpath[2], sep=',', header=0)
        plt.scatter('rt_apex','Mass [m/z]',  data=table2, marker='o', color='orange',s=2, alpha=0.6)
        Label3 = ['Inj. 3, n = '+ str(table2.shape[0])+ ', median = '+ "{0:.2e}".format(table2['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table2['rt_apex'].mean())]
        Labels.append(Label3)

    if len(table_list_bestpath) >= 4:
        table3 = pd.read_csv(table_list_bestpath[3], sep=',', header=0)
        plt.scatter('rt_apex','Mass [m/z]',  data=table3, marker='o', color='red', s=2, alpha=0.6)
        Label4 =['Inj. 4, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table3['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table3['rt_apex'].mean())]
        Labels.append(Label4)

    if len(table_list_bestpath) >= 5:
        table4 = pd.read_csv(table_list_bestpath[4], sep=',', header=0)
        plt.scatter('rt_apex','Mass [m/z]',  data=table4, marker='o', color='red', s=2, alpha=0.6)
        Label5 =['Inj. 5, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table4['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table4['rt_apex'].mean())]
        Labels.append(Label5)

    #plt.title('Features per path: '+ str(table_list[0]))
    plt.ylabel('Mass [m/z]')
    plt.xlabel('Ret. time apex (s)')

    plt.legend(labels=Labels, fontsize =4)
    plt.savefig(output_filename[:-4]+'injection_scatter_plot_mz_rt.png', dpi=300)
    plt.close()

#Best path generate feature intensity / rt figures
def make_plot_bestpath2(table_list_bestpath, output_filename):
    Labels = []
    if len(table_list_bestpath) >= 1:
        table0 = pd.read_csv(table_list_bestpath[0], sep=',', header=0)
        plt.scatter('rt_apex','intensity', data=table0, marker='o', color='blue',s=2, alpha=0.6)
        Label1 = ['Inj. 1, n = '+ str(table0.shape[0])+ ', median = '+ "{0:.2e}".format(table0['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table0['rt_apex'].mean())]
        Labels.append(Label1)

    if len(table_list_bestpath) >= 2:
        table1 = pd.read_csv(table_list_bestpath[1], sep=',', header=0)
        plt.scatter('rt_apex','intensity',  data=table1, marker='o', color='violet',s=2, alpha=0.6)
        Label2 = ['Inj. 2, n = '+ str(table1.shape[0])+ ', median = '+ "{0:.2e}".format(table1['rt_apex'].median())  + ', mean = '+ "{0:.2e}".format(table1['rt_apex'].mean())]
        Labels.append(Label2)

    if len(table_list_bestpath) >= 3:
        table2 = pd.read_csv(table_list_bestpath[2], sep=',', header=0)
        plt.scatter('rt_apex','intensity',  data=table2, marker='o', color='orange',s=2, alpha=0.6)
        Label3 = ['Inj. 3, n = '+ str(table2.shape[0])+ ', median = '+ "{0:.2e}".format(table2['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table2['rt_apex'].mean())]
        Labels.append(Label3)

    if len(table_list_bestpath) >= 4:
        table3 = pd.read_csv(table_list_bestpath[3], sep=',', header=0)
        plt.scatter('rt_apex','intensity',  data=table3, marker='o', color='red', s=2, alpha=0.6)
        Label4 =['Inj. 4, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table3['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table3['rt_apex'].mean())]
        Labels.append(Label4)

    if len(table_list_bestpath) >= 5:
        table4 = pd.read_csv(table_list_bestpath[4], sep=',', header=0)
        plt.scatter('rt_apex','intensity',  data=table4, marker='o', color='red', s=2, alpha=0.6)
        Label5 =['Inj. 5, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table4['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table4['rt_apex'].mean())]
        Labels.append(Label5)

    plt.yscale('log')
    #plt.ylim(bottom=4)
    #plt.xlim(bottom=1)

    #plt.title('Features per path: '+ str(table_list[0]))
    plt.ylabel('Feature intensity (log scale)')
    plt.xlabel('Ret. time apex (s)')

    plt.legend(labels=Labels, fontsize =4)
    plt.savefig(output_filename[:-4]+'injection_scatter_plot_intensity_rt.png', dpi=300)
    plt.close()



#Best path generate feature intensity / duration
def make_plot_bestpath3(table_list_bestpath, output_filename):
    Labels = []
    if len(table_list_bestpath) >= 0:
        table0 = pd.read_csv(table_list_bestpath[0], sep=',', header=0)
        plt.scatter('duration','intensity', data=table0, marker='o', color='blue',s=1, alpha=0.6)
        Label = ['Inj. 1, n = '+ str(table0.shape[0])+ ', median = '+ "{0:.2e}".format(table0['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table0['rt_apex'].mean())]
        Labels.append(Label)
        plt.yscale('log')
        plt.ylim(bottom=1E5)
        plt.ylabel('Feature intensity (log.  scale)')
        plt.xlabel('Duration (s)')
        plt.legend(labels=Labels, fontsize =5)

        plt.savefig(output_filename[:-4]+'injection1_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()

    if len(table_list_bestpath) >= 1:
        table1 = pd.read_csv(table_list_bestpath[1], sep=',', header=0)
        plt.scatter('duration','intensity',  data=table1, marker='o', color='violet',s=1.5, alpha=0.6)
        Label = ['Inj. 2, n = '+ str(table1.shape[0])+ ', median = '+ "{0:.2e}".format(table1['rt_apex'].median())  + ', mean = '+ "{0:.2e}".format(table1['rt_apex'].mean())]
        Labels.clear()
        Labels.append(Label)

        plt.savefig(output_filename[:-4]+'injection2_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()

    if len(table_list_bestpath) >= 2:
        table2 = pd.read_csv(table_list_bestpath[2], sep=',', header=0)
        plt.scatter('duration','intensity',  data=table2, marker='o', color='orange',s=1, alpha=0.6)
        Label = ['Inj. 3, n = '+ str(table2.shape[0])+ ', median = '+ "{0:.2e}".format(table2['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table2['rt_apex'].mean())]
        Labels.clear()
        Labels.append(Label)

        plt.savefig(output_filename[:-4]+'injection3_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()

    if len(table_list_bestpath) >= 3:
        table3 = pd.read_csv(table_list_bestpath[3], sep=',', header=0)
        plt.scatter('duration','intensity',  data=table3, marker='o', color='red', s=0.5, alpha=0.6)
        Label =['Inj. 4, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table3['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table3['rt_apex'].mean())]
        Labels.clear()
        Labels.append(Label)
        plt.savefig(output_filename[:-4]+'injection4_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()

    if len(table_list_bestpath) >= 4:
        table4 = pd.read_csv(table_list_bestpath[4], sep=',', header=0)
        plt.scatter('duration','intensity',  data=table4, marker='o', color='red', s=0.1, alpha=0.6)
        Label =['Inj. 5, n = '+ str(table4.shape[0])+ ', median = '+ "{0:.2e}".format(table4['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table4['rt_apex'].mean())]
        Labels.clear()
        Labels.append(Label)
        plt.savefig(output_filename[:-4]+'injection5_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()


if __name__ == "__main__":
    make_targeted_list_from_mzTab(str(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
