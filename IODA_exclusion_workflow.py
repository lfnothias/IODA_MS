
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

def convert_blank_range_mzTab_to_table(input_filename: str, output_filename: str):
    """Take an mzTab containing one sample, output a table with mz, charge, rt, intensities."""
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

    logger.info('   Filename(s) in the mzTab: '+str(Filenames))

    Filename1 = Filenames[0]

    # Display error message for additional samples
    for x in Filenames:
        if x == "ms_run[2]-location":
            Filename2 = Filenames[1]
            logger.info('Warning: There is more than two samples in that mzTab file. We support only one sample in the mzTab currently')

    # Read and edit the table
    main_df = pd.read_csv(input_filename,sep='\t',index_col=0, skiprows=range(0,start_row_consensus))

    # Get the columns of interest
    feat_int = main_df.loc[:, main_df.columns.str.startswith('peptide_abundance_study_variable')]
    feat_mz = main_df.loc[:, main_df.columns.str.startswith('mass_to_charge')]
    feat_charge = main_df.loc[:, main_df.columns.str.startswith('charge')]
    feat_ret = main_df[['retention_time']]
    feat_ret_list = main_df['retention_time_window'].to_list()

    split_list =[i.split('|') for i in feat_ret_list]
    feat_ret_start_end = pd.DataFrame(split_list,columns=None)
    feat_ret_start_end.rename(columns={0:'rt_start'}, inplace=True)
    feat_ret_start_end.rename(columns={1:'rt_end'}, inplace=True)

    # Concat into a master table
    df_master = pd.concat([feat_mz,feat_ret,feat_charge,feat_int], axis=1)
    df_master['rt_start'] = feat_ret_start_end['rt_start'].to_list()
    df_master['rt_end'] = feat_ret_start_end['rt_end'].to_list()

    df_master.rename(columns={'peptide_abundance_study_variable[1]':Filename1}, inplace=True)

    #Replace the sample headers for mandatory samples
    df_master.rename(columns={'mass_to_charge':"Mass [m/z]"}, inplace=True)

    # Fill the charge with 0
    df_master = df_master.sort_values('retention_time')
    df_master.to_csv(output_filename, sep=',', index=False)
    return output_filename

def make_exclusion_list(input_filename: str, sample: str, intensity:float):
    """From a table with mz, charge, rt, intensities, make an exclusion list from a single sample, above the intensity specified."""
    df_master = pd.read_csv(input_filename, sep=',')
    output_filename = input_filename
    df_master_exclusion_list = df_master[(df_master[sample] != 0)]
    df_master_exclusion_list = df_master[(df_master[sample] > intensity)]
    df_master_exclusion_list.to_csv(output_filename[:-4]+'_EXCLUSION_BLANK.csv', sep=',', index = False)
    #df_master_exclusion_list.sort_values(by=['Mass [m/z]'])
    number_of_ions = 'Initial number of ions = '+ str(df_master.shape[0])
    logger.info(number_of_ions)
    number_of_ions_after_filtering = 'Number of ions after intensity filtering = '+str(df_master_exclusion_list.shape[0]) +', with intensity >'+ str(intensity)
    logger.info(number_of_ions_after_filtering)

    
def plot_targets_exclusion(input_filename: str, blank_samplename: str, column: str, title: str):
    """From a table, make a scatter plot of a sample"""
    Labels = []
    table0 = pd.read_csv(input_filename, sep=',', header=0)
    fig = plt.figure(figsize=(8,6))
    fig = plt.scatter(column, blank_samplename, data=table0, marker='o', color='lightskyblue',s=3, alpha=0.4)
    Label1 = ['Excluded ions, n = '+ str(table0.shape[0])+ ', median abs. int. = '+ "{0:.2e}".format(table0[blank_samplename].median()) + ', mean abs. int. = '+ "{0:.2e}".format(table0[blank_samplename].mean())]
    Labels.append(Label1)
    plt.yscale('log')
    plt.ylabel('Ion intensity (log scale)', size = 11)
    plt.legend(labels=Labels, fontsize = 8, loc='best', markerscale=5)
    if column == 'Mass [m/z]':
        plt.title(title+', in m/z range', size = 11,  wrap=True)
        plt.xlabel('m/z', size = 10)
        plt.savefig('results/plot_exclusion_scatter_MZ.png', dpi=200)
    if column == 'retention_time':
        plt.title(title+', in retention time range', size =11, wrap=True)
        plt.xlabel('Ret. time (sec)', size = 10)
        plt.savefig('results/plot_exclusion_scatter_RT.png', dpi=200)
    plt.close()


def plot_targets_exclusion_range(input_filename: str, blank_samplename: str, title: str):
    Labels = []
    table0 = pd.read_csv(input_filename, sep=',', header=0)
    rt_start = table0['retention_time']-table0['rt_start']
    rt_end = table0['rt_end']-table0['retention_time']
    rt_range = [rt_start, rt_end]
    # Normalizing
    table0[blank_samplename]=((table0[blank_samplename]-table0[blank_samplename].min())/(table0[blank_samplename].max()-table0[blank_samplename].min()))*300
    gradient = table0[blank_samplename].to_list()
    plt.figure(figsize=(9,6))
    plt.errorbar('retention_time','Mass [m/z]', data=table0, xerr=rt_range, fmt='.', elinewidth=0.7, color='lightskyblue', ecolor='grey', capsize=0, alpha=0.3)
    plt.scatter('retention_time','Mass [m/z]', data=table0, s = gradient, marker = "o", facecolors='', color='', edgecolors='red', alpha=0.85, linewidth=0.45)

    Label1 = ['Red circle = intensity, Blue dot = ion apex, Horizontal line = RT range, Ions excluded (n='+ str(table0.shape[0])+')']
    Labels.append(Label1)
    plt.title(title, size =11, wrap=True)
    plt.xlabel('Ret. time (sec)')
    plt.ylabel('m/z')
    plt.legend(labels=Labels, fontsize = 8, loc='upper left', markerscale=0.45)
    plt.savefig('results/plot_exclusion_RT_range_plot.png', dpi=200)
    plt.close()

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

# Make exclusion list from two mzTabs (design to run with the OpenMS TOPPAS workflow in IODA)
def make_exclusion_from_mzTabs(input_mzML:int, min_intensity:int, rtexclusionmargininsecs:float):
    input_dir='TOPPAS_Workflow/toppas_output'
    output_dir = 'results'

    os.system('rm -r results')
    os.system('rm download_results/IODA_exclusion_results.zip')
    os.system('mkdir results')
    os.system('mkdir download_results')
    now = datetime.datetime.now()
    logger.info(now)
    os.system('rm results/logfile.txt')
    logfile('results/logfile.txt')

    # Convert the mzTabs into Tables to generate exclusion list
    logger.info('======')
    logger.info('Starting the IODA-exclusion workflow')
    logger.info('The source of the mzML file is: '+input_mzML)
    logger.info('Intermediate file input: '+input_dir+'/TOPPAS_out/mzTab_Narrow/Blank.mzTab')
    logger.info('Intermediate file input: '+input_dir+'/TOPPAS_out/mzTab_Large/Blank.mzTab')
    logger.info('======')
    logger.info('Converting mzTab to table format')
    logger.info('For narrow features')
    convert_blank_range_mzTab_to_table(input_dir+'/TOPPAS_out/mzTab_Narrow/Blank.mzTab', output_dir+'/table_narrow.csv')
    logger.info('For large features')
    convert_blank_range_mzTab_to_table(input_dir+'/TOPPAS_out/mzTab_Large/Blank.mzTab', output_dir+'/table_large.csv')
    logger.info('======')

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_dir+'/table_large.csv')
    blank_samplename = feature_table.columns[-3]

    output_filename = output_dir+'/'+blank_samplename+'.csv'
    logger.info('Blank sample name: was internally renamed as '+ blank_samplename)

    # User-defined parameters
    logger.info('User-defined parameters')
    logger.info('   Minimum ion intensity treshold (count) = '+ str(min_intensity))
    logger.info('   Additional margin for retention time range exclusion (seconds) = '+ str(rtexclusionmargininsecs))

    # Concatenating the tables from narrow and large features:
    df_narrow = pd.read_csv(output_dir+'/table_narrow.csv',sep=',')
    df_large = pd.read_csv(output_dir+'/table_large.csv',sep=',')
    df_concat = pd.concat([df_narrow,df_large])

    #We expand the margin of the retention time exclusion
    df_concat['rt_start'] = df_concat['rt_start'] - rtexclusionmargininsecs
    df_concat['rt_end'] = df_concat['rt_end'] + rtexclusionmargininsecs

    #Concatening the tables
    df_concat.to_csv(output_filename, sep=',', index=False)

    # Running the table processing
    logger.info('Running the table processing')
    make_exclusion_list(output_filename, blank_samplename, min_intensity)
    logger.info('======')

    # Convert to XCalibur format
    logger.info('Preparing list of excluded ions in XCalibur format')
    generate_QE_list_rt_range(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename,output_filename[:-4]+'_EXCLUSION_LIST_XCalibur.csv')
    logger.info('======')
    logger.info('Preparing list of excluded ions in MaxQuant.Live format')
    generate_MQL_exclusion(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, output_filename[:-4]+'_EXCLUSION_LIST_MaxQuantLive.txt')
    logger.info('======')

    # === Plot the features  ====
    logger.info('Preparing scatter plots of the excluded ions/features')
    plot_targets_exclusion_range(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, 'Distribution of the excluded ions in the blank sample')
    plot_targets_exclusion(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, 'retention_time', 'Distribution of intensities for the excluded ions in '+ blank_samplename)
    plot_targets_exclusion(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, 'Mass [m/z]', 'Distribution of intensities for the excluded ions in '+ blank_samplename)

    logger.info('======')
    logger.info('Zipping workflow results files')

    # Cleaning the files first
    os.system('mkdir results/intermediate_files')
    os.system('mv '+output_filename[:-4]+'_EXCLUSION_BLANK.csv intermediate_files/')
    os.system('mv '+output_filename+' results/intermediate_files/')
    os.system('mv results/table* results/intermediate_files/')
    get_all_file_paths('results','download_results/IODA_exclusion_results.zip')

    logger.info('======')
    logger.info('End the IODA-exclusion workflow processing')
    logger.info('======')
    logger.info('Proceed below with the results visualization')


# Make exclusion list from one mzTab
def make_exclusion_from_mzTab(input_filename:str, min_intensity:int, rtexclusionmargininsecs:float):
    #Example source filenames
    #input_filename = 'https://drive.google.com/file/d/1LYk-PKsBWl4Pv7c1TlhQwaqwkF2T6sux/view?usp=sharing'
    #input_filename = 'tests/Euphorbia/exclusion/ioda_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzTab'
    os.system('rm -r results')
    os.system('rm download_results/IODA_exclusion_results.zip')
    os.system('mkdir results')
    os.system('mkdir download_results')
    os.system('rm results/logfile.txt')
    logfile('results/logfile.txt')
    now = datetime.datetime.now()
    logger.info(now)

    logger.info('Starting the IODA exclusion-from-mzTab workflow')
    output_dir = 'results'
    logger.info('======')
    logger.info('Getting the mzTab')
    logger.info('This is the input: '+input_filename)
    if input_filename.startswith('http'):
        if 'google' in input_filename:
            logger.info('This is the Google Drive download link:'+str(input_filename))
            url_id = input_filename.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_filename = prefixe_google_download+url_id
            output_filename = output_dir+'/Exclusion_sample.csv'
        else:
            output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
            logger.info('This is the output file path: '+str(output_filename))
    else:
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))
    logger.info('======')
    logger.info('Converting mzTab to table format')
    convert_blank_range_mzTab_to_table(input_filename,output_filename)

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_filename)
    blank_samplename = feature_table.columns[-3]
    logger.info('Blank sample name for exclusion: '+ blank_samplename)

    # User-defined parameters
    logger.info('User-defined parameters')
    logger.info('   Minimum ion intensity treshold (count) = '+ str(min_intensity))
    logger.info('   Additional margin for retention time range exclusion (seconds) = '+ str(rtexclusionmargininsecs))

    # Concatenating the tables from narrow and large features:
    df_narrow = pd.read_csv(output_filename,sep=',')

    #We arbitrarly expand the exclusion range from +/- X seconds
    df_narrow['rt_start'] = df_narrow['rt_start'] - rtexclusionmargininsecs
    df_narrow['rt_end'] = df_narrow['rt_end'] + rtexclusionmargininsecs

    #Concatening the tables
    df_narrow.to_csv(output_filename, sep=',', index=False)

    # Running the table processing
    logger.info('Processing the table')
    make_exclusion_list(output_filename, blank_samplename, min_intensity)
    logger.info('======')

    # Convert to XCalibur format
    logger.info('Preparing list of excluded ions in XCalibur format')
    generate_QE_list_rt_range(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, output_filename[:-4]+'_EXCLUSION_LIST_XCalibur.csv')
    logger.info('======')
    logger.info('Preparing list of excluded ions in MaxQuant.Live format')
    generate_MQL_exclusion(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, output_filename[:-4]+'_EXCLUSION_LIST_MaxQuantLive.txt')
    logger.info('======')


    # === Plot the features  ====
    logger.info('Preparing scatter plots of the excluded ions/features')
    plot_targets_exclusion_range(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, 'Distribution of the excluded ions in the blank sample')
    plot_targets_exclusion(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, 'retention_time', 'Distribution of intensities for the excluded ions in '+ blank_samplename)
    plot_targets_exclusion(output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, 'Mass [m/z]', 'Distribution of intensities for the excluded ions in '+ blank_samplename)

    logger.info('=======================')
    logger.info('Zipping workflow results files')

    # Cleaning files first
    os.system('mkdir results/intermediate_files')
    os.system('mv results/'+output_filename[:-4]+'_EXCLUSION_BLANK.csv intermediate_files/')
    os.system('mv results/'+output_filename+' intermediate_files/')
    os.system('mv results/table* results/intermediate_files/')

    get_all_file_paths('results','download_results/IODA_exclusion_results.zip')
    logger.info('=======================')
    logger.info('End the IODA-exclusion workflow')
    logger.info('=======================')
    logger.info('Proceed below with the results visualization')


if __name__ == "__main__":
    make_exclusion_from_mzTab(str(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]))
