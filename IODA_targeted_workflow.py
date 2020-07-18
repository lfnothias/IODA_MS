
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
    logger.info('- For sample '+Filename1+' the sum of feature intensities is = '+str(column1_sum))
    column2_sum = df_master['peptide_abundance_study_variable[2]'].sum()
    logger.info('- For sample '+Filename2+' the sum of feature intensities = '+str(column2_sum))
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

def make_exclusion_list_blank(input_filename: str, sample: str):
    """From a table with mz, charge, rt, intensities, keep only features found in the sample specified"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_exclusion_list = df_master[(df_master[sample] != 0)]
    df_master_exclusion_list.to_csv(input_filename[:-4]+'_EXCLUSION_BLANK.csv', sep=',', index = False)
    #df_master_exclusion_list.sort_values(by=['Mass [m/z]'])
    logger.info('Initial number of ions ' + str(df_master.shape[0]))
    logger.info('Number of ions in the blank = ' + str(df_master_exclusion_list.shape[0]) +', with int. != 0 ')

def make_exclusion_list_shared(input_filename: str, blank: str, sample: str):
    """From a table with mz, charge, rt, intensities, keep only features shared amongst the two samples specified"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_exclusion_list = df_master[(df_master[blank] != 0) & (df_master[sample] != 0)]
    df_master_exclusion_list.to_csv(input_filename[:-4]+'_EXCLUSION_SHARED.csv', sep=',', index = False)
    #df_master_exclusion_list.sort_values(by=['Mass [m/z]'])
    #logger.info('Initial number of ions ' + str(df_master.shape[0]))
    logger.info('Number of ions shared between blank and sample = ' + str(df_master_exclusion_list.shape[0]) +', with int. != 0 ')

def make_shotgun_targeted_list(input_filename: str, sample: str):
    """From a table with mz, charge, rt, intensities, keep only features found in the sample specified"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_shotgun_list = df_master[(df_master[sample] != 0)]
    df_master_shotgun_list.to_csv(input_filename[:-4]+'_SHOTGUN.csv', sep=',', index = False)
    #logger.info('Initial number of ions ' + str(df_master.shape[0]))
    logger.info('Number of ions in the sample (shotgun) = '+ str(df_master_shotgun_list.shape[0])+', with int.!= 0 ')

def make_targeted_list_ratio(input_filename: str, blank: str, sample: str, ratio:float):
    """From a table with mz, charge, rt, intensities, keep only features that have an intensity above the specified ratio between the sample/blank"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_targeted_list_ratio = df_master[(df_master[sample] > 0) & (df_master[sample]/df_master[blank] > ratio) & (df_master[blank] == 0)]
    df_master_targeted_list_ratio.to_csv(input_filename[:-4]+'_TARGETED_RATIO.csv', sep=',', index = False,)
    #logger.info('Initial number of ions = ' + str(df_master.shape[0]))
    logger.info('Number of target ions = '+ str(df_master_targeted_list_ratio.shape[0])\
          +', with a ratio of sample/blank ratio of '+str(ratio))

def make_targeted_list_intensity(input_filename: str, blank: str, sample: str, intensity:float):
    """From a table with mz, charge, rt, intensities, keep only features that have an intensity above specified intensity"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_targeted_list_intensity = df_master[(df_master[sample] > intensity) & (df_master[blank] < intensity) & (df_master[blank] == 0)]
    df_master_targeted_list_intensity.to_csv(input_filename[:-4]+'_TARGETED_INTENSITY.csv', sep=',', index = False,)
    #logger.info('Initial number of features ' + str(df_master.shape[0]))
    logger.info('Number of target ions = ' + str(df_master_targeted_list_intensity.shape[0]) + ', with minimum intensity = '+ str(intensity))

def plot_targets_exclusion(input_filename: str, blank_samplename: str, column: str, title: str):
    """From a table, make a scatter plot of a sample"""
    Labels = []
    table0 = pd.read_csv(input_filename, sep=',', header=0)
    fig = plt.figure(figsize=(8,6))
    fig = plt.scatter(column, blank_samplename, data=table0, marker='o', color='blue',s=4, alpha=0.4)
    Label1 = ['n = '+ str(table0.shape[0])+ ', median abs. int. = '+ "{0:.2e}".format(table0[blank_samplename].median()) + ', mean abs. int. = '+ "{0:.2e}".format(table0[blank_samplename].mean())]
    Labels.append(Label1)
    plt.yscale('log')
    if column == 'Mass [m/z]':
        plt.title(title+', in m/z range', size = 13)
        plt.xlabel('m/z', size = 12)
    if column == 'retention_time':
        plt.title(title+', in retention time range range', size =13)
        plt.xlabel('Ret. time (sec)', size = 11)
    plt.ylabel('Ion intensity (log scale)', size = 11)
    plt.legend(labels=Labels, fontsize =10)
    if column == 'Mass [m/z]':
        plt.savefig(input_filename[:-4]+'_excluded_MZ_scatter_plot.png', dpi=200)
    if column == 'retention_time':
        plt.savefig(input_filename[:-4]+'_excluded_RT_scatter_plot.png', dpi=200)
    plt.close()

def plot_targets_per_groups(output_filename:str, table_list: str, output_string:str, sample: str, experiments: int):
    """From a table, make a scatter plot of up to 4 samples"""
    #Plot
    Labels = []
    if experiments >= 1:
        table0 = pd.read_csv(table_list[0], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table0, marker='o', color='blue',s=3, alpha=0.6)
        Label1 = ['Exp. 1, n = '+ str(table0.shape[0])+ ', median = '+ "{0:.2e}".format(table0[sample].median()) + ', mean = '+ "{0:.2e}".format(table0[sample].mean())]
        Labels.append(Label1)

    if experiments >= 2:
        table1 = pd.read_csv(table_list[1], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table1, marker='o', color='violet',s=3, alpha=0.6)
        Label2 = ['Exp. 2, n = '+ str(table1.shape[0])+ ', median = '+ "{0:.2e}".format(table1[sample].median())  + ', mean = '+ "{0:.2e}".format(table1[sample].mean())]
        Labels.append(Label2)

    if experiments >= 3:
        table2 = pd.read_csv(table_list[2], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table2, marker='o', color='orange',s=3, alpha=0.6)
        Label3 = ['Exp. 3, n = '+ str(table2.shape[0])+ ', median = '+ "{0:.2e}".format(table2[sample].median()) + ', mean = '+ "{0:.2e}".format(table2[sample].mean())]
        Labels.append(Label3)

    if experiments >= 4:
        table3 = pd.read_csv(table_list[3], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table3, marker='o', color='red', s=3, alpha=0.6)
        Label4 =['Exp. 4, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table3[sample].median()) + ', mean = '+ "{0:.2e}".format(table3[sample].mean())]
        Labels.append(Label4)

    plt.yscale('log')
    plt.ylim(bottom=2)

    plt.title('Ions per experiment: '+ output_string)
    plt.xlabel('Ret. time (sec)')
    plt.ylabel('Ion intensity (log scale)')

    plt.legend(labels=Labels, fontsize =8)
    plt.savefig(output_filename[:-4]+'_experiment_'+output_string+'_scatter_plot.png', dpi=300)
    plt.close()

def plot_targets_per_groups_w_shared(output_filename:str, table_list: str, output_string: str, input_filename_blank: str, sample: str, blank: str, experiments: int):
    """From a table, make a scatter plot of up to 4 samples, and plot the blank too"""
    #Plot
    Labels = []
    if experiments >= 1:
        table0 = pd.read_csv(table_list[0], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table0, marker='o', color='blue',s=1, alpha=0.6)
        Label1 = ['Exp. 1, n = '+ str(table0.shape[0])+ ', median = '+ "{0:.2e}".format(table0[sample].median()) + ', mean = '+ "{0:.2e}".format(table0[sample].mean())]
        Labels.append(Label1)

    if experiments >= 2:
        table1 = pd.read_csv(table_list[1], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table1, marker='o', color='violet',s=1, alpha=0.6)
        Label2 = ['Exp. 2, n = '+ str(table1.shape[0])+ ', median = '+ "{0:.2e}".format(table1[sample].median())  + ', mean = '+ "{0:.2e}".format(table1[sample].mean())]
        Labels.append(Label2)

    if experiments >= 3:
        table2 = pd.read_csv(table_list[2], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table2, marker='o', color='orange',s=1, alpha=0.6)
        Label3 = ['Exp. 3, n = '+ str(table2.shape[0])+ ', median = '+ "{0:.2e}".format(table2[sample].median()) + ', mean = '+ "{0:.2e}".format(table2[sample].mean())]
        Labels.append(Label3)

    if experiments >= 4:
        table3 = pd.read_csv(table_list[3], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table3, marker='o', color='red', s=1, alpha=0.6)
        Label4 =['Exp. 4, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table3[sample].median()) + ', mean = '+ "{0:.2e}".format(table3[sample].mean())]
        Labels.append(Label4)

    # Show shared features between blank and sample
    table_blank = pd.read_csv(input_filename_blank, sep=',', header=0)
    plt.scatter('Mass [m/z]', blank, data=table_blank, marker='o', color='black',s=1, alpha=0.5)
    Label2 = ['Blank, n = '+ str(table_blank.shape[0])+ ', median = '+ "{0:.2e}".format(table_blank[blank].median())  + ', mean = '+ "{0:.2e}".format(table_blank[blank].mean())]
    Labels.append(Label2)

    plt.yscale('log')
    plt.ylim(bottom=2)
    plt.title('Ions per experiment (w. blank): '+ output_string, size =10)
    plt.xlabel('Ret. time (sec)')
    plt.ylabel('Ion intensity (log scale)')

    plt.legend(labels=Labels, fontsize =8, loc='best', markerscale=5)
    plt.savefig(output_filename[:-4]+'_experiment_blank_shared_'+output_string+'_scatter_plot.png', dpi=300)
    plt.savefig('experiment_blank_shared_'+output_string+'_scatter_view.png', dpi=300)
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

# Make targeted list from mzTab
def make_targeted_list_from_mzTab(input_filename:int, experiment_number:int, ratio:float, min_intensity_value:int, pretarget_rt_exclusion_time:float, posttarget_rt_exclusion_time:float, window_bin:int):
    os.system('rm -r results_targeted')
    os.system('rm download_results/IODA_targeted_results.zip')
    os.system('mkdir results_targeted')
    os.system('mkdir download_results')
    os.system('rm results/logfile.txt')
    logfile('results_targeted/logfile.txt')

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

    output_dir = 'results_targeted'
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
    logger.info('   Ratio between sample/blank intensity for ion filtering = ' + str(ratio))
    logger.info('   Minimum intensity for ion filtering in the sample = '+ str(min_intensity_value))
    logger.info('   Number of iterative experiment(s) = ' + str(experiment_number))
    logger.info('======')

    # Hard coded parameters
    logger.info('Retention time range parameters')
    rt_window_excluded_ion = pretarget_rt_exclusion_time + posttarget_rt_exclusion_time
    logger.info('   Excluded ion retention time (sec.) = ' + str(rt_window_excluded_ion))
    logger.info('   Pre-target ion retention time range (sec.) = ' + str(pretarget_rt_exclusion_time))
    logger.info('   Post-target ion retention time range (sec.) = ' + str(posttarget_rt_exclusion_time))

    #Parameter for split features
    logger.info('Retention time window (sec.) for binning target ions = ' +str(window_bin))
    logger.info('======')

    # Running the table processing
    logger.info('Running the table processing ...')
    make_exclusion_list_blank(output_filename, blank_samplename)
    logger.info('======')
    make_exclusion_list_shared(output_filename, blank_samplename, samplename)
    logger.info('======')
    make_shotgun_targeted_list(output_filename, samplename)
    logger.info('======')
    make_targeted_list_ratio(output_filename, blank_samplename, samplename, ratio)
    logger.info('======')
    make_targeted_list_intensity(output_filename, blank_samplename, samplename, min_intensity_value)
    logger.info('======')

    # Split the tables for multiple experiment_number
    logger.info('Splitting the tables')
    from IODA_split_features import split_features
    split_features(output_filename[:-4]+'_SHOTGUN.csv', output_filename[:-4]+'_SHOTGUN.csv', samplename, window_bin, experiment_number)
    split_features(output_filename[:-4]+'_TARGETED_RATIO.csv', output_filename[:-4]+'_TARGETED_RATIO.csv', samplename, window_bin, experiment_number)
    split_features(output_filename[:-4]+'_TARGETED_INTENSITY.csv', output_filename[:-4]+'_TARGETED_INTENSITY.csv', samplename, window_bin, experiment_number)
    logger.info('======')

    # Generate the filename list
    table_list = []
    for x in range(1,experiment_number+1):
            table_list.append(output_filename[:-4]+'_TARGETED_INTENSITY_'+str(x)+'.csv')

    table_list_ratio = []
    for x in range(1,experiment_number+1):
            table_list_ratio.append(output_filename[:-4]+'_TARGETED_RATIO_'+str(x)+'.csv')

    table_list_shotgun = []
    for x in range(1,experiment_number+1):
            table_list_shotgun.append(output_filename[:-4]+'_SHOTGUN_'+str(x)+'.csv')

    # === OUTPUT FILES BELOW + LOG ====
    logger.info('Plotting the ions ... please wait ...')
    plot_targets_exclusion(output_filename[:-4]+'_EXCLUSION_SHARED.csv', blank_samplename, 'retention_time', 'Intensity distribution of ions excluded')
    plot_targets_exclusion(output_filename[:-4]+'_EXCLUSION_SHARED.csv', blank_samplename, 'Mass [m/z]', 'Intensity distribution of ions excluded')
    plot_targets_per_groups(output_filename, table_list_ratio, 'SHOTGUN', samplename, experiment_number)
    plot_targets_per_groups(output_filename, table_list, 'TARGETED_INTENSITY', samplename, experiment_number)
    plot_targets_per_groups(output_filename, table_list_ratio, 'TARGETED_RATIO', samplename, experiment_number)

    logger.info('======')

    # Convert to XCalibur format
    logger.info('Converting tables to XCalibur format ...')
    for x in range(1,experiment_number+1):
            generate_QE_list(output_filename[:-4]+'_EXCLUSION_BLANK.csv', output_filename[:-4]+'_EXCLUSION_BLANK_XCalibur_exp_'+str(x)+'.csv', rt_window_excluded_ion/2, rt_window_excluded_ion/2)
            generate_QE_list(output_filename[:-4]+'_EXCLUSION_SHARED.csv', output_filename[:-4]+'_EXCLUSION_SHARED_XCalibur_exp_'+str(x)+'.csv', rt_window_excluded_ion/2, rt_window_excluded_ion/2)
            generate_QE_list(output_filename[:-4]+'_SHOTGUN_'+str(x)+'.csv', output_filename[:-4]+'_SHOTGUN_XCalibur_exp_'+str(x)+'.csv', pretarget_rt_exclusion_time, posttarget_rt_exclusion_time)
            generate_QE_list(output_filename[:-4]+'_TARGETED_RATIO_'+str(x)+'.csv', output_filename[:-4]+'_TARGETED_RATIO_XCalibur_exp_'+str(x)+'.csv', pretarget_rt_exclusion_time, posttarget_rt_exclusion_time)
            generate_QE_list(output_filename[:-4]+'_TARGETED_INTENSITY_'+str(x)+'.csv', output_filename[:-4]+'_TARGETED_INTENSITY_XCalibur_exp_'+str(x)+'.csv', pretarget_rt_exclusion_time, posttarget_rt_exclusion_time)
    logger.info('======')

        # Convert the MaxQuant.Live format
    logger.info('Converting tables to MaxQuant.Live format ...')
    for x in range(1,experiment_number+1):
            generate_MQL_list(output_filename[:-4]+'_EXCLUSION_BLANK.csv', output_filename[:-4]+'_EXCLUSION_BLANK_MaxQuantLive_exp_'+str(x)+'.csv', 0, rt_window_excluded_ion)
            generate_MQL_list(output_filename[:-4]+'_EXCLUSION_SHARED.csv', output_filename[:-4]+'_EXCLUSION_SHARED_MaxQuantLive_exp_'+str(x)+'.csv', 0, rt_window_excluded_ion)
            generate_MQL_list(output_filename[:-4]+'_SHOTGUN_'+str(x)+'.csv', output_filename[:-4]+'_SHOTGUN_MaxQuantLive_exp_'+str(x)+'.csv', pretarget_rt_exclusion_time, posttarget_rt_exclusion_time)
            generate_MQL_list(output_filename[:-4]+'_TARGETED_RATIO_'+str(x)+'.csv', output_filename[:-4]+'_TARGETED_RATIO_MaxQuantLive_exp_'+str(x)+'.csv', pretarget_rt_exclusion_time, posttarget_rt_exclusion_time)
            generate_MQL_list(output_filename[:-4]+'_TARGETED_INTENSITY_'+str(x)+'.csv', output_filename[:-4]+'_TARGETED_INTENSITY_MaxQuantLive_exp_'+str(x)+'.csv', pretarget_rt_exclusion_time , posttarget_rt_exclusion_time))

    logger.info('======')
    logger.info('Cleaning and zipping workflow results files ...')

    # Cleaning files first

    #mkdir XCalibur
    os.system('mkdir results_targeted/XCalibur')
    os.system('mkdir results_targeted/XCalibur/exclusion')
    os.system('mkdir results_targeted/XCalibur/shotgun')
    os.system('mkdir results_targeted/XCalibur/targeted_ratio')
    os.system('mkdir results_targeted/XCalibur/targeted_intensity')
    # mv files XCalibur
    os.system('mv results_targeted/*EXCLUSION_BLANK_XCalibur* results_targeted/XCalibur/exclusion')
    os.system('mv results_targeted/*EXCLUSION_SHARED_XCalibur* results_targeted/XCalibur/exclusion')
    os.system('mv results_targeted/*SHOTGUN_XCalibur* results_targeted/XCalibur/shotgun')
    os.system('mv results_targeted/*TARGETED_INTENSITY_XCalibur* results_targeted/XCalibur/targeted_intensity')
    os.system('mv results_targeted/*TARGETED_RATIO_XCalibur* results_targeted/XCalibur/targeted_ratio')

    #mkdir XCalibur
    os.system('mkdir results_targeted/MaxQuantLive')
    os.system('mkdir results_targeted/MaxQuantLive/exclusion')
    os.system('mkdir results_targeted/MaxQuantLive/shotgun')
    os.system('mkdir results_targeted/MaxQuantLive/targeted_ratio')
    os.system('mkdir results_targeted/MaxQuantLive/targeted_intensity')
    # mv files XCalibur
    os.system('mv results_targeted/*EXCLUSION_BLANK_MaxQuantLive* results_targeted/MaxQuantLive/exclusion')
    os.system('mv results_targeted/*EXCLUSION_SHARED_MaxQuantLive* results_targeted/MaxQuantLive/exclusion')
    os.system('mv results_targeted/*SHOTGUN_MaxQuantLive* results_targeted/MaxQuantLive/shotgun')
    os.system('mv results_targeted/*TARGETED_INTENSITY_MaxQuantLive* results_targeted/MaxQuantLive/targeted_intensity')
    os.system('mv results_targeted/*TARGETED_RATIO_MaxQuantLive* results_targeted/MaxQuantLive/targeted_ratio')

    # mkdir intermediate files
    os.system('mkdir results_targeted/intermediate_files')
    os.system('mkdir results_targeted/intermediate_files/converted')
    os.system('mkdir results_targeted/intermediate_files/exclusion')
    os.system('mkdir results_targeted/intermediate_files/shotgun')
    os.system('mkdir results_targeted/intermediate_files/targeted_ratio')
    os.system('mkdir results_targeted/intermediate_files/targeted_intensity')
    os.system('mkdir results_targeted/plots')
    # mv plots
    os.system('mv results_targeted/*SHOTGUN_scatter_plot* results_targeted/plots')
    os.system('mv results_targeted/scatter_plot* results_targeted/plots')
    os.system('mv results_targeted/*TARGETED_INTENSITY_scatter_plot* results_targeted/plots')
    os.system('mv results_targeted/*TARGETED_RATIO_scatter_plot* results_targeted/plots')
    os.system('mv experiment_blank_shared_TARGETED_RATIO_scatter_view.png results_targeted/intermediate_files/')
    os.system('mv experiment_blank_shared_TARGETED_INTENSITY_scatter_view.png results_targeted/intermediate_files/')
    # mv intermediate files
    os.system('mv results_targeted/*EXCLUSION_BLANK* results_targeted/intermediate_files/exclusion')
    os.system('mv results_targeted/*EXCLUSION_SHARED* results_targeted/intermediate_files/exclusion')
    os.system('mv results_targeted/*SHOTGUN* results_targeted/intermediate_files/shotgun')
    os.system('mv results_targeted/*TARGETED_INTENSITY* results_targeted/intermediate_files/targeted_intensity')
    os.system('mv results_targeted/*TARGETED_RATIO* results_targeted/intermediate_files/targeted_ratio')

    # mv plots
    os.system('rm shotgun')
    os.system('mv '+output_filename+' results_targeted/intermediate_files/converted')
    os.system('mv results_targeted/logfile.txt results_targeted/intermediate_files/')

    get_all_file_paths('results_targeted','download_results/IODA_targeted_results.zip')

    logger.info('======')
    logger.info('END OF THE IODA-targeted-from-mzTab WORKFLOW')
    logger.info('======')
    print(' ')

if __name__ == "__main__":
    make_targeted_list_from_mzTab(str(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6]))
