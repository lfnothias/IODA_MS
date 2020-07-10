
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

    print('Filenames in the mzTab')
    print(Filenames)

    Filename1 = Filenames[0]
    Filename2 = Filenames[1]

    # Display error message for additional samples
    for x in Filenames:
        if x == "ms_run[2]-location":
            Filename3 = Filenames[2]
            print['Warning: There is more than two samples in that mzTab file. We support only two samples currently']
        if x == "ms_run[4]-location":
            Filename4 = Filenames[3]
            print['Warning: There is more than three samples in that mzTab file. We support only two samples currently.']

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
    print('#Deducing the blank sample by comparing the sum of feature intensity between samples')
    column1_sum = df_master['peptide_abundance_study_variable[1]'].sum()
    print('- For sample '+Filename1+' the sum of feature intensities is = '+str(column1_sum))
    column2_sum = df_master['peptide_abundance_study_variable[2]'].sum()
    print('- For sample '+Filename2+' the sum of feature intensities = '+str(column2_sum))
    if column1_sum > column2_sum:
        print('- The blank sample is assumed to be '+str(Filename2)+' in the mzTab-M')
        print('- The samples is assumed to be '+str(Filename1)+' in the mzTab-M')
        df_master.rename(columns={'peptide_abundance_study_variable[1]':Filename2}, inplace=True)
        df_master.rename(columns={'peptide_abundance_study_variable[2]':Filename1}, inplace=True)
    if column1_sum < column2_sum:
        print('- The blank sample is assumed to be '+str(Filename1)+' in the mzTab-M')
        print('- The samples is assumed to be '+str(Filename2)+' in the mzTab-M')
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

def make_exclusion_list_blank(input_filename: str, sample: str, window: float):
    """From a table with mz, charge, rt, intensities, keep only features found in the sample specified"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_exclusion_list = df_master[(df_master[sample] != 0)]
    df_master_exclusion_list.to_csv(output_filename[:-4]+'_EXCLUSION_BLANK.csv', sep=',', index = False)
    #df_master_exclusion_list.sort_values(by=['Mass [m/z]'])
    print('Initial number of features ' + str(df_master.shape[0]))
    print('Number of features in the blank sample = ' + str(df_master_exclusion_list.shape[0]) +', with intensity != 0 ')

def make_exclusion_list_shared(input_filename: str, blank: str, sample: str, window: float):
    """From a table with mz, charge, rt, intensities, keep only features shared amongst the two samples specified"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_exclusion_list = df_master[(df_master[blank] != 0) & (df_master[sample] != 0)]
    df_master_exclusion_list.to_csv(output_filename[:-4]+'_EXCLUSION_SHARED.csv', sep=',', index = False)
    #df_master_exclusion_list.sort_values(by=['Mass [m/z]'])
    print('Initial number of features ' + str(df_master.shape[0]))
    print('Number of features shared between the blank and the sample = ' + str(df_master_exclusion_list.shape[0]) +', with intensity != 0 ')

def make_shotgun_targeted_list(input_filename: str, sample: str, window: float):
    """From a table with mz, charge, rt, intensities, keep only features found in the sample specified"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_shotgun_list = df_master[(df_master[sample] != 0)]
    df_master_shotgun_list.to_csv(output_filename[:-4]+'_SHOTGUN.csv', sep=',', index = False)
    print('Initial number of features ' + str(df_master.shape[0]))
    print('Number of features in the "shotgun" list = '+ str(df_master_shotgun_list.shape[0])+', with intensity != 0 ')

def make_targeted_list_ratio(input_filename: str, blank: str, sample: str, window:float, ratio:float):
    """From a table with mz, charge, rt, intensities, keep only features that have an intensity above the specified ratio between the sample/blank"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_targeted_list_ratio = df_master[(df_master[sample] > 0) & (df_master[sample]/df_master[blank] > ratio) & (df_master[blank] == 0)]
    df_master_targeted_list_ratio.to_csv(output_filename[:-4]+'_TARGETED_RATIO.csv', sep=',', index = False,)
    print('Initial number of features ' + str(df_master.shape[0]))
    print('Number of features in the targeted list = '+ str(df_master_targeted_list_ratio.shape[0])\
          +', with a ratio of Sample/Blank ratio of '+str(ratio))

def make_targeted_list_intensity(input_filename: str, blank: str, sample: str, window: str, intensity:float):
    """From a table with mz, charge, rt, intensities, keep only features that have an intensity above specified intensity"""
    df_master = pd.read_csv(input_filename, sep=',')
    df_master_targeted_list_intensity = df_master[(df_master[sample] > intensity) & (df_master[blank] < intensity) & (df_master[blank] == 0)]
    df_master_targeted_list_intensity.to_csv(output_filename[:-4]+'_TARGETED_INTENSITY.csv', sep=',', index = False,)
    print('Initial number of features ' + str(df_master.shape[0]))
    print('Number of features in the targeted list = ' + str(df_master_targeted_list_intensity.shape[0]) + ', with minimum intensity = '+ str(intensity))

def plot_targets_exclusion(input_filename: str,output_string: str, sample: str,title: str):
    """From a table, make a scatter plot of a sample"""
    #Plot
    Labels = []
    table0 = pd.read_csv(input_filename, sep=',', header=0)
    plt.scatter('Mass [m/z]', sample, data=table0, marker='o', color='blue',s=1.5, alpha=0.5)
    Label1 = ['Inj. 1, n = '+ str(table0.shape[0])+ ', median = '+ "{0:.2e}".format(table0[sample].median()) + ', mean = '+ "{0:.2e}".format(table0[sample].mean())]
    Labels.append(Label1)

    plt.yscale('log')

    plt.title(title)
    plt.xlabel('Ret. time (sec)')
    plt.ylabel('Feature intensity (log scale)')

    plt.legend(labels=Labels, fontsize =8)
    plt.savefig(output_filename[:-4]+'_EXCLUSION_'+output_string+'_scatter_plot.png', dpi=200)
    plt.clf()

def plot_targets_per_groups(table_list: str, output_string:str, sample: str, injections: int):
    """From a table, make a scatter plot of up to 4 samples"""
    #Plot
    Labels = []
    if injections >= 1:
        table0 = pd.read_csv(table_list[0], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table0, marker='o', color='blue',s=3, alpha=0.6)
        Label1 = ['Inj. 1, n = '+ str(table0.shape[0])+ ', median = '+ "{0:.2e}".format(table0[sample].median()) + ', mean = '+ "{0:.2e}".format(table0[sample].mean())]
        Labels.append(Label1)

    if injections >= 2:
        table1 = pd.read_csv(table_list[1], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table1, marker='o', color='violet',s=3, alpha=0.6)
        Label2 = ['Inj. 2, n = '+ str(table1.shape[0])+ ', median = '+ "{0:.2e}".format(table1[sample].median())  + ', mean = '+ "{0:.2e}".format(table1[sample].mean())]
        Labels.append(Label2)

    if injections >= 3:
        table2 = pd.read_csv(table_list[2], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table2, marker='o', color='orange',s=3, alpha=0.6)
        Label3 = ['Inj. 3, n = '+ str(table2.shape[0])+ ', median = '+ "{0:.2e}".format(table2[sample].median()) + ', mean = '+ "{0:.2e}".format(table2[sample].mean())]
        Labels.append(Label3)

    if injections >= 4:
        table3 = pd.read_csv(table_list[3], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table3, marker='o', color='red', s=3, alpha=0.6)
        Label4 =['Inj. 4, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table3[sample].median()) + ', mean = '+ "{0:.2e}".format(table3[sample].mean())]
        Labels.append(Label4)

    plt.yscale('log')
    plt.ylim(bottom=2)

    plt.title('Features per injection: '+ output_string)
    plt.xlabel('Ret. time (sec)')
    plt.ylabel('Feature intensity (log scale)')

    plt.legend(labels=Labels, fontsize =8)
    plt.savefig(output_filename[:-4]+'_injection_'+output_string+'_scatter_plot.png', dpi=300)
    plt.clf()

def plot_targets_per_groups_w_shared(table_list: str, output_string: str, input_filename_blank: str, sample: str, blank: str, injections: int):
    """From a table, make a scatter plot of up to 4 samples, and plot the blank too"""
    #Plot
    Labels = []
    if injections >= 1:
        table0 = pd.read_csv(table_list[0], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table0, marker='o', color='blue',s=1, alpha=0.6)
        Label1 = ['Inj. 1, n = '+ str(table0.shape[0])+ ', median = '+ "{0:.2e}".format(table0[sample].median()) + ', mean = '+ "{0:.2e}".format(table0[sample].mean())]
        Labels.append(Label1)

    if injections >= 2:
        table1 = pd.read_csv(table_list[1], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table1, marker='o', color='violet',s=1, alpha=0.6)
        Label2 = ['Inj. 2, n = '+ str(table1.shape[0])+ ', median = '+ "{0:.2e}".format(table1[sample].median())  + ', mean = '+ "{0:.2e}".format(table1[sample].mean())]
        Labels.append(Label2)

    if injections >= 3:
        table2 = pd.read_csv(table_list[2], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table2, marker='o', color='orange',s=1, alpha=0.6)
        Label3 = ['Inj. 3, n = '+ str(table2.shape[0])+ ', median = '+ "{0:.2e}".format(table2[sample].median()) + ', mean = '+ "{0:.2e}".format(table2[sample].mean())]
        Labels.append(Label3)

    if injections >= 4:
        table3 = pd.read_csv(table_list[3], sep=',', header=0)
        plt.scatter('Mass [m/z]', sample, data=table3, marker='o', color='red', s=1, alpha=0.6)
        Label4 =['Inj. 4, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table3[sample].median()) + ', mean = '+ "{0:.2e}".format(table3[sample].mean())]
        Labels.append(Label4)

    # Show shared features between blank and sample
    table_blank = pd.read_csv(input_filename_blank, sep=',', header=0)
    plt.scatter('Mass [m/z]', blank, data=table_blank, marker='o', color='black',s=1, alpha=0.5)
    Label2 = ['Blank, n = '+ str(table_blank.shape[0])+ ', median = '+ "{0:.2e}".format(table_blank[blank].median())  + ', mean = '+ "{0:.2e}".format(table_blank[blank].mean())]
    Labels.append(Label2)

    plt.yscale('log')
    plt.ylim(bottom=2)

    plt.title('Features per injection and shared with blank: '+ output_string)
    plt.xlabel('Ret. time (sec)')
    plt.ylabel('Feature intensity (log scale)')

    plt.legend(labels=Labels, fontsize =8)
    plt.savefig(output_filename[:-4]+'_injection_blank_shared_'+output_string+'_scatter_plot.png', dpi=300)
    plt.clf()


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

    print('All files zipped successfully!')

# Make targeted list from mzTab
def make_targeted_list_from_mzTab(input_filename:str, output_filename:str, injection_number:int, ratio_value:float, min_intensity_value:int):
    # Convert the mzTab into a Table
    print('Starting the workflow')
    print('======')
    print('Converting mzTab to table format')
    convert_mzTab_to_table(input_filename,output_filename)
    print('======')

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_filename)
    samplename = feature_table.columns[-1]
    print('Assumed sample name: '+samplename)
    blank_samplename = feature_table.columns[-2]
    print('Assumed blank sample name: ' +blank_samplename)
    print('======')

    # User-defined parameters
    print('User-defined parameters')
    ratio = ratio_value
    print('Ratio = ' + str(ratio))
    min_intensity = min_intensity_value
    print('Minimum intensity = '+ str(min_intensity))
    injections = injection_number
    print('Injections = ' + str(injections))
    print('======')

    # Hard coded parameters
    print('Hard coded parameters')
    window_exclusion = 30
    print('Window for exclusion = '+ str(window_exclusion))
    window_targeted = 5
    print('Window for targets = ' + str(window_targeted))
    window_bin = 30
    print('Window for binning of targets = ' +str(window_bin))
    print('======')

    # Running the table processing
    print('Running the table processing')
    make_exclusion_list_blank(output_filename, blank_samplename, window_exclusion)
    print('======')
    make_exclusion_list_shared(output_filename, blank_samplename, samplename, window_exclusion)
    print('======')
    make_shotgun_targeted_list(output_filename, samplename, window_targeted)
    print('======')
    make_targeted_list_ratio(output_filename, blank_samplename, samplename, window_targeted, ratio)
    print('======')
    make_targeted_list_intensity(output_filename, blank_samplename, samplename, window_targeted, min_intensity)
    print('======')

    # Split the tables for multiple injections
    print('Splitting the tables')
    from split_features import split_features
    #split_features(output_filename[:-4]+'_EXCLUSION_BLANK.csv', 'Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_BLANK.csv',samplename, window_bin, injections)
    #split_features(output_filename[:-4]+'_EXCLUSION_SHARED.csv','Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_SHARED.csv', samplename, window_bin, injections)
    split_features(output_filename[:-4]+'_SHOTGUN.csv', 'Intermediate_files/'+output_filename[:-4]+'_SHOTGUN.csv', samplename, window_bin, injections)
    split_features(output_filename[:-4]+'_TARGETED_RATIO.csv', 'Intermediate_files/'+output_filename[:-4]+'_TARGETED_RATIO.csv', samplename, window_bin, injections)
    split_features(output_filename[:-4]+'_TARGETED_INTENSITY.csv', 'Intermediate_files/'+output_filename[:-4]+'_TARGETED_INTENSITY.csv', samplename, window_bin, injections)
    print('======')

    print('Generating filename list')
    # Generate the filename list
    table_list = []
    for x in range(1,injections+1):
            table_list.append('Intermediate_files/'+output_filename[:-4]+'_TARGETED_INTENSITY_'+str(x)+'.csv')

    table_list_ratio = []
    for x in range(1,injections+1):
            table_list_ratio.append('Intermediate_files/'+output_filename[:-4]+'_TARGETED_RATIO_'+str(x)+'.csv')

    table_list_shotgun = []
    for x in range(1,injections+1):
            table_list_shotgun.append('Intermediate_files/'+output_filename[:-4]+'_SHOTGUN_'+str(x)+'.csv')
    print('======')

    # === OUTPUT FILES BELOW + LOG ====
    print('Plotting the features')
    plot_targets_exclusion('Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_BLANK.csv', blank_samplename, 'Plots/'+output_filename[:-4]+'_EXCLUSION_BLANK.csv','Distribution of features on the exclusion list')
    plot_targets_exclusion('Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_SHARED.csv', blank_samplename, 'Plots/'+output_filename[:-4]+'_EXCLUSION_SHARED.csv','Distribution of shared features between blank and sample')
    plot_targets_per_groups(table_list_ratio, 'SHOTGUN', samplename, injections)
    plot_targets_per_groups(table_list, 'TARGETED_INTENSITY', samplename, injections)
    plot_targets_per_groups(table_list_ratio, 'TARGETED_RATIO', samplename, injections)
    plot_targets_per_groups_w_shared(table_list,'TARGETED_INTENSITY', 'Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_SHARED.csv', samplename, blank_samplename,injections)
    plot_targets_per_groups_w_shared(table_list_ratio,'TARGETED_RATIO', 'Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_SHARED.csv', samplename, blank_samplename,injections)
    print('======')

    # Convert to XCalibur format
    print('Converting tables to XCalibur format')
    for x in range(1,injections+1):
            generate_QE_list('Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_BLANK.csv', 'XCalibur/Exclusion/'+output_filename[:-4]+'_EXCLUSION_BLANK_XCalibur_'+str(x)+'.csv', window_targeted)
            generate_QE_list('Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_SHARED.csv', 'XCalibur/Exclusion/'+output_filename[:-4]+'_EXCLUSION_SHARED_XCalibur_'+str(x)+'.csv', window_targeted)
            generate_QE_list('Intermediate_files/'+output_filename[:-4]+'_SHOTGUN_'+str(x)+'.csv', 'XCalibur/Targeted/'+output_filename[:-4]+'_SHOTGUN_XCalibur_'+str(x)+'.csv', window_targeted)
            generate_QE_list('Intermediate_files/'+output_filename[:-4]+'_TARGETED_RATIO_'+str(x)+'.csv', 'XCalibur/Targeted/'+output_filename[:-4]+'_TARGETED_RATIO_XCalibur_'+str(x)+'.csv', window_targeted)
            generate_QE_list('Intermediate_files/'+output_filename[:-4]+'_TARGETED_INTENSITY_'+str(x)+'.csv', 'XCalibur/Targeted/'+output_filename[:-4]+'_TARGETED_INTENSITY_XCalibur_'+str(x)+'.csv', window_targeted)
    print('======')

        # Convert the MaxQuant.Live format
    print('Converting tables to MaxQuant.Live format')
    for x in range(1,injections+1):
            generate_MQL_list('Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_BLANK.csv', 'MaxQuantLive/Exclusion/'+output_filename[:-4]+'_EXCLUSION_BLANK_MQL_'+str(x)+'.csv', window_targeted)
            generate_MQL_list('Intermediate_files/'+output_filename[:-4]+'_EXCLUSION_SHARED.csv', 'MaxQuantLive/Exclusion/'+output_filename[:-4]+'_EXCLUSION_SHARED_MQL_'+str(x)+'.csv', window_targeted)
            generate_MQL_list('Intermediate_files/'+output_filename[:-4]+'_SHOTGUN_'+str(x)+'.csv', 'MaxQuantLive/Targeted/'+output_filename[:-4]+'_SHOTGUN_MQL_'+str(x)+'.csv', window_targeted)
            generate_MQL_list('Intermediate_files/'+output_filename[:-4]+'_TARGETED_RATIO_'+str(x)+'.csv', 'MaxQuantLive/Targeted/'+output_filename[:-4]+'_TARGETED_RATIO_MQL_'+str(x)+'.csv', window_targeted)
            generate_MQL_list('Intermediate_files/'+output_filename[:-4]+'_TARGETED_INTENSITY_'+str(x)+'.csv', 'MaxQuantLive/Targeted/'+output_filename[:-4]+'_TARGETED_INTENSITY_MQL_'+str(x)+'.csv', window_targeted)

    print('======')
    print('Workflow terminated')
