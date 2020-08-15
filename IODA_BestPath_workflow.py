
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
def make_BP_baseline_list_from_mzTab(input_filename:int, num_path:int, intensity_ratio:float, intensity_threshold:float, win_len:float, isolation:float):
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
    
    print(output_filename)
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
    logger.info('Minimum intensity for ion filtering in sample = '+ str(min_intensity))
    logger.info('Retention time window (sec.) for binning target ions = ' +str(win_len))
    logger.info('isolation window (m/z) = ' +str(isolation))
    experiements = num_path
    logger.info('Number of iterative experiment(s) = ' + str(experiements))
    logger.info('======')

    # Running the table processing
    logger.info('Running Path Finder  ...')
    print(output_filename[:-4]+'_PathFinder.csv')
    run_pathfinder_baseline(output_filename, output_filename[:-4]+'_PathFinder.csv', intensity_threshold, intensity_ratio, num_path, win_len, isolation)
    logger.info('======')
    logger.info('Running the table processing ...')
    #make_bestpath_targeted_lists_from_table(output_filename[:-4]+'_PathFinder.csv')
    logger.info('======')

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


    
### PathFinder    
# This parse one line from the BestPath output and create a output table per path. The rows to skip define which line/path is parsed.
def bestpath_format(input_filename: str, output_filename: str, rows_to_skip:int):
    df_path = pd.read_csv(input_filename, sep=' ', header=None, skiprows=rows_to_skip)

    #Make a list for the first row
    df_path_list = df_path.iloc[0].values.tolist()
    print('df_path_list')
    print(df_path_list)
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
    print('target_list')
    print(target_list)
    #Make a dataframe
    target_table = pd.DataFrame(target_list)
    target_table = target_table.rename(columns={0: 'Mass [m/z]',1: 'mz_isolation',2: 'duration',3: 'rt_start',4: 'rt_end', 5: 'intensity', 6: 'rt_apex',7: 'charge'})
    print('target_table.columns')
    print(target_table.columns)
    target_table = target_table[target_table['intensity'] > 0]

    print('For '+input_filename+', this path'+str(rows_to_skip+1)+' has number of valid targets = '+str(target_table.shape[0]))

    target_table.to_csv(output_filename, sep=',', index=False)

# This parse BestPath output file and create output tables formatted for XCalibur and MaxQuant.live
def make_bestpath_targeted_lists_from_table(input_filename:str):
    with open(input_filename) as file:
        counter = -1
        for line in file:
            try:
                counter += 1
                print("Processing path"+str(counter)+' will be renamed path'+str(counter+1))
                output_filename = input_filename[:-4]+"_"+str(counter+1)+'_formatted.txt'
                print(input_filename)
                print(output_filename)
                #Take the list and make a table
                logger.info('Formatting tables ...')
                bestpath_format(input_filename,output_filename, counter)
                logger.info('Converting tables to XCalibur format ...')
                generate_QE_list_from_BestPath(output_filename, 'XCalibur/'+output_filename[:-4]+'_QE_'+str(counter+1)+'.csv')
                #Format for MaxQuant.Live targeted experiment
                logger.info('Converting tables for MaxQuant.Live ...')
                generate_MQL_list_from_BestPath(output_filename, 'MQL/'+output_filename[:-4]+'_MQL_'+str(counter+1)+'.txt')
            except:
                raise
    
    logger.info('List of files ...')    
    table_list_bestpath = []
    for x in range(0,counter+1):
        table_list_bestpath.append(input_filename[:-4]+"_"+str(x+1)+'_formatted.txt')
    print(table_list_bestpath)
        
    logger.info('Plotting results ...')    
    make_plot_bestpath1(table_list_bestpath)
    make_plot_bestpath2(table_list_bestpath)
    make_plot_bestpath3(table_list_bestpath)

    
def run_pathfinder_baseline(input_filename:str, output_filename:str, intensity_threshold:float, intensity_ratio:float, num_path:int, win_len:float, isolation:float):
    os.system("sed -i 's/\t/ /g' "+input_filename)
    cmd_baseline = ('python3 path_finder.py baseline '+input_filename+' '+output_filename+' '+str(intensity_threshold)+' '+str(intensity_ratio)+' '+str(num_path)+' -win_len '+str(win_len)+' -isolation '+str(isolation))
    print(cmd_baseline)
    os.system(cmd_baseline)
    

#Best path generate mz / rt figures
def make_plot_bestpath1(table_list_bestpath):
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
    plt.savefig('Plots/'+output_filename[:-4]+'injection_scatter_plot_mz_rt.png', dpi=300)
    plt.close()
    
#Best path generate feature intensity / rt figures
def make_plot_bestpath2(table_list_bestpath):
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
    plt.savefig('Plots/'+output_filename[:-4]+'injection_scatter_plot_intensity_rt.png', dpi=300)
    plt.close()

    
    
#Best path generate feature intensity / duration
def make_plot_bestpath3(table_list_bestpath):
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

        plt.savefig('Plots/'+output_filename[:-4]+'injection1_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()

    if len(table_list_bestpath) >= 1:
        table1 = pd.read_csv(table_list_bestpath[1], sep=',', header=0)
        plt.scatter('duration','intensity',  data=table1, marker='o', color='violet',s=1.5, alpha=0.6)
        Label = ['Inj. 2, n = '+ str(table1.shape[0])+ ', median = '+ "{0:.2e}".format(table1['rt_apex'].median())  + ', mean = '+ "{0:.2e}".format(table1['rt_apex'].mean())]
        Labels.clear()
        Labels.append(Label)

        plt.savefig('Plots/'+output_filename[:-4]+'injection2_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()

    if len(table_list_bestpath) >= 2: 
        table2 = pd.read_csv(table_list_bestpath[2], sep=',', header=0)
        plt.scatter('duration','intensity',  data=table2, marker='o', color='orange',s=1, alpha=0.6)
        Label = ['Inj. 3, n = '+ str(table2.shape[0])+ ', median = '+ "{0:.2e}".format(table2['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table2['rt_apex'].mean())]
        Labels.clear()
        Labels.append(Label)

        plt.savefig('Plots/'+output_filename[:-4]+'injection3_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()

    if len(table_list_bestpath) >= 3: 
        table3 = pd.read_csv(table_list_bestpath[3], sep=',', header=0)
        plt.scatter('duration','intensity',  data=table3, marker='o', color='red', s=0.5, alpha=0.6)
        Label =['Inj. 4, n = '+ str(table3.shape[0])+ ', median = '+ "{0:.2e}".format(table3['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table3['rt_apex'].mean())]
        Labels.clear()
        Labels.append(Label)
        plt.savefig('Plots/'+output_filename[:-4]+'injection4_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()

    if len(table_list_bestpath) >= 4: 
        table4 = pd.read_csv(table_list_bestpath[4], sep=',', header=0)
        plt.scatter('duration','intensity',  data=table4, marker='o', color='red', s=0.1, alpha=0.6)
        Label =['Inj. 5, n = '+ str(table4.shape[0])+ ', median = '+ "{0:.2e}".format(table4['rt_apex'].median()) + ', mean = '+ "{0:.2e}".format(table4['rt_apex'].mean())]
        Labels.clear()
        Labels.append(Label)
        plt.savefig('Plots/'+output_filename[:-4]+'injection5_scatter_plot_intensity_duration.png', dpi=300)
        plt.close
        plt.clf()
        
    
if __name__ == "__main__":
    make_targeted_list_from_mzTab(str(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
