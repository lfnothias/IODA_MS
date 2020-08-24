import pandas as pd
import numpy as np
import sys
import os
from io import StringIO
import warnings
from pandas.core.common import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

# Make a list of XCalibur with retention time range for exclusion
def generate_QE_list_rt_range(input_table: str, blank_samplename:str, output_filename:str):
    """Format a table with mz, charge, rt_start, rt_end, intensities into a standard QExactive inclusion/exclusion list"""
    # Prepare the columns
    df_master = pd.read_csv(input_table, sep=',', header=0)
    df_master['Start [min]']=(df_master['rt_start'])/60
    df_master['End [min]']=(df_master['rt_end'])/60
    #Format to scientific notation the intensity
    df_master[blank_samplename] = df_master[blank_samplename].astype(float).map(lambda n: '{:.2E}'.format(n))

    #Build a comment (optional)
    df_master['block1'] = round(df_master['retention_time']*1/60,3)
    df_master['block2'] = round(df_master['retention_time'],2)
    df_master['block1'] = df_master['block1'].astype(str)
    df_master['block2'] = df_master['block2'].astype(str)
    df_master['for_comments'] = 'Apex = '+df_master['block1']+' (min) or '+df_master['block2']+' (sec), int. = '+ df_master[blank_samplename].astype(str)

    #Make the output table
    df = pd.DataFrame(data=None)
    df['Mass [m/z]'] = df_master["Mass [m/z]"].round(decimals=4)
    df['Formula [M]'] = ''
    df['Formula type'] = ''
    df['Species'] = '' # '=-H'
    df['CS [z]'] = df_master['charge']
    df['Polarity'] = 'Positive' # 'Negative'
    df["Start [min]"] = df_master["Start [min]"].round(decimals=3)
    df["End [min]"] = df_master["End [min]"].round(decimals=3)
    df['(N)CE'] = '' #Can be empty ONLY INCLUSION
    df['(N)CE type'] = '' #Can be empty ONLY INCLUSION
    df['MSX ID'] = '' #Can be empty ONLY INCLUSION
    df['Comment'] = df_master['for_comments']
    df.to_csv(output_filename, index = False, sep=',')

# Make an EXCLUSION list for MaxQuant.Live with apex retention time value
def generate_MQL_exclusion(input_table:str, blank_samplename:str, output_filename:str):
    """Format a table with mz, charge, rt, intensities into a standard MaxQuantLive list"""

    #Mass [m/z]',1: 'mz_isolation',2: 'duration',3: 'rt_start',4: 'rt_end',5: 'intensity'})
    df = pd.read_csv(input_table)
    df2 = df[['Mass [m/z]','retention_time','charge',blank_samplename]]
    df2.rename(columns={'Mass [m/z]':'Mass'}, inplace=True)
    df2['Mass']=df2['Mass'].round(decimals=5)
    df2['Retention time']= (df2['retention_time']/60)
    df2.rename(columns={blank_samplename:'Apex intensity'}, inplace=True)
    df2['Apex intensity']=df2['Apex intensity'].round(decimals=0)
    df2['placeholder'] = np.arange(len(df2)) + 1 #Mandatory for import
    df2['Modified sequence'] = np.arange(len(df2)) + 1 #Mandatory for import. Arbitrary string.
    df2['Charge'] = df2['charge']  #field is mandatory and cannot be 0
    df2['Charge'] = df2['Charge'].replace([0], 1) #MQL target list bugs if 0
    df2['Retention time'] = df2['Retention time'].round(decimals=4)
    df2['MaxIt'] = ''
    df2["Colission Energies"] = ''
    df2['RealtimeCorrection'] = 'TRUE' #Maybe a good idea to keep TRUE only for most intense ions
    df2['TargetedMs2'] = 'FALSE'
    df2['Targetedlabeled'] = 'FALSE'
    df2['TargetedMultiInjection'] = 'FALSE'
    df2['TopNExclusion'] = 'TRUE'
    df2['Fragments mz'] = ''
    df2['NCE Factors'] = ''

    df_out = df2[['placeholder','Modified sequence','Mass', 'Charge',\
                  'Retention time','Apex intensity','Fragments mz', 'MaxIt',\
                  'NCE Factors', 'Colission Energies','RealtimeCorrection','TargetedMs2',\
                  'Targetedlabeled','TargetedMultiInjection','TopNExclusion']]
    df_out.rename(columns={'placeholder':''}, inplace=True)

    df_out.to_csv(output_filename, index=None, sep='\t')

# For targeted experiment (no range, just apex)
def generate_QE_list(input_table: str, output_filename:str, pretarget_rt_exclusion_time:float, posttarget_rt_exclusion_time:float):
    """Format a table with mz, charge, rt, intensities into a standard QExactive inclusion/exclusion list"""
    # Prepare the columns
    df_master = pd.read_csv(input_table)
    df_master['Start [min]']=(df_master['retention_time']-(pretarget_rt_exclusion_time))/60
    df_master['End [min]']=(df_master['retention_time']+(posttarget_rt_exclusion_time))/60
    #Build a comment (optional) -> deactivate for lighter output
    #df_master['block1'] = round(df_master['retention_time']*1/60,3)
    #df_master['block2'] = round(df_master['retention_time'],2)
    #df_master['block1'] = df_master['block1'].astype(str)
    #df_master['block2'] = df_master['block2'].astype(str)
    #df_master['for_comments'] = 'Apex = '+df_master['block1']+' (min) '+df_master['block2']+' (sec)'

    #Make the output table
    df = pd.DataFrame(data=None)
    df['Mass [m/z]'] = df_master["Mass [m/z]"].round(decimals=4)
    df['Formula [M]'] = ''
    df['Formula type'] = ''
    df['Species'] = '' # '=-H'
    df['CS [z]'] = df_master['charge']
    df['Polarity'] = 'Positive' # 'Negative'
    df["Start [min]"] = df_master["Start [min]"].round(decimals=3)
    df["End [min]"] = df_master["End [min]"].round(decimals=3)
    df['(N)CE'] = '' #Can be empty ONLY INCLUSION
    df['(N)CE type'] = '' #Can be empty ONLY INCLUSION
    df['MSX ID'] = '' #Can be empty ONLY INCLUSION
    df['Comment'] = '' #df_master['for_comments'] # To add comments in the targeted table. \
    #Will increase the size of the table and slow down import in XCalibur

    df.to_csv(output_filename, index = False, sep=',')

# For targeted experiment
def generate_MQL_list(input_table:str, output_filename:str, pretarget_rt_exclusion_time:float, posttarget_rt_exclusion_time:float):
    """Format a table with mz, charge, rt, intensities into a standard MaxQuantLive list"""
    df = pd.read_csv(input_table)
    df2 = df[['Mass [m/z]', 'charge','retention_time']]
    df2.rename(columns={'Mass [m/z]':'Mass'}, inplace=True)
    df2.rename(columns={'charge':'Charge'}, inplace=True)
    df2.rename(columns={'retention_time':'Retention time'}, inplace=True)
    df2['Charge'] = df2['Charge'].replace([0], 1) #MQL target list bugs if 0
    # We are shifting the peak apex to allow asymetric retention time range
    df2['Retention time']= ((df2['Retention time']+pretarget_rt_exclusion_time)/60)
    df2['Retention time'] = df2['Retention time'].round(decimals=3)
    # We are setting the retention time range
    df2['Retention length']= ((pretarget_rt_exclusion_time+posttarget_rt_exclusion_time)/60)
    df2['Retention length'] = df2['Retention length'].round(decimals=3)
    df2["Intensity"] = '1'
    df2["Colission energies"] = ''
    df2['RealtimeCorrection'] = 'TRUE'
    df2['TargetedMs2'] = 'TRUE'
    df2['Targetedlabeled'] = 'FALSE'
    df2['SIM'] = 'FALSE'
    df2['ExcludeFromTopN'] = 'TRUE'
    df2['id'] = ''
    df2['placeholder'] = np.arange(len(df2)) + 1
    df2['Modified sequence'] = ''

    df_out = df2[['placeholder','id','Modified sequence','Mass', 'Charge',\
                  'Retention time','Retention length','Intensity',\
                  'Colission energies','RealtimeCorrection','TargetedMs2',\
                  'Targetedlabeled','SIM','ExcludeFromTopN']]
    df_out.rename(columns={'placeholder':''}, inplace=True)

    df_out.to_csv(output_filename, index=None, sep=',')

    

## PathFinder below
def generate_QE_list_from_MS2Planner(input_table: str, output_filename:str, rt_margin:float):
    """Format a table with mz, charge, rt_start, rt_end, intensities into a standard QExactive inclusion/exclusion list"""
    # Prepare the columns
    df_master = pd.read_csv(input_table, sep=',', header=0)
    # We extend the retention time range for target from X percent of the peak duration. Percent is user provided.
    # We extend by retention range start by X*0.25 percent and increase the retention end by X percent.
    df_master['Start [min]']=(df_master['rt_apex']-rt_margin)/60
    df_master['End [min]']=(df_master['rt_apex']+rt_margin)/60
    #Format to scientific notation the intensity
    df_master['intensity'] = df_master['intensity'].astype(float).map(lambda n: '{:.2E}'.format(n))

    #Build a comment (optional)
    df_master['block1'] = round(df_master['rt_apex']*1/60,3)
    df_master['block2'] = round(df_master['rt_apex'],2)
    df_master['block1'] = df_master['block1'].astype(str)
    df_master['block2'] = df_master['block2'].astype(str)
    df_master['for_comments'] = 'Apex = '+df_master['block1']+' (min) or '+df_master['block2']+' (sec), int. = '+ df_master['intensity'].astype(str)

    #Make the output table
    df = pd.DataFrame(data=None)
    df['Mass [m/z]'] = df_master["Mass [m/z]"].round(decimals=4)
    df['Formula [M]'] = ''
    df['Formula type'] = ''
    df['Species'] = '' # '=-H'
    df['CS [z]'] = df_master['charge']
    df['Polarity'] = 'Positive' # 'Negative'
    df["Start [min]"] = df_master["Start [min]"].round(decimals=3)
    df["End [min]"] = df_master["End [min]"].round(decimals=3)
    df['(N)CE'] = '' #Can be empty ONLY INCLUSION
    df['(N)CE type'] = '' #Can be empty ONLY INCLUSION
    df['MSX ID'] = '' #Can be empty ONLY INCLUSION
    df['Comment'] = df_master['for_comments']
    df.to_csv(output_filename, index = False, sep=',')


def generate_MQL_list_from_MS2Planner(input_table:str, output_filename:str, rt_margin:float):
    """Format a table with mz, charge, rt, intensities into a standard MaxQuantLive list"""

    #Mass [m/z]',1: 'mz_isolation',2: 'duration',3: 'rt_start',4: 'rt_end',5: 'intensity'})
    df = pd.read_csv(input_table)
    df2 = df[['Mass [m/z]','duration','rt_start','rt_end','intensity','rt_apex','charge']]
    df2.rename(columns={'Mass [m/z]':'Mass'}, inplace=True)
    df2['Mass']=df2['Mass'].round(decimals=5)
    # This is to set the RT acquisition window
    df2['Retention time']= (df2['rt_apex']-rt_margin)/60  #MQL requires in minutes
    df2['Apex intensity']= df2['intensity'].round(decimals=0)
    df2['placeholder'] = np.arange(len(df2)) + 1 #Mandatory for import
    df2['Modified sequence'] = np.arange(len(df2)) + 1 #Mandatory for import. Arbitrary string.
    df2['Charge'] = df2['charge'] #field is mandatory and cannot be 0
    df2['Charge'] = df2['Charge'].replace([0], 1) #MQL target list bugs if 0
    df2['Retention time'] = df2['Retention time'].round(decimals=4)
    # We can't customised the MaxIT from the baseline or apex mode.
    df2['MaxIt'] = ''
    df2["Colission Energies"] = ''
    df2['RealtimeCorrection'] = 'TRUE' #Maybe a good idea to keep TRUE only for most intense ions
    df2['TargetedMs2'] = 'TRUE'
    df2['Targetedlabeled'] = 'FALSE'
    df2['TargetedMultiInjection'] = 'FALSE'
    df2['TopNExclusion'] = 'TRUE'
    df2['Fragments mz'] = ''
    df2['NCE Factors'] = ''

    df_out = df2[['placeholder','Modified sequence','Mass', 'Charge',\
                  'Retention time','Apex intensity','Fragments mz', 'MaxIt',\
                  'NCE Factors', 'Colission Energies','RealtimeCorrection','TargetedMs2',\
                  'Targetedlabeled','TargetedMultiInjection','TopNExclusion']]
    df_out.rename(columns={'placeholder':''}, inplace=True)
    
    df_out.to_csv(output_filename, index=None, sep='\t')
    
    
def generate_MQL_list_from_MS2Planner_MaxIT(input_table:str, output_filename:str, transient_time:float):
    """Format a table with mz, charge, rt, intensities into a standard MaxQuantLive list"""

    #Mass [m/z]',1: 'mz_isolation',2: 'duration',3: 'rt_start',4: 'rt_end',5: 'intensity'})
    df = pd.read_csv(input_table)
    df2 = df[['Mass [m/z]','duration','rt_start','rt_end','intensity','rt_apex','charge']]
    df2.rename(columns={'Mass [m/z]':'Mass'}, inplace=True)
    df2['Mass']=df2['Mass'].round(decimals=5)
    df2['Retention time']= df2['rt_start']/60  #MQL requires in minutes
    df2['Apex intensity']= df2['intensity'].round(decimals=0)
    df2['placeholder'] = np.arange(len(df2)) + 1 #Mandatory for import
    df2['Modified sequence'] = np.arange(len(df2)) + 1 #Mandatory for import. Arbitrary string.
    df2['Charge'] = df2['charge'] #field is mandatory and cannot be 0
    df2['Charge'] = df2['Charge'].replace([0], 1) #MQL target list bugs if 0
    df2['Retention time'] = df2['Retention time'].round(decimals=4)
    #The MaxIT is derived from the 'duration' column. We take only 95% of it. We substract the transient time (MS2 and overhead time).
    df2['MaxIt'] = (df2['duration']*1000*0.95)-transient_time #MQL MaxIT in ms
    df2["Colission Energies"] = ''
    df2['RealtimeCorrection'] = 'TRUE' #Maybe a good idea to keep TRUE only for most intense ions
    df2['TargetedMs2'] = 'TRUE'
    df2['Targetedlabeled'] = 'FALSE'
    df2['TargetedMultiInjection'] = 'FALSE'
    df2['TopNExclusion'] = 'TRUE'
    df2['Fragments mz'] = ''
    df2['NCE Factors'] = ''

    df_out = df2[['placeholder','Modified sequence','Mass', 'Charge',\
                  'Retention time','Apex intensity','Fragments mz', 'MaxIt',\
                  'NCE Factors', 'Colission Energies','RealtimeCorrection','TargetedMs2',\
                  'Targetedlabeled','TargetedMultiInjection','TopNExclusion']]
    df_out.rename(columns={'placeholder':''}, inplace=True)
    
    df_out.to_csv(output_filename, index=None, sep='\t') 