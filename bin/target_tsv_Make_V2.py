#!/usr/bin/env python

import sys
import pandas as pd 
import argparse
import requests
import re

def fetch_lims_molngID(fc_id, molng):
    """
    Get the run info from LIMS. Returns the Molng-id associated to the fcid
    """
    # Get run info from lims
    NGS_LIMS = 'http://lims.stowers.org/zanmodules/molecular-biology/ngs'
    API_TOKEN = 'ca7952666a03dd4e59d0cd59e39fecc7' # Should get a new one for each pipeline, or even each user

    header = {'x-zan-apitoken': f'{API_TOKEN}', 
        'Accept': 'application/json'}
    run_info = requests.get(f'{NGS_LIMS}/flowcells/{fc_id}/samples', headers=header, verify=False)
    if not run_info.ok: 
        sys.exit('Malformed API request. Please double check your flowcell ID')
    
    lims_data = run_info.json() # Close request
    run_info.close() 
    df_run_info = pd.DataFrame.from_dict(lims_data)
    df_sample_info = df_run_info['samples'].apply(pd.Series)
    if molng == False:
        if len([v for v in df_sample_info["prnOrderNo"].unique()]) > 1:
            sys.exit(f"{fc_id} is associated with TWO or more MOLNG-IDs, please speicify a MOLNG-ID in the input.")
        molng_id = [v for v in df_sample_info["prnOrderNo"].unique()][0]
    else: 
        molng_id = molng 

    return molng_id 

def secundo3_fetch_lims_target_tsv(molng):
    NGS_LIMS = 'http://lims.stowers.org/zanmodules/molecular-biology/ngs'
    API_TOKEN = 'ca7952666a03dd4e59d0cd59e39fecc7' 
    header = {'x-zan-apitoken': f'{API_TOKEN}', 
        'Accept': 'application/json'}
    
    run_info = requests.get(f'{NGS_LIMS}/requests/{molng}/flowcells', headers=header, verify=False) # This is just for molngid as key 
    if not run_info.ok: 
        sys.exit('Malformed API request. Please double check your MOLNG ID') 
    lims_data = run_info.json() # Close request
    run_info.close()
    
    df_run_info = pd.json_normalize(lims_data) 
    return df_run_info

def fetch_sample_table_from_LIMS(fc_id, molng):
    """
    Get the run info from LIMS. Returns the /n/analysis/ folder where the primary analsysis pipeline copied all the data
    """
    # Get run info from lims
    NGS_LIMS = 'http://lims.stowers.org/zanmodules/molecular-biology/ngs'
    API_TOKEN = 'ca7952666a03dd4e59d0cd59e39fecc7' # Should get a new one for each pipeline, or even each user

    header = {'x-zan-apitoken': f'{API_TOKEN}', 
        'Accept': 'application/json'}
    run_info = requests.get(f'{NGS_LIMS}/flowcells/{fc_id}/samples', headers=header, verify=False)
    if not run_info.ok: 
        sys.exit('Malformed API request. Please double check your flowcell ID')
    
    lims_data = run_info.json() # Close request
    run_info.close() 
    df_run_info = pd.DataFrame.from_dict(lims_data)
    df_sample_info = df_run_info['samples'].apply(pd.Series)
    if molng != False:
        df_sample_info = df_sample_info[df_sample_info["prnOrderNo"]==molng]

    return df_sample_info

def replace_special_characters(input_string):
    # Define a regex pattern to match special characters
    pattern = r'[^a-zA-Z0-9\-\+]'  # Match any character that is not a letter, digit, "-", or "+"
    
    # Replace special characters with "_"
    result = re.sub(pattern, '_', input_string)
    
    # Replace "+" with "_pos_"
    result = result.replace('+', '_pos_')
    
    # Convert the result to lowercase
    result = "s_" + result.lower()
    
    return result

parser = argparse.ArgumentParser()
# parser.add_argument('-s', '--sample_report', default=None, help= "Input Samplereport.csv associated to order, generated by primary analysis pipeline")
parser.add_argument('-o', '--output_dir', default="./", help="Output directory, default=cwd, default='./'")
parser.add_argument('-l','--lims', help='provide fcid')
parser.add_argument("-m", '--molng', default=False, help='provide molng id, this can be just empty if the FCID has only 1 MOLNG-ID')
args=parser.parse_args()

df_run_info = secundo3_fetch_lims_target_tsv(fetch_lims_molngID(args.lims, args.molng))
fcid_list = [v for v in df_run_info["ngsCellId"].unique()]

df_target_tsv_comb = pd.DataFrame()

for fcid in fcid_list:
    df_sample_info = fetch_sample_table_from_LIMS(fcid, args.molng)

    order_list = df_sample_info["prnOrderNo"].to_list()
    sampleName = df_sample_info["sampleName"].to_list()
    secundoName = [replace_special_characters(v) for v in sampleName]
    flowCell = [fcid] * len(order_list)
    libraryID = df_sample_info["libID"].to_list()
    laneID = df_sample_info["laneId"].to_list()
    bardcodes = ["-".join(s['sequence'] for s in v) for v in df_sample_info["indexes"].to_list()] # This is to consider the situations where two indexes are available

    df_target_tsv = pd.DataFrame()
    df_target_tsv["Order"] = order_list
    df_target_tsv["Sample Name"] = sampleName
    df_target_tsv["Secundo Name"] = secundoName
    df_target_tsv["Flowcell"] = flowCell
    df_target_tsv["Library"] = libraryID
    df_target_tsv["Lane"] = laneID
    df_target_tsv["Barcode"] = bardcodes

    df_target_tsv_comb = pd.concat([df_target_tsv,df_target_tsv_comb], ignore_index=True)

df_target_tsv_comb.to_csv("targets.tsv",sep="\t",index=False)
