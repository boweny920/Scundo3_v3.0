#!/usr/bin/env python

import sys
import pandas as pd 
import os 
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
    df_sample_info.to_csv("lims_FCID_Info.csv",index=False)
    if molng == False:
        if len([v for v in df_sample_info["prnOrderNo"].unique()]) > 1:
            sys.exit(f"{fc_id} is associated with TWO or more MOLNG-IDs, please speicify a MOLNG-ID in the input.")
        molng_id = [v for v in df_sample_info["prnOrderNo"].unique()][0]
    else: 
        molng_id = molng 

    return molng_id 

def secundo3_fetch_lims(molng):
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
    df_run_info.to_csv("lims_MOLNG_Info.csv",index=False)
    df_nanalysis_result_path = pd.DataFrame()
    df_nanalysis_result_path["resultsPath"] = df_run_info["resultsPaths.unix"]
    df_nanalysis_result_path["pi_name"] =[v.replace("/n/analysis/","").split("/")[0] for v in df_run_info["resultsPaths.unix"]]
    df_nanalysis_result_path["requester_name"] =[v.replace("/n/analysis/","").split("/")[1] for v in df_run_info["resultsPaths.unix"]]
    
    return df_nanalysis_result_path

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
parser.add_argument('-o', '--output_dir', default="./", help="Output directory, default=cwd, default='./'")
parser.add_argument('-f','--fcid', help='provide fcid')
parser.add_argument("-m", '--molng', default=False, help='provide molng id, this can be just empty if the FCID has only 1 MOLNG-ID')
parser.add_argument("-i", '--RoboIndex_sampleSheet', default="/n/analysis/genomes/sampleSheet_ROBOINDEX_2023.csv")
args=parser.parse_args()

df_lims_info = secundo3_fetch_lims(fetch_lims_molngID(args.fcid, args.molng))
if len(df_lims_info) == 0: 
    sys.exit(f"Something is wrong with fetching {args.fcid} {args.molng} from LIMs")

df_sample_report = pd.DataFrame()
for index, row in df_lims_info.iterrows(): # Go through each flowcell to collect all data associated with on MOLNG-ID
    
    report_folder_dir = row["resultsPath"]
    pi_name = row["pi_name"]
    requester_name = row["requester_name"]
    
    if os.path.exists(os.path.join(report_folder_dir,"Sample_Report.csv")) == False:
        sys.exit(f'File {os.path.join(report_folder_dir,"Sample_Report.csv")} does not exist! Check primary analysis and or FCID/MOLNG-ID')
    df_sample_report_preCombine=pd.read_csv(os.path.join(report_folder_dir,"Sample_Report.csv"),dtype=str)
    df_sample_report = pd.concat([df_sample_report, df_sample_report_preCombine], axis=0)
    
library_ids = [v for v in df_sample_report["LibraryID"].unique()]

# Make a a dictionary of the species 
available_genomesAnnotations = {}
df_samplesheet = pd.read_csv(args.RoboIndex_sampleSheet)
for index, row in df_samplesheet.iterrows():
    # available_genomes.setdefault( str(row["name"]).replace(" ","_") + "_" + str(row["id"]), row["annotation_version"] )
    available_genomesAnnotations.setdefault( str(row["id"]), row["annotation_version"] )

sample_report_dic = {"libID":[],"molngID":[],'samplename':[],"secundoname":[], "species":[],"genome_ver":[],"fastqs":[], "pi_name":[], "requester_name":[], "readType":[], "annotation":[]}
for lib_ID in library_ids:
    sample_report_dic["libID"].append(lib_ID)
    
    ##Fetching info regarding fastqs, single reads and pair reads will be treated different##
    if len([v for v in df_sample_report[df_sample_report["LibraryID"]==lib_ID]["Type"].unique()]) > 1:
        sys.exit(f'Check LibraryID {lib_ID} sequencing type! More than one type detected')
    
    if [v for v in df_sample_report[df_sample_report["LibraryID"]==lib_ID]["Type"].unique()][0] != "Single Read":
        fastq_1 = ",".join([v for v in df_sample_report[(df_sample_report["LibraryID"]==lib_ID) & (df_sample_report["Read"]=="1")]["FastqPath"].to_list()]) #Read 1
        fastq_2 = ",".join([v for v in df_sample_report[(df_sample_report["LibraryID"]==lib_ID) & (df_sample_report["Read"]=="2")]["FastqPath"].to_list()]) #Read 2
        fastqs = " ".join([fastq_1, fastq_2]) # use space to join read 1 and read 2 fastqs
        readType = "PairEnd"
    else:
        fastqs = ",".join([v for v in df_sample_report[df_sample_report["LibraryID"]==lib_ID]["FastqPath"].to_list()]) 
        readType = "SingleEnd"
    sample_report_dic["fastqs"].append(fastqs)
    sample_report_dic["readType"].append(readType)

    refernce = [v for v in df_sample_report[df_sample_report["LibraryID"]==lib_ID]["Reference"].unique()]
    if len(refernce) > 1 :
        sys.exit(f"More than one genome version reference detected for library {lib_ID}") # IF more than 1 reference genome is specified for each library, then manual attention is required!! 
    refernce = refernce[0]
    sample_report_dic["genome_ver"].append(refernce)
    if refernce in available_genomesAnnotations:
        sample_report_dic["annotation"].append(available_genomesAnnotations[refernce])
    else : 
        sys.exit(f"The Genome version of this order {refernce} is not in current Reference Index collection.")

    species = [v for v in df_sample_report[df_sample_report["LibraryID"]==lib_ID]["Species"].unique()]
    if len(species) > 1 :
        sys.exit(f"More than one genome species detected for library {lib_ID}") # IF more than 1 reference genome is specified for each library, then manual attention is required!! 
    species = species[0]
    sample_report_dic["species"].append(species)

    sampleID = [v for v in df_sample_report[df_sample_report["LibraryID"]==lib_ID]["SampleName"].unique()]
    if len(sampleID) > 1: 
        sys.exit(f"More than one sample name: {sampleID} found associating with one lib ID: {lib_ID}")
    sample_report_dic["samplename"].append(sampleID[0])
        
    Secundo_ID = replace_special_characters(sampleID[0])
    sample_report_dic["secundoname"].append(Secundo_ID)
        
    molngID = [v for v in df_sample_report[df_sample_report["LibraryID"]==lib_ID]["Order"].unique()]
    molngID = "+".join(molngID)     # IF more than 1 MOLNG-ID is specified for each library, then join them with +, this should not be an issue 
    sample_report_dic["molngID"].append(molngID)

    sample_report_dic["pi_name"].append(pi_name) # As long as there is only one PI and requester for one MOLNG-ID, this will work file 
    sample_report_dic["requester_name"].append(requester_name)  # As long as there is only one PI and requester for one MOLNG-ID, this will work file 
        
df_info = pd.DataFrame(sample_report_dic)
outfile = os.path.join(args.output_dir, "samplesheet.tsv")
df_info.to_csv(outfile, index=False, sep="\t")