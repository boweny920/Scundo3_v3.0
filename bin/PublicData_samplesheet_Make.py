#!/usr/bin/env python

import sys
import pandas as pd 
import os 
import argparse
import re

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
parser.add_argument('-s','--samplesheet', help='provide samplesheet output from sra_download.dry_run.r')
parser.add_argument("-i" , '--RoboIndex_samplesheet')
parser.add_argument("-l" , '--lab')
parser.add_argument("-r" , '--requester')
args=parser.parse_args()

df_sample_report = pd.read_csv(args.samplesheet, dtype=str)

available_species = {}
available_genomesAnnotations = {}
df_robosamplesheet = pd.read_csv(args.RoboIndex_samplesheet)
for index, row in df_robosamplesheet.iterrows():
    available_species.setdefault( str(row["id"]), str(row["name"]))
    available_genomesAnnotations.setdefault( str(row["id"]), str(row["annotation_version"]) )

library_ids = [v for v in df_sample_report["srr_id"].unique()] #LibraryID here is SRA ID
sample_report_dic = {"libID":[],"molngID":[],'samplename':[],"secundoname":[], "species":[],"genome_ver":[],"fastqs":[], "pi_name":[], "requester_name":[], "readType":[], "annotation":[]}

for lib_ID in library_ids:
    sample_report_dic["libID"].append(lib_ID)
    
    ##Fetching info regarding fastqs, single reads and pair reads will be treated different##
    if len([v for v in df_sample_report[df_sample_report["srr_id"]==lib_ID]["read_type"].unique()]) > 1:
        sys.exit(f'Check LibraryID {lib_ID} sequencing type! More than one read type detected')
    
    if [v for v in df_sample_report[df_sample_report["srr_id"]==lib_ID]["read_type"].unique()][0] == "Paired":
        fastq_1 = ",".join([v for v in df_sample_report[(df_sample_report["srr_id"]==lib_ID) & (df_sample_report["read"]=="1")]["fq_path"].to_list()]) #Read 1
        fastq_2 = ",".join([v for v in df_sample_report[(df_sample_report["srr_id"]==lib_ID) & (df_sample_report["read"]=="2")]["fq_path"].to_list()]) #Read 2
        fastqs = " ".join([fastq_1, fastq_2]) # use space to join read 1 and read 2 fastqs
        readType = "PairEnd"
    else:
        fastqs = ",".join([v for v in df_sample_report[df_sample_report["srr_id"]==lib_ID]["fq_path"].to_list()]) 
        readType = "SingleEnd"
    sample_report_dic["fastqs"].append(fastqs)
    sample_report_dic["readType"].append(readType)
    
    refernce = [v for v in df_sample_report[df_sample_report["srr_id"]==lib_ID]["ref_genome"].unique()]
    if len(refernce) > 1 :
        sys.exit(f"More than one genome version reference detected for library {lib_ID}") # IF more than 1 reference genome is specified for each library, then manual attention is required!! 
    refernce = refernce[0]
    sample_report_dic["genome_ver"].append(refernce)
    
    if refernce in available_species:
        sample_report_dic["species"].append(available_species[refernce])
        sample_report_dic["annotation"].append(available_genomesAnnotations[refernce])
    else:
        sys.exit(f"ref_genome {refernce} not in Reference index collection")

    sampleID = [v for v in df_sample_report[df_sample_report["srr_id"]==lib_ID]["sample_name"].unique()]
    if len(sampleID) > 1: 
        sys.exit(f"More than one sample name: {sampleID} found associating with one lib ID: {lib_ID}")
    sample_report_dic["samplename"].append(sampleID[0])
    
    Secundo_ID = replace_special_characters(sampleID[0])
    sample_report_dic["secundoname"].append(Secundo_ID)
    
    molngID = [v for v in df_sample_report[df_sample_report["srr_id"]==lib_ID]["request_id"].unique()]
    molngID = "+".join(molngID)     # IF more than 1 MOLNG-ID is specified for each library, then join them with +, this should not be an issue 
    sample_report_dic["molngID"].append(molngID)
    
    sample_report_dic["pi_name"].append(str(args.lab)) # As long as there is only one PI and requester for one MOLNG-ID, this will work file 
    sample_report_dic["requester_name"].append(str(args.requester)) 

df_info = pd.DataFrame(sample_report_dic)
outfile = os.path.join(args.output_dir, "Secundo3_SampleSheet_PublicData.tsv")
df_info.to_csv(outfile, index=False, sep="\t")