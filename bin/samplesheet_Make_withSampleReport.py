#!/usr/bin/env python

import sys
import pandas as pd 
import os 
import argparse
import requests
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
parser.add_argument('-s','--Sample_Reports', type=str, help='provide a list of Sample_Report.csv tables seperated by "," ')
parser.add_argument('-r', '--requester')
parser.add_argument('-l', '--lab')
parser.add_argument("-i", '--RoboIndex_sampleSheet', default="/n/analysis/genomes/sampleSheet_ROBOINDEX_2023.csv")
args=parser.parse_args()

df_sample_report = pd.DataFrame()
for SampleReport in str(args.Sample_Reports).split(","): # Go through each Provided Sample report csv tables
    if os.path.exists(SampleReport) == False:
        sys.exit(f'File {SampleReport} does not exist!')
    df_sample_report_preCombine=pd.read_csv(SampleReport, dtype=str)
    df_sample_report = pd.concat([df_sample_report, df_sample_report_preCombine], axis=0)
    
library_ids = [v for v in df_sample_report["LibraryID"].unique()]
pi_name = args.lab
requester_name = args.requester

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
        # sys.exit(f"The Genome version of this order {refernce} is not in current Reference Index collection.")
        sample_report_dic["annotation"].append("UnKnown")

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