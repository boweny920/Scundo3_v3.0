#!/usr/bin/env python

import sys
import pandas as pd 
import os 
import argparse
import requests
import json

parser = argparse.ArgumentParser()
# parser.add_argument('-s', '--sample_report', default=None, help= "Input Samplereport.csv associated to order, generated by primary analysis pipeline")
parser.add_argument('-o', '--output_dir', default="./", help="Output directory, default=cwd, default='./'")
parser.add_argument('-s', '--samplesheet')
args=parser.parse_args()

df_sample_info = pd.read_csv(args.samplesheet, sep="\t")

lims_info_dic = {"prnOrderNo":f"{df_sample_info['molngID'].tolist()[0]}",
                 "orderType":f"Public_RNAseq_data",
                 "genome":f"{df_sample_info['genome_ver'].tolist()[0]}",
                 "analysisGoals":f"Public_RNAseqData_analysis",
                 "readLength":f"CheckSRATable",
                 "readType":f"{df_sample_info['readType'].tolist()[0]}"}

with open("lims_order.json", "w") as newfile:
    json.dump(lims_info_dic, newfile)