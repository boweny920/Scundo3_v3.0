#!/usr/bin/env python

import shutil
import argparse
import pandas as pd
import sys
import os 

parser = argparse.ArgumentParser()
parser.add_argument("-s",'--samplesheet')
parser.add_argument('-o','--outdir_of_secundo')
parser.add_argument('-a', '--annotation_Version')
parser.add_argument('-r', '--reference_genome', default=None)
parser.add_argument('-i', '--indexDir')
args = parser.parse_args()

df_samplesheet = pd.read_csv(args.samplesheet,sep="\t", dtype=str)

pi_names = list(set(v for v in df_samplesheet["pi_name"].to_list()))
if len(pi_names) > 1:
    sys.exit(f"ERROR: More than 1 PI name in the order, PI Names {pi_names}")
pi_name = pi_names[0]

requester_names = list(set(v for v in df_samplesheet["requester_name"].to_list()))
if len(requester_names) > 1:
    sys.exit(f"ERROR: More than 1 PI name in the order, PI Names {requester_names}")
requester_name = requester_names[0]

molng_ID = list(set(v for v in df_samplesheet["molngID"].to_list()))
molng_ID = molng_ID[0]

if args.reference_genome != None:
    ref_genome = args.reference_genome
    species = str(ref_genome).split("/")[0]
    genome_ver = str(ref_genome).split("/")[1]

else: 
    species = list(set(v for v in df_samplesheet["species"].to_list()))
    if len(species) > 1:
        sys.exit(f"ERROR: More than 1 species name in the order, species Names {species}")
    species = species[0]
    
    genome_ver = list(set(v for v in df_samplesheet["genome_ver"].to_list()))
    if len(genome_ver) > 1:
        sys.exit(f"ERROR: More than 1 genome_ver name in the order, genome_ver Names {genome_ver}")
    genome_ver = genome_ver[0]
    
    ref_genome = f"{species}/{genome_ver}"

output_dir = f"{args.outdir_of_secundo}/{pi_name}/{requester_name}/{molng_ID}.{genome_ver}.{args.annotation_Version}"
if os.path.exists(output_dir) == False:
    os.makedirs(output_dir)
    
### Now that the output dir is set, lets copy the necessary index files under it! ###

gene_data_txt = f"{args.indexDir}/{ref_genome}/annotation/{args.annotation_Version}/tables/{genome_ver}.{args.annotation_Version}.gene_data.txt"
if os.path.isfile(gene_data_txt) == False:
    sys.exit(f"ERROR: gene_data.txt file: {gene_data_txt}, does not exist")
    
GenomeBpTypes_txt = f"{args.indexDir}/{ref_genome}/annotation/{args.annotation_Version}/extras/{genome_ver}.{args.annotation_Version}.GenomeBpTypes.txt"
if os.path.isfile(GenomeBpTypes_txt) == False:
    sys.exit(f"ERROR: gene_data.txt file: {GenomeBpTypes_txt}, does not exist")

shutil.copy2(gene_data_txt, output_dir)
shutil.copy2(GenomeBpTypes_txt, output_dir)

print(output_dir.rstrip())