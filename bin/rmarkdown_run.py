#!/usr/bin/env python
import argparse
import subprocess
import os

parser=argparse.ArgumentParser()
parser.add_argument('-d','--secundo_dir')
parser.add_argument('-g','--Reference_Genome')
parser.add_argument('-a','--Annotation_Version')

args = parser.parse_args()

work_dir = str(args.secundo_dir).rstrip()

rmarkdown_file = os.path.join(work_dir,"analysis_rnaseq.rmd")

subprocess.run([f"R -e \"rmarkdown::render('{rmarkdown_file}', params = list(Reference_Genome = '{args.Reference_Genome}', Annotation_Version = '{args.Annotation_Version}'))\""],shell=True) # R -e "rmarkdown::render('sec_dir_path/analysis_rnaseq.rmd')"

