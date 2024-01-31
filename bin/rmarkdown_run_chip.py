#!/usr/bin/env python
import argparse
import subprocess
import os

parser=argparse.ArgumentParser()
parser.add_argument('-d','--secundo_dir')

args = parser.parse_args()

work_dir = str(args.secundo_dir).rstrip()
os.chdir(work_dir)

rmarkdown_file = os.path.join(work_dir,"analysis_chipseq.rmd")

subprocess.run([f"R -e \"rmarkdown::render('{rmarkdown_file}')\""], shell=True) 
# R -e "rmarkdown::render('sec_dir_path/analysis_rnaseq.rmd')"

