#!/bin/bash
#SBATCH -c 4 #Number of cores
#SBATCH -o %j.txt
#SBATCH --mem=80G
#SBATCH -J upload_data

gsutil -m cp -r analysis/ gs://09152021_stanford_atac_gd2car_tommy/

echo "done"
