#!/bin/bash
#SBATCH -J 01
#SBATCH -A uoa00029         # Project Account
#SBATCH --time=48:00:00     # Walltime
#SBATCH --mem-per-cpu=1024  # memory/cpu (in MB)
#SBATCH --cpus-per-task=1
###SBATCH --gres=gpu ## I could not get this part working in Apr-2019, and commented it out
#SBATCH --ntasks=1
###SBATCH --mail-type=ALL
# #SBATCH --mail-user=sknu003@aucklanduni.ac.nz

#SBATCH -o stdout_01.txt
#SBATCH -e stderr_01.txt

#load modules required
module purge

#module load Python/3.9.5-gimkl-2020a
#module load Python/3.8.2-gimkl-2020a

##module load R/4.1.0-gimkl-2020a
module load R/4.3.2-foss-2023a

#change directory to where the R code is stored
RCLIB=$(echo "Rcode_scripts")
cd $PWD
cd ../

#start the bash script that iterates over compressed gz files
#./Rcode01_get_spc_from_artspriotering.R
./${RCLIB}/Rcode11_get_spc_from_artspriotering_v02.R
