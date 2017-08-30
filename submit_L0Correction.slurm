#!/bin/bash
#Create a temporary job folder

#SBATCH -e localhost:/home/hep/matzeni/gangadir/Analysis/DoubleTuple/HDF/LogsDF_SLURM/SLURMLog3DF_${magnet}_${year}_${isplit}.err

#SBATCH -o localhost:/home/hep/matzeni/gangadir/Analysis/DoubleTuple/HDF/LogsDF_SLURM/SLURMLog3DF_${magnet}_${year}_${isplit}.out

echo "=================================================="
echo "Starting SLURM job for the creation of a joint tuple"
echo "Job ID: "$SLURM_JOB_ID


mkdir -p /scratch/$USER/$SLURM_JOB_ID

#Copy all the necessary files to the temporary job folder

cp $HOME/gangadir/Analysis/MCDataCorr/TriggerCorrection/Singlefile_L0TriggerCorrector.py /scratch/$USER/$SLURM_JOB_ID/


cd /scratch/$USER/$SLURM_JOB_ID/

echo "Python script used : Singlefile_L0TriggerCorrector.py"
echo "File: $ifile"
echo "Year: $year"
echo "Tables: $Tag_name"



echo "python Singlefile_L0TriggerCorrector.py -f $ifile -Tag $Tag_name -y $year > log3DFGIUSTO_$SLURM_JOB_ID.log"

echo "==========================================="
echo "SLURM log files saved in /home/hep/matzeni/gangadir/Analysis/DoubleTuple/HDF/LogsDF_SLURM"


# Clean up
# Removing all python  files from the scratch folder
rm -rf /scratch/$USER/$SLURM_JOB_ID/*.py*

mkdir -p $HOME/gangadir/Analysis/MCDataCorr/TriggerCorrection/Logs_L0
mv /scratch/$USER/$SLURM_JOB_ID/*.log $HOME/gangadir/Analysis/MCDataCorr/TriggerCorrection/Logs_L0

# Remove scratch folder
rm -rf /scratch/$USER/$SLURM_JOB_ID

echo "========================================="
echo "Ending SLURM job"
END_TIME=`date`
echo "End time: $END_TIME"
echo "Done."
