#!/bin/bash
#
#  Name this job.
#
#MOAB -N dome
#
#  Run this job in the "backfill" queue.
#
#MOAB -q backfill
#
#  Request 4 processors.
#
#MOAB -l nodes=1:ppn=4
#
#  Maximum wallclock time :( Hours : Minutes : Seconds )
#
#MOAB -l walltime=00:02:00
#
#  Join the two output files (standard output and standard error) into one.
#  Based on the job name, this means the single output file will be called
#  "dome.oJOBID" where JOBID is a job id number assigned when the job is
#  submitted.
#
#MOAB -j oe
#
#  This command is required to set up the Gnu version of OpenMPI.
#
module load gnu-openmpi
#
#  This command moves from your home directory to the directory
#  from which this script was submitted.
#
cd $PBS_O_WORKDIR 
#
#  Run the program
#
mpirun -np 4 ./Albany dome_input.xml
#
echo "Program terminated normally."

