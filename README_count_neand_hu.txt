#This is the README document for the scripts: 

#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'count_neand_hu.py' - cleaned .py script
*'count_neand_hu_expl.py' - .py script with more explanations

-these scripts mostly rely on 'pandas'
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'count_ext_files.sh' - .sh script used to run samtools to get the mpileup files, in interactive mode, then run .py script on those mpileup files and output results in a log file
*'count_ext_files_loop.sh' - same as above, but with bash loop
*'count_jar_files_nosam.sh' - .sh script used to run the .py script on mpileup files, in interactive mode, output results in a log file

-didn't use these scripts with 'slurm' because they were short enough to be implemented in interactive mode
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'ext_NDresults' - the results from the 'count_ext_files.sh' script
*'jar_NDresults' - the results from the 'count_jar_files_nosam.sh' script
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'results_count_neand_hu.txt' - extracted results, in a readable format
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'LJ_count_neand_hu.md' - the learning journal of the script
#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'sbatch_bam_to_mpileup.sh' - script for 'slurm' to get .mpileup files from .bam files
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'cropped_files' folder - contains the first 5000 lines of each file to create or test the scripts
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'Py_Neand_comp_hu.ipynb' - a jupyterlab file that contains details about the approach and some code
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
*'Where_is_everything' - a txt file with the paths to the files needed for the script (from Uppmax)
