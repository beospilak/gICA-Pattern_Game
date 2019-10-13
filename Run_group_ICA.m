% Input directory for the GIFT toolbox and directory for the Batch file
% containing all the parameters to run gICA

clear
clear all;
addpath('directory\gift_toolbox')
rehash path;
icatb_check_path;
icatb_batch_file_run('directory\Batch_for_group_ICA.m');