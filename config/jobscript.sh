#!/bin/sh
# properties = {properties}

# Make sure the conda install is on the path
#__conda_setup="$(CONDA_REPORT_ERRORS=false '/home/camw/miniconda3/bin/conda' shell.bash hook 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    \eval "$__conda_setup"
#else
#    if [ -f "/home/camw/miniconda3/etc/profile.d/conda.sh" ]; then
#        . "/home/camw/miniconda3/etc/profile.d/conda.sh"
#        CONDA_CHANGEPS1=false conda activate base
#    else
#        \export PATH="/home/camw/miniconda3/bin:$PATH"
#    fi
#fi
#unset __conda_setup



# unload the system python to avoid wonkiness
# module unload python

# activate the boink environment
export PATH="/home/camw/miniconda3/bin:$PATH"
source activate goetia

{exec_job}

cat $1
