#!/bin/bash

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";
workflow_loc="$SCRIPT_DIR/workflow/"

#>>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/grid/compbio/software/mambaforge/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/grid/compbio/software/mambaforge/etc/profile.d/conda.sh" ]; then
        . "/grid/compbio/software/mambaforge/etc/profile.d/conda.sh"
    else
        export PATH="/grid/compbio/software/mambaforge/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/grid/compbio/software/mambaforge/etc/profile.d/mamba.sh" ]; then
    . "/grid/compbio/software/mambaforge/etc/profile.d/mamba.sh"
fi
#<<< conda initialize <<<

conda activate snakemake-miseq

if [ "$#" -lt 1 ]; then
    echo "Please specify the number of cores."
fi

snakemake --cores $1 --rerun-incomplete --use-conda --conda-prefix ${workflow_loc}/.snakemake/envs --snakefile ${workflow_loc}Snakefile.smk ${@:2}
