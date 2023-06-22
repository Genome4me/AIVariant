#!/bin/bash

while getopts i:e:t:n:r:o:g:d: flag
do
    case "${flag}" in
        i) input_env=${OPTARG};;
        e) eval_env=${OPTARG};;
        t) tumor=${OPTARG};;
        n) normal=${OPTARG};;
        r) reference=${OPTARG};;
        o) outdir=${OPTARG};;
        g) hg=${OPTARG};;
        d) depth=${OPTARG};;
    esac 
done 

source activate base;

binpath=$CONDA_PREFIX/bin;

echo $binpath

source $binpath/deactivate;
source $binpath/activate $input_env;

python generate_input.py -t $tumor -n $normal -r $reference -o $outdir -g $hg -d $depth;

source $binpath/deactivate;
source $binpath/activate $eval_env;

python evaluate.py -t $tumor -n $normal -r $reference -o $outdir -g $hg -d $depth;

source $binpath/deactivate;