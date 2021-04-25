#!/bin/bash

_nclus_="2 3 5 10 15 25 35 45"

_resample_sz_="50 75 100 125 150"

_load_="normalized raw" 


for load in $_load_; do
    echo $load
    for nclus in $_nclus_; do
        echo -e '\t' $nclus
        for resample_sz in $_resample_sz_; do
            echo -e '\t\t' $resample_sz
            python HER2_classifier.py --data ./data/HER2/ --drug neratinib --sensitive_line WT --resistant_line T798I --load $load --nclus $nclus --out ./output/ --resample_sz $resample_sz
        done
    done 
done 