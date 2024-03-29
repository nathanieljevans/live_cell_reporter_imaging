#!/bin/bash

output=../output/

# C:\Users\natha\local\live_cell_reporter_imaging\data\HER2_SKBR3_data_6-7-21\HER2_SKBR3_data_6-7-21
skbr3_data=../data/HER2_SKBR3_data_04-07-22/

# C:\Users\natha\local\live_cell_reporter_imaging\data\HER2_EFM192A_data_11-1-21\HER2_EFM192A_data_11-1-21
efm192A_data=../data/HER2_EFM192A_data_v2/

# Sensitivity analysis params -----------------------------
_batch_correct_='True False'            # option for PC-space batch correction (centering controls)
_repeats_="1"                           # option for replicates, write as "1 2 3"
_nclus_="5 15 25"                       # 2-50
_resample_sz_="50 100 150"              # 50 - 150 
_load_="raw normalized"                 # normalized raw
_burnin_="0 8"                          # these are in units of measurement time points so multiply by measurement 
                                            # frequency (usually 15 min) to get burnin time - example intended to be 0h, 1h, 2hr
# ----------------------------------------------------------

echo 'overwriting output directories...'
rm -r $output
mkdir $output
mkdir $output/SKBR3_NERATINIB/
mkdir $output/SKBR3_TRASTUZUMAB/
mkdir $output/EFM192A_NERATINIB/
mkdir $output/EFM192A_TRASTUZUMAB/


echo '#############################################################################'
echo 'SKBR3 + NERATINIB'
echo 'sensitive line = WT'
echo 'resistant line = T798I'
echo '#############################################################################'

for _ in $_repeats_; do
    for load in $_load_; do
        echo $load
        for nclus in $_nclus_; do
            echo -e '\t' $nclus
            for burnin in $_burnin_; do
                echo -e '\t' $burnin
                for resample_sz in $_resample_sz_; do
                    echo -e '\t\t' $resample_sz
                    for correct in $_batch_correct_; do
                        echo -e '\t\t' $correct
                        python HER2_classifier.py --data $skbr3_data --drug neratinib --sensitive_line WT --resistant_line T798I --load $load --nclus $nclus --out $output/SKBR3_NERATINIB/ --resample_sz $resample_sz --burnin $burnin --batch_correct $correct
                    done
                done
            done
        done 
    done
done 

echo '#############################################################################'
echo 'SKBR3 + TRASTUZUMAB'
echo 'sensitive line = WT'
echo 'resistant line = ND611'
echo '#############################################################################'

for _ in $_repeats_; do
    for load in $_load_; do
        echo $load
        for nclus in $_nclus_; do
            echo -e '\t' $nclus
            for burnin in $_burnin_; do
                echo -e '\t' $burnin
                for resample_sz in $_resample_sz_; do
                    echo -e '\t\t' $resample_sz
                    for correct in $_batch_correct_; do
                        echo -e '\t\t' $correct
                        python HER2_classifier.py --data $skbr3_data --drug trastuzumab --sensitive_line WT --resistant_line ND611 --load $load --nclus $nclus --out $output/SKBR3_TRASTUZUMAB/ --resample_sz $resample_sz --burnin $burnin --batch_correct $correct
                    done
                done
            done
        done 
    done
done

echo '#############################################################################'
echo 'EFM192A + NERATINIB'
echo 'sensitive line = WT'
echo 'resistant line = T798I'
echo '#############################################################################'

for _ in $_repeats_; do
    for load in $_load_; do
        echo $load
        for nclus in $_nclus_; do
            echo -e '\t' $nclus
            for burnin in $_burnin_; do
                echo -e '\t' $burnin
                for resample_sz in $_resample_sz_; do
                    echo -e '\t\t' $resample_sz
                    for correct in $_batch_correct_; do
                        echo -e '\t\t' $correct
                        python HER2_classifier.py --data $efm192A_data --drug neratinib --sensitive_line WT --resistant_line T798I --load $load --nclus $nclus --out $output/EFM192A_NERATINIB/ --resample_sz $resample_sz --burnin $burnin --batch_correct $correct
                    done
                done
            done
        done 
    done
done

echo '#############################################################################'
echo 'EFM192A + TRASTUZUMAB'
echo 'sensitive line = WT'
echo 'resistant line = ND611'
echo '#############################################################################'

for _ in $_repeats_; do
    for load in $_load_; do
        echo $load
        for nclus in $_nclus_; do
            echo -e '\t' $nclus
            for burnin in $_burnin_; do
                echo -e '\t' $burnin
                for resample_sz in $_resample_sz_; do
                    echo -e '\t\t' $resample_sz
                    for correct in $_batch_correct_; do
                        echo -e '\t\t' $correct
                        python HER2_classifier.py --data $efm192A_data --drug trastuzumab --sensitive_line WT --resistant_line ND611 --load $load --nclus $nclus --out $output/EFM192A_TRASTUZUMAB/ --resample_sz $resample_sz --burnin $burnin --batch_correct $correct
                    done
                done
            done
        done 
    done
done
