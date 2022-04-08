# Live Cell Reporter Imaging
---
[![DOI](https://zenodo.org/badge/309575585.svg)](https://zenodo.org/badge/latestdoi/309575585)


> Project in Dr. Gordon [Mills Lab](https://www.ohsu.edu/school-of-medicine/mills-lab/people) at `Oregon Health and Science University` led by `Dr. Samuel Tsang`, data analysis by `Nate Evans` (evansna@ohsu.edu). 
---

To run analysis... 

0. Download data from [one drive](https://ohsuitg-my.sharepoint.com/:f:/g/personal/evansna_ohsu_edu/EnjH6gPtrgNOnFm_jayHw-cBPrOAeFQQUSfXTwAinVl_fg?email=tsangsa%40ohsu.edu&e=JlR2PQ). Email `evansna@ohsu.edu` for access.
    - if you are extracting the relevant files from Samuel's file structure, then use `./src/HER2_extract_data2.sh`
  
1. Unpack zip files `./data/` 

2. build and activate the conda environment

```bash 
# build conda env if necessary 
conda env create -f environment.yml 

# activate env 
conda activate lc_reporter 
```

2. Run analysis bash script (Note: this can take quite a while; 8+ hours): 

```bash 
$ cd ./src/
$ ./HER2_sensitivity_runs.sh
```

See `HER2_sensitivity_runs.sh` for input/output paths/ 


3. Aggregate the results by: 

```bash 

$ python agg_results.py --input ../output/ --out ../output/

```

3. Use the `sensitivity_analysis.ipynb` notebook to visualize sensitivity analysisresults.

4. Finalize results (aggregate across sensitivity analysis runs) by: 

```bash 
$ python finalize_results.py --input ../output/ --output ../output/
```

This script aggregates all results by (cell line, treatment, mutant, batch) and merges the batch effect flag.

Will save 5 files: 

    - final_results.csv -- (all results)
    - final_results-sorted-EFM192A;NERATINIB.csv
    - final_results-sorted-EFM192A;TRASTUZUMAB.csv
    - final_results-sorted-SKBR3;NERATINIB.csv
    - final_results-sorted-SKBR3;TRASTUZUMAB.csv
    
 Results are sorted by the mean probability of resistance across all sensitivity analysis runs. 

---


# Data Dictionary 

## `final_results....csv`

applies to: 

    - final_results.csv 
    - final_results-sorted-EFM192A;NERATINIB.csv
    - final_results-sorted-EFM192A;TRASTUZUMAB.csv
    - final_results-sorted-SKBR3;NERATINIB.csv
    - final_results-sorted-SKBR3;TRASTUZUMAB.csv

> **cell_line**: SKBR3 or EFM192A 	

> **treatment**: neratinib or trastuzumab 

> **mutant**: specific mutation induced in the HER2 gene 

> **batch**: dataset source 

> **mean**: average probability of resistance, aggregated across all sensitivity analysis runs 

> **std**: standard deviation of probability of resistance, aggregated across all sensitivity analysis runs 

> **q05**: 5th percentile of probability of resistance, aggregated across all sensitivity analysis runs 

> **q95**: 95th percentile of probability of resistance, aggregated across all sensitivity analysis runs 

> **min**: minimum value of probability of resistance, aggregated across all sensitivity analysis runs 

> **max**: maximum value of probability of resistance, aggregated across all sensitivity analysis runs 

> **any_flag**: proportion of batch effect flags across all sensitivity analysis runs 


## `mutant_resistance_results.csv` 

> **pc1**: First principle component, latent representation of a cell line's response to drug (many cell's reporter values over time)

> **pc2**: Second principle component

> **treatment**: The drug treatment used in a given condition. 

> **mutant**: The cell lines induced HER2 mutation used in a given condition

> **batch**: The experimental batch identifier. Each batch has it's own set of sensitive and resistant controls. 

> **cell_count**: The number of cells measured and recorded for each mutant in a given batch and treatment. Note, there will be a unique cell count for each unique (batch,mutant,treatment) observation.

> **prob_res**: Predicted probability of resistance [as defined by the resitant controls]. 

> **prob_sens**: Predicted probability of sensitivty [as defined by the sensitive controls (WT, treated)]. 

> **call**: [sens, res] Most probable call (binary)

> **run_id**: Unique analysis identifier; necessary for distinguishing results from different analysis with different hyperparameters in the sensitivity analysis. 

> **cell_line**: The cell line used in a given condition (SKBR3, EFM193A)

> **drug_check**: QC check that merge performed properly. 

> **pc1_coef**: Batch PCx adjustment; If batch correction is used, then $pcx = pcx\_uncor - pcx\_coef$. Nan if batch correction not used.  

> **pc1_pval**: Batch effect significance, as measured by ANOVA [Fit regression using controls: PCx ~ batch].

> **pc2_coef**: Batch PCx adjustment; If batch correction is used, then $pcx = pcx\_uncor - pcx\_coef$. Nan if batch correction not used.  

> **pc2_pval**: Batch effect significance, as measured by ANOVA [Fit regression using controls: PCx ~ batch].

> **pc1_uncor**: The uncorrected PCx values. Nan if batch correction not used. 

> **pc2_uncor**: The uncorrected PCx values. Nan if batch correction not used. 

## `experiment_run_results.csv`

> **accuracy(train)**: The classifier accuracy on training data (sensitive and negative controls, treated); Describes how separable the controls are. 

> **pc1_var**: Amount of variance explained by PCx

> **pc2_var**: Amount of variance explained by PCx

> **kmeans_inertia**: The inertia (within cluster sum of squares) of the time-series k-means fit. 

> **res_line**: Resistant `mutant` used in a given analysis; Note: "line" is a misnomer, it does not refer to SKBR3 or EFM192A.

> **sens_line**: Sensitive `mutant` used in a given analysis; Note: "line" is a misnomer, it does not refer to SKBR3 or EFM192A.

> **drug**: The drug used, synonymous with `treatment` in other results `.csv`s. 

> **nclus**: The number of clusters used in the time-series k-means fit. 

> **resample_sz**: The time series resample size used. 

> **load**: [normalized, raw] The data form used in the analysis. TODO: Link to description of normalization. 

> **burnin**: The burn-in used in the analysis - this is the number of measured points removed from the beginning of each cell lines time series. Intended to alieve batch effects and remove the common high-to-low reporter behavior observed at the beggining of an experiment. 

> **batch_corrected**: [True, False] whether PC-space batch correction was applied.

> **run_id**: Unique analysis identifier; necessary for distinguishing results from different analysis with different hyperparameters in the sensitivity analysis. 

> **cell_line**: [SKBR3, EFM192A] The cell line used for a given condition. 

> **drug_check**: QC check that merge performed properly. 

## `batch_effect_results.csv`

> **batch**: The experimental batch identifier. Each batch has it's own set of sensitive and resistant controls. 

> **pc1_coef**: Batch PCx adjustment; Fit regression using controls: PCx ~ batch. 

> **pc1_pval**: Batch PCx significance via ANOVA; Fit regression using controls: PCx ~ batch. 

> **pc2_coef**: Batch PCx adjustment; Fit regression using controls: PCx ~ batch. 

> **pc2_pval**: Batch PCx significance via ANOVA; Fit regression using controls: PCx ~ batch. 

> **run_id**: Unique analysis identifier; necessary for distinguishing results from different analysis with different hyperparameters in the sensitivity analysis. 
