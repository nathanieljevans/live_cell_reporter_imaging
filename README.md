# Live Cell Reporter Imaging
---

> Project in Dr. Gordon Mills Lab at Oregon Health and Science University. Led by `Dr. Samuel Tsang`, data analysis by `Nate Evans` (evansna@ohsu.edu). 
---

# **TODO ITEM LIST**

11/10/21
Meeting with Samuel ... 
- add `Burn-in` hyper-parameter to sensitivity assignment 
- write necessary `batch effect checks` and `QC` 
- Switch from `SVM linear kernel classifier` to `gaussian process` to allow for non-binary classifications (sens/res/neither)
- need final results by **Late December or early Jan** 

Later analysis steps ... 
- compare cell viability ~ predicted resistance values 
- can we predict cell viability using cluster proportions? 
- EFM192A and SKBR3 behave differently in their resitance mechanisms .. why is this? 
- what is different? (targeted descriptive analysis)
- use CCLE to compare mutation or expression patterns between EFM/SKB... 
- Focus on HER2/ERB2 and family proteins / EGFR pathway 
- Use all cell lines for dimensionality reduction and clustering to identify cell lines with similar resistance patterns? 
- Use the achilles `gene dependency` data to look for differences in protein dependencies? 
- *Samuel action item:* Come up with list of important genes that may be involved in the difference between SKBR3 and EFM192A resistance patterns. 

NOTE: 
Very challenging to be statistically powered to compare expression patterns between EFM192A and SKBR3 ... HOWEVER ... we can use the common expression values of a gene across ALL cell lines to categorize each genes into up-regulated, down-regulated, or neither. Then we could simply compare ~diffferntial expression~ between SKBR3 and EFM192A by simply boolean logic... it would still be a descriptive analysis - but results would be much more meaningful... for instance: 

We would still want to use a targeted search (e.g., just interested in genes that are expected to be relevant) but now we could create a list of genes that have statistically significantly different expression (as measured across all cell lines) between EFM192A and SKBR3. 

** ... could use this same approach with `gene dependency` ... **

--- 

To run analysis... 

0. Download data from [box archive folder](https://ohsu.app.box.com/folder/149265669941) (email `evansna@ohsu.edu` for access). 

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

The aggregated results will be saved in the `--out` directory. There will be two files: 
- `experiment_run_results.csv` : contains the general results of each experiment; accuracy, inertia, etc ... 
- `mutant_resistance_results.csv` : contains the aggregated resistance probability assignments from each experiment. 




3. Use the `final_calls.ipynb` notebook to visualize results and make aggregate sensitivity calls. (not yet implemented). 
