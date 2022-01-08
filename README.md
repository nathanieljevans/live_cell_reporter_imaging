# Live Cell Reporter Imaging
---

> Project in Dr. Gordon [Mills Lab](https://www.ohsu.edu/school-of-medicine/mills-lab/people) at `Oregon Health and Science University` led by `Dr. Samuel Tsang`, data analysis by `Nate Evans` (evansna@ohsu.edu). 
---

To run analysis... 

0. Download data from [one drive](https://ohsuitg-my.sharepoint.com/:f:/g/personal/evansna_ohsu_edu/EnjH6gPtrgNOnFm_jayHw-cBPrOAeFQQUSfXTwAinVl_fg?email=tsangsa%40ohsu.edu&e=JlR2PQ). Email `evansna@ohsu.edu` for access.
        
  
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
