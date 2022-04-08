'''
This script aggregates all results by (cell line, treatment, mutant, batch) and merges the batch effect flag.

(lc_reporter) $ python finalize_results.py --input ../output/ --output ../output/

Will save 5 files: 
    - final_results.csv 
    - final_results-sorted-EFM192A;NERATINIB
    - final_results-sorted-EFM192A;TRASTUZUMAB
    - final_results-sorted-SKBR3;NERATINIB
    - final_results-sorted-SKBR3;TRASTUZUMAB

Data Dictionary: 

cell_line: SKBR3 or EFM192A 	
treatment: neratinib or trastuzumab 
mutant: specific mutation induced in the HER2 gene 
batch: dataset source 
mean: average probability of resistance, aggregated across all sensitivity analysis runs 
std: standard deviation of probability of resistance, aggregated across all sensitivity analysis runs 
q05: 5th percentile of probability of resistance, aggregated across all sensitivity analysis runs 
q95: 95th percentile of probability of resistance, aggregated across all sensitivity analysis runs 
min: minimum value of probability of resistance, aggregated across all sensitivity analysis runs 
max: maximum value of probability of resistance, aggregated across all sensitivity analysis runs 
any_flag: proportion of batch effect flags across all sensitivity analysis runs 
'''

import argparse 
import pandas as pd 
import os

def get_args(): 
    parser = argparse.ArgumentParser(description='Finalize results of live cell reporter analysis.')
    parser.add_argument('--input', help='path to output directory.')
    parser.add_argument('--output', help='output directory to save aggregated results to.')
    return parser.parse_args()

# 5th Percentile
def q05(x):
    return x.quantile(0.05)

# 95th Percentile
def q95(x):
    return x.quantile(0.95)

if __name__ == '__main__': 
    args = get_args() 
    path = args.input

    mut_res = pd.read_csv(f'{args.input}/mutant_resistance_results.csv')
    batch_res = pd.read_csv(f'{args.input}/batch_effect_results.csv')

    # assign flags 
    batch_res2 = batch_res.assign(PC1_flag = lambda x: x.pc1_pval < 0.05,
                                  PC2_flag = lambda x: x.pc2_pval < 0.05)
    batch_res2 = batch_res2.assign(any_flag=lambda x: (x.PC1_flag | x.PC2_flag))

    # get cell_line and treatment 
    x = mut_res[['cell_line', 'treatment', 'run_id']].drop_duplicates()
    batch_res2 = batch_res2.merge(x, on='run_id', how='left')

    # aggregate across sens analysis 
    batch_agg = batch_res2[['batch', 'cell_line', 'treatment', 'any_flag']].groupby(['batch','cell_line', 'treatment']).mean().reset_index()

    # aggregate across sens analyisis 
    df = mut_res.groupby(['cell_line', 'treatment', 'mutant', 'batch'])['prob_res'].agg(['mean', 'std', q05, q95, 'min', 'max']).sort_values('mean', ascending=False).reset_index()

    # merge with batch effect flag 
    df = df.merge(batch_agg, on=['batch', 'cell_line', 'treatment'])

    df.to_csv(f'{args.output}/final_results.csv')
    df[lambda x: (x.cell_line == 'EFM192A') & (x.treatment == 'neratinib')].sort_values('mean', ascending=False).to_csv(f'{args.output}/final_results-sorted-EFM192A;NERATINIB.csv')
    df[lambda x: (x.cell_line == 'EFM192A') & (x.treatment == 'trastuzumab')].sort_values('mean', ascending=False).to_csv(f'{args.output}/final_results.csv-sorted-EFM192A;TRASTUZUMAB.csv')
    df[lambda x: (x.cell_line == 'SKBR3') & (x.treatment == 'neratinib')].sort_values('mean', ascending=False).to_csv(f'{args.output}/final_results.csv-sorted-SKBR3;NERATINIB.csv')
    df[lambda x: (x.cell_line == 'SKBR3') & (x.treatment == 'trastuzumab')].sort_values('mean', ascending=False).to_csv(f'{args.output}/final_results.csv-sorted-SKBR3;TRASTUZUMAB.csv')







