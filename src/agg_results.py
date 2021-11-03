'''
$ python agg_results.py --input ../output/ --out ../output/
'''

import argparse 
import pandas as pd 
import os

def get_args(): 
    parser = argparse.ArgumentParser(description='Aggregate results of live cell reporter analysis.')
    parser.add_argument('--input', help='path to output directory.')
    parser.add_argument('--out', help='output directory to save aggregated results to.')
    return parser.parse_args()

if __name__ == '__main__': 
    args = get_args() 
    path = args.input

    run_res = []
    prob_res = []

    for dir in os.listdir(path): 
        for exp in os.listdir(path + '/' + dir): 
            cell_line, drug = dir.split('_')
            run_res.append(pd.read_csv(path + '/' + dir + '/' + exp + '/run_results.csv').assign(cell_line = cell_line, drug_check = drug))
            prob_res.append(pd.read_csv(path + '/' + dir + '/' + exp + '/unlabeled_lines_results.csv').assign(run_id = exp, cell_line = cell_line, drug_check = drug))

    run_res = pd.concat(run_res, axis=0, ignore_index=True).drop('Unnamed: 0', axis=1)
    prob_res = pd.concat(prob_res, axis=0, ignore_index=True).drop('Unnamed: 0', axis=1)

    run_res.to_csv(args.out + '/experiment_run_results.csv', index=False)
    prob_res.to_csv(args.out + '/mutant_resistance_results.csv', index=False)


