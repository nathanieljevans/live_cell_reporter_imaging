'''
HER2 project, led by Dr. Samuel Tsang

Summary: Goal is to identify sensitive and resistant mutant lines. 

author: Nathaniel Evans
email: evansna@ohsu.edu

---

Use: 

```{bash} 

$ python HER2_classifier.py --data ./data/HER2/ --drug neratinib --sensitive_line WT --resistant_line T798I --load normalized --nclus 25 --out ./output/ --resample_sz 100
```

`--data` 
directory to data files, should be organized as: 
/data/ 
    /HER2/ 
        /dataset_name/ 
            /normalized
                -> clover_all_cell.csv
                -> mscarlet_all_cell.csv
            /raw
                -> clover_all_cell.csv
                -> mscarlet_all_cell.csv
            
`--drug`
Can be trastuzumab or neratinib

`--sensitive_line` 
The cell line to use as sensitive labels 

`--resistant_line` 
The cell line to use as resistant labels

`--load` 
Whether to use the `normalized` or `raw` data. 

`--nclus`
The number of clusters to use. 

`--out`
Directory path to save results to

`--resample_sz` 
length of time series to resample to 
---
'''

import numpy as np 
import pandas as pd
from matplotlib import pyplot as plt 
import seaborn as sbn

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC
from sklearn.preprocessing import LabelEncoder

from sklearn.decomposition import PCA
from sklearn.svm import SVC

import argparse
from datetime import datetime

import sys 
import uuid
import os

import lib 

'''
version history
1.0 - first functional version 
1.1 - moved most of the code to `lib.py` & minor changes to parsing to accomodate name convention changes in EFM192A dataset. 
'''
__VERSION_ID__ = '1.1'

def get_args(): 
    parser = argparse.ArgumentParser(description='Identify sensitive and resistant mutant lines. ')
    
    parser.add_argument('--data', type=str, nargs=1,
                        help='data directory path')
    
    parser.add_argument('--out', type=str, nargs=1,
                        help='output directory path')
    
    parser.add_argument('--drug', type=str, nargs=1,
                        help='drug to use (can be Neratinib or Trastuzumab)')
    
    parser.add_argument('--sensitive_line', type=str, nargs=1,
                        help='cell line to use for sensitive labels')
    
    parser.add_argument('--resistant_line', type=str, nargs=1,
                        help='cell line to use for resistant labels')
    
    parser.add_argument('--load', type=str, nargs=1,
                        help='whether to use `normalized` or `raw` data.')
    
    parser.add_argument('--nclus', type=int, nargs=1,
                        help='number of clusters to use')
    
    parser.add_argument('--resample_sz', type=int, nargs=1,
                        help='length of time series to resample to')
    
    args = parser.parse_args()
    
    assert args.drug[0].lower() in ['neratinib', 'trastuzumab'], '`--drug` must be either "trastuzumab" or "neratinib"' 
    assert args.load[0].lower() in ['normalized', 'raw'], '`--load` must be either "normalized" or "raw"'
    assert (args.nclus[0] <= 50) & (args.nclus[0] >= 2), '`--nclus` should be an integer between 2 and 50'    
    assert (args.resample_sz[0] <= 150) & (args.resample_sz[0] >= 25), '`--resample_sz` should be an integer between 25 and 150' 
    
    return args


if __name__ == '__main__': 
          
    args = get_args() 
          
    run_id = str(uuid.uuid4())    
          
    output_dir = args.out[0] + '/' + run_id
    
    os.mkdir(output_dir)
    out_log = output_dir + '/console_output.log'
    print('-'*25)
    print(args)
    print('console output will be logged to:', out_log)
    print('-'*25)
    
    with open(out_log, 'w') as sys.stdout: 
        
        print('config...')
        print('-'*25)
        print('script version:', __VERSION_ID__)
        print(args)
        print()
        print(datetime.now())
        print('-'*25)
        
        data, clover_sel, mscarl_sel = lib.load_data(args)

        data = lib.filter_data(args, data, clover_sel, mscarl_sel)

        X_train = lib.resample(args, data, clover_sel, mscarl_sel)

        y_pred, km = lib.fit_kmeans(args, X_train, save=output_dir)
        
        cm, lb = lib.calc_cluster_proportions(args, y_pred, data)
        
        lib.plot_cluster_corr(cm, save=output_dir)
        
        pca, res, _sens, _res, _drug = lib.dimensionality_reduction(args, cm, lb, save=output_dir)
        
        model, accuracy = lib.train_classifier(res, _sens, _res, _drug, save=output_dir)
        
        prob_res = lib.predict_new(args, res, model)

        ########### SAVE RESULTS ##############
        print('\nsaving results...')
        
        prob_res.to_csv(output_dir + '/unlabeled_lines_results.csv') 
        
        run_res = pd.DataFrame({'accuracy(train)': accuracy, 
                                'pc1_var':pca.explained_variance_ratio_[0], 
                                'pc2_var':pca.explained_variance_ratio_[1], 
                                'kmeans_inertia': km.inertia_,
                                'res_line':_res, 
                                'sens_line': _sens,
                                'drug':_drug, 
                                'nclus':args.nclus[0],
                                'resample_sz': args.resample_sz[0],
                                'load': args.load[0], 
                                'run_id':run_id}, index=[0])
        
        run_res.to_csv(output_dir + '/run_results.csv')
                         
                         
                         















