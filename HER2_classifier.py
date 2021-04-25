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

from tslearn.clustering import TimeSeriesKMeans, KernelKMeans
from tslearn.metrics import dtw

from tslearn.preprocessing import TimeSeriesScalerMeanVariance, TimeSeriesResampler

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import LabelEncoder

from sklearn.decomposition import PCA
from sklearn.svm import SVC

import argparse
from datetime import datetime

import sys 
import uuid
import os

__VERSION_ID__ = '1.0'

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
        
        ########### LOAD DATA ##############
        
        print('\nsloading data...')
        load = args.load[0].lower()
        data_dir = args.data[0]
        datasets = [x for x in os.listdir(data_dir) if os.path.isdir(data_dir + x)]
        
        print('# of datasets to load:',len(datasets))
        
        series_sel = pd.read_csv(data_dir + datasets[0] + '/' + load + '/clover_all_cell.csv').columns[1:-3]

        _datas = []
        for dataset in datasets: 
            cl_path = data_dir + dataset + '/' + load + '/clover_all_cell.csv'
            ms_path = data_dir + dataset + '/' + load + '/mscarlet_all_cell.csv'
            _clover = pd.read_csv(cl_path)
            _mscarl = pd.read_csv(ms_path)
            _data = _clover.merge(_mscarl, on=['track_index', 'cell__treatment'], how='inner')
            _data = _data.assign(dataset=dataset)
            _datas.append(_data)

        data = pd.concat(_datas, axis=0)

        clover_sel = [f'{x}_x' for x in series_sel]
        mscarl_sel = [f'{x}_y' for x in series_sel]

        data = data.assign(drug = [x.split('_', maxsplit=5)[-1] for x in data.cell__treatment])
        data = data.assign(cell_line = [x.split('_', maxsplit=5)[0] for x in data.cell__treatment])
        data = data.assign(mutant = [x.split('_', maxsplit=5)[-2] for x in data.cell__treatment])

        
        ########### FILTER TO DRUG & REMOVE NA ##############
        print('\nfiltering to drug and removing NAs...')
        
        if args.drug[0].lower() == 'neratinib': 
            drug_ = '10nm_neratinib'
        else: 
            drug_ = '10ug_ml_trastuzumab'
              
        data = data[lambda x: x.drug.isin(['untreated', drug_])]
        print('Data shape (untreated + drug):', data.shape)
        
        print('length of time series BEFORE removing time points with NA', len(clover_sel))
        clover_sel = np.array(clover_sel)[~data[clover_sel].isna().any()]
        mscarl_sel = np.array(mscarl_sel)[~data[mscarl_sel].isna().any()]
        assert len(clover_sel) == len(mscarl_sel), 'clover timeseries is different length than mscarlet time series'
        print('length of time series AFTER removing time points with NA', len(clover_sel))
        
        ########### RESAMPLE ##############
        print('\nresampling time series...')
        
        X_train = np.stack([data[clover_sel], data[mscarl_sel]], axis=2)
        print('Training data shape BEFORE resampling:', X_train.shape)

        # Make time series shorter
        X_train = TimeSeriesResampler(sz=args.resample_sz[0]).fit_transform(X_train)
        print('Training data shape AFTER resampling:', X_train.shape)
        

        ########### FIT K-MEANS ##############
        print('\nperforming time-series kmeans clustering...')
        print()
        km = TimeSeriesKMeans(n_clusters=args.nclus[0], verbose=True, random_state=0, metric='euclidean', n_jobs=8)
        y_pred = km.fit_predict(X_train)
        print()
        
        F = plt.figure(figsize=(20,10))
        for yi in range(args.nclus[0]):
            if args.nclus[0] % 5 == 0: 
                nrows = int(args.nclus[0] / 5) 
            else: 
                nrows = int(args.nclus[0] / 5)  + 1
                
            plt.subplot(nrows, 5, yi + 1)
            for xx in X_train[y_pred == yi][0:250]:
                plt.plot(xx[:,0], "r-", alpha=.05)
                plt.plot(xx[:,1], "b-", alpha=.05)

            plt.title(f'cluster sz: {len(X_train[y_pred == yi])}')
            plt.plot(km.cluster_centers_[yi][:,0], "r-", label='clover')
            plt.plot(km.cluster_centers_[yi][:,1], "b-", label='mscarlet')

            plt.xlim(0, args.resample_sz[0])
            plt.ylim(0, 1)
            plt.text(0.55, 0.85,'Cluster %d' % (yi + 1),
                     transform=plt.gca().transAxes)

        plt.tight_layout()
        plt.savefig(output_dir + '/cluster_plots.png')
        plt.close('all')
        ########### CALCULATE CLUST PROPORTIONS ##############
        print('\nquantifying experiment by cluster proportions...')
        lb = LabelEncoder()
        y_trt = lb.fit_transform([f'{x}--{y}' for x,y in zip(data.cell__treatment.values, data.dataset.values)])

        cm_cnts = {c:np.zeros(args.nclus[0]) for c in lb.classes_} 

        for i, clus, grp in zip(range(len(y_pred)), y_pred, y_trt) :
            cm_cnts[lb.classes_[grp]][clus] += 1

        cm_prob = {k:v/np.sum(v) for k,v in cm_cnts.items()}

        labels = [k for k,v in cm_prob.items()]
        cm = np.stack([v for k,v in cm_prob.items()], axis=0)
        
        ########### PLOT CLUSTER CORR ##############
        
        corr = np.corrcoef(cm, rowvar=False)

        f = plt.figure(figsize=(7,7))
        ax = sbn.clustermap(
            corr, 
            vmin=-1, vmax=1, center=0,
        )
        plt.savefig(output_dir + '/cluster_corr_plot.png')
        plt.close('all')
        
        ########### DIM. REDUCTION ##############
        _sens = args.sensitive_line[0].upper()
        _res = args.resistant_line[0].upper()
        _drug = args.drug[0].lower()
        
        print('\nperforming dim. reduction (pca)...')
        pca = PCA(n_components=2)
        PCs = pca.fit_transform(cm)

        print('PCA explained variance ratio:', pca.explained_variance_ratio_)
        print('PC shape:', PCs.shape)
        
        res = pd.DataFrame({'pc1': PCs[:,0], 'pc2':PCs[:,1], 'treatment':[x.split('--')[0].split('_')[-1] for x in lb.classes_], 'cell_line':[x.split('_')[4] for x in lb.classes_]})
        
        plt.figure(figsize=(7,7))
        sbn.scatterplot(x='pc1', y='pc2', data=res, hue='cell_line', style='treatment', s=300)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(output_dir + '/PCA_all.png', bbox_inches='tight')
        
        plt.figure(figsize=(7,7))
        sbn.scatterplot(x='pc1', y='pc2', data=res[lambda x: (x.cell_line.isin([_sens, _res]))], hue='cell_line', style='treatment', s=300)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(output_dir + '/PCA_labeled.png',bbox_inches='tight')
        plt.close('all')
        
        ########### TRAIN CLASSIFIER ##############
        print('\ntraining classifier...')
        
        print('sensitive line:', _sens)
        print('resistant line:', _res)
        print('drug:', _drug)
        
        res_drug = res[lambda x: (x.cell_line.isin([_sens, _res])) & (x.treatment == _drug)]
        print('drug + WT df size: ', res_drug.shape)

        X = res_drug[['pc1', 'pc2']].values
        y_res = ((res_drug.cell_line == _res).values)
        y_sens = ((res_drug.cell_line == _sens).values)

        assert (y_res == ~y_sens).all(), 'y class label assignment has more than 2 classes...'
        
        y = 1*y_sens
        
        print('X train shape:', X.shape)
        print('# neg class (resistant):', np.sum(y_res))
        print('# pos class (sensitive):', np.sum(y_sens))
 
        model = SVC(kernel='linear', C=10, probability=True, random_state=0)
        model.fit(X,y) 
        y_pred = model.predict(X)
        accuracy = accuracy_score(y, y_pred)
        
        plt.figure(figsize=(10, 5))
        plt.subplots_adjust(bottom=.2, top=.95)
        
        #xx = np.linspace(-1, 1, 100)
        #yy = np.linspace(-1, 1, 100).T
        xx = np.linspace(min(X[:,0]), max(X[:,0]), 100)
        yy = np.linspace(min(X[:,1]), max(X[:,1]), 100).T
        xx, yy = np.meshgrid(xx, yy)
        Xfull = np.c_[xx.ravel(), yy.ravel()]
        
        # View probabilities:
        probas = model.predict_proba(Xfull)
        n_classes = np.unique(y_pred).size
        class_names = ['resistant', 'sensitive']
        name = 'Support Vector Classifier'
        for k in range(n_classes):
            plt.subplot(1, n_classes, 0 * n_classes + k + 1)
            plt.title("%s class" % class_names[k])
            if k == 0:
                plt.ylabel(name)
            imshow_handle = plt.imshow(probas[:, k].reshape((100, 100)),
                                       extent=(-1, 1, -1, 1), origin='lower')
            plt.xticks(())
            plt.yticks(())
            idx = (y_pred == k)
            if idx.any():
                plt.scatter(X[idx, 0], X[idx, 1], marker='o', c='w', edgecolor='k')
        
        ax = plt.axes([0.15, 0.04, 0.7, 0.05])
        plt.title("Probability")
        plt.colorbar(imshow_handle, cax=ax, orientation='horizontal')

        plt.savefig(output_dir + '/classifier_results.png',bbox_inches='tight')
        plt.close('all')
        
        ########### APPLY TO UNLABELED CELL LINES ##############
        print('\npredicting unlabeled sensitivities...')
        _other = res[lambda x: ~(x.cell_line.isin([args.sensitive_line[0].lower(), args.resistant_line[0].lower()])) & (x.treatment == args.drug[0].lower())].reset_index(drop=True)

        X_all = _other[['pc1', 'pc2']].values

        y_hat = model.predict_proba(X_all)

        pres = pd.DataFrame({'prob_res':y_hat[:,0], 'prob_sens':y_hat[:,1]})
        prob_res = pd.concat([_other, pres], axis=1) 
        prob_res = prob_res.assign(call=[['res','sens'][np.argmax([x,y])] for x,y in zip(prob_res.prob_res, prob_res.prob_sens)])
        
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
                         
                         
                         















