import os
import pandas as pd 
import numpy as np
from tslearn.preprocessing import TimeSeriesResampler
from tslearn.clustering import TimeSeriesKMeans
from matplotlib import pyplot as plt
from sklearn.preprocessing import LabelEncoder
import seaborn as sbn 
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC


def load_data(args): 
    '''
    load data 
    '''
    
    print('\nloading data...')
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

        # ----------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------
        if dataset == 'H210618_Set2': 
            print('fixing minor naming convention issue in dataset "H210618_Set2 - renaming "...L755_T759del..." to "L755T759del"')
            ids = []
            i = 0 
            for id in _data.cell__treatment.values: 
                if 'L755_T759del' in id: 
                    #example: efm192a_erk_akt__L755_T759del_untreated
                    id2 = id[:21] + id[22:]
                    ids.append(id2)
                    i+=1
                else:
                    ids.append(id)
            print('number of `cell__treatment` values modified:', i)
            _data = _data.assign(cell__treatment=ids)
        # ----------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------
        
        _datas.append(_data)

    data = pd.concat(_datas, axis=0)

    clover_sel = [f'{x}_x' for x in series_sel]
    mscarl_sel = [f'{x}_y' for x in series_sel]

    data = data.assign(drug = [x.split('_')[-1].lower() for x in data.cell__treatment])
    data = data.assign(cell_line = [x.split('_')[0].upper() for x in data.cell__treatment])
    data = data.assign(mutant = [x.split('_', maxsplit=5)[-2].upper() for x in data.cell__treatment])

    return data, clover_sel, mscarl_sel


def filter_data(args, data, clover_sel, mscarl_sel): 
    '''
    filter NAs and select 2/3 treatments
    '''

    print('\nfiltering to drug and removing NAs...')

    assert data.drug.unique().shape[0] == 3, f'expected three unique treatments before drug filter, got: {data.drug.unique()}'

    '''old approach
    if args.drug[0].lower() == 'neratinib': 
        drug_ = '10nm_neratinib'
    else: 
        drug_ = '10ug_ml_trastuzumab'
    '''
    drug_ = None
    for option in data.drug.unique(): 
        if (args.drug[0].lower() in option.lower()) or (args.drug[0].lower() == option.lower()): 
            drug_ = option.lower()
    assert drug_ is not None, f'could not assign drug based on user drug choice: {args.drug[0]}'
    print('drug id:', drug_)

    data = data[lambda x: x.drug.isin(['untreated', drug_])]
    print('Data shape (untreated + drug):', data.shape)

    assert data.drug.unique().shape[0] == 2, f'expected two unique treatments, got: {data.drug.unique()}'

    print('length of time series BEFORE removing time points with NA', len(clover_sel))
    clover_sel = np.array(clover_sel)[~data[clover_sel].isna().any()]
    mscarl_sel = np.array(mscarl_sel)[~data[mscarl_sel].isna().any()]
    assert len(clover_sel) == len(mscarl_sel), 'clover timeseries is different length than mscarlet time series'
    print('length of time series AFTER removing time points with NA', len(clover_sel))

    return data

def resample(args, data, clover_sel, mscarl_sel): 
    ########### RESAMPLE ##############
    print('\nresampling time series...')
    
    X_train = np.stack([data[clover_sel], data[mscarl_sel]], axis=2)
    print('Training data shape BEFORE resampling:', X_train.shape)

    # Make time series shorter
    X_train = TimeSeriesResampler(sz=args.resample_sz[0]).fit_transform(X_train)
    print('Training data shape AFTER resampling:', X_train.shape)

    return X_train


def fit_kmeans(args, X_train, save=None): 
    '''
    fit time series k-mean clustering
    '''

    print('\nperforming time-series kmeans clustering...')
    print()
    #                                                      random_state=0, 
    km = TimeSeriesKMeans(n_clusters=args.nclus[0], verbose=True, metric='euclidean', n_jobs=8)
    y_pred = km.fit_predict(X_train)

    print('plotting...')

    F = plt.figure(figsize=(20,10))
    for yi in range(args.nclus[0]):

        if args.nclus[0] % 5 == 0: 
            nrows = int(args.nclus[0] / 5) 
        else: 
            nrows = int(args.nclus[0] / 5)  + 1

        plt.subplot(nrows, 5, yi + 1)

        for xx in X_train[y_pred == yi][0:100]:
            plt.plot(xx[:,0], "r-", alpha=.1)
            plt.plot(xx[:,1], "b-", alpha=.1)

        plt.title(f'cluster sz: {len(X_train[y_pred == yi])}')
        plt.plot(km.cluster_centers_[yi][:,0], "r-", label='clover')

        plt.plot(km.cluster_centers_[yi][:,1], "b-", label='mscarlet')


        plt.xlim(0, args.resample_sz[0])
        plt.ylim(0, 1)
        plt.text(0.55, 0.85,'Cluster %d' % (yi + 1),
                    transform=plt.gca().transAxes)


    plt.tight_layout()

    if save is not None: 
        plt.savefig(save + '/cluster_plots.png')
        plt.close('all')
    else: 
        plt.show()
    return y_pred, km

def calc_cluster_proportions(args, y_pred, data): 
    '''
    '''

    print('\nquantifying experiment by cluster proportions...')
    lb = LabelEncoder()
    y_trt = lb.fit_transform([f'{x}--{y}' for x,y in zip(data.cell__treatment.values, data.dataset.values)])

    cm_cnts = {c:np.zeros(args.nclus[0]) for c in lb.classes_} 

    for i, clus, grp in zip(range(len(y_pred)), y_pred, y_trt) :
        cm_cnts[lb.classes_[grp]][clus] += 1

    cm_prob = {k:v/np.sum(v) for k,v in cm_cnts.items()}

    labels = [k for k,v in cm_prob.items()]
    cm = np.stack([v for k,v in cm_prob.items()], axis=0)

    return cm, lb


def plot_cluster_corr(cm, save=None):
    '''
    '''
    corr = np.corrcoef(cm, rowvar=False)

    f = plt.figure(figsize=(7,7))
    ax = sbn.clustermap(
        corr, 
        vmin=-1, vmax=1, center=0,
    )
    if save is not None: 
        plt.savefig(save + '/cluster_corr_plot.png')
        plt.close('all')
    else: 
        plt.show()


def dimensionality_reduction(args, cm, lb, save=None):
    '''
    '''
    ########### DIM. REDUCTION ##############
    _sens = args.sensitive_line[0].upper()
    _res = args.resistant_line[0].upper()
    _drug = args.drug[0].lower()

    print('\nperforming dim. reduction (pca)...')
    pca = PCA(n_components=2)
    PCs = pca.fit_transform(cm)

    print('PCA explained variance ratio:', pca.explained_variance_ratio_)
    print('PC shape:', PCs.shape)

    res = pd.DataFrame({'pc1': PCs[:,0], 'pc2':PCs[:,1], 'treatment':[x.split('--')[0].split('_')[-1].lower() for x in lb.classes_], 'mutant':[x.split('__')[1].split('_')[0].upper() for x in lb.classes_]})

    plt.figure(figsize=(7,7))
    sbn.scatterplot(x='pc1', y='pc2', data=res, hue='mutant', style='treatment', s=300)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if save is not None: 
        plt.savefig(save + '/PCA_all.png', bbox_inches='tight')
    else: 
        plt.show()

    plt.figure(figsize=(7,7))
    sbn.scatterplot(x='pc1', y='pc2', data=res[lambda x: (x.mutant.isin([_sens, _res]))], hue='mutant', style='treatment', s=300)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if save is not None: 
        plt.savefig(save + '/PCA_labeled.png', bbox_inches='tight')
        plt.close('all')
    else: 
        plt.show()

    return pca, res, _sens, _res, _drug

def train_classifier(res, _sens, _res, _drug, save=None): 
    ''''''
    ########### TRAIN CLASSIFIER ##############
    print('\ntraining classifier...')

    print('sensitive line:', _sens)
    print('resistant line:', _res)
    print('drug:', _drug)

    res_drug = res[lambda x: (x.mutant.isin([_sens, _res])) & (x.treatment == _drug)]
    print('drug + WT df size: ', res_drug.shape)

    X = res_drug[['pc1', 'pc2']].values
    y_res = ((res_drug.mutant == _res).values)
    y_sens = ((res_drug.mutant == _sens).values)

    assert (y_res == ~y_sens).all(), 'y class label assignment has more than 2 classes...'

    y = 1*y_sens

    print('X train shape:', X.shape)
    print('# neg class (resistant):', np.sum(y_res))
    print('# pos class (sensitive):', np.sum(y_sens))
        
    #                                                 , random_state=0
    # rbf or linear
    model = SVC(kernel='linear', C=10, probability=True)
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

    if save is not None: 
        plt.savefig(save + '/classifier_results.png',bbox_inches='tight')
        plt.close('all')
    else: 
        plt.show()

    return model, accuracy


def predict_new(args, res, model): 
    ''''''
    ########### APPLY TO UNLABELED CELL LINES ##############
    print('\npredicting unlabeled sensitivities...')
    _other = res[lambda x: ~(x.mutant.isin([args.sensitive_line[0].lower(), args.resistant_line[0].lower()])) & (x.treatment == args.drug[0].lower())].reset_index(drop=True)

    X_all = _other[['pc1', 'pc2']].values

    y_hat = model.predict_proba(X_all)

    pres = pd.DataFrame({'prob_res':y_hat[:,0], 'prob_sens':y_hat[:,1]})
    prob_res = pd.concat([_other, pres], axis=1) 
    prob_res = prob_res.assign(call=[['res','sens'][np.argmax([x,y])] for x,y in zip(prob_res.prob_res, prob_res.prob_sens)])

    return prob_res 