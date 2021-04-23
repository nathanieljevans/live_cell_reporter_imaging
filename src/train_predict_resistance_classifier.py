import argparse
import os 
import pandas as pd

def get_args(): 

    parser = argparse.ArgumentParser(description='train classifier and predict sensitivities')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')
    args = parser.parse_args()


def load_data(load='normalized'): 

    assert load in ['normalized', 'raw'], f'`load` should be `normalized` or `raw`, you passed: {load}'

    series_sel = pd.read_csv('../data/HER2/H210122_SKBR3/normalized/clover_all_cell.csv').columns[1:-3]

    _datas = []
    for dataset in os.listdir('../data/HER2'): 
        cl_path = '../data/HER2/' + dataset + '/' + load + '/clover_all_cell.csv'
        ms_path = '../data/HER2/' + dataset + '/' + load + '/mscarlet_all_cell.csv'
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

    return data


if __name__ == '__main__': 

    args = get_args() 

    if args.load_norm: 
        load = 'normalized'
    else: 
        load = 'raw'

    data = load_data(load=load) 

    if args.treatment == 'neratinib': 
        drug = '10nm_neratinib'
    elif args.treatment == 'trastuzumab': 
        drug = '10uL_ml_trastuzumab'

    data = data[lambda x: x.drug.isin(['untreated', drug])]