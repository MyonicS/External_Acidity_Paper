
import pandas as pd
import numpy as np
path = 'Catalyst_database.csv'

def write_db(catalyst,property,value,db_path=path):
    df = pd.read_csv(db_path, sep=',',header=0,index_col=0)
    df.loc[catalyst,property] = value
    df.to_csv(db_path, sep=',',header=True,index=True)

def read_db(catalyst,property,db_path=path):
    df = pd.read_csv(db_path, sep=',',header=0,index_col=0)
    return df.loc[catalyst,property]

def read_unit(property,db_path=path):
    df = pd.read_csv(db_path, sep=',',header=0,index_col=0)
    return df.loc['Unit',property]

def get_db(db_path=path):
    return pd.read_csv(db_path, sep=',',header=0,index_col=0)
    
def get_catname(catcode, db_path=path):
    if catcode == 'Z13':
        return r'ZY$_{{{}}}$'.format(''.join(str(int(np.round(read_db(catcode,'SAR_ICP'))))+'\mathrm{st}'))
    elif catcode == 'Z15':
        return 'FCC'
    elif catcode == 'Z16':
        return 'FCCcr'
    elif catcode == 'Z14':
        return r'ZY$_{{{}}}$'.format(''.join('5')+'\mathrm{B}')
    else:
        return r'ZY$_{{{}}}$'.format(''.join(str(int(np.round(read_db(catcode,'SAR_ICP'))))))