from src.sp_score.__init__ import __version__
import src.utils.sascore as sascore
from rdkit import Chem
import pandas as pd
import statistics
from multiprocessing.pool import ThreadPool as Pool
import time
from multiprocessing import  Pool
import numpy as np


def ncores():
    ncore = 1 #28
    return ncore

def test_version():
    assert __version__ == '0.1.0'

if __name__ == "__main__":

    df= pd.read_pickle("SFHP_pscores_forward_synthesis.pkl")
    #>>> df.columns
    #Index(['SFHP_ID', 'pathway', 'pscore', 'class', 'smiles', 'pscore',
    #       'reactant1', 'zinc_id1', 'reactant2', 'zinc_id2'],
    #      dtype='object')
    #>>>

    ctr=0
    

    def calc_hsum(df):
        """The split df to calculate harmonic mean of SAscores for monomers"""
        for ind in df.index:
            _id = df['SFHP_ID'][ind]
            _pw = df['pathway'][ind]
            #_ps = df['pscore'][ind]   #error
            _cls= df['class'][ind]
            _sm = df['smiles'][ind]
            _r1 = df['reactant1'][ind]
            _id1= df['zinc_id1'][ind]
            _r2 = df['reactant2'][ind]
            _id2= df['zinc_id2'][ind]
        
            mp=Chem.MolFromSmiles(_sm)
            mr1=Chem.MolFromSmiles(_r1)
            mr2=Chem.MolFromSmiles(_r2)
            print("Test {},{},{},{}".format(_id, _sm, _id1, _id2))
        
            if mp and mr1 and mr2:
                #works
                #[sascore.sa_score(smiles=dfk['reactant1'][ind]) for ind in dfk.index]

                [sascore.sa_score(smiles=dfk['reactant1'][ind]) for ind in dfk.index]

                #not working...
                _sar1=sascore.sa_score(smiles=Chem.MolToSmiles(mr1))
                _sar2=sascore.sa_score(smiles=Chem.MolToSmiles(mr2))
        
                _sa = statistics.harmonic_mean([_sar1, _sar2])
                print("{},{},{},{}".format(_id, _sm, _id1, _id2))
            else:
                ctr+=1
        return None

    #ncores
    n = ncores()
    df_split = np.array_split(df, n)
    with Pool(processes=n) as pool:
        results = pool.map(calc_hsum, df_split)
    pool.close()
    pool.join()
    
    end   = time.time()
    print("Excluded", ctr)    
    print("{} seconds".format({time.time() - start}))
    print("Done")            
    
       # _sa=round((10.0 - sa)/9.0, 4)
    
    #Given two smiles calulate harmonic_mean of SAscores for monomer
    
    
