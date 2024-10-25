from rdkit import Chem
from pathos.multiprocessing import ProcessingPool as Pool
import time
import pandas as pd
import numpy as np
from multiprocessing import  Pool
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from src.sp_score.sp_score import PScore

# org_df = pd.read_csv("data/ps_known.csv")
fout   = open("graphs.dat", "w+")

#G2G data
# org_df = pd.read_csv('data/g2g_800k_smiles.txt', sep=" ", header=None)
# org_df.columns = ["smiles"]

# start = time.time()

# graphs = []
# index  = 0
# for smi in org_df[100:150]["smiles"].values:
#     index += 1
#     print()
#     print(index, smi)
#     x = PScore(smiles_string=smi, known=True)
#     p = x[0]
#     g = x[1]
#     g2= x[2]
#     if g:
#         graphs.append(g)

#     print(x[0])
#     print(x[1])
#     print(x[2])
#     # print(x[3])

#     fout.write("\n{}".format(x[0]))
#     fout.write("\n{}".format(x[1]))
#     fout.write("\n{}".format(x[2]))
#     # fout.write("\n{}".format(x[3])) 
#     fout.write("\n")   

# df           = pd.DataFrame()
# df["graphs"] = list(set(graphs))
# df.to_csv("data/graph/known_graphs.csv")      #save to data/graph/
# print("#graphs", len(df))
# print("Done")

#Parallelize this operation




start = time.time()

# org_df = pd.read_csv('g2g_800k_smiles.txt', sep=" ", header=None)
# org_df.columns = ["smiles"]
# known = False


known = True 
org_df = pd.read_csv("data/ps_known.csv")

def add_features(df):
    array_smiles = []
    array_p = []
    array_g = []
    array_g2 = []
    array_delta = []
    for ind in df.index:
        smiles = df['smiles'][ind] 

        if "[e]" in smiles:
            smiles = smiles.replace('[e]', '*')
            smiles = smiles.replace('[d]', '*')
            smiles = smiles.replace('[t]', '*')
            smiles = smiles.replace('[g]', '*')
            smiles= Chem.CanonSmiles(smiles) 
        else:
            smiles= Chem.CanonSmiles(smiles)  

        res = PScore(smiles_string=smiles, known=known)
        p_s, g, g2, delta = res[0], res[1], res[2], res[3]

        array_p.append(p_s)
        array_smiles.append(smiles)
        array_g.append(g)
        array_g2.append(g2)
        array_delta.append(delta)
    df['smiles'] = list(array_smiles)
    df['p_score'] = list(array_p)
    df['graph'] = list(array_g)
    df['graph2'] = list(array_g2)
    df["delta"] = array_delta
    #writing should be done here to avoid any later data-loss
    #df.to_csv("name"+number+".csv")    #file names should be different to avoid any issues
    fname = time.strftime("%Y%m%d_%H%M%S.csv")
    # df.to_csv("result/"+fname, mode="a")
    df.to_csv("data/graph/"+fname, mode="a")

    return df

def parallelize_dataframe(df, func, n_cores=30):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

train = parallelize_dataframe(org_df, add_features)


print("type:", type(train))
print(train.head())

# train.to_csv("result/g2g_800k_pscore.csv")
train.to_csv("data/graph/known_graphs.csv")      #save to data/graph/


end   = time.time()
print("{} seconds".format({time.time() - start}))
