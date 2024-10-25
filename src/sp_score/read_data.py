import pandas as pd
import ast

#read functional group based classified CSV files
#provide only pkl later
acidchlorides          = pd.read_csv('src/data/fingerprint/acidchlorides_and_FPs.csv',          usecols=["smiles", "fragment_fp", "fragment_fp2"])
acids                  = pd.read_csv('src/data/fingerprint/acids_and_FPs.csv',                  usecols=["smiles", "fragment_fp", "fragment_fp2"]) 
alcohols               = pd.read_csv('src/data/fingerprint/alcohols_and_FPs.csv',               usecols=["smiles", "fragment_fp", "fragment_fp2"]) 
aldehydes              = pd.read_csv('src/data/fingerprint/aldehydes_and_FPs.csv',              usecols=["smiles", "fragment_fp", "fragment_fp2"]) 
amides                 = pd.read_csv('src/data/fingerprint/amides_and_FPs.csv',                 usecols=["smiles", "fragment_fp", "fragment_fp2"]) 
amines                 = pd.read_csv('src/data/fingerprint/amines_and_FPs.csv',                 usecols=["smiles", "fragment_fp", "fragment_fp2"])
esters                 = pd.read_csv('src/data/fingerprint/esters_and_FPs.csv',                 usecols=["smiles", "fragment_fp", "fragment_fp2"])
ketones                = pd.read_csv('src/data/fingerprint/ketones_and_FPs.csv',                usecols=["smiles", "fragment_fp", "fragment_fp2"])
heteroatomics          = pd.read_csv('src/data/fingerprint/heteroatomics_and_FPs.csv',          usecols=["smiles", "fragment_fp", "fragment_fp2"])
hydrocarbons           = pd.read_csv('src/data/fingerprint/hydrocarbons_and_FPs.csv',           usecols=["smiles", "fragment_fp", "fragment_fp2"])
other_groups           = pd.read_csv('src/data/fingerprint/other_groups_and_FPs.csv',                usecols=["smiles", "fragment_fp", "fragment_fp2"])
phosphonyl_acids       = pd.read_csv('src/data/fingerprint/phosphonyl_acids_and_FPs.csv',            usecols=["smiles", "fragment_fp", "fragment_fp2"])
phosphonyl_PO_PO2_PO3_PO4= pd.read_csv('src/data/fingerprint/phosphonyl_PO_PO2_PO3_PO4_and_FPs.csv', usecols=["smiles", "fragment_fp", "fragment_fp2"])
PN                     = pd.read_csv('src/data/fingerprint/PN_and_FPs.csv',                          usecols=["smiles", "fragment_fp", "fragment_fp2"])
sulfides               = pd.read_csv('src/data/fingerprint/sulfides_and_FPs.csv',                    usecols=["smiles", "fragment_fp", "fragment_fp2"])
sulfonic_acids         = pd.read_csv('src/data/fingerprint/sulfonic_acid_and_FPs.csv',               usecols=["smiles", "fragment_fp", "fragment_fp2"])
sulfonyl_SO2_or_SO3    = pd.read_csv('src/data/fingerprint/sulfonyl_SO2_or_SO3_and_FPs.csv',         usecols=["smiles", "fragment_fp", "fragment_fp2"])


group_info             = {
    'acidchloride'          :  acidchlorides, 
    'acid'                  :  acids, 
    'alcohol'               :  alcohols, 
    'aldehyde'              :  aldehydes, 
    'amide'                 :  amides, 
    'amine'                 :  amines,
    'ester'                 :  esters, 
    'ketone'                :  ketones, 
    'heteroatomic'          :  heteroatomics, 
    'other_groups'          :  other_groups, 
    'phosphonyl_acids'      :  phosphonyl_acids,
    'phosphonyl_PO_PO2_PO3_PO4' : phosphonyl_PO_PO2_PO3_PO4,
    'PN'                    :  PN,
    'sulfonic_acid'         :  sulfonic_acids, 
    'sulfide'               :  sulfides,
    'sulfonyl_SO2_or_SO3'   :  sulfonyl_SO2_or_SO3,
    'hydrocarbon'           :  hydrocarbons}  


#read fragments and pair  statistics
df_ufrag       =  pd.read_csv("src/data/statistics/unique_fragments2.csv")
df_upair       =  pd.read_csv("src/data/statistics/unique_pairs2.csv")

ufrag_smiles   = df_ufrag["fragment_smiles"]
ufrag_freq     = df_ufrag["fragment_frequency"]  
ufrag_fp       = df_ufrag["fragment_fp"]         
ufrag_fp2      = df_ufrag["fragment_fp2"]        

upair_smiles   = df_upair["fragment_pair_smiles"]
upair_ismiles  = df_upair["fragment_ismiles"]
upair_jsmiles  = df_upair["fragment_jsmiles"]                    
upair_freq     = df_upair["fragment_pair_frequency"]
upair_fp       = df_upair["fragment_pair_fp"]
upair_fp2      = df_upair["fragment_pair_fp2"]    

stat_info ={
    'ufrag_smiles' : ufrag_smiles,
    'ufrag_freq'   : ufrag_freq,
    'ufrag_fp'     : ufrag_fp,
    'ufrag_fp2'    : ufrag_fp2,
    'upair_smiles' : upair_smiles,
    'upair_ismiles': upair_ismiles,
    'upair_jsmiles': upair_jsmiles,
    'upair_freq'   : upair_freq,
    'upair_fp'     : upair_fp,
    'upair_fp2'    : upair_fp2
}


master_graphs = pd.read_csv("src/data/graph/known_graphs.csv")

master_graphs = master_graphs["graph"].values
