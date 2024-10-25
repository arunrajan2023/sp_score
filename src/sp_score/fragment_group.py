from rdkit import Chem
import ast
from src.sp_score.read_data import *

def fragment_group(frag):
    '''
    Identify functional group associated with each fragment.

    Though the best way is to straight away find the functional group of the fragment, we try to look at the list of available fragments and their known functional goups.

    The reason is, sometimes, otherwise the functional groups would be wrongly reported by rdkit.

    This issue can be resolved in the future. 
    '''
    
    #find the fingerprint of the given fragment
    jnstar             = frag.count("*")
    jfmol              = Chem.MolFromSmiles(frag)
    jnfatm             = jfmol.GetNumAtoms()   #- jnstar    #actual number of atoms
    jncycle            = jfmol.GetRingInfo().NumRings()
    jfp, jfp_additional= Chem.RDKFingerprint(jfmol).ToBitString(), [jnstar, jnfatm, jncycle]     

    fragments_in_class = group_info                     #invokes group_info from read_csv

    #comparison: given frag SMILES and/or its FP, with stored values
    for key, values in fragments_in_class.items():
        fp_found        = False
        # fragments       = values                           #fragment SMILES per class
        fingerprints    = fragments_in_class[key]            #FP array per class

        # print(key, fingerprints.columns)

        frag_smiles     = fingerprints["smiles"]
        fingerprints_fp = fingerprints["fragment_fp"]
        fingerprints_fp2= fingerprints["fragment_fp2"]
        xdata           = len(fingerprints)
        
        for idata in range(xdata):
            _fp         =  fingerprints_fp[idata]  
            _fp2        =  fingerprints_fp2[idata]
            _smi        =  frag_smiles[idata]  

            _fp2 = ast.literal_eval(_fp2)                   #imp: convert to string list
            if _fp == jfp and _fp2 == jfp_additional:  
                # print("match", key, frag, _smi)
                fp_found = True
                break
        if fp_found:
            break

    if fp_found:
        frag_group = key
        #include FP for t-SNE in this part 
    else:
        frag_group = "unknown"
    return frag_group 
