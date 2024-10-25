import re
import os
import sys
import ast
import time
import numpy as np
import pandas as pd
import networkx as nx

from rdkit import Chem
from rdkit.Chem import Fragments, MolStandardize
from rdkit.Chem.BRICS import BRICSDecompose, BRICSBuild, FindBRICSBonds, BreakBRICSBonds

from src.sp_score.fragment_group import fragment_group

unknowns_array = []
def generate_graph(smiles, verbose=False):
        """
        1. Frags
        2. Pairs
        3. Graph (two types)              NOTE: Turned off the unknown frags
        """
        if "[e]" in smiles:
            smiles = smiles.replace('[e]', '*')
            smiles = smiles.replace('[d]', '*')
            smiles = smiles.replace('[t]', '*')
            smiles = smiles.replace('[g]', '*')
            csmiles= Chem.CanonSmiles(smiles) 
        else:
            csmiles= Chem.CanonSmiles(smiles) 

        graph             = nx.Graph() 
        graph2            = nx.Graph()

        # print(csmiles)
        mol        = Chem.MolFromSmiles(csmiles) 

        if not mol:
            print("-----------------------------not rdkit mol--------------------------")
            return None

        #BRICS fragmentation of the SMILES
        # start = time.time()
        natm       = mol.GetNumAtoms()
        bonds      = list(FindBRICSBonds(mol))
        nbond      = len(list(bonds))
        npair      = nbond 

        mol_bbb    = BreakBRICSBonds(mol)
        frags_mol  = Chem.GetMolFrags(mol_bbb, asMols=True)  #frags mol
        frags_idx  = Chem.GetMolFrags(mol_bbb, asMols=False) #frags atoms (indices)

        frags      = [Chem.MolToSmiles(x, False) for x in frags_mol] 
        nfrag      = len(frags) 
        sys.stdout.flush() 
        # time.sleep(1)
        # end   = time.time()
        # print("nfrag= {}, BRICS: {} seconds".format(nfrag, time.time() - start))        

        if verbose:       
            print(frags) 

        #find SMILES pairs
        pairs      = []

        def count_bondtype(mol):  
            #count bond types       
            ubond_rdkit = []
            ubond_element=[]
            for ijbond, ijlabel in FindBRICSBonds(mol):
            #     print(ijlabel)
                if ijlabel not in ubond_rdkit:
                    ubond_rdkit.append(ijlabel)
                    iatm, jatm= ijbond[0], ijbond[1]
                    ijelement = mol.GetAtomWithIdx(iatm).GetSymbol() + mol.GetAtomWithIdx(jatm).GetSymbol()
                    ubond_element.append(ijelement)

            ubond_freq = {}
            for upos, ubond in enumerate(ubond_rdkit):
                ubond_freq[upos] = 0 
                for ijbond, ijlabel in FindBRICSBonds(mol):
                        if ijlabel == ubond:
                                ubond_freq[upos] += 1

            #     print("Unique RDKit bonds  :", ubond_rdkit)
            #     print("Unique bond elements:", ubond_element)
            #     if len(ubond_element) != len(np.unique(ubond_element)):
            #         print("Warning: Same elements of different chemical classes found")

            return ubond_freq   
        ubond_freq =  count_bondtype(mol)

        #frag atoms,  nneighbors
        frags_element_index = []
        for ictr, frag in enumerate(frags_idx):
            frag_atoms  =  [x for x in frag if x in range(natm)]          #atoms
            frag_dummy  =  [x for x in frag if x not in range(natm)]      #dummy atoms     
            frags_element_index.append(frag_atoms)      

        for ictr, ifrag in enumerate(frags):
            ielements = frags_element_index[ictr]
            ineb      = []
            out_ineb  = []
            for iatm in ielements:
                neb = [x.GetIdx() for x in mol.GetAtoms()[iatm].GetNeighbors()]   
                ineb.append(iatm)
                for _ in neb:
                    if _ in ielements:
                            ineb.append(_)
                    else:
                            out_ineb.append(_)
            for jctr, jfrag in enumerate(frags):
                if ictr >= jctr: continue
                jelements = frags_element_index[jctr]
                jneb      = []
                out_jneb  = []
                for jatm in jelements:
                    neb = [x.GetIdx() for x in mol.GetAtoms()[jatm].GetNeighbors()]   
                    jneb.append(jatm)
                    for _ in neb:
                        if _ in jelements:
                            jneb.append(_)
                        else:
                            out_jneb.append(_)
                if list(set(ineb) & set(out_jneb)) and list(set(jneb) & set(out_ineb)):
                    pair = [ifrag, jfrag]
                    if verbose:
                        print("pair", *pair, sep='\t')
                    pairs.append(pair)
                    iclass = fragment_group(ifrag)                                       #calling fragment_group 
                    jclass = fragment_group(jfrag)
                    # print("classes:", iclass, jclass)
                    if iclass == "unknown":
                        # print(iclass, "frag smiles:",ifrag, "smiles:", smiles)
                        unknowns_array.append(ifrag)
                        break
                    if jclass == "unknown":
                        # print(jclass, "frag smiles:", jfrag, "smiles:", smiles)
                        unknowns_array.append(jfrag)
                        break

                    # print(ifrag, jfrag)
                    # if ifrag.count("*") != 2: ifrag.replace('*', '[Bi]')
                    # if jfrag.count("*") != 2: jfrag.replace('*', '[Bi]')
                    # fp1 = PolymerSmiles(ifrag).fingerprint
                    # fp2 = PolymerSmiles(jfrag).fingerprint

                    # accuracy=5
                    # fp1 = {key : round(fp1[key], accuracy) for key in fp1} 
                    # fp2 = {key : round(fp2[key], accuracy) for key in fp2}
                    # print(fp1, fp2)
                    graph.add_edge(ictr, jctr, iclass=iclass, jclass=jclass)             #adding Graph's edge and functional groups 
                    graph2.add_edge(ifrag, jfrag, iclass=iclass, jclass=jclass)  

        # print("Graph:",graph.edges(data=True))    #graph as list 
        # if unknowns_array:
        #     frags, pairs, graph, graph2, graph3 = [], [], [], [], []
        # else:
        #     if npair != len(pairs):
        #         print("Warning! Unequal pairs.")
        #         npair = len(pairs)
        #         print("npair", npair)

        #     if verbose:
        #         print("pairs:", pairs) 

        graph3 = graph
        graph = graph.edges(data=True)
        graph2= graph2.edges(data=True)
        return frags, pairs, graph, graph2, graph3
