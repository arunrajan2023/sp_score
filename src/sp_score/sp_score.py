import re
import os
import sys
import ast
import networkx as nx
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Fragments, MolStandardize
from rdkit.Chem.BRICS import BRICSDecompose, BRICSBuild, FindBRICSBonds, BreakBRICSBonds

from src.sp_score.read_data import *  
from src.sp_score.generate_graph import generate_graph as generate_graph
import src.utils.sascore as sascore

class PScore:
    """
    Polymer synthetic-accessibility score (sp-score) 

    Remaining: master_graphs
    """
    def __init__(self, **kwargs):                        
        return None

    def __new__(self, **kwargs):
        """
        Read the input SMILES from the user
        """

        def delta_score():
            """
            Check if graph is isomorphic wrt the known polymer graphs
            """

            def regenerate_graph(g):
                """
                Regenerate graph using the edges and attributes read from file
                """
                _graph = nx.Graph()

                if g == []: 
                    pass     
                else:
                    # g = ast.literal_eval(g)
                    for _g in g:
                        ictr, jctr, attr = _g
                        ikey, jkey = list(attr)
                        iclass = attr[ikey]     
                        jclass = attr[jkey]     
                        _graph.add_edge(ictr, jctr, iclass=iclass, jclass=jclass)
                return _graph  
            
            #default value   
            delta = 0                            

            #known graph versus novel (isomorphic or not) graph
            if known:
                delta = 1 

            else:
                
                # source = regenerate_graph(list(graph))
                source = graph3
                for igraph in master_graphs:
                    igraph = ast.literal_eval(igraph)
                    # if nx.is_empty(igraph): continue
                    igraph = regenerate_graph(igraph)
                    output = nx.is_isomorphic(igraph, source, edge_match=lambda x, y: x==y)
                    if output == True:
                        delta =1 #isomorphic to known graph
                        break
            return delta

        def fragment_score():
            """
            Fragments' contribution and score per polymer
            """

            #use the statistics of fragments and fragment-pairs
            ufrag_smiles = stat_info['ufrag_smiles']
            ufrag_freq = stat_info['ufrag_freq']
            ufrag_fp = stat_info['ufrag_fp']
            ufrag_fp2 = stat_info['ufrag_fp2']

            #maximum value of fragment frequency
            fi_max      = max(list(ufrag_freq))

            if verbose:
                print("fi_max,  fip_max   :", fi_max,  fip_max)

            ifrequency_sum      =  0

            #fragment FPs from the master dataset
            fp_array            =  ufrag_fp 
            fp_additional_array =  ufrag_fp2 

            fragment_frequencies=  []

            nfrag = len(frags)           
            if verbose:        
                    print(nfrag, "frags")
                    print(frags)

            for j, fsmiles in enumerate(frags):
                    if verbose:
                            print(fsmiles)
                    #jfrag FP
                    jnstar       = fsmiles.count("*")
                    jfmol        = Chem.MolFromSmiles(fsmiles)
                    jnfatm       = jfmol.GetNumAtoms()
                    jncycle      = jfmol.GetRingInfo().NumRings()
                    jfp, jfp_additional= Chem.RDKFingerprint(jfmol).ToBitString(), [jnstar, jnfatm, jncycle]            

                    for index, fp in enumerate(fp_array):
                            fp_additional = fp_additional_array[index]
                            fp_additional = ast.literal_eval(fp_additional) #imp: convert to list

                            if fp == jfp and fp_additional == jfp_additional:
                                    if verbose:
                                        print("frag index", index)
                                        print("frag freq ", ufrag_freq[index])
                                    fi                = ufrag_freq[index] 
                                    ifrequency_sum   += np.log10(fi) 
                                    fragment_frequencies.append(fi)
            frag_isolated_contribution  = round(ifrequency_sum  / (nfrag * np.log10(fi_max)),4)             
            return frag_isolated_contribution 
                                                        
        def fragment_pair_score():
            """
            Calculate fragment-pair contribution per polymer
            """

            upair_smiles = stat_info['upair_smiles']
            upair_ismiles = stat_info['upair_ismiles']
            upair_jsmiles = stat_info['upair_jsmiles']
            upair_freq = stat_info['upair_freq']
            upair_fp = stat_info['upair_fp']
            upair_fp2 = stat_info['upair_fp2'] 

            #maximum value of pair frequency
            fip_max     = max(list(upair_freq))                 

            ijfrequency_sum          =  0
            
            pair_smiles              =  upair_smiles           #fragment-pair smiles       
            pair_ismiles             =  upair_ismiles 
            pair_jsmiles             =  upair_jsmiles          # check again?          
            
            fp_pair_array            =  upair_fp 
            fp_pair_additional_array =  upair_fp2                       
            fragment_pair_frequencies=  []

            npair = len(pairs)           
            
            if npair == 0:
                frag_pair_contribution = 0.0
                if verbose:
                    print("frag. pair contribution:", frag_pair_contribution)                    
            else: 
                for j in range(npair):
                    pair_ismiles, pair_jsmiles = pairs[j]    #ith and jth frags in a pair 

                    #first FP
                    j1nstar              = pair_ismiles.count("*") 
                    j1fmol               = Chem.MolFromSmiles(pair_ismiles)
                    j1nfatm              = j1fmol.GetNumAtoms()
                    j1ncycle             = j1fmol.GetRingInfo().NumRings()
                    j1fp, j1fp_additional= Chem.RDKFingerprint(j1fmol).ToBitString(), [j1nstar, j1nfatm, j1ncycle] 

                    #second FP 
                    j2nstar              = pair_jsmiles.count("*") 
                    j2fmol               = Chem.MolFromSmiles(pair_jsmiles)
                    j2nfatm              = j2fmol.GetNumAtoms()
                    j2ncycle             = j2fmol.GetRingInfo().NumRings()
                    j2fp, j2fp_additional= Chem.RDKFingerprint(j2fmol).ToBitString(), [j2nstar, j2nfatm, j2ncycle] 

                    #pair (i,j) FP
                    j12fp                = [j1fp,            j2fp]
                    j12fp_additional     = [j1fp_additional, j2fp_additional]   

                    #pair (j,i) FP
                    j21fp                = [j2fp,            j1fp] 
                    j21fp_additional     = [j2fp_additional, j1fp_additional]   

                    condition12, condition21 = False, False 
                    for index, fp in enumerate(fp_pair_array):
                        fp_pair_additional = fp_pair_additional_array[index]

                        fp = ast.literal_eval(fp)                                   #imp: convert to list
                        fp_pair_additional = ast.literal_eval(fp_pair_additional)   #imp: convert to list

                        condition12= (fp == j12fp) and (fp_pair_additional == j12fp_additional)
                        condition21= (fp == j21fp) and (fp_pair_additional == j21fp_additional)

                        if True in [condition12, condition21]:                         
                            if verbose:
                                print("pair index",  index)
                                print("pair freq ",  upair_freq[index])
                            fip                = upair_freq[index] 
                            ijfrequency_sum   += np.log10(fip) 
                            fragment_pair_frequencies.append(fip)  
                            if verbose:
                                print("fip", fip)  

                frag_pair_contribution  = round((ijfrequency_sum) / (npair * np.log10(fip_max)),4)    
            return frag_pair_contribution

        def mod_sascore():
            """
            Calculate the modified SA score of monomer molecule in a range 0-1
            """
            mol   = Chem.MolFromSmiles(self.smiles)

            if not mol:
                print(self.smiles, "SMILES error")
                mod_sa_score = None

            else:
                try:
                    new_mol    = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts('[#0]'))   #chop-off the *s
                except:
                    new_mol    = mol  

                #synthetic accessibility score 
                _sa              = sascore.sa_score(smiles=Chem.MolToSmiles(new_mol)) 

                #convert to 0-1 range (0: impossible, 1:ideal)
                mod_sa_score     = round((10.0 - _sa)/9.0, 4) 
            return mod_sa_score

        self.verbose = False
        self.smiles = False
        self.known = False

        for key, value in kwargs.items(): 
            if key == "smiles_string":
                self.smiles  = value  
            elif key == "verbose":
                self.verbose = value  
            elif key == "known":
                self.known = value
            else:
                    print("*"*30)
                    print("Key error", key)

        verbose =  self.verbose
        smiles =  self.smiles 
        known = self.known       

        #generate BRICS fragments, pairs and the graph
        #here, graph3 is in nx.graph format
        frags, pairs, graph, graph2, graph3 =  generate_graph(smiles)

        if not graph:
            delta, sp_score = 0, 'None'
        else:
            #scores
            delta = delta_score() 
            if delta == 1:
                frag_isolated_contribution = fragment_score()
                frag_pair_contribution = fragment_pair_score()

                try:
                    mod_sa_score = mod_sascore() 
                except:
                    sp_score = 'None'
                    graph = []
                    graph2 = [] 
                    delta = 0
                    return sp_score, graph, graph2, delta  

                # print("scores:",frag_isolated_contribution, frag_pair_contribution, mod_sa_score)

                scores = [frag_isolated_contribution, frag_pair_contribution, mod_sa_score]
                factor = 1.0 / 3.0
                sp_score = round(delta * factor * sum([num for num in scores if num]), 4)
                delta = 'one'
            else:
                sp_score = 'None'
                scores = ['None']
                delta  = 'zero'

        return sp_score, graph, graph2, delta                                         

       
