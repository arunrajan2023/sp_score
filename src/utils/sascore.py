#
# calculation of synthetic accessibility score as described in:
#
# Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
# Peter Ertl and Ansgar Schuffenhauer
# Journal of Cheminformatics 1:8 (2009)
# http://www.jcheminf.com/content/1/1/8
#
# several small modifications to the original paper are included
# particularly slightly different formula for marocyclic penalty
# and taking into account also molecule symmetry (fingerprint density)
#
# for a set of 10k diverse molecules the agreement between the original method
# as implemented in PipelinePilot and this implementation is r2 = 0.97
#
# peter ertl & greg landrum, september 2013
#

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import math
from collections import defaultdict
import pickle as cPickle
import sys,gzip,time
import os

#================================================================================================================
#Rewritten script on July 16th, 2020 by Arun   (excluded functions, used class)
#================================================================================================================

class sa_score:
        def __init__(self, **kwargs):
                return None

        def __new__(self, **kwargs):
                self.smiles          = None
                self._fscores = None
                
                for key, value in kwargs.items(): 
                    if key == "smiles":
                            self.smiles_tag = key
                            self.smiles     = value
                    else:
                            print("*"*30)
                            print("Key error", key)

                outf     = None 
                ismiles         = self.smiles
                if '[e]' in ismiles:
                        ismiles = ismiles.replace('[e]', '[*]')
                        ismiles = ismiles.replace('[d]', '[*]')
                        ismiles = ismiles.replace('[t]', '[*]')
                        ismiles = ismiles.replace('[g]', '[*]')        
                ismiles = Chem.CanonSmiles(ismiles)           #canonicalize            
                mols    = [Chem.MolFromSmiles(ismiles)]

                sasmiles_array, sascore_array = [], []             
                count = {}
                for i,m in enumerate(mols):
                    if m is None:
                        continue

                    #calculate SA score
                    if self._fscores is None:    
                            name          = 'fpscores'
                            self._fscores = cPickle.load(gzip.open('src/utils/%s.pkl.gz'%name))

                            outDict  = {}
                            for i in self._fscores:
                                for j in range(1,len(i)):
                                    outDict[i[j]] = float(i[0])
                            self._fscores = outDict                

                    # fragment score
                    fp = rdMolDescriptors.GetMorganFingerprint(m,2)  #<- 2 is the *radius* of the circular fingerprint
                    fps = fp.GetNonzeroElements()
                    score1 = 0.
                    nf = 0
                    for bitId,v in fps.items():
                        nf += v
                        sfp = bitId
                        score1 += self._fscores.get(sfp,-4)*v
                    score1 /= nf

                    # features score
                    nAtoms             = m.GetNumAtoms()
                    nChiralCenters     = len(Chem.FindMolChiralCenters(m,includeUnassigned=True))
                    ri                 = m.GetRingInfo()
                    mol                = m
                    if ri is None:
                        ri=mol.GetRingInfo()
                    arings = [set(x) for x in ri.AtomRings()]
                    spiros=set()
                    for i,ari in enumerate(arings):
                        for j in range(i+1,len(arings)):
                            shared=ari&arings[j]
                            if len(shared)==1:
                                spiros.update(shared)
                    nSpiro=len(spiros)

                    # find bonds that are shared between rings that share at least 2 bonds:
                    nBridge=0
                    brings = [set(x) for x in ri.BondRings()]
                    bridges=set()
                    for i,bri in enumerate(brings):
                        for j in range(i+1,len(brings)):
                          shared=bri&brings[j]
                          if len(shared)>1:
                            atomCounts=defaultdict(int)
                            for bi in shared:
                              bond = mol.GetBondWithIdx(bi)
                              atomCounts[bond.GetBeginAtomIdx()]+=1
                              atomCounts[bond.GetEndAtomIdx()]+=1
                            tmp=0
                            for ai,cnt in atomCounts.items():
                              if cnt==1:
                                tmp+=1
                                bridges.add(ai)
                            #if tmp!=2: # no need to stress the users
                              #print 'huh:',tmp
                    nBridgeheads,nSpiro= len(bridges),nSpiro  
                    nMacrocycles       = 0
                    for x in ri.AtomRings():
                        if len(x)>8: nMacrocycles+=1

                    sizePenalty       = nAtoms**1.005 - nAtoms
                    stereoPenalty     = math.log10(nChiralCenters+1)
                    spiroPenalty      = math.log10(nSpiro+1)
                    bridgePenalty     = math.log10(nBridgeheads+1)
                    macrocyclePenalty = 0.
                    # ---------------------------------------
                    # This differs from the paper, which defines:
                    #  macrocyclePenalty = math.log10(nMacrocycles+1)
                    # This form generates better results when 2 or more macrocycles are present
                    if nMacrocycles > 0: macrocyclePenalty = math.log10(2)

                    score2 = 0. -sizePenalty -stereoPenalty -spiroPenalty -bridgePenalty -macrocyclePenalty

                    # correction for the fingerprint density
                    # not in the original publication, added in version 1.1
                    # to make highly symmetrical molecules easier to synthetise
                    score3 = 0.
                    if nAtoms > len(fps):
                        score3 = math.log(float(nAtoms) / len(fps)) * .5

                    sascore = score1 + score2 + score3

                    # need to transform "raw" value into scale between 1 and 10
                    min = -4.0
                    max = 2.5
                    sascore = 11. - (sascore - min + 1) / (max - min) * 9.
                    # smooth the 10-end
                    if sascore > 8.: sascore = 8. + math.log(sascore+1.-9.)
                    if sascore > 10.: sascore = 10.0
                    elif sascore < 1.: sascore = 1.0                     
                    return sascore

#================================================================================================================
#================================================================================================================

#usage: 
# import sa_score
# sa_score(smiles="*C*")

#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#                  
