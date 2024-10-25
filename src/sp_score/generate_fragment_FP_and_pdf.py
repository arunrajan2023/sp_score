import pandas as pd
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
from rdkit import Chem
import os
# from pscore.read_data import *


def generate_fragment_FP_and_pdf():
    """
    Steps to be completed if you add new fragments and/or their functional groups:
        1. Edit the following entries (filenames) to include or exclude them 
        2. Check read_data.py to ensure the generated fingerprint files (in /data/fingerprint) are read
    """
    #Names of CSV files inside data/input (without extension) folder  
    gs  = [
        'acidchlorides',
        'acids',
        'alcohols',
        'aldehydes',
        'amides',
        'amines',
        'esters',
        'heteroatomics',
        'hydrocarbons',
        'ketones',
        'phosphonyl_acids',
        'phosphonyl_PO_PO2_PO3_PO4',
        'PN',
        'sulfides',
        'sulfonic_acid',
        'sulfonyl_SO2_or_SO3'
        ]     

    def fingerprint_and_pdf(gs):
        """
        Step: Read CSV files and store FPs in data/fingerprint folder 
        """

        def pdf_image(df, labels, fname):
            """
            Render PDF files corresponding to each list of fragment SMILES and their chemical functional group.
            """
            #---------------------------------------------------------------------------------------
            pdf_fname             = fname
            png_fname             = fname
            factor                = 16
            tags, group           = list(labels), list(smiles_strings)
            verbose, render_image = False, False
            nfrag                 = len(group)
            nimages               = int(nfrag/factor) if nfrag % factor == 0 else int(nfrag/factor) + 1
            image_list            = []

            for iimage in range(1, nimages+1):
                if verbose: print("Drawing ...",iimage, "of", nimages, " image files:")     
                image_tags = png_fname + str(iimage)

                fname   = image_tags+'.png'
                istart  = (iimage -1) * factor
                iend    = istart + factor
                ss      = group[istart:iend]
                legends = [str(int(istart + ictr + 1)) + " " + frag  for ictr,  frag in enumerate(ss)]
                if verbose: print("legends:", legends)
                mols    = [Chem.MolFromSmiles(_) for _ in ss]
                img     = Chem.Draw.MolsToGridImage(mols, molsPerRow=4, legends=legends, subImgSize=(400, 200))
                if render_image: img.save(fname)  #(m x n) grid PNGs
                if iimage > 1:
                    image_list.append(img.convert('RGB'))
                if iimage == 1:
                    img0   = img                  #first image img0

            if nfrag > 0:
                draw = ImageDraw.Draw(img)

                #write PDF
                pdf_fname  = pdf_fname + "_images.pdf"
                img0.save(pdf_fname,   save_all=True, append_images=image_list)     #append on img0 as PDF

            if os.path.exists(pdf_fname):
                print(pdf_fname, "rendered!")
            else:
                if verbose: print("PDF rendering failed.")
            return None

        #reading CSV files and generating FPs
        for file in gs:
            df             = pd.read_csv("data/input/"+file+".csv")
            smiles_strings = df["smiles"]
            labels         = df['Unnamed: 0']
            pdf_image(df, labels, file)
            
            def function(smi):
                a = smi.count("*")
                m = Chem.MolFromSmiles(smi) 
                b = m.GetNumAtoms()
                c = m.GetRingInfo().NumRings()
                return [a,b,c]
                
            df['fragment_fp']  = df.apply(lambda row: Chem.RDKFingerprint(Chem.MolFromSmiles(row["smiles"])).ToBitString(), axis = 1)
            df['fragment_fp2'] = df.apply(lambda row: function(row["smiles"]), axis = 1)
            df.to_csv(path + "data/fingerprint/" + file + "_and_FPs" + ".csv")
            print("Modified the fingerprints!")    
            print()        

    #function call to fingerprint the fragments
    fingerprint_and_pdf(gs)