#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function

import pandas as pd
import numpy as np
from rdkit import Chem

import pickle
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults
from rdkit.Chem import AllChem
import math
import pickle
from collections import OrderedDict

######### Similarity ##########
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import AllChem

############
from copy import deepcopy

import os, errno

import sys
print(sys.version)
print(sys.version_info)
pd.options.display.max_rows = 10
pd.options.display.float_format = '{:.1f}'.format


from IPython.display import Image
from IPython.display import display

from rdkit.Chem import PandasTools
import pandas as pd
import os
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

import argparse

# In[53]:


class Pep2struc:
    def __init__(self, in_df_path=None, out_df_path=None, topology=None):
        self.in_df_path  = in_df_path
        self.out_df_path = out_df_path
        self.topology    = topology


    def condensation_reaction(self, molecules, c, toplogyType, bad_monomer):
        if (bad_monomer): return None, True
        n = len(molecules)
        #print("Total number of monomers: ", n )
        bad_composition = False
        growingPeptide = {}
        #changeheadcules[0] = changehead(molecules[0])
        #changeheadcules[n-1] = changehead(molecules[n-1])
        head = molecules[0]
        tail = molecules[n-1]
        try:
            for i in range(0, len(molecules)):
                m = molecules[i]
                #print("i:", i)
                new_possibilities = []
                if (not growingPeptide):
                    m = self.protectAllAtoms(m)
                    m = self.isItAminoAcid(m)
                    s = Chem.MolToSmiles(m)
                    growingPeptide[i]= new_possibilities.append([m,(-1,-1), s])
                else:
                    possibilities = growingPeptide[len(growingPeptide)-1]
                    while (len(possibilities)!= 0):
                        p = possibilities.pop()
                        m1= p[0]
                        m2= m
                        new_possibilities.extend(self.conductingReaction(i, m1, m2, n, toplogyType))
                        #print("new_possibilities: ", len(new_possibilities))
                        #print()

                growingPeptide[i] = new_possibilities
                #print("growingPeptide_" + str(i) +"_size: ", len(growingPeptide[i]))
                #print("#################")

            #print("growingPeptideMap_keys_size: ", len(growingPeptide))

            if(toplogyType != 'linear'):
                growingPeptide = self.conductingReaction2(growingPeptide, toplogyType, head, tail)
            #print("done")

            finals= self.namingAllPossibilities2(growingPeptide, toplogyType)
            #print("Finals: ", len(finals))
            #print("#########################################")
            if (len(finals)==0): bad_composition = True
            return finals, bad_composition

        except Exception as e:
            bad_composition = True
            #raise
            print(e)
            return None, n, bad_composition

    def protecAmidbond(self, m):
        amidep = Chem.MolFromSmarts('[N;$(NC=[O,S])]')
        for match in m.GetSubstructMatches(amidep):
                m.GetAtomWithIdx(match[0]).SetProp('_protected','1')
        return m

    def protectAllAtoms(self, m):
        for a in m.GetAtoms():
            #print(a.GetSymbol())
            a.SetProp('_protected','1')
        return m

    def changehead(self, m):
        repl = Chem.MolFromSmiles('[NX3,NX4+][C(Cl)(Cl)(Cl)][CX3](=[OX1])[O,N]')
        patt = Chem.MolFromSmarts('[NX3,NX4+][C][CX3](=[OX1])[O,N]')
        rms = AllChem.ReplaceSubstructs(m,patt,repl)
        display(rms[0])
        Chem.MolToSmiles(rms[0])
        return rms[0]

    def isItAminoAcid(self, m):
        generic_amino_acid = Chem.MolFromSmarts('[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]')
        fake_amino_acid = Chem.MolFromSmarts('[OX2H][CX4H]([*])[CX3](=[OX1])[O,N]')
        gly = Chem.MolFromSmarts('[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N])]')
        generic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
        hydroxyl = Chem.MolFromSmarts('[OX2H]')
        amines = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)].[NX3;H2,H1;!$(NC=O)]')

        all_matches = m.GetSubstructMatches(generic_amino_acid)
        #print("All_matches_len: ", len(all_matches))
        if(len(all_matches)>=1):
            for match in all_matches:
                for i in range (len(match)):
                    a = m.GetAtomWithIdx(match[i])
                    a.ClearProp('_protected')
            return m

        all_matches_gly = m.GetSubstructMatches(gly)
        #print("All_matches_len_gly: ", len(all_matches_gly))
        if(len(all_matches_gly)>=1):
            for a in m.GetAtoms():
                a.ClearProp('_protected')
            return m

        all_matches_fake_aa = m.GetSubstructMatches(fake_amino_acid)
        #print("All_matches_len_fake: ", len(all_matches_fake_aa))
        if(len(all_matches_fake_aa)>=1):
            for match in all_matches_fake_aa:
                for i in range (len(match)):
                    a = m.GetAtomWithIdx(match[i])
                    a.ClearProp('_protected')
            return m

        all_matches_acid = m.GetSubstructMatches(generic_acid)
        all_matches_hydroxyl = m.GetSubstructMatches(hydroxyl)
        all_matches_amines = m.GetSubstructMatches(amines)
        if(len(all_matches_acid)>=1):
            for match in all_matches_acid:
                for a in m.GetAtoms():
                    a.ClearProp('_protected')

            for match in all_matches_hydroxyl:
                for a in m.GetAtoms():
                    a.ClearProp('_protected')

            for match in all_matches_amines:
                for a in m.GetAtoms():
                    a.ClearProp('_protected')
            return m

        elif(len(all_matches_hydroxyl)>=1):
            for match in all_matches_hydroxyl:
                for a in m.GetAtoms():
                    a.ClearProp('_protected')

            for match in all_matches_amines:
                for a in m.GetAtoms():
                    a.ClearProp('_protected')
            return m

        elif(len(all_matches_amines)>=1):
            for match in all_matches_amines:
                for a in m.GetAtoms():
                    a.ClearProp('_protected')
            return m
        else:
            for a in m.GetAtoms():
                    a.ClearProp('_protected')
            return m

        return m

    def conductingReaction(self, reaction_number, m1, m2, n, toplogyType):
        #print("m1: ", Chem.MolToSmiles(m1))
        #print("m2 ", Chem.MolToSmiles(m2))
        #print(display(m1), display(m2))
        m2 = self.protectAllAtoms(m2)
        m2 = self.isItAminoAcid(m2)

        possibilities = []
        if ((reaction_number == n-1) and (toplogyType == "linear")):
            ps, reaction_code = self.possibleReactions(m1, m2)
            if(ps):
                for i in range(len(ps)):
                    for j in range(len(ps[0])):
                        s=Chem.MolToSmiles(ps[i][j])
                        peptide=ps[i][j]
                        infos=[peptide, (i,j), s, "s_" + str(reaction_number) +"_"+ reaction_code]
                        name = str(toplogyType) + '_' + str(infos[-1]) + "_" + str(infos[1])
                        peptide.SetProp("_Name", name)
                        possibilities.append(infos)
            return possibilities

        else:
            ps, reaction_code = self.possibleReactions(m1, m2)
            if(ps):
                for i in range(len(ps)):
                    for j in range(len(ps[0])):
                        s=Chem.MolToSmiles(ps[i][j])
                        infos=[ps[i][j] , (i,j), s, "s_" + str(reaction_number) + "_" + reaction_code]
                        possibilities.append(infos)
            return possibilities

    def conductingReaction2(self, growingPeptide, toplogyType, head, tail):
        last_key = len(growingPeptide)-1
        possibilities=[]
        for p in growingPeptide[last_key]:
            m = Chem.MolFromSmiles(p[2])
            self.protecAmidbond(m)
            #print("Cyclic_start: ", p[2])
            #display(m)
            ps, cyclic_code = self.cycliczation(m)
            infos=[]
            if(ps):
                for i in range(len(ps)):
                    for j in range(len(ps[0])):
                        s=Chem.MolToSmiles(ps[i][j])
                        #print("Cyclic: ", s)
                        peptide = ps[i][j]
                        infos = [peptide, (i,j), s, p[-1] + "_" + cyclic_code]
                        name = str(toplogyType) + '_' + str(infos[-1]) + "_" + str(infos[1])
                        peptide.SetProp("_Name", name)
                        possibilities.append(infos)

        if((toplogyType=='double cyclic') or (toplogyType=='other')):
            temp_possibilities = deepcopy(possibilities)
            for t in temp_possibilities :
                m=Chem.MolFromSmiles(p[2])
                ps, cyclic_code = cycliczation(m)
                infos=[]
                if(ps):
                    for i in range(len(ps)):
                        for j in range(len(ps[0])):
                            s=Chem.MolToSmiles(ps[i][j])
                            peptide=ps[i][j]
                            infos=[peptide,(i,j), s, p[-1] + "_" + cyclic_code]
                            name = str(toplogyType) + '_' + str(infos[-1]) + "_" + str(infos[1])
                            peptide.SetProp("_Name", name)
                            possibilities.append(infos)

        #print("total_cyclization: ", len(possibilities))
        growingPeptide[toplogyType]= possibilities
        return growingPeptide

    def possibleReactions(self, m1, m2):
        #amide_inter1 = AllChem.ReactionFromSmarts('[N:1][C:2][C:3](=[O:4])-[OD1].[N!H0:5]>>[N:1][C:2][C:3](=[O:4])[N:5]')
        amide_inter2 = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]')
        #amide_inter3 = AllChem.ReactionFromSmarts('[C:1][C:2](=[O:3])-[O:4].[N!H0:5]>>[C:1][C:2](=[O:3])[N:5]')
        #amide_inter4 = AllChem.ReactionFromSmarts('[C:1][O:2].[N!H0:3]>>[C:1][N:3]')
        #ester_inter1 = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[C:3](=[O:4])-[OD1]>>[C:1](=[O:2])-[OD1][C:3](=[O:4])')
        ester_inter2 = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[C:3]-[OD1]>>[C:1](=[O:2])-[OD1][C:3]')
        aromatic_condensation_inter = AllChem.ReactionFromSmarts('[c:1][OH:2].[N!H0:4]>>[c:1][N:4]')
        amine_inter1 = AllChem.ReactionFromSmarts('[N!H0:1].[C:2](=[O:3])-[OD1]>>[N:1][C:2](=[O:3])')
        amine_inter2 = AllChem.ReactionFromSmarts('[N!H0:1].[C:2][O:3]>>[N:1][C:2]')


        #ps = amide_inter1.RunReactants((m1, m2))
        #if(ps):
            #print("pass1")
            #return ps,"r_1"

        ps = amide_inter2.RunReactants((m1, m2))
        if(ps):
            #print("pass2")
            return ps,"r_2"

        #ps = amide_inter3.RunReactants((m1, m2))
        #if(ps):
            #print("pass3")
            #return ps,"r_3"

        #ps = amide_inter4.RunReactants((m1, m2))
        #if(ps):
            #print("pass4")
            #return ps,"r_4"

        #ps = ester_inter1.RunReactants((m1, m2))
        #if (ps):
            #print("pass5")
            #return ps,"r_5"

        ps = ester_inter2.RunReactants((m1, m2))
        if (ps):
            #print("pass6")
            return ps,"r_6"

        ps = aromatic_condensation_inter.RunReactants((m1, m2))
        if(ps):
            #print("pass7")
            return ps,"r_7"

        ps = amine_inter1.RunReactants((m1, m2))
        if(ps):
            #print("pass8")
            return ps,"r_8"

        ps = amine_inter2.RunReactants((m1, m2))
        if(ps):
            #print("pass9")
            return ps,"r_9"

        return None, "no_reaction"

    def cycliczation(self, m):
        ##
        #amide_intra1 = AllChem.ReactionFromSmarts('([C:1][C:2](=[O:3])-[OH].[N!H0:4])>>[C:1][C:2](=[O:3])[N:4]')
        ##
        amide_intra1 = AllChem.ReactionFromSmarts('([C:1][C:2](=[O:3])-[OH].[N!H0:4])>>[C:1][C:2](=[O:3])[N:4]')
        ester_intra1 = AllChem.ReactionFromSmarts('([C:1](=[O:2])-[OD1].[C:3]-[OD1])>>[C:1](=[O:2])-[OD1][C:3]')
        ester_intra2 = AllChem.ReactionFromSmarts('([C:1](=[O:2])-[OD1].[C:3]-[O:4])>>[C:1](=[O:2])-[OD1][C:3]')

        ps = amide_intra1.RunReactant(m,0)
        if(ps):
            #print("pass_C1")
            return ps,"c_1"

        #ps = ester_intra1.RunReactant(m,0)
        #if(ps):
            #return ps,"c_2"

        #ps = ester_intra2.RunReactant(m,0)
        #if(ps):
            #return ps,"c_3"

        return None, "no_cyclization"


    def namingAllPossibilities2(self, growingPeptide, toplogyType):
        finals=[]
        if(toplogyType!='linear'):
            last_key = toplogyType
            name_key = toplogyType
        else:
            last_key = len(growingPeptide)-1
            name_key = toplogyType

        return growingPeptide[last_key]


    def runConversion(self, in_df_path, topology):
        ##### You main dataframe with all the information and sequences ########
        main_path = "C:/Users/Sheri/Dropbox/Projects/NRP_QSAR/final_data/entire_data_as_Jan_2019"

        #### A file with all the monomers, in your case it is just 20 amino acids, please make sure to have smiles for each monomer##
        monomers= pd.read_csv(main_path + "/" + "all_monomers_in_Norine_DB.csv")

        #### a Map of monomers and their corresponding smiles ##########
        m2smilesMap={}
        for index, row in monomers.iterrows():
            m2smilesMap[str(row[7].strip())]=str(row[9].strip())


        data = pd.read_csv(in_df_path)
        print("input_File_shape: ", data.shape)
        data_positive = data[data["M1"]==1]
        print("data_positive_shape: ", data_positive.shape)

        bad_seqs = []

        list_df = []
        dict_f= {}
        counter = 0
        for index, row in data_positive.iterrows():
            if(index == counter):
                print(index)
                counter = index+5000
            c = row[1]
            #print(c)
            monos = c.split(",")
            if all((((i.capitalize()).strip()) in m2smilesMap) for i in monos):
                monos=[Chem.MolFromSmiles(m2smilesMap[(i.capitalize()).strip()]) for i in monos ]
                finals, bad_composition = self.condensation_reaction(monos, c, topology, False)
                #print("Finals_size: ", len(finals), "bad_composition: ", bad_composition)
                dict_f["size"] = len(finals)
                #dict_f["finals"] = finals
                dict_f["smiles"] = self.getAllSmiles(finals)
                dict_f["composition"] = c
                df_temp = pd.DataFrame.from_dict(dict_f)
                list_df.append(df_temp)
            else:
                bad_seqs.append(c)

        results = pd.concat(list_df)
        return results, bad_seqs

    def getAllSmiles(self, finals):
        list_s=[]
        for final in finals:
            list_s.append(final[2])
        return list_s

    def writeResult(self, results, out_df_path, topology):
        print(results.shape)
        results.to_csv(out_df_path + "peptides2smiles_2_" + topology + ".csv")
        display(results.head(n=5))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-in','--in_df_path',
                        help='Path of input file path.',  required=True, default=None)
    parser.add_argument('-o', '--out_df_path',
                        help='Path to new .csv file for saving the potential peptide dataframe.', required=True, default=None)
    parser.add_argument('-t', '--topology',
                        help='cyclic or linear.', required=True, default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    pep2struc = Pep2struc( in_df_path  = args.in_df_path,
                           out_df_path = args.out_df_path,
                           topology    = args.topology)

    results_l, bad_seqs_l = pep2struc.runConversion(args.in_df_path, args.topology)
    pep2struc.writeResult(results_l, args.out_df_path, args.topology)
