import numpy as np
import pandas as pd
import os
from pathlib import Path
import csv

def isint(value):
    '''Check if a value (line startwith) an expression of an integer'''
    try:
        j = int(value)
        return True
    except ValueError:
        return False

'''reading the translative databases'''
# database of human pfam entries from UniProt
df_pfam_uniprot = pd.read_csv('C:\\D_drive\\Schuldiner\\Human All Pfam Uniprot.csv')
print(df_pfam_uniprot.head(), df_pfam_uniprot.shape)
# database of PDB entries from UniProt
df_pdb = pd.read_csv('C:\\D_drive\\Schuldiner\\pdb_chain_uniprot.csv')
print(df_pdb.head(), df_pdb.shape)
# database of pfam descriptions from UniProt
df_mapping_pfam = pd.read_csv('C:\\D_drive\\Schuldiner\\Pfam_mapping_descriptions.txt', sep='\t')
print(df_mapping_pfam.head(), df_mapping_pfam.shape)

i=0
semi_dict = {}

directory_in_str = 'C:\\D_drive\\Schuldiner\HHbits Yeast\\'
pathlist = Path(directory_in_str).glob('**/*.hhr')
for path in pathlist:
    path_in_str = str(path)
    f = open(path_in_str)
    print('######################' + path_in_str)
    line = f.readline()
    if 'dubious' in str(line) or 'Dubious' in str(line):
        continue
    line = str(line).split()
    Yeast_Systematic_Name = line[1]
    Yeast_Symbol = line[2]
    Yeast_ID = line[3][:-1]
    for sample in f:
        line = str(sample).split()
        if line != []:
            if isint(line[0]):
                name = line[1]
                real_line = str(sample[35::]).split()
                prob = real_line[0]
                if not float(prob) >= 95.0:
                    continue
                E_value = real_line[1]
                P_value = real_line[2]
                Score = real_line[3]
                human_counter = 0
                g = open(path_in_str)
                for row in g:
                    phrase = str(row).split()
                    if phrase != []:
                        if phrase[0].startswith('>') and phrase[0][1::] == name:
                            print(str(row).rstrip('\n'))
                            human_counter = 1
                            full_row_g = str(row)
                            print('1')
                            break
                if name[:2] == 'PF':# and human_counter == 1: # If the protein is Pfam scribed:
                    print(name)
                    name = str(name.split('.')[0])
                    try:
                        df_slice = df_mapping_pfam[df_mapping_pfam['PFAM AC'] == name]
                        pfam_id = df_slice['PFAM ID'].iloc[0]
                        pfam_clan_id = df_slice['PFAM Clan ID'].iloc[0]
                        pfam_clan_name = df_slice['PFAM Clan Name'].iloc[0]
                        pfam_description = df_slice['PFAM Description'].iloc[0]
                    except:
                        pfam_id = ''
                        pfam_clan_id = ''
                        pfam_clan_name = ''
                        pfam_description = ''
                    semi_dict[i] = []
                    semi_dict[i].append(Yeast_Systematic_Name)
                    semi_dict[i].append(Yeast_Symbol)
                    semi_dict[i].append(Yeast_ID)
                    semi_dict[i].append(str(name))
                    semi_dict[i].append('')
                    semi_dict[i].append('')
                    semi_dict[i].append('')
                    semi_dict[i].append('')
                    semi_dict[i].append('')
                    semi_dict[i].append(pfam_id)
                    semi_dict[i].append(pfam_clan_id)
                    semi_dict[i].append(pfam_clan_name)
                    semi_dict[i].append(pfam_description)
                    semi_dict[i].append(prob)
                    semi_dict[i].append(E_value)
                    semi_dict[i].append(P_value)
                    semi_dict[i].append(Score)
                    semi_dict[i].append(full_row_g)
                    i+=1
                    break
                elif human_counter == 1:
                    if name[0] == 'd' and len(name) == 7: #If the protein is d_ prefixed:
                        print(name)
                        name = name.lstrip('d')
                        real_name = name[1:5].lower()
                        chain = name[5::]
                        if chain.endswith('_'):
                            chain = chain[:-1]
                    elif name[0] == 'g' and '.' in name:
                        print(name)
                        name = name.lstrip('g')
                        real_name = name[1:5].lower()
                        chain = 'pass'
                    else: #If protein has a regular PDB name:
                        print(name)
                        wait_name = name.split('_')
                        real_name = str(wait_name[0]).lower()
                        chain = str(wait_name[1])
                        name = str(name).lower()
                    df_temp3 = df_pdb[df_pdb['PDB'] == real_name]
                    id = ''
                    symbol = ''
                    names = ''
                    brief = ''
                    description = ''
                    uni = ''
                    if not df_temp3.empty:
                        for entry in range(df_temp3.shape[0]):
                            if str(df_temp3['CHAIN'].iloc[entry]) == chain or str(df_temp3['CHAIN'].iloc[entry]) == chain.lower()\
                                    or str(df_temp3['CHAIN'].iloc[entry]) == chain.upper() or chain == 'pass':
                                uni = df_temp3['SP_PRIMARY'].iloc[entry]
                                df_temp2 = df_pfam_uniprot[df_pfam_uniprot['Uniprot ID'] == uni]
                                if not df_temp2.empty:
                                    id = str(df_temp2['ID'].iloc[0])
                                    symbol = str(df_temp2['Symbol'].iloc[0])
                                    names = str(df_temp2['Names'].iloc[0])
                                    brief = str(df_temp2['Brief Description'].iloc[0])
                                    description = str(df_temp2['Description'].iloc[0])
                    semi_dict[i] = []
                    semi_dict[i].append(Yeast_Systematic_Name)
                    semi_dict[i].append(Yeast_Symbol)
                    semi_dict[i].append(Yeast_ID)
                    semi_dict[i].append(str(name))
                    semi_dict[i].append(id)
                    semi_dict[i].append(symbol)
                    semi_dict[i].append(uni)
                    semi_dict[i].append(names)
                    semi_dict[i].append(brief)
                    semi_dict[i].append('')
                    semi_dict[i].append('')
                    semi_dict[i].append('')
                    semi_dict[i].append(description)
                    semi_dict[i].append(prob)
                    semi_dict[i].append(E_value)
                    semi_dict[i].append(P_value)
                    semi_dict[i].append(Score)
                    semi_dict[i].append(full_row_g)
                    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    i+=1
            elif line[0][0] == '>':
                human_counter = 0
                print('0')
                break

df_output = pd.DataFrame.from_dict(semi_dict, orient='index',) #columns=['Yeast Systematic Name', 'Yeast Symbol',
                                                                        #'Yeast ID', 'Human Database Name', 'Human ID',
                                                                        #'Human Symbol', 'Human Name',
                                                                        #'Human Brief Description', 'Human Description',
                                                                        #'Human Pfam ID', 'Human Uniprot ID',
                                                                        #'HHpred prob', 'HHpred E_value',
                                                                        #'HHpred P_value', 'HHpred Score'])
export_csv = df_output.to_csv(r'C:\\D_drive\\Schuldiner\\try222222.csv',
                              index = False, header=True) #Don't forget to add '.csv' at the end of the path
