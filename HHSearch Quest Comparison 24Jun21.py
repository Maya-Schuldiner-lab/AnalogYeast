import numpy as np
import pandas as pd
import os
from pathlib import Path
import csv


'''Loading databases as DataFrames'''
# the Quest database
df_QUEST = pd.read_csv('C:\\D_drive\\Schuldiner\\ORTHOLOGY-ALLIANCE_COMBINED_yeast human homologs 2020 modified 12Jun21.csv')
print(df_QUEST.head(), df_QUEST.shape)
# HHsearch results that include only human hits
df_HHpred = pd.read_csv('C:\\D_drive\\Schuldiner\\DBv7 only human.csv')
print(df_HHpred.head(), df_HHpred.shape)
# SGD database of all Yeast proteins
df_all_yeast_orfs = pd.read_csv("C:\\D_drive\\Schuldiner\\Yeast All ORFs.csv")
print(df_all_yeast_orfs.head(), df_all_yeast_orfs.shape)

# remove merged proteins that were discounted as actual ORFs
removed_orfs = ['S000000082', 'S000000279','S000000304','S000000509','S000000518','S000000565',
                'S000000567','S000000625','S000000651','S000000654','S000000658','S000000676',
                'S000001682','S000001683','S000001851','S000001888','S000002196','S000002882',
                'S000003014','S000003504','S000003554','S000003555','S000003558','S000005614',
                'S000005766','S000006294','S000007609','S000028531','S000028549','S000028559',
                'S000028570','S000028577','S000028807','S000028842','S000082441','S000113724',
                'S000115154','S000117700','S000178054']
all_orfs_list = [x for x in list(df_all_yeast_orfs["Yeast ID"]) if x not in removed_orfs]


quest_list = []
for entry in range(df_QUEST.shape[0]):
    id_yeast_list = df_QUEST['Gene1ID'].iloc[entry]
    id_yeast_list = id_yeast_list.split(':')
    s_yeast_id = id_yeast_list[1]
    id_yeast = str(id_yeast_list[0]) + 'ID:' + id_yeast_list[1]
    if s_yeast_id in all_orfs_list:
        quest_list.append(id_yeast)

hhsearch_list = []
for entry in range(df_HHpred.shape[0]):
    id_yeast = str(df_HHpred.loc[entry, "SGDID"])
    s_yeast_id = id_yeast.split(':')[1]
    if s_yeast_id in all_orfs_list:
        hhsearch_list.append(id_yeast)

quest_list = list(set(quest_list))
hhsearch_list = list(set(hhsearch_list))
print(quest_list, hhsearch_list)
unique_quest = [x for x in quest_list if x not in hhsearch_list]
unique_hhsearch = [x for x in hhsearch_list if x not in quest_list]
cross = [x for x in hhsearch_list if x in quest_list]
union = list(set(quest_list + hhsearch_list))

print("Length of genome: " + str(len(all_orfs_list)))
print("Unique to Quest: " + str(len(unique_quest)))
print("Unique to HHsearch: " + str(len(unique_hhsearch)))
print("Union: " + str(len(union)))
print("Cross: " + str(len(cross)))

