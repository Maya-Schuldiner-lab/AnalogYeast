import pandas as pd

'''Opening databases'''
# database of HHsearch results
df_HHbit = pd.read_csv('C:\\D_drive\\Schuldiner\\DBv7 ec.csv')
print(df_HHbit.head(), df_HHbit.shape)
# database of SGD yeast proteins
df_all_ORFs = pd.read_csv('C:\\D_drive\\Schuldiner\\Yeast All ORFs + EC 20Apr20.csv')
print(df_all_ORFs.head(), df_all_ORFs.shape)
# list of uncharacterized proteins
df_unknown_genes = pd.read_csv('C:\\D_drive\\Schuldiner\\Unknown proteins.csv')
print(df_unknown_genes.head(), df_unknown_genes.shape)

'''Humans'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Homo sapiens' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Human Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Human Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Nematodes'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Caenorhabditis elegans' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Nematodes Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Nematodes Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''E. Coli'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Escherichia coli' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of E.coli Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of E.coli Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Drosophila'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Drosophila melanogaster' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Drosophila Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Drosophila Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Arabidopsis'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Arabidopsis thaliana' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Arabidopsis Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Arabidopsis Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Zebrafish'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Danio rerio' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Zebrafish Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Zebrafish Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Mouse'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Mus musculus' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Mouse Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Mouse Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Stacking'''
'''Human and Mouse'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Mus musculus' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Homo sapiens' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Human+Mouse Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Human+Mouse Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Human and Mouse and Zebrafish'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Mus musculus' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Homo sapiens' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Danio rerio' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Human+Mouse+Zebrafish Homologs of all ORFs: ' + str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))


human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Human+Mouse+Zebrafish Homologs of Unknown ORFs: ' + str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Human and Mouse and Zebrafish and Drosophila'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Mus musculus' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Homo sapiens' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Danio rerio' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Drosophila melanogaster' in df_HHbit.loc[i, 'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Human+Mouse+Zebrafish+Drosophila Homologs of all ORFs: ' +
      str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Human+Mouse+Zebrafish+Drosophila Homologs of Unknown ORFs: ' +
      str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Human and Mouse and Zebrafish and Drosophila and Nematode'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Mus musculus' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Homo sapiens' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Danio rerio' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Drosophila melanogaster' in df_HHbit.loc[i, 'Organism']:
        human_homologs_idxs.append(i)
    elif 'Caenorhabditis elegans' in df_HHbit.loc[i, 'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Human+Mouse+Zebrafish+Drosophila+Nematode Homologs of all ORFs: ' +
      str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Human+Mouse+Zebrafish+Drosophila+Nematode Homologs of Unknown ORFs: ' +
      str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Human and Mouse and Zebrafish and Drosophila and Nematode and Arabidopsis'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Mus musculus' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Homo sapiens' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Danio rerio' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Drosophila melanogaster' in df_HHbit.loc[i, 'Organism']:
        human_homologs_idxs.append(i)
    elif 'Caenorhabditis elegans' in df_HHbit.loc[i, 'Organism']:
        human_homologs_idxs.append(i)
    elif 'Arabidopsis thaliana' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Human+Mouse+Zebrafish+Drosophila+Nematode+Arabidopsis Homologs of all ORFs: ' +
      str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Human+Mouse+Zebrafish+Drosophila+Nematode+Arabidopsis Homologs of Unknown ORFs: ' +
      str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))


'''Human and Mouse and Zebrafish and Drosophila and Nematode and Arabidopsis and E.coli'''
human_homologs_idxs = []
for i in range(df_HHbit.shape[0]):
    if type(df_HHbit.loc[i,'Organism']) == float:
        continue
    if 'Mus musculus' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Homo sapiens' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Danio rerio' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Drosophila melanogaster' in df_HHbit.loc[i, 'Organism']:
        human_homologs_idxs.append(i)
    elif 'Caenorhabditis elegans' in df_HHbit.loc[i, 'Organism']:
        human_homologs_idxs.append(i)
    elif 'Arabidopsis thaliana' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
    elif 'Escherichia coli' in df_HHbit.loc[i,'Organism']:
        human_homologs_idxs.append(i)
df_human_HHbit_homologs = df_HHbit.iloc[human_homologs_idxs,:]
all_yeast_ORFs = list(set(list(df_all_ORFs['Yeast Systematic Name'])))
all_yeast_unknown_ORFs = list(set(list(df_unknown_genes['ORF'])))

human_homologs_counter = 0
for gene in all_yeast_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_homologs_counter += 1
print('\nPercentage of Human+Mouse+Zebrafish+Drosophila+Nematode+Arabidopsis+E.coli Homologs of all ORFs: ' +
      str(human_homologs_counter/len(all_yeast_ORFs)))
print(human_homologs_counter, len(all_yeast_ORFs))

human_unknown_homologs_counter = 0
for gene in all_yeast_unknown_ORFs:
    if gene in list(df_human_HHbit_homologs['Yeast Systematic Name']):
        human_unknown_homologs_counter += 1
print('Percentage of Human+Mouse+Zebrafish+Drosophila+Nematode+Arabidopsis+E.coli Homologs of Unknown ORFs: ' +
      str(human_unknown_homologs_counter/len(all_yeast_unknown_ORFs)))
print(human_unknown_homologs_counter, len(all_yeast_unknown_ORFs))