import pandas as pd

'''Opening databases'''
# database of HHsearch results
df_HHbit = pd.read_csv('C:\\D_drive\\Schuldiner\\DBv7.csv')
print(df_HHbit.head(), df_HHbit.shape)
# database of SGD yeast proteins
df_all_ORFs = pd.read_csv('C:\\D_drive\\Schuldiner\\Yeast All ORFs + EC 20Apr20.csv')
print(df_all_ORFs.head(), df_all_ORFs.shape)
# list of uncharacterized proteins
df_unknown_genes = pd.read_csv('C:\\D_drive\\Schuldiner\\Unknown proteins.csv')
print(df_unknown_genes.head(), df_unknown_genes.shape)
# database of MalaCards
df_malacards = pd.read_csv('C:\\D_drive\\Schuldiner\\MalaCards_Genes.txt', delimiter='\t')
print(df_malacards.head(), df_malacards.shape)

df_HHbit['Diseases'] = ''
for entry in range(df_HHbit.shape[0]):
    gene = df_HHbit.loc[entry,'Yeast Systematic Name']
    print(gene)
    organism = df_HHbit.loc[entry,'Organism']
    if type(organism) == float:
        continue
    if 'homo sapiens' in organism:
        symbol = df_HHbit.loc[entry,'Human Symbol']
        df_malacards_temp = df_malacards[df_malacards['GeneSymbol'] == symbol]
        diseases = str(list(set(list(df_malacards_temp['DiseaseName'])))).rstrip(']').lstrip('[').replace("'","")
        df_HHbit.loc[entry,'Diseases'] = diseases

export_csv = df_HHbit.to_csv(r'C:\D_drive\Schuldiner\Yeast_All_ORFs HHbit +95threshold all_organisms 27Feb21.csv',
                             index = False, header=True)


