import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

'''Opening databases'''
# database of HHsearch results
df_HHbit = pd.read_csv('C:\\D_drive\\Schuldiner\\DBv7.csv')
print(df_HHbit.head(), df_HHbit.shape)
# database of metabolic pathway holes
df_pathways = pd.read_csv('C:\\D_drive\\Schuldiner\\Yeast Metabolic Pathway Holes.csv',encoding='latin1')
print(df_pathways.head(), df_pathways.shape)
# database of SGD Yeast proteins that are functionally uncharacterized
df_curated_unknown = pd.read_csv('C:\\D_drive\\Schuldiner\\Yeast Curated Unknown Genes.csv')
print(df_curated_unknown.head(), df_curated_unknown.shape)


big_dict = {}
for entry in range(df_HHbit.shape[0]):
    orf = str(df_HHbit['Yeast Systematic Name'].iloc[entry])
    if not orf in list(df_curated_unknown['Yeast Systematic Name']):
        continue
    human_id = str(df_HHbit['HumanID'].iloc[entry])
    if not 'ENSG' in human_id:
        continue
    if orf in big_dict:
        if human_id in big_dict[orf]:
            continue
    else:
        big_dict[orf] = []
    big_dict[orf].append(human_id)
    ec = str(df_HHbit['EC'].iloc[entry])
    if ec == 'nan' or ec == '':
        continue
    if ec.count(',') >= 1:
        ec = ec.split(',')
        ec_set = set(ec)
        ec = list(ec_set)
        if '' in ec:
            ec.remove('')
    else:
        ec = str(ec)
    print(orf, ec)
    for specify in range(df_pathways.shape[0]):
        incomplete = 0
        multiple = 0
        target = str(df_pathways['EC#'].iloc[specify])
        if target == '' or target == 'nan':
            continue
        elif target[-1] == '-':
            incomplete = 1
        elif '/' in target:
            multiple = 1
        if incomplete == 0 and multiple == 0:
            if type(ec) == str:
                shoot = target.split('.')
                real_ec = ec.split('.')
                for i in [3,2,1]:
                    if shoot[:i+1] == real_ec[:i+1]:
                        before = str(df_pathways[str(i+1)+'-num fit prot'].iloc[specify])
                        if before == 'nan':
                            df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = orf
                        else:
                            quick = before + ',' + orf
                            df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = str(quick)
                        break
            elif type(ec) == list:
                for instance in ec:
                    shoot = target.split('.')
                    occurrence = instance.split('.')
                    for i in [3,2,1]:
                        if shoot[:i+1] == occurrence[:i+1]:
                            before = str(df_pathways[str(i+1)+'-num fit prot'].iloc[specify])
                            if before == 'nan':
                                df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = orf
                            else:
                                quick = before + ',' + orf
                                df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = str(quick)
                            break
        elif incomplete == 1:
            if type(ec) == str:
                shoot = target.split('.')
                real_ec = ec.split('.')
                for i in [2,1]:
                    if shoot[:i+1] == real_ec[:i+1]:
                        before = str(df_pathways[str(i+1)+'-num fit prot'].iloc[specify])
                        if before == 'nan':
                            df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = orf
                        else:
                            quick = before + ',' + orf
                            df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = str(quick)
                        break
            elif type(ec) == list:
                for instance in ec:
                    shoot = target.split('.')
                    occurrence = instance.split('.')
                    for i in [2,1]:
                        if shoot[:i+1] == occurrence[:i+1]:
                            before = str(df_pathways[str(i+1)+'-num fit prot'].iloc[specify])
                            if before == 'nan':
                                df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = orf
                            else:
                                quick = before + ',' + orf
                                df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = str(quick)
                            break
        elif multiple == 1:
            if type(ec) == str:
                shoot = target.split('/')
                for aim in shoot:
                    aim = aim.split('.')
                    real_ec = ec.split('.')
                    for i in [3,2,1]:
                        if aim[:i+1] == real_ec[:i+1]:
                            before = str(df_pathways[str(i+1)+'-num fit prot'].iloc[specify])
                            if before == 'nan':
                                df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = orf
                            else:
                                quick = before + ',' + orf
                                df_pathways[str(i+1)+'-num fit prot'].iloc[specify] = str(quick)
                            break
            elif type(ec) == list:
                for instance in ec:
                    shoot = target.split('/')
                    for aim in shoot:
                        aim = aim.split('.')
                        occurrence = instance.split('.')
                        for i in [3,2,1]:
                            if aim[:i+1] == occurrence[:i+1]:
                                before = str(df_pathways[str(i+1) + '-num fit prot'].iloc[specify])
                                if before == 'nan':
                                    df_pathways[str(i+1) + '-num fit prot'].iloc[specify] = orf
                                else:
                                    quick = before + ',' + orf
                                    df_pathways[str(i+1) + '-num fit prot'].iloc[specify] = str(quick)
                                break
for entry in df_pathways.index:
    complete_id = str(df_pathways.loc[entry, "4-num fit prot"])
    if str(complete_id) != str(np.nan):
        complete_id = list(set(str(df_pathways.loc[entry, "4-num fit prot"]).split(",")))
        df_pathways.loc[entry,'4-num fit prot (name, #hits)'] = len(complete_id)
    else:
        df_pathways.loc[entry, '4-num fit prot (name, #hits)'] = 0
    partly_complete_id = df_pathways.loc[entry, "3-num fit prot"]
    if str(partly_complete_id) != str(np.nan):
        partly_complete_id = list(set(str(df_pathways.loc[entry,"3-num fit prot"]).split(",")))
        df_pathways.loc[entry,"3-num fit prot (name, #hits)"] = len(partly_complete_id)
    else:
        df_pathways.loc[entry, "3-num fit prot (name, #hits)"] = 0
export_csv = df_pathways.to_csv(r'C:\\D_drive\\Schuldiner\\Yeast Metabolic Pathway Holes + All Organisms Candidates 11Jun21.csv',index=False, header=True)  # Don't forget to add '.csv' at the end of the path

df_pathways = df_pathways.sort_values(["4-num fit prot (name, #hits)","3-num fit prot (name, #hits)"], ascending=False).reset_index()
fig, ax1 = plt.subplots(figsize=(9, 6))
ax1.plot(df_pathways.index, df_pathways["4-num fit prot (name, #hits)"],
         'o-', color='forestgreen', markersize=2, alpha=0.65)
# plt.ylim(-1,)
plt.xticks(df_pathways.index,rotation=90, fontsize=7)
plt.yticks(range(int(min(df_pathways["4-num fit prot (name, #hits)"])),
                 int(max(df_pathways["4-num fit prot (name, #hits)"]) + 1)))
ax1.set_ylabel('# of Proteins with 4-digits EC match', color='forestgreen', fontsize=12)
ax1.tick_params(axis='y', labelcolor='forestgreen')

ax2 = ax1.twinx()
ax2.plot(df_pathways.index, df_pathways["3-num fit prot (name, #hits)"],
         'o-', color='tab:blue', markersize=2, alpha=0.65)
ax1.set_xticklabels(df_pathways['Pathway(s) needing this reaction'])
# plt.ylim(-1,18.5)
plt.yticks(range(int(min(df_pathways["3-num fit prot (name, #hits)"])),
                 int(max(df_pathways["3-num fit prot (name, #hits)"]) + 2), 3))
ax2.set_ylabel('# of Proteins with 3-digits EC match', color='tab:blue', fontsize=12)
ax2.tick_params(axis='y', labelcolor='tab:blue')
plt.tight_layout()
jerry = plt.savefig("C:\\D_drive\\Schuldiner\\EC.pdf")



