import pandas as pd
import numpy as np

def index_list1_list2(list1, list2): #list 2 longer than list 1
    idx_list1_in_list2 = []
    idx_list2_in_list1 = []
    for i in range(len(list2)):
        for j in range(len(list1)):
            if list2[i] == list1[j]:
                idx_list1_in_list2.append(i)
                idx_list2_in_list1.append(j)
    return idx_list1_in_list2, idx_list2_in_list1

def Diff_list(li1, li2):
    """
    Inform on the difference between two lists

    Parameters
    ----------
    li1: list numero 1
    li2: list numero 2

    Returns
    -------
    Elements than are in list 1 but not in list 2
    """

    return (list(set(li1) - set(li2)))


df_speculoos_michael = pd.read_csv('speculoos.txt',delimiter=' ')

Filter = ['I+z'] * len(df_speculoos_michael)
texp_spc = [0] * len(df_speculoos_michael)
nb_hours_treshold = np.zeros(len(df_speculoos_michael))
nb_hours_surved = np.zeros(len(df_speculoos_michael))
telescope = [[]]  * len(df_speculoos_michael)
df_speculoos_michael.insert(2, "Filter",Filter,True)
df_speculoos_michael.insert(3, "texp_spc",texp_spc,True)
df_speculoos_michael.insert(4, "nb_hours_threshold",nb_hours_treshold,True)
df_speculoos_michael.insert(5, "nb_hours_surved",nb_hours_surved,True)
df_speculoos_michael.insert(len(df_speculoos_michael.columns)-1, "telescope",telescope,True)

df_speculoos_elsa = pd.read_csv('SPECULOOS_target_list_v2.txt',delimiter=' ')
index_michael_in_elsa,index_elsa_in_michael = index_list1_list2(df_speculoos_michael['Sp_ID'], df_speculoos_elsa['Sp_ID'])

df_speculoos_michael['Sp_ID'][index_elsa_in_michael] = df_speculoos_elsa['Sp_ID'][index_michael_in_elsa]
df_speculoos_michael['Filter'][index_elsa_in_michael] = df_speculoos_elsa['Filter'][index_michael_in_elsa]
df_speculoos_michael['texp_spc'][index_elsa_in_michael] = df_speculoos_elsa['texp_spc'][index_michael_in_elsa]
df_speculoos_michael['telescope'][index_elsa_in_michael] = df_speculoos_elsa['telescope'][index_michael_in_elsa]
df_speculoos_michael['nb_hours_threshold'][index_elsa_in_michael] = df_speculoos_elsa['nb_hours_threshold'][index_michael_in_elsa]
df_speculoos_michael['nb_hours_surved'][index_elsa_in_michael] = df_speculoos_elsa['nb_hours_surved'][index_michael_in_elsa]
df_speculoos_michael['telescope'][index_elsa_in_michael] = df_speculoos_elsa['telescope'][index_michael_in_elsa]

program1 = np.where((df_speculoos_michael['Program'] ==  1))[0]
program2 = np.where((df_speculoos_michael['Program'] ==  2))[0]
program3 = np.where((df_speculoos_michael['Program'] ==  3))[0]
#nb_hours_treshold = [100] * len(df_speculoos_michael)
nb_hours_treshold[program1] = 200
nb_hours_treshold[program2] = 100
nb_hours_treshold[program3] = 100


print()