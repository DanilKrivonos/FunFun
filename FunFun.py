import os
import argparse
from pandas import read_csv
from pandas.core.dtypes.missing import na_value_for_dtype
from src.model import get_functionality

parser = argparse.ArgumentParser()

group1 = parser.add_argument_group("Base arguments")
group2 = parser.add_argument_group("Algorithm Parameter")
# If ypu want to find its inner genome
group1.add_argument('-ITS1', # Should be concatinate file
                    type=str,
                    default=None,
                    help='Sequence of ITS1 in fasta format.')
group1.add_argument('-ITS2', 
                    type=str,
                    default=None,
                    help='Sequence of ITS2 in fasta format.')
group1.add_argument('-CONCAT', 
                    type=str,
                    default=None,
                    help='Sequence of ITS1, 5.8rRNA, ITS2 in fasta format.')
group2.add_argument('-K', 
                    type=str,
                    default=10,
                    help='K nearest neighbors.')
group2.add_argument('-e', 
                    type=str,
                    default=0.5,
                    help='Epsilon distance for nearest neighbor. ε ⊆ [0; 1]')
group1.add_argument('-out', 
                    type=str,
                    default='./FunFun_output',
                    help='Get output of your file.')
args = parser.parse_args()

# Get arguments
ITS1 = args.ITS1
ITS2 = args.ITS2
CONCAT = args.CONCAT
K = args.K
e = args.e
out = args.out
# Making directory
try:

    os.mkdir(out)

except:

    print('Directory is exist!')

# Fasta availability check
if  ITS1 is None and ITS2 is None and CONCAT is None:

    print('Give an ITS!')
    import sys
    sys.exit()

from Bio.SeqIO import parse

def get_names(ITS_fasta):

    names = []

    with parse(ITS_fasta, 'fasta') as fasta:
        for line in fasta:

            names.append(line.id)

    return names

names = get_names(ITS1)

# Make realignment on base
from subprocess import call

if ITS1 is not None:

    call('muscul -profile -in1 ./data/all_ITS1.msa -in2 {} -clwstrictout ./{}/realigned_{}.clw'.format(ITS1, out, ITS1[: -3]), shell=True)
    # Get distance matrix 
    call('clustalo --p1 ./{}/realigned_{}.clw --distmat-out=ITS1_dist_mat.txt --full --force'.format(ITS1[: -3]))

if ITS2 is not None:
    
    call('muscul -profile -in1 ./data/all_ITS2.msa -in2 {} -clwstrictout ./{}/realigned_{}.clw'.format(ITS2, out, ITS1[: -3]), shell=True)
    call('clustalo --p1 ./{}/realigned_{}.clw --distmat-out=ITS2_dist_mat.txt --full --force'.format(ITS1[: -3]))

if CONCAT is not None:
    
    call('muscul -profile -in1 ./data/all_CONCAT.msa -in2 {} -clwstrictout ./{}/realigned_{}.clw'.format(CONCAT, out, ITS1[: -3]), shell=True)
    call('clustalo --p1 ./{}/realigned_{}.clw --distmat-out=CONCAT_Kimura_dist_mat.txt --full --force'.format(ITS1[: -3]))

# Parsing matrix

from pandas import DataFrame
import numpy as np

realigned_seq = './{}/realigned_{}.clw'.format(out, ITS1[: -3])

def get_matrix(realigned_seq):

    identity_persent = []
    columns = []
    dist_matrix = read_csv(realigned_seq, sep='\t', comment='#', header=None)

    for i in dist_matrix.index:
        try:
            if len(dist_matrix.iloc[i].values[0].split()) == 1:
                continue

            identity_persent.append(np.array(eval(', '.join(dist_matrix.iloc[i].values[0].split()[1:]))))
            columns.append('_'.join(dist_matrix.iloc[i].values[0].split()[0].split('_')[1:4]).replace('\n', ''))
            
        except:
            continue

    identity_persent = np.array(identity_persent)
    data_dict = {}

    for i in range(len(identity_persent)):
        for j in range(len(identity_persent[i])):
            
            if columns[i] not in  data_dict:
                
                data_dict[columns[i]]= {}
                
            data_dict[columns[i]][columns[j]] = identity_persent[i][j]
            
    Distanece_DF = DataFrame(data=data_dict)

    return Distanece_DF

Distanece_DF = get_matrix(realigned_seq)
# Importing functions data

functions = read_csv('./data/kofam_ontology_NEW_DATA.csv', sep='\t', index_col=[0])

# Make prediction and generating report file 

metafunctional = []

def make_fun_report(names):
    
    for fun in names:

        # Get FUNction for FUNgi

        dict_of_methabolic_function = get_functionality(fun, Distanece_DF, functions, n_neigbors=K, epsilon=e)
        metafunctional.append(list(dict_of_methabolic_function.values()))
        # Function ...
        make_report = lambda function: res_file.write('{}\t{}\n'.format(function, dict_of_methabolic_function[function]))

        with open('{}/Results_{}.tsv'.format(out, fun), 'w') as res_file:

            res_file.write('KEGG function\tShare of all functions\n')        
            list(map(make_report, list(dict_of_methabolic_function.keys())))
    
    function_list = list(dict_of_methabolic_function.keys())

    return function_list

# Activate function

function_list = make_fun_report(names)

# Methafungal functional 
metafunctional = np.array(metafunctional)
metafunctional = np.mean(metafunctional, axis=0)

def get_metafunctional(metafunctional, function_list):

    with open('{}/Results_METAFUNCTIONAL.tsv'.format(out), 'w') as res_file:

        for idx in range(function_list):

            res_file.write('{}\t{}\n'.format(function_list[idx], metafunctional[idx]))

    print('Methafunctional annotated')

# Activate function 

get_metafunctional(metafunctional, function_list)

    


print('Job is done!')