import os
import argparse
from pandas import read_csv
from src.Construct_data import get_aglighnment
from src.KNN import add_in_matix, get_functionality

parser = argparse.ArgumentParser()

group1 = parser.add_argument_group("Base arguments")
# If ypu want to find its inner genome
group1.add_argument('-ITS1', 
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
group1.add_argument('-out', 
                    type=str,
                    default='./FunFun_output',
                    help='Get output of your file.')

args = parser.parse_args()

# Get arguments
ITS1 = args.ITS1
ITS2 = args.ITS2
out = args.out
# Making directory
try:

    os.mkdir(out)

except:

    print('Directory is exist!')


if  ITS1 is None and ITS2 is None:

    print('Give an ITS!')
    import sys
    sys.exit()

# Find ITS
#get_itss(primerA, primerB, MIN_ITS_length, MAX_ITS_length, mismatches, genome, out)
# get_aglighnment
its1_path_out, its2_path_out, fungi = get_aglighnment(ITS1, ITS2, out)
# get new matrix
matrix = add_in_matix(its1_path_out, its2_path_out, fungi, out)
# get unsupervised knn 
kofam_ontology = read_csv('./data/kofam_ontology.tsv', sep='\t', index_col=[0])
dict_of_methabolic_function = get_functionality(fungi, matrix, kofam_ontology, n_neigbors=2)

#return  answer
with open('{}/Results.tsv'.format(out), 'w') as res_file:
    res_file.write('KEGG function\tShare of all functions\n')
    for key in dict_of_methabolic_function.keys():

        res_file.write('{}\t{}\n'.format(key, dict_of_methabolic_function[key]))

print('Job is done!')