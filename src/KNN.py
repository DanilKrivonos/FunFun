import numpy as np
from os import listdir
from numpy import intersect1d
from pandas import read_csv
from pandas import DataFrame
from src.Get_score import get_score
from Bio.SeqIO import parse

def add_in_matix(path_to_its1, path_to_its2, fungi, out):
    """

    """
    if not(path_to_its1 is None):
        
        matrix = read_csv('./data/matrix_possitions_ITS1.csv', index_col=[0])
        data_vector = {}
        path_aln = listdir('{}/Fasta_alignment/ITS1/'.format(out))
        target = ''

        for aln in path_aln:
            
            parsed_aln = parse('{}/Fasta_alignment/ITS1/{}'.format(out, aln), 'fasta')
            data_vector[aln.split('vs')[1][: -len('.aln')]] = [get_score(parsed_aln, 1, 1, 0, 1)]

        # Add vector in matrix
        
        vector = DataFrame(data=data_vector, index=[fungi])
        matrix = matrix.append(vector).T[: -1].append(vector)
        
    if not(path_to_its2 is None):
        
        matrix = read_csv('./data/matrix_possitions_ITS2.csv', index_col=[0])
        data = {}
        path_aln = listdir('{}/Fasta_alignment/ITS2/'.format(out))
        target = ''

        for aln in path_aln:

            parsed_aln = parse('{}/Fasta_alignment/ITS2/{}'.format(out, aln), 'fasta')
            data[aln.split('vs')[1]: -len('.aln')] = get_score(parsed_aln, 1, 1, 0, 1)

        # Add vector in matrix
        vector = DataFrame(data)
        matrix = matrix.append(vector).T[:-1].append(vector)

    
    return matrix

def get_functionality(fungi, matrix_possitions, kofam_ontology, n_neigbors):

    min_vals = np.sort(matrix_possitions.loc[fungi].values)[: n_neigbors + 1][1: ]
    neighbors = []
    neighbors_func = []
    
    for i in min_vals:
        
        neighbor = matrix_possitions.loc[fungi][matrix_possitions.loc[fungi] == i].index[0]
        neighbors.append(neighbor)
        neighbors_func.append(list(kofam_ontology[neighbor].index))

    neighbor_func_info = kofam_ontology[neighbors].mean(axis=1)
    dict_of_methabolic_function = dict(zip(neighbor_func_info.index, neighbor_func_info.values))
    
    return dict_of_methabolic_function
