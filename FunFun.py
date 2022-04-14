import os
import sys
import logging
import argparse
import shutil
import numpy as np
from pandas import read_csv, DataFrame
from Bio.SeqIO import parse
from itertools import product 
from subprocess import call
from scipy.spatial.distance import pdist, squareform
import plotly.express as px
#from pandas.core.dtypes.missing import na_value_for_dtype

def main():
    
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Base arguments")
    group2 = parser.add_argument_group("Algorithm Parameter")

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
    import sys
    if len(sys.argv)==1:

        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Get arguments
    ITS1 = args.ITS1
    ITS2 = args.ITS2
    CONCAT = args.CONCAT
    K = args.K
    e = args.e
    out = args.out
    functionality = read_csv('./data/functionality.tsv', sep='\t', index_col=[0])
    # Making directory
    a_logger = logging.getLogger()
    a_logger.setLevel(logging.DEBUG)

    if out is None:
        
        out = './FunFun_output'
        
    if os.path.exists(out):

        shutil.rmtree(out)
    
    os.mkdir(out)
    
    output_file_handler = logging.FileHandler(f"{out}/FunFun.log")
    stdout_handler = logging.StreamHandler(sys.stdout)

    # Fasta availability check
    if  ITS1 is None and ITS2 is None and CONCAT is None:

        a_logger.debug('Give an ITS!')
        sys.exit()

    def get_names(ITS_fasta):

        names = []

        fasta = parse(ITS_fasta, 'fasta')

        for line in fasta:

            names.append(line.id)

        return names

    # Make realignment on base

    if ITS1 is not None:

        base = read_csv('./data/ITS1_base.tsv', sep='\t', index_col=[0])
        marker_seq = ITS1
    if ITS2 is not None:
        
        base = read_csv('./data/ITS2_base.tsv', sep='\t', index_col=[0])
        marker_seq = ITS2

    if CONCAT is not None:
        
        base = read_csv('./data/CONCAT_base.tsv', sep='\t', index_col=[0])
        marker_seq = CONCAT
    # Parsing matrix

    def normalize(kmers):
        
        norm = sum(list(kmers.values()))
        
        for kmer in kmers.keys():
            
            kmers[kmer] = kmers[kmer]/ norm
        
        return kmers

    print('uibibui')
    def get_kmers_fereq(seq, k=2):
        
        kmers = {"".join(kmer) : 0 for kmer in list(product("AGTC", repeat=k))}
        
        step = 1
        start = 0
        end = k
        cken = []
        while end != len(seq) - 1:
            
            kmers[str(seq[start: end])] += 1
            start, end = start + step, end + step
            
        step = 1
        start = 0
        end = k

        while end != len(seq.reverse_complement()) - 1:
            
            kmers[str(seq.reverse_complement()[start: end])] += 1        
            start, end = start + step, end + step
            
        kmers = normalize(kmers)
        
        return kmers

    def get_functionality(fungi, matrix_possitions, kofam_ontology, n_neigbors, epsilont=0.5):

        neighbor = matrix_possitions.loc[fungi].sort_values()
        neighbor = neighbor.drop(fungi)
        neighbor = neighbor[neighbor <= epsilont]

        if len(neighbor) > n_neigbors:
            
            neighbor = neighbor[: n_neigbors]
            
        if len(neighbor[neighbor == 0].index) > 0:
            
            dict_of_methabolic_function = kofam_ontology[neighbor[neighbor == 0].index].mean(axis=1)
            
        else:
            
            dict_of_methabolic_function = kofam_ontology[neighbor.index].mean(axis=1)

        return dict_of_methabolic_function

    def get_functionality(fungi, matrix_possitions, kofam_ontology, n_neigbors, epsilont=0.5):

        neighbor = matrix_possitions.loc[fungi].sort_values()
        neighbor = neighbor.drop(fungi)
        neighbor = neighbor[neighbor <= epsilont]

        if len(neighbor) > n_neigbors:
            
            neighbor = neighbor[: n_neigbors]
            
        if len(neighbor[neighbor == 0].index) > 0:
            
            dict_of_methabolic_function = kofam_ontology[neighbor[neighbor == 0].index].mean(axis=1)
            
        else:
            
            dict_of_methabolic_function = kofam_ontology[neighbor.index].mean(axis=1)

        return dict_of_methabolic_function

    Meta_micom = {}
    marker_seq = parse(marker_seq, 'fasta')
    c = 0

    for its in marker_seq:

        print(f'Processing {c}', end='')
        fungi_sample = its.id
        marker_frequence = get_kmers_fereq(its.seq, k=5)
        marker_frequence['Fungi'] = fungi_sample
        base_subset= base.append(DataFrame(data=marker_frequence, index=[fungi_sample])[base.columns])
        distance_matrix = squareform(pdist(base_subset.values, 'cosine'))
        distance_matrix = DataFrame(data=distance_matrix, index=base_subset.index, columns=base_subset.index)
        its_function = get_functionality(fungi_sample, distance_matrix, functionality, n_neigbors=K, epsilont=e)
        
        if "Ortology group" not in Meta_micom:

            Meta_micom["Ortology group"] = list(its_function.keys())
        
        Meta_micom[f'Fraction score {fungi_sample}'] = its_function.values
        print('\r', end='')

    Meta_micom = DataFrame(Meta_micom)
    print(Meta_micom.T)
    fig = px.bar(Meta_micom.T, y="Ortology group", x=Meta_micom.columns[1: ], width=1500, height=1000)
    fig.write_html(f"{out}/Functional_barplot.html")
    Meta_micom.to_csv(f'{out}/Results.tsv', sep='\t', index=False)
    a_logger.debug('Job is done!')

if __name__ == "__main__":
    main()
