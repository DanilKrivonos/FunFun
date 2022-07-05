import os
import sys
import logging
import argparse
import shutil
import warnings
import numpy as np 
import plotly.express as px
from FunFun.src.Get_distance import get_matrix
from pandas import read_csv, DataFrame
from Bio.SeqIO import parse
from subprocess import call

warnings.filterwarnings('ignore')

def main():
    
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Base arguments")
    group2 = parser.add_argument_group("Algorithm Parameter")

    group1.add_argument('-its', # Should be concatinate file
                        type=str,
                        default=None,
                        help='Sequence of ITS1 in fasta format.')
    group1.add_argument('-type', 
                        type=str,
                        default='concatenate',
                        help='Type of its region: its1, its2, concatenate.')
    group1.add_argument('-out', 
                        type=str,
                        default='./FunFun_output',
                        help='Get output of your file.')

    group2.add_argument('-K', 
                        type=int,
                        default=10,
                        help='K nearest neighbors. Default k=10.')
    group2.add_argument('-e', 
                        type=float,
                        default=0.5,
                        help='Epsilon distance for nearest neighbor. 0 <= ε <= 1. Default ε=0.5.')
    if len(sys.argv)==1:

        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    # Get arguments
    ITS = args.its
    TYPE = args.type.upper()
    K = args.K
    e = args.e
    out = args.out
    functionality = read_csv(os.path.dirname(os.path.abspath(__file__)) + '/data/functionality.tsv', sep='\t', index_col=[0])
    # Making directory
    if out is None:
        
        out = './FunFun_output'
        
    if os.path.exists(out):

        shutil.rmtree(out)
    
    os.mkdir(out)
    a_logger = logging.getLogger()
    a_logger.setLevel(logging.DEBUG)
    output_file_handler = logging.FileHandler(f"{out}/FunFun.log")
    stdout_handler = logging.StreamHandler(sys.stdout)
    a_logger.addHandler(output_file_handler)
    a_logger.addHandler(stdout_handler)
    # Fasta availability check
    if  ITS is None:

        a_logger.debug('Give an ITS!')
        sys.exit()

    # Make realignment on base
    base = read_csv(os.path.dirname(os.path.abspath(__file__)) + f'/data/{TYPE}_base.tsv', sep='\t', index_col=[0])

    # Parsing matrix
    a_logger.debug('Calculation functionality ...')
    marker_seq = parse(ITS, 'fasta')
    Meta_micom = {}
    non_pred = 0
    Ortology_group = []
    c = 0

    for its in marker_seq:
        
        Ortology_group, Meta_micom, non_pred = get_matrix(its, base, Ortology_group, functionality, Meta_micom, non_pred, K, e)
        c += 1

    predicted_persentage = ((c-non_pred)/c) * 100
    a_logger.debug(f'Functional was predicted for {predicted_persentage}% of all fungies.')
    
    # Generation output 
    Meta_micom['Function'] = Ortology_group
    Meta_micom = DataFrame(Meta_micom).set_index('Function')
    Meta_micom = Meta_micom.assign(m=Meta_micom.mean(axis=1)).sort_values('m').drop('m', axis=1)
    fig = px.bar(Meta_micom.T, x=Meta_micom.columns, y=Meta_micom.index, width=1500, height=1000)
    fig.update_layout(legend_traceorder="reversed")
    fig.write_html(f"{out}/Functional_community.html")
    Meta_micom.to_csv(f'{out}/Results.tsv', sep='\t')
    a_logger.debug('Job is done!')

if __name__ == "__main__":
    main()
