import argparse
from pandas import read_csv
import numpy as np 


parser = argparse.ArgumentParser()

group1 = parser.add_argument_group("Base arguments")
group2 = parser.add_argument_group("Primer arguments")
# If ypu want to find its inner genome
group1.add_argument('-genome', 
                    type=str,
                    default=None,
                    help='Sequence of genome assemble in faa format.')
group1.add_argument('-ITS1', 
                    type=str,
                    default=None,
                    help='Sequence of ITS1 in fasta format.')
group1.add_argument('-ITS2', 
                    type=str,
                    default=None,
                    help='Sequence of ITS2 in fasta format.')
group1.add_argument('-out', 
                    type=str,
                    default='./FunFun_output'
                    help='Get output of your file.')

group2.add_argument('-primerA', 
                    type=str,
                    default=None,
                    help='Primers to extract ITS region.')
group2.add_argument('-MIN_ITS_length',
                    type=int,
                    default=100,
                    help='If you know minimal lengtht of potential ITS.')
group2.add_argument('-MAX_ITS_length',
                    type=int,
                    default=1000,
                    help='If you know maximum lengtht of potential ITS.')
group2.add_argument('-mismatches',
                    type=int,
                    type=3,
                    default='Number of mismatches for primer.')
group2.add_argument('-primerB', 
                    type=str,
                    default=None,
                    help='Primers to extract ITS region.')

args = parser.parse_args()

# Get arguments
genome = args.genome
primerA = args.primerA
primerB = args.primerB
ITS1 = args.ITS1
ITS2 = args.ITS2
out = args.out
MIN_ITS_length = args.MIN_ITS_length
MAX_ITS_length = args.MAX_ITS_length
mismatches = args.mismatches

if  ITS1 is None and ITS2 is None and genome is None:

    print('Give a genome or one of IST!')
    import sys
    sys.exit()

# Find ITS 
get_itss(primerA, primerB, MIN_ITS_length, MAX_ITS_length, mismatches, genome, out)