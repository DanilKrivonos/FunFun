import os
from os import listdir
from pandas import read_csv
from subprocess import call 

def get_name(ITS):

    with open(ITS, 'r') as its:
        for line in its:
            if '>' in line:

                name = line.split(' ')[0][1: ].replace('\n', '')

    return name

def make_concat(ITS1, ITS2, out):
    """
    Making concatenate fasta.
    

    """
    if not(ITS1 is None):

        ITS1_BASE = './data/ITS1_BASE/'
        try:

            os.mkdir('{}/ITS1_concat/'.format(out))

        except:

            print('Directory is exist!') 

        name = get_name(ITS1)

        for fungi in listdir(ITS1_BASE):
            with open('{}/ITS1_concat/{}vs{}.fasta'.format(out, name, fungi[: -len('_genomic.fasta')]), 'w') as concat:
                with open(ITS1, 'r') as its1:
                    for line in its1:
                        
                        concat.write(line)
                
                with open('{}/{}'.format(ITS1_BASE, fungi), 'r') as vs_part:
                    for line in vs_part:

                        concat.write(line)
        with open('{}/ITS1_concat/{}vs{}.fasta'.format(out, name, name), 'w') as concat:
            with open(ITS1, 'r') as its1:
                for line in its1:
                    
                    concat.write(line)
            with open(ITS1, 'r') as its1:
                for line in its1:
                    
                    concat.write(line)
            

    if not(ITS2 is None):

        ITS2_BASE = './data/ITS2_BASE'

        try:

            os.mkdir('{}/ITS2_concat/'.format(out))

        except:

            print('Directory is exist!')

        with open(ITS2, 'r') as its2:

            name = get_name(its2)

            for fungi in listdir(ITS2_BASE):
                with open('{}/ITS2_concat/{}vs{}.fasta'.format(out, name, fungi[: -len('_genomic.fasta')]), 'w') as concat:
                    for line in its2:

                        concat.write(line)
                    
                    with open('{}/{}'.format(ITS2_BASE, fungi), 'r') as vs_part:
                        for line in vs_part:

                            concat.write(line)
    return name

def get_aglighnment(ITS1, ITS2, out):
    """
    Making alighnment with muscle.

    """
    its1_path_out = None
    its2_path_out = None
    name = make_concat(ITS1, ITS2, out)

    try:
        
        os.mkdir('{}/Fasta_alignment/'.format(out))

    except:

        print('Directory is exist!')

    if not(ITS1 is None):

        its1_path_out = '{}/Fasta_alignment/ITS1/'.format(out)

        try:

            os.mkdir(its1_path_out)

        except:

            print('Directory is exist!')

        path_concat = listdir('{}/ITS1_concat'.format(out))

        for concat in  path_concat:

            call('muscle  -quiet -in {}/ITS1_concat/{} -out {}/{}.aln'.format(out, concat, its1_path_out, concat[:- 6]), shell=True)
            
    if not(ITS2 is None):

        its2_path_out = '{}/Fasta_alignment/ITS2/'.format(out) 

        try:
            
            os.mkdir(its2_path_out)

        except:

            print('Directory is exist!')

        path_concat = listdir('{}/ITS2_concat'.format(out))

        for concat in  path_concat:
            
            call('muscle  -quiet -in {}/ITS2_concat/{} -out {}/{}.aln'.format(out, concat, its2_path_out, concat[: -6]), shell=True)
      
    return its1_path_out, its2_path_out, name