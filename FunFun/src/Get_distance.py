import logging
import numpy as np
from pandas import DataFrame
from itertools import product 
from scipy.spatial.distance import cosine
from FunFun.src.Get_functional import get_functionality

a_logger = logging.getLogger()
a_logger.setLevel(logging.DEBUG)


def normalize(kmers):
    """
    The functuion noramlize kmers counts.
    ----------
    kmers : kmers dict
        Kmers counts dictionary.
    Returns
    -------
    kmers : dict
        Normalized khmer dictionary.
    """
    norm = sum(list(kmers.values()))
    
    for kmer in kmers.keys():
        
        kmers[kmer] = kmers[kmer]/ norm
    
    return kmers

def get_kmers_fereq(seq, k=2):
    """"
    Calculating k mers abundance.
    ----------
    kmers : Bio.Seq.Seq
        ITS sequence.
    Returns
    -------
    kmers : dict
        Normalized khmer dictionary.
    """
    kmers_combination = np.sort([''.join(kmer) for kmer in list(product("AGTC", repeat=k))])
    kmers = {kmer : 0 for kmer in kmers_combination}
    step = 1
    start = 0
    end = k

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

def get_matrix(its, *args):
    """
    Calculating cosine distance vector for target fungi.
    ----------
    its : Bio.Seq.Seq
        ITS sequence.
    base : DataFrame
        Base with 5-mers counts.
    Ortology_group : list
        list with kofam ortology groups.
    functionality : DataFrame
        Representation of group ortology for every fungi in base.
    Meta_micom : dict
        Functional content for input fungi.
    non_pred : int
        Counts of unpredicted samples.
    K : int
        K nearest neighbours.
    e : float
        ε neighborhood.
    Returns
    -------
    Ortology_group : list
        list with kofam ortology groups.
    Meta_micom : dict
        Functional content for input fungi.
    non_pred : int
        Counts of unpredicted samples.
    """
    base, Ortology_group, functionality, Meta_micom, non_pred, K, e = args
    fungi_sample = its.id
    marker_frequence = list(get_kmers_fereq(its.seq, k=5).values())
    distance_vector = list(map(lambda x: cosine(x, marker_frequence), base.values))
    base_subset = DataFrame(data={fungi_sample : distance_vector}, index=list(base.index))

    try:

        its_function = get_functionality(fungi_sample, base_subset, functionality, n_neigbors=K, epsilont=e)
        a_logger.debug(f'Functional for {fungi_sample} was succesfully predicted!')
        if Ortology_group == []:

            Ortology_group = list(its_function.keys())
        
        predicted_vector = its_function.values

        if str(predicted_vector[0]) == 'nan':

            a_logger.debug(f'Functional for {its.id} was not predicted!\nTry to change -e ...')    
            non_pred += 1

        Meta_micom[fungi_sample] = predicted_vector
    except:

        a_logger.debug(f'Fuctional for {fungi_sample} was not predicted! Try to change e and K params.')
        pass

    return Ortology_group, Meta_micom, non_pred
