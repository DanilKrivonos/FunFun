from subprocess import call 
from Bio.SeqIO import parse

def make_primer_file(primerA, primerB, MIN_ITS_length, MAX_ITS_length, out):
    """
    Recorfing file with primers for ipcress
    Parameters
    ----------
    primerA : str
        Sequence of forvard primer.
    primerB : str
        Sequence of reverse primer.
    MIN_ITS_length : int
        Minimal length of insilico PCR product
    MAX_ITS_length: int 
        Maximum length of insilico PCR product.
    out : str
        Path to output directory.
    """
    with open('{}/primer_file.txt'.format(out)) as primer_file:

        primer_file.write('YOURE_PRIMER\t{}\t{}\t{}\t{}'.format(primerA,
                                                                primerB,
                                                                MIN_ITS_length,
                                                                MAX_ITS_length))

def call_ipcress(primer_path, genome, mismatches, out):
    """
    """

    call('ipcress -i {} -s > {} -m {} -p False -P True > ./IPCR_{}.prod'.foramt(primer_path, 
                                                                                genome, 
                                                                                mismatches, 
                                                                                out))

def make_fasta():



def get_itss(primerA, primerB, MIN_ITS_length, MAX_ITS_length, mismatches, genome, out):
    """
    """
    # Gnerating file with primers for ipcress
    if primerA not is None and primerB not is None:
    
        make_primer_file(primerA, primerB, MIN_ITS_length, MAX_ITS_length, out)
        primer_path = './primer_file.txt'
    # PCR stage
    call_ipcress(primer_path, genome, mismatches, out)
