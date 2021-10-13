import numpy as np

def calcuate_penalty(seq1_, seq2_, gap_penalty, gap_extension, mismatch_penalty):
        
        PENALTY = 0
        open_gap = False
        length_gap = 0

        for idx in range(len(seq1_)):
            if seq1_[idx] == '-' or seq2_[idx] == '-':
                if open_gap == True:
                    
                    length_gap += 1
                    
                    
                else:

                    length_gap += 1
                    open_gap = True

            else:
                if length_gap != 0:
                    
                    PENALTY += - gap_penalty - gap_extension * (length_gap - 1)
                    
                if seq2_[idx] != seq1_[idx]:

                    PENALTY -= mismatch_penalty

                open_gap = False
                length_gap = 0
                    
        return PENALTY
    
def get_score(aln, match, gap_penalty, gap_extension, mismatch_penalty):
    """
    Affin function of scorich 
    Penatly_of_gap = N∑i=1 (Gap_penalty[i]+ Gap_extension_penalty[i] * i)
    The function normalize score on length of aligned sequences.
    """
    SCORE = 0
    sequences = []
    
    
    for seq in aln:
        
        sequences.append(seq.seq)
    
    seq1, seq2 = np.array([i for i in sequences[0]]), np.array([i for i in sequences[1]])
    SCORE += len(seq1[seq1 == seq2]) * match # Adding matches
    seq1_ = seq1[seq1 != seq2] # Parts only with gaps and mismatches
    seq2_ = seq2[seq1 != seq2]
    PENALTY = calcuate_penalty(seq1_, seq2_, gap_penalty, gap_extension, mismatch_penalty)
    SCORE += PENALTY
    SCORE = SCORE / (len(seq1) * match)
    
    return np.round(SCORE, 3)
