import numpy as np
import random
import sys
import itertools as it
from sequence_alignment import minimum_penalty
########## from http://www.martinbroadhurst.com/greedy-set-cover-in-python.html #######
def set_cover(universe, subsets):
    """Find a family of subsets that covers the universal set"""
    elements = set(e for s in subsets for e in s)
    # Check the subsets cover the universe
    if elements != universe:
        return None
    covered = set()
    cover = []
    # Greedily add the subsets with the most uncovered points
    while covered != elements:
        subset = max(subsets, key=lambda s: len(s - covered))
        cover.append(subset)
        covered |= subset
 
    return cover
##########################################################################################

####### from https://stackoverflow.com/a/20072009 #######
def find_shortest_superstring(sequences):
    LONGEST_SUPERSTRING = ''.join(sequences)
    current_shortest = LONGEST_SUPERSTRING
    trim = len(current_shortest)-1
    seen_prefixes = set()
    for perm in it.permutations(sequences):
        candidate_string = ''.join(perm)[:trim]
        if candidate_string in seen_prefixes:
            continue
        seen_prefixes.add(candidate_string)
        while is_superstring(candidate_string, sequences):
            current_shortest = candidate_string
            candidate_string = candidate_string[:-1]
            trim = len(current_shortest)-1
    # print(len(current_shortest))
    return current_shortest

def is_superstring(s, sequences):
    return all(seq in s for seq in sequences)
##########################################################

##### from https://stackoverflow.com/a/20074602 #####
from collections import defaultdict
def dijkSuperstring(originalSeqs):
  paths = defaultdict(set)
  paths[0] =  { '' }
  while paths:
    minLength = min(paths.keys())
    while paths[minLength]:
      candidate = paths[minLength].pop()
      seqAdded = False
      for seq in originalSeqs:
        if seq in candidate:
          continue
        seqAdded = True
        for i in reversed(range(len(seq)+1)):
          if candidate.endswith(seq[:i]):
            newCandidate = candidate + seq[i:]
            paths[len(newCandidate)].add(newCandidate)
      if not seqAdded:  # nothing added, so all present?
        return candidate
    del paths[minLength]
########################################################


def check_match(gene1, gene2, len1, len2, _max, result):
    for i in range(1, min(len1, len2)+1):
        if (gene1[len1-i:] == gene2[:i]):
            if(_max < i):
                _max = i
                result=gene1+gene2[i:]    
    return result, _max

def overlap(gene1, gene2):
    result = ''
    _max=-sys.maxsize-1
    len1, len2 = len(gene1), len(gene2)
    # primeiro checa sufixo de gene1 com prefixo de gene2
    result, _max = check_match(gene1, gene2, len1, len2, _max, result)
    # primeiro checa prefixo de gene1 com sufixo de gene2
    result, _max = check_match(gene2, gene1, len2, len1, _max, result)
    return result

def shotgun(gene):
    _sum = len(gene)
    n = 12
    rnd_array = np.random.multinomial(_sum, np.ones(n)/n, size=1)[0]
    partial_genes = []
    j = 0
    for i in range(0, len(rnd_array)):
        partial_genes.append(gene[j:j+rnd_array[i]])
        j+=rnd_array[i]
    # print(partial_genes)
    random.shuffle(partial_genes)
    return partial_genes

def merge(genes):
    universe=set(g for g in genes)
    subsets=[set([i]) for i in universe]

    for i in range(len(genes)):
        for j in range(len(genes)):
            if(i!=j):
                ovij = overlap(genes[i], genes[j])
                if(len(ovij) != 0):
                    subsets.append(set([genes[i], genes[j]]))

    # print(subsets)
    # x=sorted(subsets,key=len)

    cover = set_cover(universe, subsets)
    overs=[]
    for c in cover:
        aux=c.pop()
        for l in c:
            aux=overlap(aux, l)    
        overs.append(aux)
        # break
    # print(overs)
    random.shuffle(overs)
    res=overs.pop(0)    
    for o in overs:
        ov=overlap(res, o)
        res=ov if len(ov)!=0 else res+o

    return res

gene = '''
GCAAGGAGTCGACCGTCGCAAACAAACGCCATGTGTCTCATGTTCTCGTTCGGTCTATCC
CAAAGCGGATCAGCTGGAAGCGGGATTCGAAAACAGCGACATCCCTGTCTTGGGCAAACG
CATTGCACACATTCCTTCGGGTGAAACAATGAAGTTACACACTCAGTAAGGCTTCTTGGC
TTATCACTTGTCATTGGGTGCTCAGTGGCCGATAGATGCAGTCTGTATATAAAGATGATC
AACACAACCTATAATGACGGCCATGTAACACCCCTAAGCAAACACCCTGAACCAAAACAT
CCTCTACACTCTATCTAGTCAATCTTCGGAATGGTGAACACCTGCACCTATCTCCCCCTC
AGCGGCAAGGTCGCCCTCGTGACCGGCGGCGGCCGAGGCATCGGGGCCGGCATTGCCCTC
GAACTGGCCCGTCGCGGGGCCTCCGTGGCCATCAACTACGGCCACAGCGCCAAGTCGGCC
CAAGAGGTCGTCGAAGCCATCCAGGCCATCGGCCGCCAGGCGGTCGCCATCCAGGCGGAT
CTGACGTGCGTGCCCAACATCGAGTCGCTCATCCAGGAGGTCGTCCGCCACTTCGGCCGG
CTGGACATTGTGGTGTCCAACTCCGGCATGGAAAAGTTCAAGCCCCTGGAGGAGACCACG
CTGGAGGACTTCAACGAGGTGTTCAATCTCAACACGCGGGCGCAGATGTTCGTCGCGCGC
TACGCCTACGACCACATCCAGCCGGGGGGGCGGGTGATCCTCATGTCCTCAATTGCCGCG
GGGCTGGGGGTGCCGGGCCACGCGCTGTACGCGGGCAGCAAGGCGGCCATCGAGGGGTTC
ACGCGCTGCTTGGCGGCCGATTTTGGGCGCAAGGGGTGCACGGTCAATGCGATCGCGCCG
GCGGGGGTCAAGAGCGACATGTGGCGGGAGAACGCATGGCGGTACGCTCCCGGCTGCGAC
AAGAGCTCGTCGCTGGAGGAGATCGAGACGGCGCTGGCGAGCGGGAGTCCCTTGAAGCGG
TGTGGTGTGCCCGAGGACATTGGCAAGGTGGTTTCATTCCTGGCGAGTCCGGATGCGGAA
TGGGTGAACGGTGAGTTCTTCTCCCCCCCTCCATGCAAGGGTCCCCTCCCCGGTGGGGAT
TTGGAATGCTGATCATTGCGCATCTGAACGAACAGGCCAGATCATCCCTGTCAATGGAGG
CGCAAACATCTGATTGAGACTGCGTGTACGGAGTGTATGGTTGCATTGAGTTGATTTGGA
ATTGTCTCGCCGAGTTTGTTCGGTCTGCTGTATCCTGATGATGTTTGACTTTCTTCAATA
CTACTGCGCCTAAAGTATAATCTTGCGGAGGCAGGGCAGCTAGGATAAATACGCCCATCA
GAATCCAATTCCCTTCATCCTGGAACTTTCTTGGTTGCTGTCGTGTGTGGCAACGTCAAT
TGACTTTGTATATGGCCTGATGCACGACGCTATGAGGGTTGCGCTCGAAGAGACCGATCC
ATTTCCAACTCATCAGACCTGAGTTCCAGAATGAATACCTCTATGGAGGCTGTCATGCCA
TGCTGCTCGGGGCAGATCACCATGTATGAGTCAAACATTGCACATCCTCTATGTCATGGA
CTCG
'''
gene = gene.replace('\n', '')
# print(len(gene))
# print(gene)
fragments = shotgun(gene)

# res1 = find_shortest_superstring(fragments) # muito lento
# res2 = dijkSuperstring(fragments) # muito lento
res3 = merge(fragments) # rapido

# print(res1, len(res1))
# print(res2, len(res2))
# print(res3, len(res3))

minimum_penalty(gene, res3, 1, 1)
# referencias
# 
# https://fileadmin.cs.lth.se/cs/Personal/Andrzej_Lingas/superstring.pdf
# https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_scs.pdf
# https://www.researchgate.net/profile/Sunil_Kumar_Sahu/post/Does_anyone_have_an_e-book_about_bioinformatics_How_can_i_find_it/attachment/59d6234ac49f478072e995f8/AS%3A272115389927425%401441888772501/download/Bioinformatics+Sequence+and+Genome+Analysis+-+David+W.+Mount.pdf
# https://www8.cs.umu.se/kurser/TDBAfl/VT06/algorithms/BOOK/BOOK5/NODE209.HTM
# https://math.mit.edu/~goemans/18434S06/superstring-lele.pdf
# https://www.geeksforgeeks.org/shortest-superstring-problem-set-2-using-set-cover/
#
# http://sunny.today/generate-random-integers-with-fixed-sum/
# http://www.vision.ime.usp.br/~jmena/mgwt/datasets/