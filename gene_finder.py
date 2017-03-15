# -*- coding: utf-8 -*-
"""
Gene Finder

@author: Lakhvinder Jordan

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'
    else:
        print('Failure in get_complement. Input did not match options.')
        return


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    result = ''
    n = len(dna)
    for i in range(n):
            index = n-1-i
            x = get_complement(dna[index])
            result += x
    return result


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGTAAGTAAGA")
    'ATGAGTAAG'
    >>> rest_of_ORF("ATGAGTAAGTAAGACTAGATG")
    'ATGAGTAAG'
    """
    n = len(dna)/3
    n = int(n)
    for i in range(1, n):
        s = 3*i
        e = s+3
        if dna[s:e] == 'TGA' or dna[s:e] == 'TAA' or dna[s:e] == 'TAG':
            return dna[0:s]
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("AAAATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    ORFs = []
    count = 0
    n = 0
    k = len(dna)
    while k-n >= 3:
        while k-n >= 3:
            if dna[n:n+3] == 'ATG':
                break
            else:
                n = n + 3
        if not k-n >= 3:
            break
        ORFs.insert(count, rest_of_ORF(dna[n:k]))
        n = len(ORFs[count])+3+n
        count += 1
    return ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    list0 = find_all_ORFs_oneframe(dna)
    list1 = find_all_ORFs_oneframe(dna[1:len(dna)])
    list2 = find_all_ORFs_oneframe(dna[2:len(dna)])
    list3 = list0 + list1 + list2
    return list3


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    lista = find_all_ORFs(dna)
    listb = find_all_ORFs(get_reverse_complement(dna))
    listc = lista + listb
    return listc


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORFs = find_all_ORFs_both_strands(dna)
    return max(ORFs, key=len)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    results = []
    for i in range(num_trials):
        shuffled = shuffle_string(dna)
        longest = longest_ORF(shuffled)
        results.append(longest)
    return len(max(results, key=len))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    protein = ''
    n = len(dna)
    n += -(n % 3)
    n = int(n)
    for i in range(0, n, 3):
        amino_acid = aa_table[dna[i:i+3]]
        protein += amino_acid
    return protein


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    results = []
    threshold = longest_ORF_noncoding(dna, 1500)
    potentials = find_all_ORFs_both_strands(dna)
    for i in potentials:
        if len(i) > threshold:
            results.append(coding_strand_to_AA(i))
    return results

# if __name__ == "__main__":
    # import doctest
    # doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)


from load import load_seq
dna = load_seq("./data/X73525.fa")
print(gene_finder(dna))
