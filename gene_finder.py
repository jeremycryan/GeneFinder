# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: JEREMY RYAN

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
    # Replaces character "nucleotide" with corresplonding complement
    if nucleotide == "T":
        return "A"
    if nucleotide == "A":
        return "T"
    if nucleotide == "G":
        return "C"
    if nucleotide == "C":
        return "G"
    # TODO: implement this
    pass


def get_reverse_complement(dna):
    reverse = "".join([get_complement(i) for i in dna])
    # reverse = dna[::-1]
    # reverse = reverse.lower()
    # #   Replaces all characters in dna with lowercase so that replace function
    # #   does not interfere with itself.
    # reverse = reverse.replace("a", get_complement("A"))
    # reverse = reverse.replace("t", get_complement("T"))
    # reverse = reverse.replace("c", get_complement("C"))
    # reverse = reverse.replace("g", get_complement("G"))
    return reverse[::-1]
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    pass


def rest_of_ORF(dna):
    index = 0
    for letter in dna:
        if index%3 == 0 and letter == "T":  #   Scans for frames divisible by 3 that begin with T
            frame = dna[index:index + 3]
            if frame == "TAG" or frame == "TAA" or frame == "TGA":
                return dna[:index]
                #   Returns dna up to the index of first in frame stop codon
        index = index + 1
    return dna
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
    Stop codons include TAG, TAA, and TGA.
    """
    # TODO: implement this
    pass


def find_all_ORFs_oneframe(dna):
    index = 0
    listnum = 0
    list1 = []
    for letter in dna:
        if dna[index:index + 3] == "ATG":
            restofdna = dna[index:]
            orf = (rest_of_ORF(restofdna))
            list1.append(orf)
            listnum = listnum + 1
            index = index + len(orf)
        index = index + 3   #   Searches only every 3 characters
    return list1
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
    """
    # TODO: implement this
    pass


def find_all_ORFs(dna):
    list1 = []
    frame = 0
    while frame < 3:
        frameorfs = find_all_ORFs_oneframe(dna[frame:])
        list1.append(frameorfs)
        frame = frame + 1
    list1 = list1[0] + list1[1] + list1[2]
    return list1

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
    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    orfs = find_all_ORFs(dna)
    revdna = get_reverse_complement(dna)
    revorfs = find_all_ORFs(revdna)
    return orfs + revorfs

    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass


def longest_ORF(dna):
    dnalist = find_all_ORFs_both_strands(dna)
    greatestlen = 0
    greatestorf = "N/A"
    for orf in dnalist:
        if len(orf) > greatestlen:
            greatestlen = len(orf)
            greatestorf = orf
    return greatestorf
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    dnalist = []
    greatestlen = 0
    for char in dna:
        dnalist.append(char)    #   Convert string to list to randomize
    trial = 0
    while trial < num_trials:
        random.shuffle(dnalist)     #   Randomize!
        dnaliststr = "".join(dnalist)
        strlength = len(longest_ORF(dnaliststr))
        if strlength > greatestlen:
            greatestlen = strlength
            greatestorf = longest_ORF(dnaliststr)
            print(greatestlen)
        trial = trial + 1
    return greatestorf
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    index = 0
    aminochain = []
    while index < len(dna) - 2:
        frame = dna[index:index + 3]
        aminochain.append(aa_table[frame])
        index = index + 3
    amino = "".join(aminochain)
    return aminochain
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents a protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    orf = longest_ORF(dna)
    return coding_strand_to_AA(orf)
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
