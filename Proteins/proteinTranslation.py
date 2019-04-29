
# Table mapping mRNA to protein. '*' corresponds to STOP peptide.
RNAtoPeptideTable = {"AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N", "ACA": "T", "ACC": "T",
                     "ACG": "T", "ACU": "T", "AGA": "R", "AGC": "S", "AGG": "R", "AGU": "S",
                     "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I", "CAA": "Q", "CAC": "H",
                     "CAG": "Q", "CAU": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
                     "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R", "CUA": "L", "CUC": "L",
                     "CUG": "L", "CUU": "L", "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D",
                     "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A", "GGA": "G", "GGC": "G",
                     "GGG": "G", "GGU": "G", "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
                     "UAA": "*", "UAC": "Y", "UAG": "*", "UAU": "Y", "UCA": "S", "UCC": "S",
                     "UCG": "S", "UCU": "S", "UGA": "*", "UGC": "C", "UGG": "W", "UGU": "C",
                     "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F"}


def splitStringLength(string, n):
    """
    Split string into blocks of size n.

    :param string: Single string.
    :param n: Integer, size of block
    :return: List containing splits of size n.
    """
    return [string[i:i + n] for i in range(0, len(string), n)]


def proteinTranslate(sequence):
    """

    :param sequence: String of RNA sequence
    :return: Peptide string of protein translation.
    """
    split = splitStringLength(sequence, 3)
    ignore_stop = ["UAA", "UGA", "UAG"]
    translation = []

    for codon in split:
        if codon not in ignore_stop:
            translation.append(RNAtoPeptideTable.get(codon))

    return ''.join(translation)


def dnaComplement(sequence):
    # Used to obtain complement for DNA string.
    dna_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'}
    return "".join([dna_complement[base] for base in sequence])

def rnaTranscription(sequence):
    # Used to turn 'T' to 'U' in sequence
    rna_transcript = {'T': 'U', 'A': 'A', 'C': 'C', 'G': 'G'}
    return "".join([rna_transcript[base] for base in sequence])


def RNAtoDNA(sequence):
    # Used to turn 'U' to 'T' in sequence
    translation = {'U': 'T', 'A': 'A', 'C': 'C', 'G': 'G'}
    return "".join([translation[base] for base in sequence])


def rotations(peptide):
    # Obtain all circular rotations of peptide string.
    rotates = [peptide[i:] + peptide[:i] for i in range(len(peptide))]
    return rotates


def peptideEncoding(DNA, peptide):
    """

    :param DNA: String of DNA sequence.
    :param peptide: String of peptide sequence.
    :return: List of all DNA sub-strings which translates to peptide string.
    """
    # Obtain total codon size of the peptide chain.
    codon_len = len(peptide) * 3

    # Transcribing coding strand
    transcript = rnaTranscription(DNA)

    # Retrieve reverse complement of DNA and transcribe
    reverse_complement = dnaComplement(DNA[::-1])
    reverse_transcript = rnaTranscription(reverse_complement)

    # Split the sequences relative to the codon size of the peptide chain, consolidate both lists.
    splits = [transcript[i:i+codon_len] for i in range(len(transcript)-codon_len+1)]
    reverse_splits = [reverse_transcript[i:i+codon_len] for i in range(len(reverse_transcript)-codon_len+1)]

    # Translation of total splits.
    trans_splits = [proteinTranslate(split) for split in splits]
    reverse_trans_splits = [proteinTranslate(split) for split in reverse_splits]

    # Obtain indexes if peptide matches a split
    locations = [index for index, tran in enumerate(trans_splits) if peptide in tran]
    reverse_locations = [index for index, tran in enumerate(reverse_trans_splits) if peptide in tran]

    # Change 'U' nucleotides to 'T'
    codons = [RNAtoDNA(splits[index]) for index in locations]

    # reverse_codons needs to be reconverted using dnaComplement method.
    reverse_codons = [dnaComplement(reverse_splits[index][::-1]) for index in reverse_locations]

    # Consolidate results
    result = codons + reverse_codons

    return result


