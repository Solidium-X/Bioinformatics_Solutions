

aminoMasses = {'G': '57', 'A': '71', 'S': '87',
               'P': '97', 'V': '99', 'T': '101',
               'C': '103', 'I': '113', 'L': '113',
               'N': '114', 'D': '115', 'K': '128',
               'Q': '128', 'E': '129', 'M': '131',
               'H': '137', 'F': '147', 'R': '156',
               'Y': '163', 'W': '186'}


def countSubPeptides(N):
    # Returns integer of subpeptides of a cyclic peptide of integer length N.
    return N * (N-1)


def rotations(peptide):
    # Obtain all circular rotations of peptide string.
    rotates = [peptide[i:] + peptide[:i] for i in range(len(peptide))]
    return rotates


def circularPeptideSubs(peptide):
    # return all of the peptide substrings within the circular peptide string.
    rotes = rotations(peptide)
    combos = []
    for rote in rotes:
        for i in range(len(rote)-1):
            combos.append(rote[:i+1])
    combos.append(peptide)

    return combos


def generateThoereticalSpectrum(peptide):
    """

    :param peptide: String of peptide
    :return: List with masses of all substrings of peptide.
    """
    all_masses = []
    combos = circularPeptideSubs(peptide)

    for subPeptide in combos:
        mass = 0
        for amino in subPeptide:
            mass += int(aminoMasses[amino])
        all_masses.append(mass)

    all_masses.append(0)

    return sorted(all_masses)

