import timeit
from collections import Counter

# Edited to remove L and Q as I/L and K/Q are considered the same.
aminoMasses = {'G': 57, 'A': 71, 'S': 87,
               'P': 97, 'V': 99, 'T': 101,
               'C': 103, 'I': 113, 'N': 114,
               'D': 115, 'K': 128, 'E': 129,
               'M': 131, 'H': 137, 'F': 147,
               'R': 156, 'Y': 163, 'W': 186}

# aminoMasses dictionary with key/values reversed
massesAmino = {v: k for k, v in aminoMasses.items()}

# list of amino and masses
aminos = list(aminoMasses.keys())
masses = list(aminoMasses.values())


def getInput(file):
    """
    Read input file and extract relevant variables.

    :param file: string with file name
    :return: List of int
    """

    f = open(file, "r")
    line = f.readline()
    line.strip()
    strings = line.split()
    ints = [int(i) for i in strings]
    f.close()

    return ints


def peptide_count_from_mass(m):
    # Calculate number of peptide products that can formulate given spectrometer (m)ass
    m += 1
    peptides = [0 for i in range(m)]
    for i in range(m):

        for mass in masses:
            if i - mass > 0:
                peptides[i] += peptides[i - mass]
            elif i - mass == 0:
                peptides[i] += 1

    return peptides[-1]


def subpeptide_count_from_length(n):
    # Calculate amount of sub peptides available from length of peptide.
    return int(n*(n+1)/2+1)


def rotations(peptide):
    # Obtain all circular rotations of peptide string.
    rotates = [peptide[i:] + peptide[:i] for i in range(len(peptide))]
    return rotates


def circularPeptideSubs(peptide):
    # return all of the peptide substrings within a circular peptide string.
    rotes = rotations(peptide)
    combos = []
    for rote in rotes:
        for i in range(len(rote)-1):
            combos.append(rote[:i+1])
    combos.append(peptide)

    return combos


def linearPeptideSubs(peptide):
    # return all of the peptide substrings within a linear peptide string.
    rotes = rotations(peptide)
    combos = []
    start = len(peptide)+1

    for rote in rotes:
        start -= 1
        for i in range(start):
            combos.append(rote[:i+1])

    return combos


def generateThoereticalSpectrum(peptide, linear=False):
    """
    Generate the theoretical mass spectrometry spectrum for all possible sub-peptides of
    given circular or linear peptide string.

    :param peptide: String of peptide
    :param linear: Default to circular peptide, unless selected to linear.
    :return: List with masses of all substrings of peptide.
    """
    all_masses = []

    if linear == True:
        combos = linearPeptideSubs(peptide)
    else:
        combos = circularPeptideSubs(peptide)

    for subPeptide in combos:
        mass = 0
        for amino in subPeptide:
            mass += amino
        all_masses.append(mass)

    all_masses.append(0)

    return sorted(all_masses)


def expand(peptide_lst, mass_rng):
    """
    Expand peptides with peptide mass_rng
    :param peptide_lst: List with list-elements to be expanded
    :param mass_rng: List of elements to be expanded by.
    :return: Expanded list of list.
    """
    expansion = []
    for peptides in peptide_lst:
        for mass in mass_rng:
            extend = peptides + [mass]
            expansion.append(extend)
    return expansion


def is_multiset(sub, target_list):
    """
    Check if subset of list, accounts for duplicates.
    :param sub: subset list.
    :param target_list: target list.
    :return: True or False.
    """
    c1, c2 = Counter(sub), Counter(target_list)
    return not c1-c2


def cycloPeptideSequencing(spectrum, strFormat = False):
    """

    :param spectrum: Experimental spectrum as list of int.
    :param strFormat: Default to False, if True it will return list in a print-friendly format.
    :return: List of all viable cyclo-peptides.
    """
    mass_pool = [mass for mass in spectrum if mass in masses]
    candidatePeptides = [[]]
    finalPeptides = []

    while candidatePeptides != []:
        good_peptides = []
        candidatePeptides = expand(candidatePeptides, mass_pool)

        for peptide in candidatePeptides:
            peptide_spectrum = generateThoereticalSpectrum(peptide, linear=True)
            if sum(peptide) >= spectrum[-1]:
                if is_multiset(peptide_spectrum, spectrum):
                    if peptide not in finalPeptides:
                        finalPeptides.append(peptide)

            elif is_multiset(peptide_spectrum, spectrum):
                good_peptides.append(peptide)

        candidatePeptides = good_peptides

    if strFormat == True:
        strings = []
        for i in finalPeptides:
            strings.append('-'.join(str(x) for x in i))
        return strings

    else:
        return finalPeptides


# Test script
start = timeit.default_timer()
spectrum = getInput('dataset_100_6.txt')
result = (cycloPeptideSequencing(spectrum, strFormat=True))
elapsed = timeit.default_timer()-start
print(f'Elapsed time for function in seconds: {elapsed}')
print(*result, sep=' ')


