import timeit

# Edited to remove L and Q as I/L and K/Q are considered the same.
aminoMasses = {'G': 57, 'A': 71, 'S': 87,
               'P': 97, 'V': 99, 'T': 101,
               'C': 103, 'I': 113, 'L': 113,
               'N': 114, 'D': 115, 'K': 128,
               'Q': 128, 'E': 129, 'M': 131,
               'H': 137, 'F': 147, 'R': 156,
               'Y': 163, 'W': 186}

masses = list(aminoMasses.values())


def getInput(file):
    """
    Read input file and extract relevant variables. Adjustable mode for different functions.

    :param file: string with file name
    :return: DNA as list of strings. k as int.
    """

    f = open(file, "r")
    lines = f.readlines()
    lines = [line.strip() for line in lines]
    N = int(lines[0])
    spectrum = lines[1].split(' ')
    spectrum = [int(i) for i in spectrum]
    f.close()
    return N, spectrum


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
            # Select one of the lines below if peptide is given as strings(1st line) or mass(2nd line)
            mass += aminoMasses[amino]
            #mass += amino
        all_masses.append(mass)

    all_masses.append(0)

    return sorted(all_masses)


def peptideScore(theoretical, experimental):
    count = 0
    # Create a new experimental list to avoid mutating the argument.
    temp_experimental = experimental.copy()

    for score in theoretical:
        if score in temp_experimental:
            count += 1
            temp_experimental.remove(score)
    return count


def trim(leaderboard, spectrum, N):
    """
    Trim list of peptides based on their score to given spectrum down to N amount.
    Includes ties on the N-th scores.

    :param leaderboard: List with peptide candidates.
    :param spectrum: List of experimental spectrum.
    :param N: Integer cut off-value
    :return: Trimmed list of qualifying peptides.
    """
    leaderScores = []
    qualified = []

    # Return leaderboard if can not trim by N due to insufficient length.
    if len(leaderboard) < N:
        return leaderboard

    # Generate scores for peptides.
    for candidate in leaderboard:
        score = peptideScore(generateThoereticalSpectrum(candidate, linear=True), spectrum)
        leaderScores.append((score, candidate))

    # Reverse rank order.
    rank = sorted(leaderScores, reverse=True)

    # Append to final list if N-th score matches.
    for i in range(len(rank)-1):
        if rank[i][0] >= rank[N-1][0]:
            qualified.append(rank[i][1])

    return qualified


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


def leaderboardCycloPeptideSequencing(spectrum, N, strFormat = False):
    """

    :param spectrum: Experimental spectrum as list of int.
    :param N: Int for leaderboard cut-off.
    :param strFormat: Default to False, if True it will return list in a print-friendly format.
    :return: List of all viable cyclo-peptides.
    """

    candidatePeptides = [[]]
    leaderPeptide = []
    leader_spectrum = []

    while candidatePeptides != []:
        good_peptides = []
        candidatePeptides = expand(candidatePeptides, masses)

        for peptide in candidatePeptides:
            sum_peptide = sum(peptide)
            peptide_spectrum = generateThoereticalSpectrum(peptide, linear=False)
            if sum_peptide == spectrum[-1]:
                if peptideScore(peptide_spectrum, spectrum) > peptideScore(leader_spectrum, spectrum):
                    leaderPeptide = peptide
                    leader_spectrum = generateThoereticalSpectrum(leaderPeptide, linear=False)
            elif sum_peptide < spectrum[-1]:
                good_peptides.append(peptide)

        candidatePeptides = trim(good_peptides, spectrum, N)

    if strFormat == True:
        strings = '-'.join(str(i) for i in leaderPeptide)
        return strings

    else:
        return leaderPeptide


# Test script
start = timeit.default_timer()
N, spectrum = getInput('dataset_102_8.txt')
result = leaderboardCycloPeptideSequencing(spectrum, N, strFormat=True)
elapsed = timeit.default_timer()-start
print(f'Elapsed time for function in seconds: {elapsed}')
print(result)
