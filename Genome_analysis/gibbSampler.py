
import random
from functools import reduce

def getInput(file):
    """
    Read input file and extract relevant variables.

    :param file: string with file name
    :return: DNA as list of strings. k, t and N as int.
    """
    f = open(file, "r")
    lines = f.readlines()
    elements = lines[0].split()
    k = int(elements[0])
    t = int(elements[1])
    N = int(elements[2])
    DNA = [seq.strip() for seq in lines[1:]]
    f.close()

    return DNA, k, t, N


def countMotif(motifs, pseudo=False):
    """
    Count occurrence of nucleotides in column of motif matrix.

    :param motifs: List of motifs as strings.
    :param psuedo: If true, adds psuedocounts of 1 to each count.
    :return: Count of each column in A, C, G, T format.
    """
    result = []
    for i in range(len(motifs[0])):
        column = [row[i] for row in motifs]
        a_count = column.count("A")
        c_count = column.count("C")
        g_count = column.count("G")
        t_count = column.count("T")

        if pseudo is True:
            column_counts = [a_count + 1, c_count + 1, g_count + 1, t_count + 1]
        else:
            column_counts = [a_count, c_count, g_count, t_count]
        result.append(column_counts)
    return result


def profile(motifs, pseudo=False):
    """
    Produces a profile of probabilities for each nucleotide per column.

    :param motifs: List of DNA strings/motifs.
    :param pseudo: If True, will utilise Pseudocounts, if false/default then normal count.
    :return: Dictionary with probabilities of nucleotide per column.
    """
    if pseudo is True:
        counts = countMotif(motifs, pseudo=True)
    else:
        counts = countMotif(motifs)
    result = {'A': [], 'C': [], 'G': [], 'T': []}
    for column in counts:
        percent = [element / sum(column) for element in column]
        result['A'].append(percent[0])
        result['C'].append(percent[1])
        result['G'].append(percent[2])
        result['T'].append(percent[3])
    return result


def scoreMotif(motifs):
    """
    Score motifs and returns an integer equal to the amount of nucleotide differences between the motifs.
    0 = Completely identical.

    :param motifs: List of motifs strings.
    :return: Integer equal to the nucleotide differences between all strings.
    """
    counts = countMotif(motifs)
    result = [sum(i)-max(i) for i in counts]
    return sum(result)


def computeProfile(consensus, profile):
    """
    Computes the probability of consensus string based on the given profile.

    :param consensus: Single string of nucleotide.
    :param profile: Dictionary with probabilities of nucleotide per column
    :return: Float equal to the probability of string ooccurringbased on the profile values.
    """
    # Function requires string consensus and dictionary with list as values as profile
    result = []

    # Loop through consensus string, obtain dictionary value key.
    for i in range(len(consensus)):
            probability = profile.get(consensus[i])
            result.append(probability[i])

    # Multiplies everything within result list
    return reduce(lambda x, y: x*y, result)


def randomMotif(DNA, k):
    """
    :param DNA: List of strings containing DNA.
    :param k: Integer corresponding to K-mer size.
    :return: List of randomly selected K-mers for each string in DNA.
    """
    return [randomKmer(seq, k) for seq in DNA]


def randomKmer(seq, k):
    """
    :param seq: Single string
    :param k: Integer corresponding to K-mer size.
    :return: Randomly selected slice of seq string of size K.
    """
    position = random.randint(0, len(seq)-k)
    return seq[position:position+k]


def randomSelect(probabilities):
    """
    Randomly select a K-mer by giving the corresponding index based on weighted probabilities.

    :param probabilities: List of probabilities.
    :return: Integer corresponding to index of selected K-mer.
    """
    total = float(sum(probabilities))
    incSum = 0.0
    rand = random.random()
    for i in range(len(probabilities)):
        incSum += probabilities[i]/total
        if incSum >= rand:
            return i


def gibbsSampler(DNA, k, t, N):
    """
    Primary function to find optimal motifs of list of DNA based on rolling-dice method.

    :param DNA: List of DNA strings/motifs.
    :param k: Integer corresponding to K-mer size.
    :param t: Integer corresponding to length of DNA list.
    :param N: Integer corresponding to number of iterations should occur.
    :return: Best motifs obtained from this occurrence of the function as a list as well as score as integer.
    """
    motifs = randomMotif(DNA, k)
    bestMotifs = motifs

    for iteration in range(N):
        randInt = random.randint(0, t-1)

        # Create profile excluding randomly selected string.
        minus_profile = profile(motifs[:randInt] + motifs[randInt+1:], pseudo=True)

        # Calculate probabilities of all K-mers in randomly selected string.
        prob_lst = [computeProfile(DNA[randInt][i:i+k], minus_profile) for i in range(len(DNA[randInt])-k+1)]

        # Randomly select K-mer based on weighted probabilities.
        select = randomSelect(prob_lst)

        # Replace randomly selected motif with new motif.
        motifs[randInt] = DNA[randInt][select:select+k]

        if scoreMotif(motifs) < scoreMotif(bestMotifs):
            bestMotifs = motifs[:]

    return scoreMotif(bestMotifs), bestMotifs


def repeatGibbsSampler(DNA, k, t, N, R):
    """
    Repeats the primary gibbsSampler by R amount of times. Keeps the best motifs and scores of all repeats.
    This allows new random starting points to occur, increasing chances of converging towards best motifs.

    :param DNA: List of DNA strings/motifs.
    :param k: Integer corresponding to K-mer size.
    :param t: Integer corresponding to length of DNA list.
    :param N: Integer corresponding to number of iterations should occur.
    :param R: Integer corresponding to number of repeats that gibbsSampler function should be repeated
    :return: List of strings consisting of best motifs obtained.
    """
    bestMotifs = randomMotif(DNA, k)
    bestScore = scoreMotif(bestMotifs)
    for i in range(R):
        print(i)
        cur_score, cur_motifs = gibbsSampler(DNA, k, t, N)
        if cur_score < bestScore:
            bestMotifs = cur_motifs
            bestScore = cur_score
    return bestMotifs


t