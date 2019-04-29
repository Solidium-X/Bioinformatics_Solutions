from Week_2.skew import hammingDistance, neighbors
from Week_1.count_substring import numberToPattern
from collections import Counter
import itertools
import timeit
import math
from functools import reduce


def window(seq, k):
    for i in range(1 + len(seq) - k):
        yield seq[i:i + k]


def motifEnumeration(DNA, k, d):
    patterns = set()
    temp_motifs = []
    for i in DNA:
        for j in range(len(i)-k):
            temp_patt = i[j:j+k]
            neigh = neighbors(temp_patt, d)
            temp_motifs.append(neigh)

    for i in temp_motifs:
        for combo in i:
            if all(any(hammingDistance(combo, pat) <= d for pat in window(string, k)) for string in DNA):
                patterns.add(combo)

    return list(patterns)


def countMotif(motifs):
    result = []
    for i in range(len(motifs[0])):
        column = [row[i] for row in motifs]
        a_count = column.count("A")
        c_count = column.count("C")
        g_count = column.count("G")
        t_count = column.count("T")

        column_counts = [a_count, c_count, g_count, t_count]
        result.append(column_counts)
    return result


def countMotifPseudo(motifs):
    result = []
    for i in range(len(motifs[0])):
        column = [row[i] for row in motifs]
        a_count = column.count("A")
        c_count = column.count("C")
        g_count = column.count("G")
        t_count = column.count("T")

        column_counts = [a_count + 1, c_count + 1, g_count + 1, t_count + 1]
        result.append(column_counts)
    return result


def scoreMotif(motifs):
    counts = countMotif(motifs)
    result = [sum(i)-max(i) for i in counts]
    return sum(result)


def consensus(motifs):
    result = []
    for i in range(len(motifs[0])):
        column = [row[i] for row in motifs]
        column_con = max(column, key=column.count)
        result.append(column_con)
    return "".join(result)


def profile(motifs):
    counts = countMotif(motifs)
    result = {'A': [], 'C': [], 'G': [], 'T': []}
    for column in counts:
        percent = [element / sum(column) for element in column]
        result['A'].append(percent[0])
        result['C'].append(percent[1])
        result['G'].append(percent[2])
        result['T'].append(percent[3])
    return result


def profilePseudo(motifs):
    counts = countMotifPseudo(motifs)
    result = {'A': [], 'C': [], 'G': [], 'T': []}
    for column in counts:
        percent = [element / sum(column) for element in column]
        result['A'].append(percent[0])
        result['C'].append(percent[1])
        result['G'].append(percent[2])
        result['T'].append(percent[3])
    return result


def fastScoreMotif(DNA):
    con = consensus(DNA)
    score = 0
    for i in DNA:
        score += hammingDistance(con, i)
    return score


def entropy(lst):
    result = 0
    for i in lst:
        if i != 0:
            result += i * math.log2(i)
    return result


def motifEntropy(motifs):
    prof = profile(motifs)
    result = []
    for i, j, k, d in zip(prof['A'], prof['C'], prof['G'], prof['T']):
        result.append(entropy([i, j, k, d]))
    return abs(sum(result))


def patStringDistance(Pattern, DNA):
    k = len(Pattern)
    distance = 0
    for i in DNA:
        hd = len(i)

        for j in range(len(i)-k+1):
            kmer = i[j:j+k]
            temp_hd = hammingDistance(Pattern, kmer)
            if temp_hd < hd:
                hd = temp_hd

        distance += hd

    return distance


def medianString(DNA, k):
    distance = k * len(DNA)
    for i in range(int(math.pow(4, k))):
        pattern = numberToPattern(i, k)
        if distance > patStringDistance(pattern, DNA):
            distance = patStringDistance(pattern, DNA)

            median = pattern

    return median


def computeProfile(consensus, profile):
    # Function requires string consensus and dictionary with list as values as profile
    result = []

    # Loop through consensus string, obtain dictionary value key.
    for i in range(len(consensus)):
            probability = profile.get(consensus[i])
            result.append(probability[i])

    # Multiplies everything within result list
    return reduce(lambda x, y: x*y, result)


def profileMostProbable(string, k, profile):
    """
    Function utilises computeProfile to determine probability of all k-mers of string
    and return the most probable k-mer.

    :param string: Takes in genome data as single string.
    :param k: integer for k-mer size.
    :param profile: Takes in matrix of probabilities as dictionary.
    :return: String of most probable k-mer.
    """

    result = {}
    for i in range(len(string) - k + 1):
        kmer = string[i:i+k]
        prob = computeProfile(kmer, profile)
        result[kmer] = prob
    return max(result, key=result.get)


def greedyMotif(DNA, k, t):
    bestMotifs = [string[:k] for string in DNA]
    for i in range(len(DNA[0])-k+1):
        first = DNA[0]
        motif = [first[i:i + k]]
        for string in DNA[1:t]:
            motif_profile = profile(motif)
            result = profileMostProbable(string, k, motif_profile)
            motif.append(result)

        if scoreMotif(motif) < scoreMotif(bestMotifs):
            bestMotifs = motif
    return bestMotifs


def greedyMotifPseudo(DNA, k, t):
    bestMotifs = [string[:k] for string in DNA]
    for i in range(len(DNA[0])-k+1):
        first = DNA[0]
        motif = [first[i:i + k]]
        for string in DNA[1:t]:
            motif_profile = profilePseudo(motif)
            result = profileMostProbable(string, k, motif_profile)
            motif.append(result)

        if scoreMotif(motif) < scoreMotif(bestMotifs):
            bestMotifs = motif
    return bestMotifs
