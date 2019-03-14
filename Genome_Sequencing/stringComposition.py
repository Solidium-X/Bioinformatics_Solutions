from collections import defaultdict

def getInput(file, mode=None):
    """
    Read input file and extract relevant variables. Adjustable mode for different functions.

    :param file: string with file name
    :return: DNA as list of strings. k as int.
    """
    if mode in ['comp', 'deBruijn']:
        f = open(file, "r")
        lines = f.readlines()
        k = int(lines[0])
        DNA = [seq.strip() for seq in lines[1:]]
        f.close()
        return DNA, k

    elif mode in ['ss', 'overlap', 'dbk']:
        f = open(file, "r")
        lines = f.readlines()
        DNA = [seq.strip() for seq in lines]
        f.close()
        return DNA


    else:
        Exception('Please specify a mode!')



def composition(text, k):
    '''
    String Composition Problem: Generate the k-mer composition of a string.

    :param text: String of genome snippet
    :param k: Integer corresponding to k-mer size.
    :return: Composition-k(Text), where the k-mers are arranged in lexicographic order.
    '''

    k_list = [text[i:i+k] for i in range(len(text)-k+1)]
    k_list.sort()
    return k_list


def stringSpelled(patterns):
    '''
    Reconstruct a string from its genome path.

    :param patterns: List of path k-mers such that k-1 symbols of pattern[i] are equal
    :return:
    '''

    genome = ''
    for i in range(len(patterns)):
        genome += patterns[i][:1]
    genome += patterns[-1][1:]
    return genome


def overlapGraph(patterns, mode=None):
    '''

    :param patterns: List containing collection of patters of K-mers.
    :return: Overlapping graph. Mode sets return to dictionary or print friendly format.
    '''

    results = {}
    prefix_lst = [patt[:-1] for patt in patterns]

    for patt in patterns:
        suffix = patt[1:]
        # If suffix of pattern matches anything in prefix_lst, return list of indexes.
        indices = [i for i, prefix in enumerate(prefix_lst) if prefix == suffix]

        # If indices is not empty, create key and value list with no repeats.
        if indices:
            results[patt] = set([patterns[i] for i in indices])

    # Return a print-friendly or dictionary format.
    if mode in ['Print', 'print']:
      for key, value in results.items():
          print(key + ' ->', ', '.join(value))

    else:
        return results


def deBruijn(text, k, mode=None):
    '''

    :param text: Sting consisting of sequence
    :param k: Integer corresponding to size of k
    :param mode: Decide to print values in graph format or return dictionary (default).
    :return: deBruijin graph of text and k-mer in adjacency list form.
    '''
    result = defaultdict(list)
    edges = [text[i:i+k] for i in range(len(text)-k+1)]
    nodes = [(edge[:-1], edge[1:]) for edge in edges]

    for left, right in nodes:
        if left[1:] == right[:-1]:
            result[left].append(right)

    if mode in ['Print', 'print']:
        for key in sorted(result.keys()):
            print(key, "->", ', '.join(sorted(result[key])))

    else:
        return result


def deBruijnFromKmers(patterns, mode=None):
    '''

    :param patterns: List containing k-mer patterns.
    :param mode: Decide to print values in graph format or return dictionary (default).
    :return: deBruijin graph of text and k-mer in adjacency list form.
    '''
    result = defaultdict(list)
    nodes = [(edge[:-1], edge[1:]) for edge in patterns]

    for l, r in nodes:
        if l[1:] == r[:-1]:
            result[l].append(r)

    if mode in ['Print', 'print']:
        for key in sorted(result.keys()):
            print(key, "->", ', '.join(sorted(result[key])))

    else:
        return result


def kmers(string, k):
    '''

    :param string: Single string constituting sequence.
    :param k: Int corresponding to desired size of k-mer.
    :return: List of k-mers
    '''
    result = []
    for i in range(len(string)-k):
        result.append(string[i:i+k])
    if len(set(result)) == len(result):
        return True
    else:
        return False

