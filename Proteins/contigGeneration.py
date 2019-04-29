from collections import defaultdict


def getInput(file):
    """
    Read input file and extract relevant variables. Adjustable mode for different functions.

    :param file: string with file name
    :return:
    """

    f = open(file, "r")
    lines = f.readlines()
    lines = [seq.strip() for seq in lines[:]]
    f.close()
    return lines


def contigGeneration(patterns):
    """

    :param patterns: List of k-mer patterns.
    :return: List of all contigs derived from the patterns deBruijn graph.
    """
    graph = deBruijnFromKmers(patterns)
    print(graph)
    result = maximalNonBranchingPaths(graph)
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


def maximalNonBranchingPaths(dB):
    # Initialise empty lists and dictionary
    paths = []
    cycles = []
    node_in_out = {}
    values = dB.values()
    keys = list(dB.keys())

    # Populate dictionary with in and out edges.
    for x in dB.keys():
        node_in_out[x] = [sum(v == x for value in values for v in value), len(dB.get(x))]

    # Iterate through each node in graph.
    for node in node_in_out:
        # Continue if node is not a 1-in-1-out node & out edge is more than 0.
        if node_in_out[node] != [1, 1] and node_in_out[node][1] > 0:
            # For each out_edge from node
            for out_edge in dB[node]:
                position = out_edge
                path = [node, out_edge]

                while position in keys and node_in_out[position] == [1, 1]:
                    path.append(dB[position][0])
                    position = dB[position]
                    print(position)
                cycles.append(path)
    print(cycles)
    # Generate paths for each isolated cycle
    for isolated_cycle in cycles:
        paths.append(stringSpelled(isolated_cycle))

    return paths


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


# Test for script
kmers = getInput('ctest1.txt')
contig_lst = contigGeneration(kmers)

for contig in contig_lst:
    print(contig)
