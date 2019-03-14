from collections import defaultdict
import random

def getInput(file):
    """
    Read input file and extract relevant variables. Adjustable mode for different functions.

    :param file: string with file name
    :return: pairComp as list of tuples.
    """

    f = open(file, "r")
    lines = f.readlines()
    numbers = [int(numb) for numb in lines[0].split()]
    lines = [seq.strip() for seq in lines[1:]]
    k = numbers[0]
    d = numbers[1]
    f.close()
    pairComp = [tuple(line.split('|')) for line in lines]
    return pairComp, k, d



def pairComposition(sequence, k, d):
    """
    Obtains pairs of k-mers separated by distance.

    :param sequence: String
    :param k: Int of k-mer size.
    :param d: Int of distance between k-mers.
    :return: List of tuples that have k-mer pairs.
    """
    result = []
    for i in range(len(sequence) - k - k - d + 1):
        result.append((sequence[i:i+k], sequence[i+k+d:i+k+k+d]))
    return sorted(result)


def stringSpelled(patterns):
    """
    Reconstruct a string from its genome path.

    :param patterns: List of path k-mers such that k-1 symbols of pattern[i] are equal.
    :return: String of k-mer path.
    """

    genome = ''
    for i in range(len(patterns)):
        genome += patterns[i][:1]
    genome += patterns[-1][1:]
    return genome


def stringSpelledPair(pathPairPatterns, k, d):
    """
    Reconstruct a string from its genome path of pair reads.

    :param pathPairPatterns: List of tuples of k-mer pairs that are in eulerian order.
    :param k: Int of k-mer size.
    :param d: Int of distance between k-mers.
    :return: String of k-mer pair path.
    """
    firstPatterns = [str(pair[0]) for pair in pathPairPatterns]
    secondPatterns = [str(pair[1]) for pair in pathPairPatterns]

    prefixString = stringSpelled(firstPatterns)
    suffixString = stringSpelled(secondPatterns)


    for i in range((k+d+1),len(prefixString)):
        if prefixString[i] != suffixString[i-k-d]:
            return "There is no string spelled by the gap patterns"

    return prefixString + suffixString[-k-d:]


def deBruijnFromPairs(pairComp):
    """

    :param pairComp: List of tuples containing pair compositions.
    :return: deBruijin graph of text and k-mer in adjacency list form.
    """
    result = defaultdict(list)
    prefix = [(pair[0][:-1], pair[1][:-1]) for pair in pairComp]
    suffix = [(pair[0][1:], pair[1][1:]) for pair in pairComp]

    for left, right in zip(prefix, suffix):
        if left[0][1:] == right[0][:-1] and left[1][1:] == right[1][:-1]:
            result[left].append(right)

    return result


def eulerianPath(graph):
    """
    Finds the path through a graph where each in and out edge or 'bridge' is used only once.

    :param graph: Directed adjacency list of Eulerian graph as dictionary.
    :return: List in order representing Eulerian path
    """
    # Initialize empty stack, circuit and in/out lists.
    stack = []
    circuit = []
    in_edges = []
    out_edges = []

    # Obtain length count of all values in graph dictionary and list version of keys.
    values = graph.values()
    count = sum(len(value) for value in values)
    listed_keys = list(graph.keys())

    # Populate individual check lists.
    for x in graph.keys():
        out_edges.append(len(graph.get(x)))
        in_edges.append(sum(v == x for value in values for v in value))

    # Populate compared check list.
    check = [x - y for x, y in zip(out_edges, in_edges)]

    # Check if Eulerian path is available and select starting location based on in-out edges.
    if in_edges == out_edges:
        location = random.choice(listed_keys)
        print('All in and out edges of vectors are equal.')

    elif sorted(set(check)) == [0, 1] or sorted(set(check)) == [-1, 0, 1]:
        location = listed_keys[check.index(1)]
        print('All but two vectors in and out edges are equal')

    else:
        raise ValueError('No Eulerian path available in this graph!')

    # While loop until length of circuit is of equal length to total length of values in graph.values()
    while len(circuit) != count:
        # Current key iteration
        key = graph.get(location)

        # Runs if key is empty.
        if not key:
            circuit.append(location)
            location = stack.pop()

        # Runs if it reaches a key-node that still has out-edges left.
        if key:
            # Append location to stack
            stack.append(location)
            # Select index and pop out of key-value list element.
            index = random.randint(0, len(key) - 1)
            select = key.pop(index)
            # Reset location to selection
            location = select

    # Finalise circuit path.
    circuit.append(location)

    # Return reversed circuit for correct order.
    return circuit[::-1]


def stringReconstructFromPairs(pairs, k, d):
    """
    Uses several of the above functions to obtain the eulerian string constructions from read pairs.

    :param pairs: List of tuples consiting of pairs
    :param k: Int of k-mer size.
    :param d: Int of distance between k-mers.
    :return: String of constructed sequence.
    """
    dB = deBruijnFromPairs(pairs)
    path = eulerianPath(dB)
    result = stringSpelledPair(path, k, d)
    return result


# Test for input and stringSpelled functions
sample, k, d = getInput('pairsample.txt')
print(stringSpelledPair(sample, k, d))


# Test for deBruijnFromPairs
dB = deBruijnFromPairs(sample)
print(dB)

# Test for eulerianPath
path = eulerianPath(dB)
print(path)

# Test for stringSpelled using path
result = stringSpelledPair(path, k, d)
assert result == 'GTGGTCGTGAGATGTTGA'
print(result)


# Large dataset test
test, k, d = getInput('dataset_204_16.txt')
print(stringReconstructFromPairs(test, k, d))