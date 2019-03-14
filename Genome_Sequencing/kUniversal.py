from collections import defaultdict
import itertools
import random
import timeit


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
    check = [x-y for x, y in zip(out_edges, in_edges)]

    # Check if Eulerian path is available and select starting location based on in-out edges.
    if in_edges == out_edges:
        location = random.choice(listed_keys)
        print('All in and out edges of vectors are equal.')

    elif sorted(set(check)) == [0, 1] or sorted(set(check)) == [-1,0,1]:
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

    # Finalise circuit path!
    circuit.append(location)

    # Return reversed circuit for correct order.
    return circuit[::-1]


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


def kUniversalCircular(k):
    """
    Function which runs the above functions to provide a string where a sliding k-mer window can only occur once,
    even as a cycle.

    :param k: Int, size of k-mer.
    :return: String, Binary universal cycle.
    """
    kmers = list(itertools.product('01', repeat=k))
    kmers = [''.join(x)for x in kmers]
    print('Kmers complete.')
    graph = deBruijnFromKmers(kmers)
    print('Graph complete.')
    cycle = eulerianPath(graph)
    print('Cycle complete.')
    result = stringSpelled(cycle[:-(k-1)])
    print('Result complete.')
    return result

# Test
start = timeit.default_timer()
print(kUniversalCircular(10))
eclipsed = timeit.default_timer() - start
print(f'Time for function: {eclipsed} seconds.')
