import re
import random


def getInput(file, mode=None):
    """
    Read input file and extract relevant variables. Adjustable mode for different functions.

    :param file: string with file name
    :return: DNA as list of strings. k as int.
    """
    if mode in ['eulerian']:

        f = open(file, "r")
        lines = f.readlines()
        node_strip = [seq.strip() for seq in lines[:]]
        node_lst = [re.findall(r'\d+', row) for row in node_strip]
        node_dct = {node[0]:node[1:] for node in node_lst}

        f.close()
        return node_dct

    else:
        Exception('Please specify a mode!')


def eulerianCycle(graph):
    """
    Finds the cycle through a graph where each in and out edge or 'bridge' is used only once.

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

    elif sorted(set(check)) == [-1, 0, 1]:
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

    # Finalise circuit loop!
    circuit.append(circuit[0])

    # Return reversed circuit for correct order.
    return circuit[::-1]


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

