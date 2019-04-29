

def skew(genome):
    genome = str(genome).upper()
    skewcount = 0
    skewlist = [0]
    for i in range(len(genome)):
        if genome[i] == "C":
            skewcount -= 1
        elif genome[i] == "G":
            skewcount += 1
        skewlist.append(skewcount)

    # returns list in string format, separated by spaces
    return " ".join(map(str, skewlist))


def minSkew(genome):
    # reconvert string to list of int
    temp_list = skew(genome)
    temp_list = [int(x) for x in temp_list.split(' ')]

    result = []
    min_value = min(temp_list)

    for i in range(len(temp_list)):
        if temp_list[i] == min_value:
            result.append(i)

    # returns list in string format, separated by spaces
    return " ".join(map(str, result))


def hammingDistance(p, q):
    count = 0
    if len(p) != len(q):
        raise Exception('Cannot compare sequences of different lengths.')
    if len(p) == len(q):

        for i, j in zip(p, q):
            if i == 'A' and j != 'A' \
            or i == 'T' and j != 'T' \
            or i == 'C' and j != 'C' \
            or i == 'G' and j != 'G':
                count += 1
    return count


def aproxPattPos(seq, patt, d):
    position = []
    for i in range(len(seq)-len(patt)+1):
        temp = seq[i:i+len(patt)]
        if hammingDistance(temp, patt) <= int(d):
            position.append(i)
    position = " ".join(map(str, position))
    return position


def aproxPattCount(seq, patt, d):
    t_count = 0
    for i in range(len(seq)-len(patt)+1):
        temp = seq[i:i+len(patt)]
        if hammingDistance(temp, patt) <= int(d):
            t_count += 1
    return t_count


def immediateNeighbors(pattern):
    neighborhood = []
    for i in range(len(pattern)):
        symbol = pattern[i]
        for x in range(len(pattern)):
            if x != symbol:
                neighbor = symbol - x
                neighborhood.append(neighbor)
    return neighborhood
            

def neighbors(pattern, d):
    char = 'ACGT'
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ["A", "C", "G", "T"]

    neighborhood = []
    firstsymbol = pattern[:1]
    suffix = pattern[1:]

    suffixNeighbors = neighbors(suffix, d)

    for text in suffixNeighbors:
        if hammingDistance(text, suffix) < d:
            for x in char:
                neighborhood.append(x+text)
        else:
            neighborhood.append(firstsymbol+text)
    return neighborhood

