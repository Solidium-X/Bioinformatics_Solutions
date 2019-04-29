import math
import timeit

def pattCount(seq, patt):
    t_count = 0
    for i in range(len(seq)-len(patt)):
        temp = seq[i:i+len(patt)]
        if temp == patt:
            t_count += 1
    return t_count


def pattPos(seq, patt):
    position = []
    for i in range(len(seq)-len(patt)):
        temp = seq[i:i+len(patt)]
        if temp == patt:
            position.append(i)
    position = " ".join(map(str, position))
    return position


def frequentWords(text, k):
    freqpatt = []
    count = []
    for i in range(len(text)-k):
        pattern = text[i:i+k]
        temp = pattCount(text, pattern)
        count.append(temp)
    max_count = max(count)
    for i in range(len(text)-k):
        if count[i] == max_count:
            freqpatt.append(text[i:i+k])
            freqpatt = list(set(freqpatt))
            str1 = ' '.join(freqpatt)
    return str1


def fastFrequentWords(text, k):
    freq_patt = []
    freq_array = computingFreq(text, k)
    max_count = max(freq_array)
    for i in range(int(math.pow(4, k))):
        if freq_array[i] == max_count:
            pattern = numberToPattern(i, k)
            freq_patt.append(pattern)
    return freq_patt


def patternToNumber(seq):
    result = 0
    nucleo = ["A", "C", "G", "T"]
    for i in range(len(seq)):
        expo = len(seq) - 1 - i
        result += (nucleo.index(seq[i]) * math.pow(4, expo))
    return int(result)


def numberToPattern(number, k):
    nucleo = ["A", "C", "G", "T"]
    result = ""
    tempnum = number

    for i in range(k):
        remainder = int(tempnum % 4)
        tempnum = int(tempnum/4)
        result += nucleo[remainder]

    return result[::-1]

def computingFreq(text, k):
    freqArray = []
    for i in range(int(math.pow(4, k))):
        freqArray.append(0)

    for i in range(len(text) - k + 1):
        value = text[i:i+k]
        num = patternToNumber(value)
        freqArray[num] += 1
    #return " ".join(map(str, freqArray))
    return freqArray



def clumpFinding(genome, k, L, t):
    freq_patt = []
    clump = []

    for i in range(int(math.pow(4, k))):
        clump.append(0)
    for i in range(len(genome) - L + 1):
        text = genome[i:i+L]
        freq_array = computingFreq(text, k)
        for index in range(int(math.pow(4, k))):
            if freq_array[index] >= t:
                clump[index] = 1

    for i in range(int(math.pow(4, k))):
        if clump[i] == 1:
            pattern = numberToPattern(i, k)
            freq_patt.append(pattern)

    return freq_patt


def betterClumpFinding(genome, k, L, t):
    freq_patt = []
    clump = []

    for i in range(int(math.pow(4, k))):
        clump.append(0)
    text = genome[:L]
    freq_array = computingFreq(text, k)
    for i in range(int(math.pow(4, k))):
        if freq_array[i] >= t:
            clump[i] = 1
    for i in range(len(genome) - L + 1):
        first_patt = genome[i:i+k]
        index = patternToNumber(first_patt)
        freq_array[index] -= 1
        last_patt = genome[i + L - k: i + L]
        index = patternToNumber(last_patt)
        freq_array[index] += 1

        if freq_array[index] >= t:
            clump[index] = 1

    for i in range(int(math.pow(4, k))):
        if clump[i] == 1:
            pattern = numberToPattern(i, k)
            freq_patt.append(pattern)

    return freq_patt


