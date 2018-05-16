import re

def match_length(cigar):
    matches = re.findall("\d+M", cigar)
    return sum([int(x[0:-1]) for x in matches])

def overlap(cigar, readstart, start, end):
    matches = re.findall("\d+[MNID]", cigar)
    matches = [[int(x[0:-1]), x[-1]] for x in matches]
    regions = []
    overlaps = 0
    position = int(readstart)
    for x in matches:
        if x[1] == "M":
            regions.append([position, position+x[0]])
            for x in range(position, position+x[0]):
                if x >=int(start) and x<=int(end):
                    overlaps += 1
        else:
            position += x[0]

    return overlaps