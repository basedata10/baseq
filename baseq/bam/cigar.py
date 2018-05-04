import re

def match_length(cigar):
    matches = re.findall("\d+M", cigar)
    return sum([int(x[0:-1]) for x in matches])