#!/usr/bin/env python3

'''
This script works on Python 3.6 or more.
'''

from struct import *
import re
import argparse
from decimal import Decimal

parser = argparse.ArgumentParser(description='convert the latest ascii psiblast checkpoint file into the pdb2vall checkpoint file.')
parser.add_argument("-p","--pdb", dest="pdb", help="Required. PDB+chain ID. i.e. 7odcA", required=True)

args = parser.parse_args()

Blosum62 = """0.0215
0.0023 0.0178
0.0019 0.0020 0.0141
0.0022 0.0016 0.0037 0.0213
0.0016 0.0004 0.0004 0.0004 0.0119
0.0019 0.0025 0.0015 0.0016 0.0003 0.0073
0.0030 0.0027 0.0022 0.0049 0.0004 0.0035 0.0161
0.0058 0.0017 0.0029 0.0025 0.0008 0.0014 0.0019 0.0378
0.0011 0.0012 0.0014 0.0010 0.0002 0.0010 0.0014 0.0010 0.0093
0.0032 0.0012 0.0010 0.0012 0.0011 0.0009 0.0012 0.0014 0.0006 0.0184
0.0044 0.0024 0.0014 0.0015 0.0016 0.0016 0.0020 0.0021 0.0010 0.0114 0.0371
0.0033 0.0062 0.0024 0.0024 0.0005 0.0031 0.0041 0.0025 0.0012 0.0016 0.0025 0.0161
0.0013 0.0008 0.0005 0.0005 0.0004 0.0007 0.0007 0.0007 0.0004 0.0025 0.0049 0.0009 0.0040
0.0016 0.0009 0.0008 0.0008 0.0005 0.0005 0.0009 0.0012 0.0008 0.0030 0.0054 0.0009 0.0012 0.0183
0.0022 0.0010 0.0009 0.0012 0.0004 0.0008 0.0014 0.0014 0.0005 0.0010 0.0014 0.0016 0.0004 0.0005 0.0191
0.0063 0.0023 0.0031 0.0028 0.0010 0.0019 0.0030 0.0038 0.0011 0.0017 0.0024 0.0031 0.0009 0.0012 0.0017 0.0126
0.0037 0.0018 0.0022 0.0019 0.0009 0.0014 0.0020 0.0022 0.0007 0.0027 0.0033 0.0023 0.0010 0.0012 0.0014 0.0047 0.0125
0.0004 0.0003 0.0002 0.0002 0.0001 0.0002 0.0003 0.0004 0.0002 0.0004 0.0007 0.0003 0.0002 0.0008 0.0001 0.0003 0.0003 0.0065
0.0013 0.0009 0.0007 0.0006 0.0003 0.0007 0.0009 0.0008 0.0015 0.0014 0.0022 0.0010 0.0006 0.0042 0.0005 0.0010 0.0009 0.0009 0.0102
0.0051 0.0016 0.0012 0.0013 0.0014 0.0012 0.0017 0.0018 0.0006 0.0120 0.0095 0.0019 0.0023 0.0026 0.0012 0.0024 0.0036 0.0004 0.0015 0.0196"""
aadict = {1:'A', 2:'B', 3:'C', 4:'D', 5:'E', 6:'F', 7:'G', 8:'H', 9:'I', 10:'K',
          11:'L', 12:'M', 13:'N', 14:'P', 15:'Q', 16:'R', 17:'S', 18:'T', 19:'V', 20:'W', 21:'X', 22:'Y'}


def parse_new_checkpoint_file(filename):
    nrow = ncol = 0
    aafilter = [-1,0,-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,-1,19,-1,-1,-1,-1,-1]

    f = open(filename,"r")
    line = f.readline().strip()
    seq = ""
    while line:
        if "seq-data ncbistdaa" in line:
            seq += line.split("'")[-1]
            line = f.readline().strip()
            while not "'" in line:
                seq += line
                line = f.readline().strip()
            seq += line.split("'")[0]

        if "numRows" in line:
            ld = re.split('[{,._-} ]', line)
            ld = [x for x in ld if len(x) > 0]
            nrow = int(ld[-1])
        if "numColumns" in line:
            ld = re.split('[{,._-} ]', line)
            ld = [x for x in ld if len(x) > 0]
            ncol = int(ld[-1])


        #if "weightedResFreqsPerPos" in line:
        if "freqRatios" in line:
            assert nrow >= 22, "numRows is insufficient"
            data = [[0.0]*20 for i in range(ncol)]

            line = f.readline().strip()
            for i in range(ncol):
                for j in range(nrow):
                    ld = re.split('[{,._-} ]', line)
                    ld = [float(x) for x in ld if len(x) > 0]
                    assert len(ld) == 3 , "premature pssm data"
                    val = ld[0] * ld[1] ** ld[2]
                    if aafilter[j] >= 0:
                        data[i][aafilter[j]] = float(val)
                    line = f.readline().strip()
        line = f.readline().strip()
    f.close()

    decode_seq = ""
    for i in range(0,len(seq),2):
        d = "0X" + seq[i:i+2]
        d = int(d,0)
        decode_seq += aadict[d]
    #print(decode_seq)
    return decode_seq, data

def parse_checkpoint_file(filename):
    aa_order = "ACDEFGHIKLMNPQRSTVWY"
    altschul = [0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18]
    f = open(filename,'rb')

    buf = f.read(4)
    seqlen = unpack('i',buf)[0]
    output = [[0.0]*20 for _ in range(seqlen)]

    buf = f.read(seqlen)
    seqstr = buf.decode()

    #いきなり全部まとめて unpack できないみたいなので
    #altschul の逆写像定義すれば良いんですけど相違点大きくなるので…
    for i in range(seqlen):
        w = [0.0]*20
        for j in range(20):
            buf = f.read(8)
            w[j] = unpack('d',buf)[0]
        for j in range(20):
            output[i][j] = w[altschul[j]]

    f.close()
    return seqstr,output

def finish_checkpoint_matrix(seq,matrix):
    assert len(seq) == len(matrix) , "Length mismatch between sequence and chekpoint file!\n"
    blos_aa = [0,14,11,2,1,13,3,5,6,7,9,8,10,4,12,15,16,18,19,17]

    aaNum = dict()
    aa_order = "ACDEFGHIKLMNPQRSTVWY"
    for i in range(20): aaNum[aa_order[i]] = i
    aaNum['X'] = 0 #"cheep fix for now"

    b62 = [[0.0]*20 for _ in range(20)]
    blosum62 = Blosum62.split("\n")
    for i in range(20):
        words = list(map(float,blosum62[i].split()))
        for j in range(len(words)):
            b62[i][j] = b62[j][i] = words[j]

    for i in range(20):
        sum = 0.0
        for j in range(20): sum += b62[i][j]
        for j in range(20): b62[i][j] /= sum
    for i in range(len(matrix)):
        sum = 0.0
        for j in range(20): sum += matrix[i][j]
        if sum == 0.0:
            for j in range(20): matrix[i][j] = b62[aaNum[seq[i]]][j]
    return matrix

def write_checkpoint_file(filename,seq,matrix):
    f = open(filename,"x")
    assert len(seq) == len(matrix) , "Length mismatch between sequence and chekpoint file!\n"

    f.write(f"{len(matrix)}\n")
    for row in range(len(matrix)):
        f.write(f"{seq[row]} ")
        for col in range(20):
            f.write("{0:6.3f} ".format(matrix[row][col]))
        f.write("\n")
    f.write("END")
    f.close()

def main():
    seq, matrix = parse_new_checkpoint_file(args.pdb+".2.chk")
    matrix = finish_checkpoint_matrix(seq, matrix)
    write_checkpoint_file(args.pdb+".checkpoint", seq, matrix)

if __name__ == '__main__':
    main()
