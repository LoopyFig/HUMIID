#!usr/bin/env python3

import re, os, sys
import classDef, signatureDef, orthologyDef, geneDef
import orthologySignatureMap

# Reads MinPath output and prior calculated
# orthology relative abundance scores to calculate
# relative abundance scores at the signature level
# as well as signature completion scores
def main():
    # reads through MinPath.py output and diamondScore.py output
    # and generates relative abundance scores for signatures
    # and signature completion scores
    for line in open(sys.argv[3], 'r'):
        result = line.split('\t')
        orthologyDef.orthologyDef[int(result[0])][-1] = float(result[-1])

    #
    details = open(sys.argv[2], 'r')
    detail = details.readline()
    orthologyScores = open(sys.argv[4], 'w')
    signatureMatches = list()
    signatureScores = open(sys.argv[5], 'w')
    for line in open(sys.argv[1], 'r'):
        result = re.split(' +', line)
        if (result[7] == '1'):
            signatureScore = 0.0
            signatureMatch = int(result[1])
            signatureDef.signatureDef[signatureMatch][-2] = float(result[11])
            signatureMatches.append(signatureMatch)
            detail = details.readline()
            while (detail != '' and detail[:4] != 'path'):
                orthologyMatch = int(re.search(r'# o([0-9]+)', detail).group(1))
                orthologyEntry = orthologyDef.orthologyDef[orthologyMatch]
                orthologyScores.write(
                    '\t'.join([str(orthologyMatch), orthologyEntry[0], orthologyEntry[1],
                    str(signatureMatch), signatureDef.signatureDef[signatureMatch][0],
                    signatureDef.signatureDef[signatureMatch][1], str(orthologyEntry[2]),
                    str(orthologyEntry[3]), str(orthologyEntry[4])]) + '\n')
                signatureDef.signatureDef[signatureMatch][-1] += orthologyEntry[-1]
                detail = details.readline()
            signatureEntry = signatureDef.signatureDef[signatureMatch]
            signatureScores.write(
                '\t'.join([str(signatureMatch), signatureEntry[0], signatureEntry[1],
                str(signatureEntry[2]), classDef.classDef[signatureEntry[2]], 
                str(signatureEntry[3]), str(signatureEntry[4])
                , str(signatureEntry[4]/signatureEntry[3]), str(signatureEntry[5])]) + '\n')
    details.close()
    orthologyScores.close()
    signatureScores.close()

if __name__ == '__main__':
    main()
