#!usr/bin/env python3

import re, os, sys
import signatureDef, orthologyDef, geneDef
import orthologySignatureMap

# Reads diamond output specified as the first argument
# and prepares MinPath input and abundance scores at
# gene and orthology level
def main():
    # reads through diamond output and uses prepared
    # lists to help generate abundance scores
    readId = ''
    alignments = list()
    identityAdjustSum = 0.0
    geneMatches = list()
    orthologyMatches = list()
    for line in open(sys.argv[1], 'r'):
        result = line.split('\t')
        geneInfo = result[1].split(':')

        # checks for new read to finish abundance calculations for last read
        # and start new alignment list
        if (result[0] != readId):
            for alignment in alignments:
                if (geneDef.geneDef[alignment[0]][-1] == 0):
                    geneMatches.append(alignment[0])
                    if (orthologyDef.orthologyDef[alignment[1]][-1] == 0):
                        orthologyMatches.append(alignment[1])
                geneDef.geneDef[alignment[0]][-1] += (alignment[2]/identityAdjustSum)
                orthologyDef.orthologyDef[alignment[1]][-1] += (
                    (alignment[2]*orthologyDef.orthologyDef[alignment[1]][-2])/
                    (identityAdjustSum*orthologyDef.orthologyDef[alignment[1]][-3]))
            readId = result[0]
            alignments = list()
            identityAdjustSum = 0.0
        # iteratively calculates total squared identity and prepares list
        # for abundance calculation
        identityAdjust = pow(float(result[2]), 2)
        alignments.append([int(geneInfo[0]), int(geneInfo[1]), identityAdjust])
        identityAdjustSum += identityAdjust

    # because the prior abundance calculation is based on detecting a new read
    # and there is no new read at the end, the abundance calculation block is repeated
    for alignment in alignments:
        if (geneDef.geneDef[alignment[0]][-1] == 0):
            geneMatches.append(alignment[0])
            if (orthologyDef.orthologyDef[alignment[1]][-1] == 0):
                orthologyMatches.append(alignment[1])
        geneDef.geneDef[alignment[0]][-1] += (alignment[2]/identityAdjustSum)
        orthologyDef.orthologyDef[alignment[1]][-1] += (
            (alignment[2]*orthologyDef.orthologyDef[alignment[1]][-2])/
            (identityAdjustSum*orthologyDef.orthologyDef[alignment[1]][-3]))

    # writes gene abundance scores and the input for MinPath
    geneScores = open(sys.argv[3], 'w')
    minPathInput = open(sys.argv[2], 'w')
    for geneMatch in geneMatches:
        geneEntry = geneDef.geneDef[geneMatch]
        geneScores.write(
            '\t'.join([str(geneMatch), geneEntry[0], geneEntry[1],
            str(geneEntry[2]), str(geneEntry[3]), orthologyDef.orthologyDef[geneEntry[3]][0],
            orthologyDef.orthologyDef[geneEntry[3]][1], str(geneEntry[4])]) + '\n')
        minPathInput.write(str(geneMatch) + '\to' + str(geneEntry[3]) + '\n')

    # writes orthology abundance scores, both uniquely and also duplicated for signatures
    orthologyScores = open(sys.argv[4], 'w')
    orthologyDupeScores = open(sys.argv[4] + '.dup', 'w')
    for orthologyMatch in orthologyMatches:
        orthologyEntry = orthologyDef.orthologyDef[orthologyMatch]
        orthologyScores.write(
            '\t'.join([str(orthologyMatch), orthologyEntry[0], orthologyEntry[1],
            str(orthologyEntry[2]), str(orthologyEntry[3]), str(orthologyEntry[4])]) + '\n')
        for signatureMatch in orthologySignatureMap.orthologySignatureMap[orthologyMatch]:
            orthologyDupeScores.write(
                '\t'.join([str(orthologyMatch), orthologyEntry[0], orthologyEntry[1],
                str(signatureMatch), signatureDef.signatureDef[signatureMatch][0],
                signatureDef.signatureDef[signatureMatch][1], str(orthologyEntry[2]),
                str(orthologyEntry[3]), str(orthologyEntry[4])]) + '\n')
    geneScores.close()
    minPathInput.close()
    orthologyScores.close()
    orthologyDupeScores.close()

if __name__ == '__main__':
    main()
