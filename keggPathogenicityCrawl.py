#!usr/bin/env python3

import urllib.request, urllib.error, urllib.parse
import re, os, sys

def main():
    # description incoming
    signatureOrthologyPattern = (
        r'<area shape="rect" id="node[^"]*" href="([^"]+)" title="([^"]+)" alt="[^"]*" coords="[^"]*"/>\\n'
        )
    orthologyGeneListPattern = (
        r'Genes</nobr></th>\\n<[^<>]*>(.*?)\\n'
        )
    genePattern = (
        r'<nobr>[^>:]*:&nbsp;</nobr></td><td><a href="([^"]+)">([^><]+)</a>([^><]*)'
        )
    definitionPattern = r'Definition</nobr></th>\\n<[^<>]*><[^<>]*><[^<>]*>(.*?)<br>\\n'
    aaPattern = (
        r'([0-9]+) aa <.*?AA seq.*?>DB search</button><br>\\n([A-Z<br>\\n]+)</td></tr>\\n'
        )

    classDefDict = dict()
    classDef = open('classDef.py', 'w')
    classDef.write('classDef = [\n')
    classDefKey = 0
    signatureDefDict = dict()
    signatureDef = open('signatureDef.py', 'w')
    signatureDef.write('signatureDef = [\n')
    signatureDefKey = 0
    orthologyDefDict = dict()
    orthologyDef = open('orthologyDef.py', 'w')
    orthologyDef.write('orthologyDef = [\n')
    orthologyDefKey = 0
    geneDefDict = dict()
    geneDef = open('geneDef.py', 'w')
    geneDef.write('geneDef = [\n')
    geneDefKey = 0
    orthologySignatureMapDict = dict()
    orthologySignatureMap = open('orthologySignatureMap.py', 'w')
    orthologySignatureMap.write('orthologySignatureMap = {\n')
    fasta = open('pathogenicityGenes.fa', 'w')

    for line in open('signature.csv', 'r'):
        signature = line.split('\t')
        if (not (signature[0] in classDefDict)):
            classDefDict[signature[0]] = classDefKey
            classDef.write('    "' + signature[0] + '",\n')
            classDefKey = classDefKey + 1

        #
        signatureDefDict[signature[1]] = signatureDefKey
        signatureDef.write('    ["' + signature[1]
            + '", "' + signature[2] + '", "'+ str(classDefDict[signature[0]]) + '", 0],\n')
        signatureDefKey = signatureDefKey + 1

        #
        signatureData = str(urllib.request.urlopen(signature[3]).read())
        print(signature)
        for signatureOrthology in re.findall(signatureOrthologyPattern, signatureData):
            print('\t' + str(signatureOrthology))
            orthologyData = str(urllib.request.urlopen(signatureOrthology[0]).read())
            #print(orthologyData)
            orthologyDefinition = re.search(definitionPattern, orthologyData).group(1)
            #print('\t\t' + orthologyDefinition)
            if (not (signatureOrthology[1] in orthologyDefDict)):
                orthologyDefDict[signatureOrthology[1]] = orthologyDefKey
                orthologySignatureMapDict[orthologyDefKey] = list()
                orthologyDefKey = orthologyDefKey + 1

                #
                orthologyGeneList = re.search(orthologyGeneListPattern, orthologyData).group(1)
                #print(orthologyGeneList)
                seqSumLen = 0
                seqCount = 0
                for orthologyGene in re.findall(genePattern, orthologyGeneList):
                    #print('\t\t\t' + str(orthologyGene))
                    geneName = orthologyGene[1] + orthologyGene[2]
                    geneData = str(urllib.request.urlopen('https://www.genome.jp' + orthologyGene[0]).read())
                    geneDefinition = re.search(definitionPattern, geneData).group(1)
                    #print('\t\t\t\t' + geneDefinition)
                    if (not (geneName in geneDefDict)):
                        aaSeq = re.search(aaPattern, geneData)
                        seqSumLen = seqSumLen + int(aaSeq.group(1))
                        seqCount = seqCount + 1
                        geneDefDict[geneName] = geneDefKey
                        geneDef.write('    ["' + geneName
                            + '", "' + geneDefinition
                            + '", ' + aaSeq.group(1)
                            + ', ' + str(orthologyDefDict[signatureOrthology[1]]) + ', 0],\n')
                        geneDefKey = geneDefKey + 1
                        aaSeq = re.search(aaPattern, geneData)
                        fasta.write(
                            '>' + str(geneDefDict[geneName]) +
                            ':' + str(orthologyDefDict[signatureOrthology[1]]) +
                            ':' + aaSeq.group(1) + '\n' +
                            aaSeq.group(2).replace('<br>\\n', '\n') + '\n'
                            )
                    else:
                        print('Orthology-Gene Key Violated:' + str(orthologyGene))
                orthologyDef.write('    ["' + signatureOrthology[1]
                    + '", "' + orthologyDefinition + '", ' + str(seqSumLen) + ', ' + str(seqCount) + ', 0],\n')
            else:
                print('Orthology with Two Signatures: ' + str(signatureOrthology))
            orthologySignatureMapDict[orthologyDefDict[signatureOrthology[1]]].append(signatureDefDict[signature[1]])

if __name__ == '__main__':
    main()
