#!usr/bin/env python3

import urllib.request, urllib.error, urllib.parse
import re, os, sys

# Reads a signature file specified as the first argument
# and generates python files with lists for:
# signature classes, signatures, orthologies, and genes
#
# Also generates the following:
# a fasta file with all reference genes
# a python file containing mappings from orthologies to signatures
# a mapping file for signatures and orthologies for use with MinPath
#
# The original signature file has the following tab-separated fields:
# signature_class, signature_id, signature_definition, signature_url
#
# Run with: python3 keggPathogenicityCrawl.py signature.csv
# where signature.csv is any signature file
def main():
    # regex patterns for searching webpages
    signatureOrthologyPattern = (
        r'<area shape="rect" id="node[^"]*" href="([^"]+)" title="([^"]+)" alt="[^"]*" coords="[^"]*"/>\\n'
        )
    orthologyGeneListPattern = (
        r'Genes</nobr></th>\\n<[^<>]*>(.*?)\\n'
        )
    genePattern = (
        r'<nobr>[^>:]*:&nbsp;</nobr></td><td><a href="([^"]+)">([^><]+)</a>([^><]*)'
        )
    definitionPattern = r'Definition</nobr></th>\\n<.*?>([^><]*?)<br>\\n'
    aaPattern = (
        r'([0-9]+) aa <.*?AA seq.*?>DB search</button><br>\\n([A-Z<br>\\n]+)</td></tr>\\n'
        )

    # data structures used to generate outputs
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
    orthologySignatureMapList = list()
    orthologySignatureMap = open('orthologySignatureMap.py', 'w')
    orthologySignatureMap.write('orthologySignatureMap = [\n')
    fasta = open('pathogenicityGenes.fa', 'w')
    signatureOrthologyMap = open('signatureOrthologyMap.sig', 'w')

    # reads through signature file and crawls through webpages
    # to extract data
    for line in open(sys.argv[1], 'r'):
        signature = line.split('\t')

        # generates class definition as a python list:
        # class_id
        # class key corresponds to position in list
        if (not (signature[0] in classDefDict)):
            classDefDict[signature[0]] = classDefKey
            classDef.write('    "' + signature[0] + '",\n')
            classDefKey = classDefKey + 1

        # sets internal data structure for signature
        signatureDefDict[signature[1]] = signatureDefKey
        signatureDefKey = signatureDefKey + 1

        # opens signature url to get list of orthologies and crawl through them
        signatureData = str(urllib.request.urlopen(signature[3]).read())
        print(signature)
        orthologyCount = 0
        for signatureOrthology in re.findall(signatureOrthologyPattern, signatureData):
            print('\t' + str(signatureOrthology))
            orthologyCount = orthologyCount + 1
            # reads orthology data on orthology url
            orthologyData = str(urllib.request.urlopen(signatureOrthology[0]).read())
            #print(orthologyData)
            orthologyDefinition = re.search(definitionPattern, orthologyData).group(1)
            #print('\t\t' + orthologyDefinition)
            if (not (signatureOrthology[1] in orthologyDefDict)):
                orthologyDefDict[signatureOrthology[1]] = orthologyDefKey
                orthologySignatureMapList.append(list())
                orthologyDefKey = orthologyDefKey + 1

                # searches for list of genes in orthology data to crawl through them
                orthologyGeneList = re.search(orthologyGeneListPattern, orthologyData).group(1)
                #print(orthologyGeneList)
                seqSumLen = 0
                seqCount = 0
                for orthologyGene in re.findall(genePattern, orthologyGeneList):
                    #print('\t\t\t' + str(orthologyGene))
                    geneName = orthologyGene[1] + orthologyGene[2]

                    # reads gene data, checks to see if page discontinued/empty before retrieving
                    # definition and amino acid sequence data
                    geneData = str(urllib.request.urlopen('https://www.genome.jp' + orthologyGene[0]).read())
                    try:
                        geneDefinition = re.search(definitionPattern, geneData).group(1)
                    except:
                        print('Empty Gene Entry: ' + str(orthologyGene))
                        continue
                    #print('\t\t\t\t' + geneDefinition)
                    if (not (geneName in geneDefDict)):
                        aaSeq = re.search(aaPattern, geneData)
                        seqSumLen = seqSumLen + int(aaSeq.group(1))
                        seqCount = seqCount + 1
                        geneDefDict[geneName] = geneDefKey
                        # writes gene information as a list:
                        # gene_id, gene_definition, sequence_length, orthology_key, relative_abundance_placeholder
                        # gene key corresponds to position in list
                        #
                        # fasta file is written at the same time with following format:
                        # gene_key:orthology_key:sequence_length
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
                        # alerts situation with violated one to many relationship
                        # for orthology to genes
                        print('Orthology-Gene Key Violated:' + str(orthologyGene))
                # orthology definition written as a python list:
                # orthology_id, orthology_definition, sequence_length, count_of_sequences, relative_abundance_placeholder
                orthologyDef.write('    ["' + signatureOrthology[1]
                    + '", "' + orthologyDefinition + '", ' + str(seqSumLen) + ', ' + str(seqCount) + ', 0],\n')
            else:
                print('Orthology with Two Signatures: ' + str(signatureOrthology))
            # keeps track of signature-orthology relationships, assumed to be many-to-many
            orthologySignatureMapList[orthologyDefDict[signatureOrthology[1]]].append(signatureDefDict[signature[1]])
        # genearates signature definition as a python list:
        # signature_id, signature_definition, class_key, orthology_count, orthology_found_placeholder, abundance_score_placeholder
        # signature key corresponds to position in list
        signatureDef.write('    ["' + signature[1]
            + '", "' + signature[2] + '", '+ str(classDefDict[signature[0]]) + ', ' + str(orthologyCount) + ', 0, 0],\n')

    # this section writes the signature-orthology relationship as a
    # python list of lists (with sublists containing signatures and
    # position in list corresponding to orthology key) and as a
    # tab-separated map of signatures to orthologies for use with MinPath
    orthologyKey = 0
    for signatureMapList in orthologySignatureMapList:
        orthologySignatureMap.write('    ' + str(signatureMapList) + ',\n')
        for signatureKey in signatureMapList:
            signatureOrthologyMap.write(str(signatureKey) + '\to' + str(orthologyKey) + '\n')
        orthologyKey = orthologyKey + 1
    classDef.write(']\n')
    classDef.close()
    signatureDef.write(']\n')
    signatureDef.close()
    orthologyDef.write(']\n')
    orthologyDef.close()
    geneDef.write(']\n')
    geneDef.close()
    orthologySignatureMap.write(']\n')
    orthologySignatureMap.close()
    fasta.close()
    signatureOrthologyMap.close()

if __name__ == '__main__':
    main()
