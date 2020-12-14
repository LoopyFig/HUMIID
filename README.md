# HUMIID
HUMIID is a tool to identify KEGG pathogenicity signatures in metagenomic data.
To use HUMIID you must build or pull the docker image, assumed to be named loopyfig/humiid,
and run the following command (replace loopyfig/humiid with actual image name if different):

sudo docker run -v {inputDirectory}:/data loopyfig/humiid pathogenicity {inputFASTA} {outputDirectory}

{inputDirecotry}: the directory containing the input file and output directory, full path
{inputFASTA}: a fasta format metagenomic dataset file in the input directory, path from input directory
{outputDirectory}: the directory that will contain the output, created if not already present, path from input directory

The output directory includes the following output files:

out.tsv: The DIAMOND blastx output

minPath.mp: The input for minPath, tab separated fields - gene_key, signature_key

geneScores.gout: The relative abundance scores for the genes, tab separated fields - gene_key, gene_id, gene_description, sequence_length, orthology_key, orthology_id, orthology_description, relative_abundance_score

orthologyScores.oout: The relative abundance scores for the orthologies, tab separated fields - orthology_key, orthology_id, orthology_description, sum_of_sequence_lengths, count_of_genes, relative_abundance_score

orthologyScores.oout.dup: The relative abundance scores for the orthologies duplicated for shared signatures, tab separated fields - orthology_key, orthology_id, orthology_description, signature_key, signature_id, signature_description, sum_of_sequence_lengths, count_of_genes, relative_abundance_score

minPath.mp.minpath: The MinPath output

minPath.mp.minpath.details: The MinPath output details

orthologyScores.oout.sig: The relative abundance scores for the orthologies duplicated for shared signatures selected by MinPath, tab separated fields - orthology_key, orthology_id, orthology_description, signature_key, signature_id, signature_description, sum_of_sequence_lengths, count_of_genes, relative_abundance_score

signatureScores.sout: The relative abundance scores for the signatures, tab separated fields - signature_key, signature_id, signature_description, orthology_count, orthologies_found, found_total_ratio, relative_abundance_score

Two other options include the help command (this one) run with:

sudo docker run loopyfig/humiid help

and the run command which will allow you to input a command that is run inside the container, run with:

sudo docker run loopyfig/humiid run "{commands}"

{commands}: the desired commands, flanked by quotes; this allows use of internal packages including diamond and bbmap
