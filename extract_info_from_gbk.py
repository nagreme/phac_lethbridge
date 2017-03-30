#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 31, 2017
# Purpose: Parse .gbk files and extract start, end, strand, product,
#          locus tag, gene name, protein ID, codon start, and transl
#          table. (If present of course) (tab-delimited output .csv)
#          And add a column with a selected gene name (4-letter > protein_ID > product)
#----------------------------------------------------------------------

import sys
import re
from Bio import SeqIO


#---------------------------------Main---------------------------------

def main():
    #Establish filenames
    inFilename = sys.argv[1]
    outFilename = inFilename[0:inFilename.find(".gbk")] + "_genes_info.csv"

    #open the output file
    outFile = open(outFilename, 'w')

    #print header for output file
    outFile.write("genome_name\tfeature_type\tstart\tend\tstrand\tproduct\tlocus_tag\tgene_name\tprotein_ID\tcodon_start\ttransl_table\tselected_gene_name\n")

    #setup to read the genbank file
    for record in SeqIO.parse(inFilename, "genbank"):
        #iterate_features = iter(record.features)

        for feature in record.features:
            feat_type = feature.type
            #we only want features we actually use
            if (feat_type == "CDS" or feat_type == "tRNA" or feat_type == "rRNA" or feat_type == "tmRNA"):

                #not all features have all fields so guard against this
                if ('gene' not in feature.qualifiers):
                    gene = ""
                else:
                    gene = "GN:{}".format(repr(feature.qualifiers['gene']).strip("[]'\""))

                if ('codon_start' not in feature.qualifiers):
                    codon_start = ""
                else:
                    codon_start = repr(feature.qualifiers['codon_start']).strip("[]'")

                if ('transl_table' not in feature.qualifiers):
                    transl_table = ""
                else:
                    transl_table = repr(feature.qualifiers['transl_table']).strip("[]'")

                if ('product' not in feature.qualifiers):
                    product = ""
                else:
                    product = repr(feature.qualifiers['product']).strip("[]'")

                if ('protein_id' not in feature.qualifiers):
                    protein_id = ""
                else:
                    protein_id = "Pid:{}".format(repr(feature.qualifiers['protein_id']).strip("[]'"))

                if ('locus_tag' not in feature.qualifiers):
                    locus_tag = ""
                else:
                    locus_tag = repr(feature.qualifiers['locus_tag']).strip("[]'")

                #select a gene name depending on what's there (4-letter > protein_ID > product)
                if (gene):
                    selected_gene_name = gene
                elif (protein_id):
                    selected_gene_name = protein_id
                else:
                    selected_gene_name = ""

                #build and write the line for that feature
                outFile.write("{genome_name}\t{feature_type}\t{start}\t{end}\t{strand}\t{product}\t{locus_tag}\t{gene_name}\t{protein_ID}\t{codon_start}\t{transl_table}\t{selected_gene_name}\n".format(
                genome_name=record.name, feature_type=feature.type, start=feature.location.start, end=feature.location.end, strand=feature.location.strand, product=product,
                locus_tag=locus_tag, gene_name=gene, protein_ID=protein_id, codon_start=codon_start, transl_table=transl_table, selected_gene_name=selected_gene_name))


    #close the file
    outFile.close()





if __name__ == '__main__':
    main()
