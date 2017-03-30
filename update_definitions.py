import json
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--alleles', required = True, help = "Allele directory")
    parser.add_argument('-j', '--jsons', required = True, help = "JSON directory")
    parser.add_argument('-t', '--test', required = True, help = "Test name (i.e. column 2 in .markers file)")

    return parser.parse_args()

def load_data(jsonpath):
    """Opens JSON file and returns it
    as Python dictionary.
    """

    with open(jsonpath, 'r') as f:
        data = json.load(f)
    return data

def get_known_alleles(allele_dir):
    """Reads each Fasta in allele_dir and
    stores the sequence definitions in a list.

    Returns a dictionary. Keys are gene names, value is list of seqs.
    """

    # strips path and extension
    get_name = lambda x: os.path.basename(os.path.splitext(x)[0])

    known_alleles = {}

    for fname in os.listdir(allele_dir):

        fpath = os.path.join(allele_dir, fname)
        with open(fpath, 'r') as f:
            alleles = [str(x.seq) for x in SeqIO.parse(f, 'fasta')]

        known_alleles[get_name(fname)] = alleles

    return known_alleles

def update_alleles(known_alleles, allele_dir):
    """Overwrites allele files with updated definitions."""

    for gene in known_alleles:

        fname = os.path.join(allele_dir, gene + ".fasta")
        with open(fname, 'w') as f:
            out = ""
            counter = 0
            for allele in known_alleles[gene]:
                out += ">{}\n".format(counter + 1)
                out += known_alleles[gene][counter] + "\n"
                counter += 1
            f.write(out)

def revcomp(seq):
    #
    # comp = {'A': 'T', 'a': 't', 'T': 'A', 't': 'a', 'G':'C', 'g': 'c', 'C': 'G', 'c': 'g', 'N':'N'}
    # return ''.join([comp[x] for x in reversed(seq)])
    return Seq(seq).reverse_complement().__str__() #Mod March 7, 2017 (NP-V)

def update(json_dir, known_alleles, test):
    """Reads JSON output from MIST.

    If a new allele is discovered, it is appended to the list of known alleles
    and the JSON is updated accordingly.

    Returns dict of known alleles and JSON objects for later writing.
    """

    for fname in os.listdir(json_dir):

        json_name = os.path.join(json_dir, fname)
        data = load_data(json_name)

        genes = data["Results"][0]["TestResults"][test]

        for g in genes:

            gene = genes[g]

            if gene['BlastResults'] is None:
                continue

            if gene['IsContigTruncation']:
                continue

            if not gene["CorrectMarkerMatch"]:
                br = gene["BlastResults"]

                sj = br["SubjAln"].replace('-', '')

                br["QueryAln"] = sj
                br["SubjAln"] = sj

                rc = revcomp(sj)

                if sj not in known_alleles[g]:

                    known_alleles[g].append(sj)

                br["Mismatches"] = 0
                br["Gaps"] = 0
                br["PercentIdentity"] = 100.0
                gene["Mismatches"] = 0
                gene["BlastPercentIdentity"] = 100.0
                gene["CorrectMarkerMatch"] = True


                gene["MarkerCall"] = str(known_alleles[g].index(sj) + 1)

                gene["AlleleMatch"] = gene["MarkerCall"]

                gene["BlastResults"] = br
                genes[g] = gene

            data["Results"][0]["TestResults"][test] = genes


        with open(json_name, 'w') as f:
            json.dump(data, f, indent = 4,
                      separators = (',', ': '), sort_keys = True)

    return known_alleles

def update_definitions(alleles, jsons, test):

    known_alleles = get_known_alleles(alleles)

    known_alleles = update(jsons, known_alleles, test)

    update_alleles(known_alleles, alleles)

def main():

    args = arguments()

    update_definitions(args.alleles, args.jsons, args.test)

if __name__ == '__main__':
    main()
