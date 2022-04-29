#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Take a list of KEGG KO identifiers, look up 
their functional annotations using the KEGG
API (internet access required), and output
a new file with the annotations.
"""
import argparse
import csv
from bioservices import KEGG, KEGGParser, logger
from pathlib import Path


k = KEGG(verbose=False)
kp = KEGGParser()
ko_info_db = {}

# suppress all the SYMBOL warnings
k.logging.level = 'ERROR'
kp.logging.level = 'ERROR'


def partition_KO_rec(ko_field):
    k_range = ['K'+str(i) for i in range(6)]
    total = []
    curr = []
    for line in ko_field:
        if line[:2] in k_range:
            curr.append(line)
            total.append(['' for _ in range(4-len(curr))] + list(curr))
            curr = list()
        else:
            curr.append(line)
    return total

def KO_info(ko_ids):
    ko_entries = {}
    to_check = set(ko_ids).difference(ko_info_db)
    for koid in list(to_check):
        try:
            info = kp.parse(k.get(koid))
        except AttributeError as ae:
            print(koid)
            continue
        except ValueError as ve:
            print(f"{koid} not found")
            continue
        if 'BRITE' in info:
            brite = [e.strip() for e in info['BRITE'].splitlines()]
            # extract all KO pathways
            in_ko = False
            ko_field = []
            for line in brite:
                if line.startswith("KEGG Orthology"):
                    in_ko = True
                    continue
                if "[BR:" in line and "KEGG Orthology" not in line:
                    break
                if in_ko:
                    ko_field.append(line)
            ko_entries[koid] = partition_KO_rec(ko_field)
    return ko_entries


def write_ko_ann(out_fp, ko_ids_ann):
    with open(out_fp, "w") as outf:
        outf.write("\t".join(['ID', 'Level 1', 'Level 2', 'Level 3', 'Level 4'])+"\n")
        for ko in ko_ids_ann:
            for entry in ko_ids_ann[ko]:
                outf.write("\t".join([ko] + entry))
                outf.write("\n")


def handle_program_options():
    """
    Parses the given options passed in at the command line.
    """
    parser = argparse.ArgumentParser(description="Take a list of KEGG KO identifiers,"
                                     " look up their functional annotations using the"
                                     " KEGG API (internet access required), and output"
                                     " a new file with the annotations.")
    parser.add_argument('ko_ids_fp',
                        help="Path to a list of KO IDs.")
    parser.add_argument('-o', '--output_fp', default="KO_annotations.tsv",
                        help="Path to the output file. Default: KO_annotations.tsv")
    parser.add_argument('-c', '--chunksize', default=20, type=int,
                        help='The number of KO annotations to look up at once.'
                             ' Defaults to 20.')

    return parser.parse_args()


def main():
    args = handle_program_options()

    all_KOs = []
    with open(args.ko_ids_fp) as inf:
        all_KOs = [entry.strip() for entry in inf]
    print(f"KOs to be annotated: {len(all_KOs)}")
    
    all_KOs_ann = {}
    chunksize = args.chunksize

    for i, chunk in enumerate(range(chunksize, len(all_KOs), chunksize)):
        print(f"Annotating chunk {i+1} ({chunk}) of {int(len(all_KOs)/chunksize)}")
        all_KOs_ann.update(KO_info(all_KOs[chunk-chunksize:chunk]))

    write_ko_ann(args.output_fp, all_KOs_ann)
    print(f"Annotation Complete.\nOutput written to: {args.output_fp}")

if __name__ == "__main__":
    main()