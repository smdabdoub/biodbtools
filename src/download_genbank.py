#!/usr/bin/env python
# coding: utf-8
"""
Given a list of GenBank accession IDs, download all complete
genomes and contigs for non-complete genomes. Write out the
results to both GB and FASTA-format files for use with 
K-SLAM and CLARK respectively.
"""
import argparse
import itertools
import io
import os
import os.path as osp
import sys
from urllib.error import HTTPError
from Bio import Entrez, SeqIO

def expand_entries(gb_ids):
    """
    Checking each row for the presence of multiple 
    accession IDs and split by a comma separator, 
    flattening the file.
    Typically the multiple IDs correspond to,
    first, the organism genome, and the rest being
    plasmids or other separate genomes contained
    within the cell.
    
    >>> expand_entries(["A,B", "D, E", "G,H, I"])
    ['A', 'B', 'D', 'E', 'G', 'H', 'I']
    """
    ids = []
    for row in gb_ids:
        ids.extend([entry.strip() for entry in row.split(",")])
    return ids


def pp_list(items):
    """
    Pretty-print a list of items.
    
    >>> pp_list(['a', 'b', 'c'])
    a, b, c
    """
    step = 10
    for i in range(0, len(items), step):
        print(", ".join(items[i:i+step]))


def esearch_range(accn, db="nuccore", retmax=100, webenv=None):
    """
    Given the first and last in a series of GenBank accessions,
    retrieve the database identifiers.
    
    retmax: By default, ESearch returns only the first 20 records
            (retmax=20) and discards the rest. Max of 100,000.
    webenv: If a WebEnv string is passed, retmax will be ignored
            (max 20 still returned), and all search results will be
            stored in your personal history with an associated
            query key ('QueryKey' in returned dictionary).
    """
    MAX_RETMAX = 100000
    kwds = {}
    
    if webenv is not None:
        kwds = {'usehistory': 'y', 'WebEnv': webenv}
    else:
        kwds['retmax'] = retmax if retmax <= MAX_RETMAX else MAX_RETMAX
        
    with Entrez.esearch(db="nuccore", term=f"{':'.join(accn)}[accn]", **kwds) as handle:
        return Entrez.read(handle)


def epost(id_list, db="nuccore"):
    """
    Given a list of NCBI identifiers (e.g. from ESearch), submit them
    to the Entrez History server for later use with EFetch.
    
    Returns: WebEnv and QueryKey values
    """
    with Entrez.epost(db, id=",".join(id_list)) as request:
        result = Entrez.read(request)
        return result["WebEnv"], result["QueryKey"]


def batch_efetch(id_list=None, batch_size=100, db="nuccore", 
                 rettype="gbwithparts", retmode="text",
                 webenv=None, query_key=None, query_len=None):
    """
    Using a pre-posted query, use the WebEnv and QueryKey identifiers to
    retrieve the actual data behind the query items. The theoretical max
    should be 500 per EFetch when using History (see: 
    https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large),
    but here the default is set to 100 as a safe margin.
    """
    kwds = {}
    if webenv is None or query_key is None:
        webenv, query_key = epost(id_list, db)
        id_count = len(id_list)
    else:
        id_count = query_len

    kwds['query_key'] = query_key
    kwds['WebEnv'] = webenv
    
    complete = []
    contigs = {}
    for start in range(0, id_count, batch_size):
        end = min(id_count, start+batch_size)
        print(f"Downloading records {start+1} to {end}...", end='')
        if id_list is not None:
            pp_list(id_list[start:end])
        else:
            print("\n")
        
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                gb_handle = Entrez.efetch(db=db, rettype=rettype, retmode=retmode,
                                          retstart=start, retmax=batch_size,
                                          **kwds)
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    time.sleep(15)
                else:
                    raise
        
        for seq_record in SeqIO.parse(gb_handle, "gb"):
            if "wgs" in seq_record.annotations:
                contigs[seq_record.id] = seq_record.annotations["wgs"]
            else:
                complete.append(seq_record)
    
    return complete, contigs, webenv


def write_gb_fasta(gb_records, out_dir, genome_id=None, write_gb=True):
    for rec in gb_records:
        filename=rec.id if genome_id is None else genome_id
        mode = 'w' if genome_id is None else 'a'

        if write_gb:
            with open(osp.join(out_dir,"gbff", f"{filename}.gbff"), mode) as gb_f:
                gb_f.write(rec.format("gb"))

        with open(osp.join(out_dir,"fasta",f"{filename}.fa"), mode) as fasta_f:
            fasta_f.write(rec.format("fasta"))


def handle_program_options():
    """Parses the given options passed in at the command line."""
    parser = argparse.ArgumentParser(description="Given a list of GenBank "
                                     "accession IDs, download all complete "
                                     "genomes and contigs for non-complete "
                                     "genomes. Write out the results to both "
                                     "GBFF and FASTA-format files.")
    parser.add_argument("genbank_ids", metavar='genbank-ids',
                        help="Path to a file containing a list of GenBank "
                        "accession IDs (one per line).")
    parser.add_argument("-b", "--batch-size", type=int, default=100, 
                        choices=range(1,501), metavar="[1-500]",
                        help="The number of records to attempt to fetch at "
                        "once. NCBI mandates an upper limit of 500, but you "
                        "may experience connection issues with too large a "
                        "number. 100 to 200 is usually a safe range. Defaults "
                        "to 100.")
    parser.add_argument("-o", "--out-dir", default=".",
                        help="Downloaded data will be written to the gbff and "
                        "fasta directories inside this directory.")
    parser.add_argument("--no-gbff", action="store_false", 
                        help="If specified, disable writing GBFF files.")

    return parser.parse_args()


def main():
    args = handle_program_options()

    # check/create output dir
    if not osp.exists(args.out_dir):
        os.mkdir(args.out_dir)

    with open(args.genbank_ids) as inf:
        gb_ids = expand_entries([row for row in inf])

    # identify this script to NCBI
    Entrez.email = ""
    Entrez.tool = ""
    
    print(f"Downloading GenBank records for {len(gb_ids)} accessions.\n")
    
    for i in range(0, len(gb_ids), args.batch_size):
        gb_chunk = gb_ids[i:i+args.batch_size]
        complete, contigs, webenv = batch_efetch(gb_chunk, batch_size=args.batch_size)
        print(f"\nDownloaded {len(complete)} complete genomes ", end='')
        print(f"and accession IDs for {len(contigs)} incomplete genomes.")
        
        # download GB sequence data for complete genomes
        print(f"\nWriting sequence files for {len(complete)} complete genomes.")
        write_gb_fasta(complete, args.out_dir, write_gb=args.no_gbff)
        
        # download the GB sequence data for all the contigs
        print("Downloading contig data for incomplete genomes.")
        for genome in contigs:
            es = esearch_range(contigs[genome], webenv=webenv)
            print(f"\nFetching {es['Count']} GenBank contig records for genome {genome}...\n")
            contig_complete, _, _ = batch_efetch(webenv=webenv, query_key=es['QueryKey'], 
                                                 query_len=int(es['Count']))
            print("\n...complete.")
        
            # write out all data to both GenBank format and FASTA format files
            print(f"\nWriting sequence files for {es['Count']} contigs.")
            
            write_gb_fasta(contig_complete, args.out_dir, genome_id=genome, 
                           write_gb=args.no_gbff)


if __name__ == '__main__':
    main()
