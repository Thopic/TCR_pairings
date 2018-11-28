#!/usr/bin/env python3

""" Transform a (sorted) file listing all sequences
    The returned file (given by the -o option) will contain all
    the unique sequences, the well they are in and the number of copy
    in each well.
    -o : name of the returned file
    -w : minimum number of wells necessary to keep a sequence (default 3)
    -p : if given, with a float, compute a table containing
         the number of wells wxy two sequences must share to be paired
         knowing that the sequences are respectively in wx and wy wells.
         The float is the log of the maximal p-value wanted.
"""

import argparse
import scipy.special
import numpy as np

def comb(n, k):
    """ Return n choose k 
    """
    return scipy.special.comb(n, k, exact=True)


def local_pvalues(wxy, wx, wy, W):
    """ return the p-value,
        The p-value (wx, wy, wxy) is the probability of randomly 
        finding two sequences (x, y) in more than wxy shared wells
        knowing that they were initially in wx and wu wells
    """
    return sum([comb(wx, x)*comb(W-wx, wy-x)
                for x in range(wxy, min(wx, wy)+1)])/comb(W, wy)


def create_p_values_table(W, logprecision):
    """ Create the table of p-values for W wells
    """
    tablogpv = np.zeros((W+1, W+1, W+1), dtype=np.float64)
    fw1 = open("pvalues.tsv", "w")
    for wa in range(W+1):
        for wb in range(W+1):
            for wab in range(1, min(wa, wb) + 1):
                pv = local_pvalues(wab, wa, wb, W)
                fw1.write("{} {} {} {}\n".format(wa, wb, wab, pv))
                tablogpv[wa, wb, wab] = np.log10(pv)
    fw1.close()
    
    tabpvbool = (tablogpv < logprecision)
    valid = [[] for _ in range(W+1)]
    wab_min = np.full((W+1, W+1), W+4)
    for wa in range(W+1):
        for wb in range(W+1):
            if any(tabpvbool[wa, wb, :]):
                valid[wa].append(wb)
                wab_min[wa, wb] = np.argmax(tabpvbool[wa, wb, :])

    fw2 = open("dvalues.tsv", "w")
    fw2.write("{}\n".format(W))
    for n1 in range(W+1):
        for n2 in range(W+1):
            fw2.write("{} {} {}\n".format(n1, n2, max(int(n1 + n2 - 2*wab_min[n1, n2]), -1)))
    fw2.close()


def write_sequence(fw, nb_sequence, seqs, wnames, min_nb_wells):
    """ Return the final string associated with a (list of) sequences
    """
    copies = [0]*len(wnames) # number of copies:
    wells = [0]*len(wnames)
    for t in seqs:
        wells[wnames[t[-1].strip()]] = 1
        copies[wnames[t[-1].strip()]] += int(t[2])
    if(sum(wells) >= min_nb_wells):
        t = seqs[0]
        line = (str(nb_sequence) + "\t" + t[0] + "\t" + t[1]
            + "\t" + str(sum(copies)) + "\t" + "\t".join(t[3:-1])
            + "\t" + str(wells) + "\t" + str(copies)
            + "\t" + str(sum(wells))
            + "\t" + "".join([str(u) for u in wells]) + "\n")
        fw.write(line)
        return nb_sequence + 1
    return nb_sequence

     

def list_sequences(filein, fileout, min_wells):
    """ Main function, browse the entry file and write the unique sequences
        Return the number of wells
    """
    colname = ["index", "sequence", "amino", "copy", "cdr3Length", "vname",
               "dname", "jname", "vdel", "d5del", "d3del", "jdel", "n2ins",
               "n1ins", "status", "wells", "copies", "nb_wells", "short_wells"]

    ## Number the wells
    wnames = []
    with open(filein, 'r') as f:
        for line in f:
            wname = line.split()[-1]
            if wname not in wnames:
                wnames.append(wname.strip())
    wnames.sort()
    wnames = {v:i for i, v in enumerate(wnames)}

    ## Write the output
    nb_sequence = 0
    with open(filein, 'r') as f, open(fileout, 'w') as fw:
        fw.write("\t".join(colname)+"\n")
        previous_sequences = []
        for l in f:
            if(l[0] not in ["A", "T", "G", "C"]): # remove unwanted headers
                continue
            t = l.split("\t")
            if(previous_sequences != [] and
               t[0] != previous_sequences[0][0]): # write the previous sequence
                if len(previous_sequences) >= min_wells:
                    nb_sequence = write_sequence(fw, nb_sequence,
                                                 previous_sequences,
                                                 wnames, min_wells)
                previous_sequences = []
            previous_sequences.append(t)
        if len(previous_sequences) >= min_wells: # write potential last seq
            write_sequence(fw, nb_sequence, previous_sequences, wnames, min_wells)
            
    return len(wnames)


def main():
    """ Parse the arguments and call list_sequences
    """
    parser = argparse.ArgumentParser(
        description="Create a list of unique sequences")
    parser.add_argument("input", type=str,
                        help="Input file, list of all sequences to be parsed")
    parser.add_argument("output", type=str, help="Output file")
    parser.add_argument("-w", "--min_number_wells", type=int, default=3,
                        help="Minimal number of wells necessary to keep a copy")
    parser.add_argument("-p", "--p_values", type=float,
                        help="Create a table of possible p-values")

    args = parser.parse_args()
    print("List all sequences... (minimum number of wells: {})".format(args.min_number_wells))
    nb_wells = list_sequences(args.input, args.output, args.min_number_wells)
    if(args.p_values):
        print("Create table of p-values... (log of the precision {})".format(args.p_values))
        create_p_values_table(nb_wells, args.p_values)


if __name__ == "__main__":
    main()
