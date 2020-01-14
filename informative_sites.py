import os
from smallBixTools import smallBixTools as st
import argparse
import sys


def condense(file_to_condense, outdir, caseSensitive):
    """
    takes a filename of a fasta file.
    removes all conserved sites
    writes out the reduced file, containing no conserved sites.
    :param d: a fasta file path.
    :param outdir: the place to write the output file to
    :param caseSensitive: bool indicates if the case of the characters in the sequences counts. If True, then we
    should respect the case of the sequences. "a" != "A".
    :return: no return
    """

    # establish filenames
    base_filename = os.path.splitext(os.path.split(file_to_condense)[1])[0]
    out_fn = os.path.join(outdir, base_filename+"_condensed.fasta")

    # read in file into a dictionary
    d = st.fasta_to_dct(file_to_condense)
    seqids, seqs = [], []
    positions_to_pop = []

    for k, v in d.items():
        seqids.append(k)
        if not caseSensitive:
            seqs.append(v.upper())
        else:
            seqs.append(v)

    # step over each site finding sites which have no variability.
    for i in range(len(seqs[0])):  # columns of the sequence.
        pos_j_chars = []
        for j in range(len(seqs)):  # sequences within the file.
            pos_j_chars.append(seqs[j][i])
        if len(set(pos_j_chars)) == 1:
            positions_to_pop.append(i)

    # reverse the list, so we pop from the back
    positions_to_pop = positions_to_pop[::-1]

    for i in positions_to_pop:
        for seqid in d.keys():
            seq = d[seqid]
            mod_seq = seq[:i] + seq[(i+1):]
            d[seqid] = mod_seq

    seqs = list(d.values())
    print("After removing conserved sites, the alignment has {} sites.".format(len(seqs[0])))
    st.dct_to_fasta(d, out_fn)


def main(infile, outdir, case_sens):
    """
    takes in an aligned fasta file, and an output dir.
    Reduces the input fasta file to only colunns with differences.
    Writes that to file in the output directory - with the same filename as the input file, except with
    "_condensed.fasta" extension.
    Can ignore case.
    :param infile: the fasta file to use as source
    :param outdir: the place to write the output to.
    :param case_sens: bool indicates if the case of the characters in the sequences counts.
    :return: no return
    """

    # confirm that there are no duplicate names in the infile.
    seqids = []
    with open(infile, "r") as fh:
        for line in fh:
            if (line[0] == ">"):
                seqids.append(line[1:])
    if len(list(set(seqids))) != len(seqids):
        print("It appears that there are duplicate sequenceIDS in your fasta format file.\nNote that we cannot "
              "currently handle duplicate sequenceIDS.\nPlease ensure unique IDs and try again.\nNow exiting.")
        sys.exit(1)

    dct = st.fasta_to_dct(infile)
    # confirm that there is at least one sequence in the infile.
    if len(dct) < 1:
        print("It appears that there is no readable fasta format data in the infile specified.\nPlease confirm the "
              "input file and try again.\nNow exiting.")
        sys.exit(1)

    seqs = list(dct.values())
    print("Input alignment is {} sites long.".format(len(seqs[0])))

    # account for case sensitivity setting from user
    if not case_sens:
        tmp = []
        for s in seqs:
            tmp.append(s.upper())
        seqs = tmp
        del tmp

    # Check the sequences in the infile are the same length
    if len(set(list(map(len, seqs)))) != 1:
        print("Not all sequences in the input file are the same length.\nPlease make sure you specify an alignment with"
              "equally lengthed sequences.\nNow exiting")
        sys.exit(1)

    # finally, call condense
    condense(infile, outdir, case_sens)

    print("Completed. \nNow exiting")
    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='takes a fasta file, finds the "informative sites" (aka the sites that'
                                                 'have some variation) and spits them out into a fasta file with '
                                                 'almost the same name as the input file, but in the specified'
                                                 'output directory'
                                                 '')
    parser.add_argument('-in', '--infasta', type=str,
                        help='input fasta file to be reduced',
                        required=True)
    parser.add_argument('-out_dir', '--out_dir', type=str,
                        help='path to where the output condensed version of the fasta file will go.'
                             'It will have the same name as the input file -- with "_condensed.fasta" extension',
                        required=True)
    parser.add_argument('--caseSensitive', dest='caseSensitivity', action='store_true')
    parser.add_argument('--caseInSensitive', dest='caseSensitivity', action='store_false')
    parser.set_defaults(caseSensitivity=False)

    args = parser.parse_args()
    infasta = args.infasta
    out_dir = args.out_dir
    case_sens = args.caseSensitivity

    main(infasta, out_dir, case_sens)
