import os
from smallBixTools import smallBixTools as st
import argparse
import sys


def condense(file_to_condense, outdir):
    """
    takes a filename of a fasta file.
    removes all conserved sites
    writes out the reduced file, containing no conserved sites.
    :param d: a fasta file path.
    :param outdir: the place to write the output file to
    :return: no return
    """

    base_filename = os.path.splitext(os.path.split(file_to_condense)[1])[0]
    out_fn = os.path.join(outdir, base_filename+"_condensed.fasta")

    d = st.fasta_to_dct(file_to_condense)
    seqids, seqs = [], []
    positions_to_pop = []
    for k, v in d.items():
        seqids.append(k)
        seqs.append(v)
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

    st.dct_to_fasta(d, out_fn)


def main(infile, outdir):
    """
    takes in a fasta file, and an output dir.
    Reduces the input fasta file to only colunns with differences.
    Writes that to file in the output directory - with the same filename as the input file, except with
    "_condensed.fasta" extension.
    :param infile: the fasta file to use as source
    :param outdir: the place to write the output to.
    :return: no return
    """
    dct = st.fasta_to_dct(infile)
    # Check the sequences in the infile are the same length
    seqs = list(dct.values())
    if len(set(list(map(len, seqs)))) != 1:
        print("Not all sequences in the input file are the same length. Please make sure you specify an alignment with"
              "equally lengthed sequences. Now exiting")
        sys.exit(1)

    condense(infile, outdir)

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

    args = parser.parse_args()
    infasta = args.infasta
    out_dir = args.out_dir

    main(infasta, out_dir)
