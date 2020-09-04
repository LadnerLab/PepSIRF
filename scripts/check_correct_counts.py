#!/usr/bin/env python3
import sys
import protein_oligo_library as oligo


def main():
    if len( sys.argv ) != 4:
        print( "Verify correctness of exact-match output.")
        print( "USAGE: check_correct_counts.py designed_fasta pepsirf_demux_output reads_fastq")
        sys.exit( 1 )

    design_lib   = sys.argv[ 1 ]
    demux_output = sys.argv[ 2 ]
    reads_fastq  = sys.argv[ 3 ]

    start_index = 43
    length      = 90


    demux_counts = get_counts( demux_output )
    d_n, d_s     = oligo.read_fasta_lists( design_lib )
    fastq_counts = get_fastq_counts( reads_fastq, start_index, length )

    seq_dict = {}
    for name, seq in zip( d_n, d_s ):
        seq_dict[ name ] = seq

    for name, count in demux_counts.items():
        if seq_dict[ name ] in fastq_counts:
            try:
                assert( count == fastq_counts[ seq_dict[ name ] ] )
            except AssertionError as e:
                print( seq_dict[ name ], count, fastq_counts[ seq_dict[ name ] ] )

def get_counts( filename ):
    out_dict = {}
    to_str = lambda x: [ int( item ) for item in x ]
    with open( filename, 'r' ) as open_file:
        for lineno, line in enumerate( open_file ):
            if lineno:
                split_str = line.strip().split( ',' )
                out_dict[ split_str[ 0 ] ] = sum( to_str( split_str[ 1::] ) )
    return out_dict

def get_fastq_counts( reads_file, start_idx, length ):
    read_dict = {}
    with open( reads_file, 'r' ) as open_file:
        for lineno, line in enumerate( open_file ):
            if not ( ( lineno+3 ) % 4 ):
                substr = line.strip()[ start_idx: start_idx + length ]
                if substr not in read_dict:
                    read_dict[ substr ] = 0
                read_dict[ substr ] += 1
    return read_dict



if __name__ == '__main__':
    main()
