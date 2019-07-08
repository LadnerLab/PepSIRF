#!/usr/bin/env python3
import argparse
import itertools
import math
import pandas as pd
from pprint import pprint



def main():
    argp = argparse.ArgumentParser( description = "Count index pairs for each sample" )

    argp.add_argument( '--r1', type = str, help = "Filename of R1 reads" )
    argp.add_argument( '--i1', type = str, help = "Filename of I1 reads" )
    argp.add_argument( '--r1_idx_loc', type = str, help = "Comma-separated location of R1 index" )
    argp.add_argument( '--i1_idx_loc', type = str, help = "Comma-separated location of I1 index" )
    argp.add_argument( '--mapping', type = str, help = "Map of Sample names to sequences." )
    argp.add_argument( '--min_score', type = int, help = "Minimum score for an item to be output", default = 0 )

    args = argp.parse_args()

    samples_df        = parse_snames( args.mapping )
    concatemer_sample, idx_seqs = df_2_sdict( samples_df )

    sample_names = set( samples_df[ 'Sample' ] )

    r1_idx = to_tuple( args.r1_idx_loc )
    i1_idx = to_tuple( args.i1_idx_loc )

    r1_file = open( args.r1, 'r' )
    i1_file = open( args.i1, 'r' )

    n = 100000

    r1_lines = read_n_lines_fasta( r1_file, n )
    i1_lines = read_n_lines_fasta( i1_file, n )

    while r1_lines and i1_lines:
        for r1, i1 in zip( r1_lines, i1_lines ):
            f_seq = get_seq( r1, r1_idx )
            r_seq = get_seq( i1, i1_idx )
            concat = f_seq + r_seq

            if concat not in concatemer_sample:
                concatemer_sample[ concat ] = [ 'No Sample', 0, idx_seqs.get( f_seq, '' ), idx_seqs.get( r_seq, '' ) ]
            concatemer_sample[ concat ][ 1 ] += 1

        r1_lines = read_n_lines_fasta( r1_file, n )
        i1_lines = read_n_lines_fasta( i1_file, n )


    r1_file.close()
    i1_file.close()

    spl = lambda x: [ x[ 0: r1_idx[ 1 ] ], x[ r1_idx[ 1 ]:: ] ]
    st  = lambda x: [ str( item ) for item in x ]
    fun = lambda x: '\t'.join( list( st( [
                                         x[ 3 ],
                                         spl( x[ 0 ] )[ 0 ], 
                                         x[ 4 ],
                                         spl( x[ 0 ] )[ 1 ],
                                         x[ 1 ],
                                         x[ 2 ]
                                          ]
                                       )
                                   )
                             )

    print( 'Forward Name\tForward Read\tReverse Name\tReverse Read\tSample\tCount' )
    for key, val in concatemer_sample.items():
        if val[ 1 ] >= args.min_score or val[ 0 ] in sample_names:
            print( fun( [ key ] + val ) )


def get_seq( seq, tup ):
    return seq[ tup[ 0 ] : tup[ 0 ] + tup[ 1 ] ]

    
def read_n_lines_fasta( fh, n ):
    count = 0
    out = list()                              

    while count < n:
        try:
            out.append( read_line_fasta( fh ) )
            count += 1
        except StopIteration:
            return out
    return out

                        

def read_line_fasta( fh ):
    for i in range( 4 ):
        line = next( fh )
        if not ( i + 3 ) % 4:
            out = line.strip()
    return out

    
def to_tuple( string ):
    l = [ int( item ) for item in string.split( ',' ) ]
    return l

def df_2_sdict( df ):
    out_dict = {}
    out_idx  = {}

    reverse_d = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' } 
    reverse = lambda y: ''.join( [ reverse_d[ i ] for i in reversed( y ) ] )

    for idx, row in df.iterrows():
        if not str( row[ 'Sample' ] ).lower() == 'nan':
            rev = reverse( row[ 'Reverse Sequence' ] )
            concat = row[ 'Forward Sequence' ] + rev
            out_dict[ concat ] = [ row[ 'Sample' ], 0,
                                   row[ 'Forward Name' ], row[ 'Reverse Name' ]
                                 ]
            out_idx[ row[ 'Forward Sequence' ] ] = row[ 'Forward Name' ]
            out_idx[ row[ 'Reverse Sequence' ] ] = row[ 'Reverse Name' ]

    return out_dict,out_idx

def parse_snames( fname ):
    csv = pd.read_csv( fname )

    wanted_cols = [ 'Sample', 'Forward Name', 'Forward Sequence',
                    'Reverse Name', 'Reverse Sequence'
                  ]
    return csv[ wanted_cols ]

if __name__ == '__main__':
    main()
