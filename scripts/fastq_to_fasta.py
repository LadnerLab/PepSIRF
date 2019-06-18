#!/usr/bin/env python3
import sys

def main():
    if len( sys.argv ) != 2:
        print( "USAGE: fastq_to_fasta.py in_file" )
        sys.exit( 1 )

    in_file = sys.argv[ 1 ]
    out_fname = '.'.join( in_file.split( '.' )[ :-1: ] ) + ".fasta"

    n, s = list(), list()

    with open( in_file, 'r' ) as open_file:
        for lineno, line in enumerate( open_file ):
            if not( lineno % 4 ):
                n.append( line.strip().split( '@' )[ 1 ] )
            elif not( ( lineno + 3 ) % 4 ):
                s.append( line.strip() )

    assert( len( n ) == len( s ) )
    assert( not( lineno + 1 ) == 0 )

    write_fasta( out_fname, n, s )

def write_fasta( fname, n, s ):
    with open( fname, 'w' ) as ofile:
        for name, seq in zip( n, s ):
            ofile.write( '>%s\n%s\n' % ( name, seq ) )

if __name__ == '__main__':
    main()
