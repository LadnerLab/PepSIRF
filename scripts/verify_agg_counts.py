#!/usr/bin/env python3

import sys
def main():
    if len( sys.argv ) != 3:
        print( "USAGE: ", sys.argv[ 0 ], "regular_infile agg_infile" )
        sys.exit( 1 )

    regular_in = sys.argv[ 1 ]
    agg_in     = sys.argv[ 2 ]

    reg_d = parse_counts( regular_in )
    agg_d = parse_counts( agg_in )

    for k, v in reg_d.items():
        assert( v == agg_d[ k ] )


def parse_counts( fname, trimchar = '-' ):
    odict = {}
    with open( fname, 'r' ) as of:
        for lineno, line in enumerate( of ):
            if lineno:
                tline = line.strip().split( '\t' )
                spline = tline[ 0 ].split( trimchar )[ 0 ]
                
                numline = [ int( item ) for item in tline[ 1:: ] ]

                if spline not in odict:
                    odict[ spline ] = [ 0 ] * len( numline )
                for idx, num in enumerate( numline ):
                    odict[ spline ][ idx ] += num
    return odict

if __name__ == '__main__':
    main()
