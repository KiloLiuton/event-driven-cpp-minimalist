import numpy as np
from tqdm import tqdm

N = 20000
K = 100

neighbor_list = []
indexes = []
kernel_sizes = []

def get_neighbors( n ):
    idx = indexes[ n ]
    res = []
    for i in range( 2*K ):
        res.append( neighbor_list[ idx + i ] )
    return res

for i in tqdm( range( N ) ):
    indexes.append( i * 2*K )
    kernel_sizes.append( 2 * K )
    for k in range( -K, K + 1 ):
        if k == 0:
            continue
        n = i + k
        if n < 0:
            n += N
        if n >= N:
            n -= N
        neighbor_list.append( n )

with open( 'topology.txt', 'w' ) as f:

    n = 0
    f.write( 'neighbor_list\n' )
    for i in tqdm( range( len( neighbor_list ) ) ):
        s = str( neighbor_list[ i ] )
        f.write( s )
        n += len( s )
        if i < len( neighbor_list ) - 1:
            f.write( ',' )
            n += 1
        if n > 80:
            n = 0
            f.write( '\n' )
    f.write( '\n\n' )
    print( 'neighbor_list done.' )

    n = 0
    f.write( 'indexes\n' )
    for i in tqdm( range( len( indexes ) ) ):
        s = str( indexes[ i ] )
        f.write( s )
        n += len( s )
        if i < len( indexes ) - 1:
            f.write( ',' )
            n += 1
        if n > 80:
            n = 0
            f.write( '\n' )
    f.write( '\n\n' )
    print( 'indexes done.' )

    n = 0
    f.write( 'kernel_sizes\n' )
    for i in tqdm( range( len( kernel_sizes ) ) ):
        s = str( kernel_sizes[ i ] )
        f.write( s )
        n += len( s )
        if i < len( kernel_sizes ) - 1:
            f.write( ',' )
            n += 1
        if n > 80:
            n = 0
            f.write( '\n' )
    f.write( '\n' )
    print( 'kernel_sizes done.' )
