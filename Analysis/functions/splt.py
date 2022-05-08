from turtle import color
import numpy as np
import matplotlib.pyplot as plt
from . import splt

_plot_params = 0
def init( numrows, numcols, axis = False ):
    plt.figure( figsize = ( numcols * 6, numrows * 4  ) )
    global _plot_params
    _plot_params = [ numrows, numcols, 0 ]
    if( axis ):
        next()
def next():
    _plot_params[2] += 1
    plt.subplot( *_plot_params )
def goto( i_row, i_col ):
    _plot_params[2] = i_row * _plot_params[1] + i_col + 1
    plt.subplot( *_plot_params )

def iterate( count, func_row, func_col, func_leg, func_x, func_y, params, each_plot ):

    rows = set()
    cols = set()
    legs = set()
    for d in range(count):
        rows.add( func_row( d ) )
        cols.add( func_col( d ) )
        legs.add( func_leg( d ) )

    splt.init( len( rows ), len( cols ) )
    color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for r in rows:
        for c in cols:
            splt.next()
            for i_l, l in enumerate( legs ):
                temp_x = []
                temp_y = []
                for d in range(count):
                    if(
                        func_col( d ) == c and
                        func_row( d ) == r and
                        func_leg( d ) == l
                    ):
                        temp_x.append( func_x( d ) )
                        temp_y.append( func_y( d ) )
                
                sort_idx = np.argsort( temp_x )
                x = np.array( temp_x )[ sort_idx ]
                y = np.array( temp_y )[ sort_idx ]
                
                if( len( y ) > 0 ):
                    for i in range( len( y[0] ) ):
                        plt.plot( x, [ yi[i] for yi in y ], color = color_list[ i_l ], **params( r, c, l )[i] )

                    plt.legend()
                
                each_plot( r, c )