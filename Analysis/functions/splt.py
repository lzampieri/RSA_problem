from turtle import color
import numpy as np
import matplotlib.pyplot as plt
from . import splt, stats

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
def hline( y, format = ':k', xl = [] ):
    if( len( xl ) == 0 ):
        xl = plt.xlim()
    plt.plot( xl, [y, y], format )
    plt.xlim( xl )
def vline( x, format = ':k', yl = [] ):
    if( len( yl ) == 0 ):
        yl = plt.ylim()
    plt.plot( [x, x], yl, format )
    plt.ylim( yl )

def iterate( count, func_row, func_col, func_leg, func_x, func_y, params, each_plot ):

    rows = np.unique( [ func_row( i ) for i in range(count) ] )
    cols = np.unique( [ func_col( i ) for i in range(count) ] )
    legs = np.unique( [ func_leg( i ) for i in range(count) ] )

    splt.init( len( rows ), len( cols ) )
    color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for r in rows:
        for c in cols:
            splt.next()
            all_x = set()
            for i_l, l in enumerate( legs ):
                temp_x = []
                temp_y = []
                for d in range(count):
                    if(
                        func_row( d ) == r and
                        func_col( d ) == c and
                        func_leg( d ) == l
                    ):
                        temp_x.append( func_x( d ) )
                        temp_y.append( func_y( d ) )
                        all_x.add( func_x( d ) )
                
                sort_idx = np.argsort( temp_x )
                x = np.array( temp_x )[ sort_idx ]
                y = np.array( temp_y )[ sort_idx ]
                
                if( len( y ) > 0 ):
                    for i in range( len( y[0] ) ):
                        plt.plot( x, [ yi[i] for yi in y ], color = color_list[ i_l ], **params( r, c, l )[i] )

                    plt.legend()
                
            plt.xticks( list( all_x ) )
            each_plot( r, c )

def oneitem_iterate( count, func_item, func_leg, func_x, func_y, params, each_plot ):

    items = np.unique( [ func_item( i ) for i in range(count) ] )
    n_items = len( items )
    legs = np.unique( [ func_leg( i ) for i in range(count) ] )

    if( n_items < 4 ):
        splt.init( n_items, 1 )
    elif( n_items % 3 == 0 ):
        splt.init( n_items / 3, 3 )
    elif( n_items % 2 == 0 ):
        splt.init( n_items / 2, 2 )
    else:
        splt.init( n_items, 3 )

    color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for it in items:
        splt.next()
        all_x = set()
        for i_l, l in enumerate( legs ):
            temp_x = []
            temp_y = []
            for d in range(count):
                if(
                    func_item( d ) == it and
                    func_leg( d ) == l
                ):
                    temp_x.append( func_x( d ) )
                    temp_y.append( func_y( d ) )
                    all_x.add( func_x( d ) )
            
            sort_idx = np.argsort( temp_x )
            x = np.array( temp_x )[ sort_idx ]
            y = np.array( temp_y )[ sort_idx ]
            
            if( len( y ) > 0 ):
                for i in range( len( y[0] ) ):
                    plt.plot( x, [ yi[i] for yi in y ], color = color_list[ i_l ], **params( it, l )[i] )

                plt.legend()
            
        plt.xticks( list( all_x ) )
        each_plot( it )

def oneitem_iterate_errorbar(
        count, func_item, func_leg, func_x, func_y, params, each_plot,
        pvals_ys = [], pvals_exp = 0
    ):

    items = np.unique( [ func_item( i ) for i in range(count) ] )
    n_items = len( items )
    legs = np.unique( [ func_leg( i ) for i in range(count) ] )

    if( n_items < 4 ):
        splt.init( n_items, 1 )
    elif( n_items % 3 == 0 ):
        splt.init( n_items / 3, 3 )
    elif( n_items % 2 == 0 ):
        splt.init( n_items / 2, 2 )
    else:
        splt.init( n_items, 3 )

    color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for it in items:
        splt.next()
        all_x = set()
        for i_l, l in enumerate( legs ):
            temp_x = []
            temp_y = []
            for d in range(count):
                if(
                    func_item( d ) == it and
                    func_leg( d ) == l
                ):
                    temp_x.append( func_x( d ) )
                    temp_y.append( func_y( d ) )
                    all_x.add( func_x( d ) )
            
            sort_idx = np.argsort( temp_x )
            x = np.array( temp_x )[ sort_idx ]
            y = np.array( temp_y )[ sort_idx ]
            
            if( len( y ) > 0 ):
                for i in range( len( y[0] ) ):
                    plt.errorbar(
                        x, [ yi[i].n for yi in y ], [ yi[i].s for yi in y ],
                        color = color_list[ i_l ], **params( it, l )[i] )

                    if( len( pvals_ys ) == len( legs ) ):
                        for the_x, the_y in zip( x, y ):
                            plt.annotate(
                                stats.get_pvalue_string( the_y[i], pvals_exp ),
                                ( the_x, pvals_ys[ i_l ] ),
                                ha='center',
                                color = color_list[ i_l ]
                            )

                plt.legend( loc='lower right' )
            
        plt.xticks( list( all_x ) )
        each_plot( it )