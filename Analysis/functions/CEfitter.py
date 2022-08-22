from cmath import log
from re import L
import numpy as np
from uncertainties import *
from uncertainties import unumpy as unp
from functions import splt

def is_straight( logx, y, B, Bmax ):
    if( B >= Bmax ):
        return np.inf
    logy = np.log2( np.array( y ) - B )
    _, residuals, _, _, _ = np.polyfit( logx, logy, 1, full=True )
    return np.log( residuals[0] )

import matplotlib.pyplot as plt

def estimate_B( x, y, plot = False ):
    logx = np.log2( x )
    Bmin = 0
    Bmax = np.min( y )
    BabsMax = np.min( y )
    cont = 0

    while( Bmax - Bmin > 0.00000001 * BabsMax ):
        Bs = np.linspace( Bmin, Bmax, 500 )
        res = [ is_straight( logx, y, B, BabsMax ) for B in Bs ]
        if( plot ):
            plt.semilogx( Bs, res, '.b' )
        i = np.argmin( res )
        Bmin = Bs[ max( 0, np.min( i ) - 2 ) ]
        Bmax = Bs[ min( len( Bs ) - 1, np.max( i ) + 2 ) ]
        cont = cont + 1

    B = ufloat( np.mean( [Bmin, Bmax ] ), ( Bmax - Bmin ) / np.sqrt( 12 ) )
    if( plot ):
        plt.plot( B.n, is_straight( logx, y, B.n, BabsMax ), 'rs' )
    logyB = unp.uarray(
        np.log2( np.array( y ) - B.n ),
        B.s / np.abs( B.n - np.array( y ) ) / np.log( 2 )
    )


    return B, logx, logyB

def estimate_B_log( x, y ):
    logx = np.log2( x )
    Bmin = 1e-15
    Bmax = np.min( y )
    BabsMax = np.min( y )
    cont = 0

    while( np.abs( np.log10( Bmax ) - np.log10( Bmin ) ) > 0.001 ):
        Bs = np.logspace( np.log10( Bmin ), np.log10( Bmax ), 500 )
        res = [ is_straight( logx, y, B, BabsMax ) for B in Bs ]
        i = np.argmin( res )
        Bmin = Bs[ max( 0, np.min( i ) - 2 ) ]
        Bmax = Bs[ min( len( Bs ) - 1, np.max( i ) + 2 ) ]
        cont = cont + 1

    B = ufloat( np.mean( [Bmin, Bmax ] ), ( Bmax - Bmin ) / np.sqrt( 12 ) )
    logyB = unp.uarray(
        np.log2( np.array( y ) - B.n ),
        B.s / np.abs( B.n - np.array( y ) ) / np.log( 2 )
    )

    return B, logx, logyB

def estimate_B_errony( x, y, plot = False ):
    logx = np.log2( x )
    y_n = unp.nominal_values( y )
    Bmin = 0
    Bmax = np.min( y_n )
    BabsMax = np.min( y_n )
    cont = 0

    while( True ):
        Bs = np.linspace( Bmin, Bmax, 500 )
        res = np.array( [ is_straight( logx, y_n, B, BabsMax ) for B in Bs ] )
        res[ res == np.inf ] = np.max( res[ res < np.inf ] ) # remove infinite

        i = np.where( res < np.max( res ) - 0.9 * np.ptp( res ) )[0]
        if( np.ptp( res ) == 0 ):
            break

        imin = np.min( i )
        imax = np.max( i )

        imin = np.max( [0, np.min( [ imin, imax - 3 ]) ] )
        imax = np.min( [len( Bs ) - 1, np.max( [ imax, imin + 3 ]) ] )

        Bmin = Bs[ imin ]
        Bmax = Bs[ imax ]

        cont = cont + 1

        if( cont > 3 ): # Make at least 3 steps
            x_simplified = np.arange( len( Bs ) )
            p = np.polyfit( x_simplified, res, 2 )
            res_expected = p[0] * ( x_simplified**2 ) + p[1] * x_simplified + p[2]
            
            if( np.std( res - res_expected ) * 2 > np.ptp( res_expected ) ):
                break
            if( cont > 30 ):
                break

        old_bs = Bs
        old_res = res


    if( plot ):
        print( cont )
        plt.plot( x_simplified, old_res, label = cont - 1 )
        plt.plot( x_simplified, res, label = cont )
        plt.plot( x_simplified, res_expected )
        plt.legend()
        plt.ylim( np.min( old_res ), np.max( old_res ) )

    B = ufloat( np.mean( [Bmin, Bmax ] ), ( Bmax - Bmin ) )
    logyB = unp.uarray(
        np.log2( np.array( y_n ) - B.n ),
        np.sqrt( B.s **2 + unp.std_devs( y ) **2 ) / np.abs( B.n - np.array( y_n ) ) / np.log( 2 )
    )

    return B, logx, logyB

def estimate_oneovernu( logx, logyB ):
    res, cov = np.polyfit(
        logx, unp.nominal_values( logyB ), 1,
        cov = True, w = 1 / unp.std_devs( logyB )
        )
    return ufloat( - res[0], np.sqrt( cov[0,0] ) ), ufloat( res[1], np.sqrt( cov[1,1] ) )

def estimate_oneovernu_bis( logx, logyB, errmin = 0 ):
    resc = np.polyfit(
        logx, unp.nominal_values( logyB ), 1
        )
    resp, cov = np.polyfit(
        logx, unp.nominal_values( logyB ) + unp.std_devs( logyB ), 1, cov = True
        )
    resm = np.polyfit(
        logx, unp.nominal_values( logyB ) - unp.std_devs( logyB ), 1
        )
    res = np.mean( [ resc, resp, resm ], axis = 0 )
    errs_ptp = np.ptp( [ resc, resp, resm ], axis = 0 ) / np.sqrt( 12 )
    errs_cov = np.sqrt( np.diag( cov ) )
    errs_errmin = errs_cov * 0 + errmin
    errs = np.max( [ errs_ptp, errs_cov, errs_errmin ], axis = 0 )
    return ufloat( - res[0], errs[0] ), ufloat( res[1], errs[1] )    