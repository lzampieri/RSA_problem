from cmath import log
from re import L
import numpy as np
from uncertainties import *
from uncertainties import unumpy as unp

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

    while( Bmax - Bmin > 0.00000001 * BabsMax ):
        Bs = np.linspace( Bmin, Bmax, 500 )
        res = [ is_straight( logx, y_n, B, BabsMax ) for B in Bs ]
        if( plot ):
            plt.semilogx( Bs, res, '.b' )
        i = np.argmin( res )
        Bmin = Bs[ max( 0, np.min( i ) - 2 ) ]
        Bmax = Bs[ min( len( Bs ) - 1, np.max( i ) + 2 ) ]
        cont = cont + 1

    B = ufloat( np.mean( [Bmin, Bmax ] ), ( Bmax - Bmin ) / np.sqrt( 12 ) )
    if( plot ):
        plt.plot( B.n, is_straight( logx, y_n, B.n, BabsMax ), 'rs' )
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

def estimate_oneovernu_bis( logx, logyB ):
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
    errs_errmin = errs_cov * 0
    errs = np.max( [ errs_ptp, errs_cov, errs_errmin ], axis = 0 )
    return ufloat( - res[0], errs[0] ), ufloat( res[1], errs[1] )    