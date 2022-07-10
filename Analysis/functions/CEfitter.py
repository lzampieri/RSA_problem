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

def estimate_B( x, y ):
    logx = np.log2( x )
    Bmin = 0
    Bmax = np.min( y )
    BabsMax = np.min( y )
    cont = 0

    while( Bmax - Bmin > 0.0000001 * BabsMax ):
        Bs = np.linspace( Bmin, Bmax, 500 )
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

def estimate_oneovernu( logx, logyB ):
    res, cov = np.polyfit(
        logx, unp.nominal_values( logyB ), 1,
        cov = True, w = 1 / unp.std_devs( logyB )
        )
    return ufloat( - res[0], np.sqrt( cov[0,0] ) ), ufloat( res[1], np.sqrt( cov[1,1] ) )

def estimate_oneovernu_bis( logx, logyB ):
    res = np.polyfit(
        logx, unp.nominal_values( logyB ), 1
        )
    resp = np.polyfit(
        logx, unp.nominal_values( logyB ) + unp.std_devs( logyB ), 1
        )
    resm = np.polyfit(
        logx, unp.nominal_values( logyB ) - unp.std_devs( logyB ), 1
        )
    errs = np.ptp( [ res, resp, resm ], axis = 0 )
    return ufloat( - res[0], errs[0] ), ufloat( res[1], errs[1] )    