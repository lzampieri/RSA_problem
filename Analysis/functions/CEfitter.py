import numpy as np
from uncertainties import *

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

    print( f"{cont} iter." )
    B = np.mean( [ Bmin, Bmax ] )
    return B, logx, np.log2( np.array( y ) - B )

def estimate_oneovernu( logx, logyB ):
    res, cov = np.polyfit( logx, logyB, 1, cov = True )
    return ufloat( - res[0], np.sqrt( cov[0,0] ) ), ufloat( res[1], np.sqrt( cov[1,1] ) )