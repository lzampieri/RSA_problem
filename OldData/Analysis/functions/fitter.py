from math import sqrt, pi, log, inf
from struct import Struct
from scipy.optimize import curve_fit, minimize
import numpy as np
from matplotlib import pyplot as plt
from . import splt
from scipy.special import gamma, polygamma
from functools import partial

def chi_squared( real_ys, exp_ys ):
    sum = 0
    count = 0
    for ( real_y, exp_y ) in zip( real_ys, exp_ys ):
        if( real_y > 0 ):
            sum = sum + ( real_y - exp_y ) ** 2 / real_y
            count = count + 1
    return sum

def Gauss(x, A, mu, s):
        y = A / s / sqrt( 2 * pi ) *np.exp(-0.5 * ( ( x - mu ) / s )**2 )
        return y
    
def Gumbel(x, A, mu, s):
    alpha = 1.283 / s
    u = mu - 0.45 * s
    z = ( x - u ) * alpha
    return A * alpha * np.exp( - z - np.exp( - z ) )

def GG(x, A, mu, s, a ):
    b = 1 / s * sqrt( polygamma(1, a) )
    u = mu - 1 / b * ( log( a ) - polygamma(0, a) )
    pref = A * ( a**a ) * b / gamma( a )
    return pref * np.exp( - a * ( b * ( x - u ) + np.exp( - b * ( x - u ) ) ) )


def classic_fit( x, h, function, p0, bounds = ( -inf, inf ) ):
    try:
        pars, errs = curve_fit( function, x, h, p0, bounds= bounds )
    except (RuntimeError, RuntimeWarning) as e:
        pars = p0
    estimated_y = function( x, *pars )

    return {
        'est_y': estimated_y,
        'chi2' : chi_squared( h, estimated_y ),
        'pars' : pars
    }

def classic_fits( x, h, item, threshold = 0, plot = True ):
    x = np.array( x )
    h = np.array( h )

    mean = np.sum( x * h ) / np.sum( h )
    std = np.sqrt( np.sum( h * ( ( x - mean )**2 ) ) )
    step = np.min( x[1:] - x[:-1] )

    functions = [
        Gauss,
        Gumbel,
        GG
    ]

    p0s = [
        [ step, mean, std ],
        [ step, mean, std ],
        [ step, mean, std, 1.1 ]
    ]

    bounds = [
        ( [ 0, 0, 0 ], [ inf, 1, inf ]),
        ( [ 0, 0, 0 ], [ inf, 1, inf ]),
        ( [ 0, 0, 0, 0.01 ], [ inf, 1, inf, 100 ])
    ]

    x_tofit = x[ h > max(h) * threshold ]
    h_tofit = h[ h > max(h) * threshold ]

    fits = {}
    for f, p0, bound in zip( functions, p0s, bounds ):
        fits[ f.__name__ ] = classic_fit( x_tofit, h_tofit, f, p0, bound )

    if( plot ):
        plt.scatter( x, h, marker='+', color='k' )
        for key, res in fits.items():
            plt.plot( x_tofit, res['est_y'], label= fr"{key} - $\chi^2\sim{res['chi2']:.0e}$" )
        plt.legend()

        plt.title( fr"Side: {item['side']} - q = {item['defects_frac']} - $\gamma$ = {item['gamma']} - {item['dep_polymers']}" )
        plt.xlim( min( x_tofit ), max( x_tofit ) )
        plt.ylim( - max( h_tofit ) * 0.1, max( h_tofit ) * 1.1 )

    return fits

def logGauss(x, mu, s):
        y = np.log( 1.0 / s / sqrt( 2 * pi ) ) - 0.5 * ( ( x - mu ) / s )**2
        return y
    
def logGumbel(x, mu, s):
    alpha = 1.283 / s
    u = mu - 0.45 * s
    z = ( x - u ) * alpha
    return np.log( alpha ) - z - np.exp( - z )

def logGG(x, mu, s, a ):
    b = 1 / s * sqrt( polygamma(1, a) )
    u = mu - 1 / b * ( log( a ) - polygamma(0, a) )
    pref = ( a**a ) * b / gamma( a )
    return np.log( pref ) - a * ( b * ( x - u ) + np.exp( - b * ( x - u ) ) )

def entropy( x, func, param ):
    probs = func( x, *param ) 
    return np.average( - probs )

def entropyc_fit( x, function, p0, bound, x_to_est ):

    score = lambda pars: entropy( x, function, pars )
    
    results = minimize( score, x0 = p0, bounds=bound, tol = 1e-15 )

    estimated_y = np.exp( function( x_to_est, *results.x ) ) * np.min( x_to_est[1:] - x_to_est[:-1] )

    return {
        'est_y': estimated_y,
        'pars' : results.x
    }

def entropyc_fits( x, x_to_est ):
    x = np.array( x )
    x_to_est = np.array( x_to_est )

    mean = np.mean( x )
    std  = np.std ( x )

    functions = [
        logGauss,
        logGumbel,
        logGG
    ]

    p0s = [
        [ mean, std ],
        [ mean, std ],
        [ mean, std, 1.1 ]
    ]

    bounds = [
        ( ( 0, 1 ), ( 1e-9, inf ) ),
        ( ( 0, 1 ), ( 1e-9, inf ) ),
        ( ( 0, 1 ), ( 1e-9, inf ), ( 0.001, 100 ) )
    ]

    fits = {}
    for f, p0, bound in zip( functions, p0s, bounds ):
        fits[ f.__name__[3:] ] = entropyc_fit( x, f, p0, bound, x_to_est )

    return fits


def fits_twocols( data, threshold_1 = 0, threshold_2 = 0.5, data_edit = lambda x, y : ( x, y ), plot = True ):
    all_fits = []
    all_fits_2 = []

    if( plot ):
        splt.init( len( data ), 2 )


    for d in data:
        x, h = data_edit( d['h_x'], d['h'] )

        if( plot ):
            splt.next()
        all_fits  .append( fits( x, h, d, threshold_1, plot ) )
        if( plot ):
            splt.next()
        all_fits_2.append( fits( x, h, d, threshold_2, plot ) )

    return all_fits, all_fits_2

def rebin( x, h, fact ):
    if( fact == 1 ):
        return x, h
    step = np.min( x[1:] - x[:-1] )
    new_step = step * fact
    
    start_x = x[0] - step/2 + new_step / 2
    current_xi = 0
    current_x = lambda : start_x + current_xi * new_step
    current_h = 0
    new_x = []
    new_h = []

    for xi, hi in zip( x, h ):
        while( xi > current_x() + new_step / 2 ):
            if( current_h > 0 ):
                new_x.append( current_x() )
                new_h.append( current_h  )
            current_xi = current_xi + 1
            current_h  = 0
        current_h = current_h + hi

    if( len( new_x ) < 5 ):
        raise Exception("Rebinning too large!")

    # print( len( new_x ), len( new_h ), np.sum( h ), np.sum( new_h ) )

    return new_x, new_h

def rebin_width( x, h, width ):
    step = np.min( x[1:] - x[:-1] )
    return rebin( x, h, max( round( width / step ), 1 ) )

def rebin_n( x, h, n ):
    max_size = np.ptp( x )
    return rebin_width( x, h, max_size / n )

def rebinner_byfact( fact ):
    return lambda x, h: rebin( x, h, fact )

def rebinner_bywidth( width ):
    return lambda x, h: rebin_width( x, h, width )

def rebinner_byn( n ):
    return lambda x, h: rebin_n( x, h, n )
