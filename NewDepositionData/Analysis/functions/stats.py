from scipy.stats import kurtosis
import numpy as np
from uncertainties import ufloat
import matplotlib.pyplot as plt

def the_kurtosis( data ):
    kurt = []
    for d in data:
        kurt.append( kurtosis( d['chunks'] ) )
    return kurt

def G( thearray ):
    thearray = np.array( thearray )
    n = len( thearray )
    avg = np.mean( thearray )
    sum3 = np.sum( ( thearray - avg ) ** 3 ) / n
    sum232 = ( np.sum( ( thearray - avg ) ** 2 ) / n ) ** 1.5
    return np.sqrt( n * ( n - 1 ) ) / ( n - 2 ) * sum3 / sum232

def compute_in_decades( thearray, func, method = 'auto', plot = False, hist = False ):
    num_sectors = 10
    tot_len = len( thearray )
    extremities = np.linspace( 0, tot_len, num_sectors, dtype = int )
    thearray = np.array( thearray )

    results = []
    for i in range( num_sectors - 1 ):
        # if( i == 0 ):
        #     data = thearray[ extremities[1]: ]
        # elif( i == num_sectors - 1 ):
        #     data = thearray[ :extremities[num_sectors - 2] ]
        # else:
        #     data = np.hstack([ thearray[ :extremities[i] ], thearray[ extremities[i+1]: ] ])
            
        data = thearray[ extremities[i] : extremities[i+1] ]
        
        results.append( func( data ) )
    
    result = func( thearray )

    if( method == 'ptp' ):
        std = np.ptp( results ) / 2
    elif( method == 'uniform' ):
        std = np.ptp( results ) / np.sqrt( 12 )
    elif( method == 'std' ):
        std = np.std( results, ddof = 1 )
    elif( method == 'auto' ):
        std = np.max( [ np.std( results, ddof = 1 ), np.ptp( results ) / np.sqrt( 12 ) ] )
    else:
        raise ValueError( f"{method} is not a recognized method. Use ptp, uniform, std or auto" )
    
    if( plot ):
        plt.plot( np.arange( len( results ) ), results, 'o' )
        plt.errorbar( [ len( results ) / 2 ], [ result ], [ std ]  )

    if( hist ):
        plt.hist( results, 20 )
        plt.plot( [ result - std, result + std], np.array( [ 1, 1 ] ) * np.mean( plt.ylim() ) )

    return ufloat( result, std )

def get_pvalue_string( value, expected ):
    z = np.abs( value.n - expected ) / value.s
    if( z > 3.7 ):
        return '⁎⁎⁎⁎'
    if( z > 3.09 ):
        return '⁎⁎⁎'
    if( z > 2.33 ):
        return '⁎⁎'
    if( z > 1.65 ):
        return '⁎'
    return 'ns'