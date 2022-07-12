from scipy.stats import kurtosis
import numpy as np
from uncertainties import ufloat

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

def compute_in_decades( thearray, func, method = 'auto' ):
    tot_len = len( thearray )
    extremities = np.linspace( 0, tot_len, 10, dtype = int )
    thearray = np.array( thearray )

    results = []
    for i in range( 10 ):
        if( i == 0 ):
            data = thearray[ extremities[1]: ]
        elif( i == 9 ):
            data = thearray[ :extremities[9] ]
        else:
            data = thearray[ np.r_[ :extremities[i], extremities[i+1]: ] ]
        results.append( func( data ) )
    
    if( method == 'ptp' ):
        return ufloat( func( thearray ), np.ptp( results ) )
    elif( method == 'uniform' ):
        return ufloat( func( thearray ), np.ptp( results ) / np.sqrt( 12 ) )
    elif( method == 'std' ):
        return ufloat( func( thearray ), np.std( results ) )
    elif( method == 'auto' ):
        return ufloat( func( thearray ), np.max( [ np.std( results ), np.ptp( results ) / np.sqrt( 12 ) ] ) )
    else:
        raise ValueError( f"{method} is not a recognized method. Use ptp, uniform, std or auto" )