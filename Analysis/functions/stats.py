import imp
from scipy.stats import kurtosis

def the_kurtosis( data ):
    kurt = []
    for d in data:
        kurt.append( kurtosis( d['chunks'] ) )
    return kurt