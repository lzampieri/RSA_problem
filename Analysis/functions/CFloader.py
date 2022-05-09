import numpy as np
import pandas as pd
import json
from pathlib import Path
import os
from glob import glob

from sqlalchemy import false

# Load single file
def load_file( filename ):
    text = filename.read_text()
    text = text.replace("nan","0")
    text = text.replace("False", "false" )
    text = text.replace("True", "true" )
    text = text.replace(",]", "]" )
    return json.loads( text )

def load_data( path ):
    data = []

    for file in glob("../" + path + "/**/details.txt", recursive=True):

        item = {}
        d = Path( os.path.dirname( file ) )

        # Verify that everything exists
        if( not (
            ( d / 'details.txt').exists() and
            ( d / 'CF_avg.txt').exists() ) ):
            continue

        # Load data
        item.update( load_file( d / 'details.txt'     ) )
        item.update( { 'CF_D': pd.read_csv( d / 'CF_avg.txt', sep='\t', names=['x', 'y'], index_col=False ) } );

        data.append( item )

    return np.array( data )

