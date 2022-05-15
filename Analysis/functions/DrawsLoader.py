import numpy as np
import pandas as pd
import json
from pathlib import Path
import os
from glob import glob

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

    for file in glob("../" + path + "/**/draw_*.txt", recursive=True):

        item = {}
        d = Path( os.path.dirname( file ) )

        # Load data
        item.update( load_file( d / 'details.txt'     ) )
        item.update( { 'draw': pd.read_csv( file, sep='\t', names=['x', 'y', 'v'], index_col = False )} );

        data.append( item )

    return np.array( data )

