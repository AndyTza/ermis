"""

- Read your .db catalog
- Define several funcitons that will:
    - for a given (ra, dec): calculate the pointings per filter
    - for a given (ra, dec): calculate the theoretical E(B-V) using __
    - Return the Delta_t cadence per specified filter
"""
from astropy.io import ascii
from astropy.time import Time
import numpy as np
from astropy.coordinates import SkyCoord
import generator as generator
import models as phot_models
import astro_tools
import dummy_functions as dummy


'''
Sampler will consist of a fetching funciton that will do a quick query to to the
epoch matching. Can also enable modes of searching
'''



class simulation:
    def __init__(self, table):
        self.table = table
    
