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


# LSST constants from https://smtn-002.lsst.io/
lsst_bands = ['u', 'g', 'r', 'i', 'z', 'y']
m5_limiting = [23.87, 24.82, 24.36, 23.93, 23.36, 22.47]
filter_zeropoint = [27.03, 28.38, 28.16, 27.85, 27.46, 26.68]
