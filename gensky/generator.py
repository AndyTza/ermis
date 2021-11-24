"""
- Generate light curves for a specified model
- Using the sampling algorithm at hand
- Specify the epochs range
- Specify the number of light curves per pointing
"""

from astropy.table import Table
from astropy.io import ascii
import numpy as np
from astropy.time import Time
import warnings
import os
import sys
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u

if not sys.warnoptions:
    warnings.simplefilter("ignore")

## This will be the LSST generator code. It will by default load the table.


class simulation:
    """ we want to create a class that will load the lsst table"""

    def load_table(filename="sim_table.ascii"):
        """ Unpack the data of the whole simulator. Default data path should be in: ../data/sim_data.

        Current output will give back the: ra, dec, filter, mjd of observation, sky background, and m5 limiting magnitude

        Input
        -----
        filename (str): Default name of the catalog should be: sim_table.ascii

        Output
        ------
        astropy.table (astropy.table): table containing all information about the survey simulation data

        """

        file_path = f"data/sim_data/{filename}"
        # Head directory path
        dir_path = os.path.dirname(os.path.realpath('.'))

        full_path = os.path.join(dir_path, file_path)

        # Simulation table
        sim_table = ascii.read(full_path)

        return sim_table

""" A series of functions that can help us fetch and query the simulation table """

def fetch_cadence_info(table, ra, dec, filter_band='all' , fov=3.5):
    """
    Fetch the table information for a given poinitng and filter. Retuns data from all epochs.

    Input
    -----
    table (astropy.table): Astropy table of simulation data: ra, dec, filter, mjd, 5sigma, skyB
    ra, dec (float, float): Coordinates of your field of interest
    filter_band (str): Supports either specific LSST-bandpass (i.e u,g,r,i,z,y) or all; that returns all observations in that band
    fov (float): Field of view of the average LSST pointing (approximately ~9.8 sqr deg)

    Output
    ------
    table_new (astropy.table): new table under the query pointing and filter conditions
    """

    # Convert pointing to astropy SkyCoords
    target_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    table_coordinates = SkyCoord(ra=table['ra']*u.deg, dec=table['dec']*u.deg)
    delta_sep = target_coords.separation(table_coordinates).degree # seperation in degrees

    # The maximum seperation our coordinates should have from the central pointing
    max_seperation = fov/2

    if filter_band=='all':
        return table[delta_sep<=max_seperation]
    else:
        return table[np.where((delta_sep<=max_seperation) & (table['filter']==filter_band))]

def simple_mjd_sampler(table, time_separation=365, mode='random'):
    """ Given a range of dates, this funciton will return the cadence for a given field under a time-seperation.

        Input
        -----
        table (astropy.table): simulation data (ideally restricted to specific field and bands)
        time_separation (int): the time you want to query data
        mode (str): How each time query will be performed:
                    random: will randomly select a range of dates that satisfied max(T)-min(T)=time_separation
                    start_end: will select a range of dates starting from the first pointing of the simulation

        Output
        ------
        table (astropy.table): simulation data at specified time mode
    """

    if mode=='start_end':
        # Choose all samples mjd times that START from the begning of the simulation (per pointing) and add time_separation

        start_time = Time(min(table['mjd']), format='mjd') # minimum time that it has been observed

        end_time_isot = Time(start_time.isot, format='isot') + time_separation # and add the seperated time

        end_time = Time(end_time_isot, format='isot').mjd # convert to mjd

        return table[np.where((table['mjd']>=start_time.mjd) & (table['mjd']<=end_time))]

    elif mode=='random':
        # All the mjd times observed for the table (ideally this should be per pointing)
        mjd_times = Time(table['mjd'], format='mjd') # survey times

        # From these times, randomly select some MJD start (assuming that it is not too close to the end)
        # Choose a random index from start time to end_time - time_request
        min_time_select = Time(max(mjd_times.isot), format='isot')
        min_time_select -= time_separation # subtract the time duration

        min_time_select_mjd = Time(Time(min_time_select, format='isot').mjd, format='mjd')

        mjd_time_mod = mjd_times[mjd_times<=min_time_select_mjd] # times that are less than that time stamp

        # now select a random index within that time_stamp
        rand_time = np.random.randint(0, len(mjd_time_mod)) # select a random index

        start_time = mjd_time_mod[rand_time]

        # Now from this random time, we need to add our time separation
        time_end_isot = Time(Time(start_time, format='mjd').isot, format='isot') + time_separation

        end_time = Time(time_end_isot, format='isot').mjd

        return table[np.where((table['mjd']>=start_time.mjd) & (table['mjd']<=end_time))]
