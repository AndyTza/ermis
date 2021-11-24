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
import numpy.ma as mask
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import models
import generator
import dummy_func
import astro_tools
import json
from astro_tools import m52snr


# LSST constants from https://smtn-002.lsst.io/
lsst_bands = ['u', 'g', 'r', 'i', 'z', 'y']
m5_limiting = [23.87, 24.82, 24.36, 23.93, 23.36, 22.47]
filter_zeropoint = [27.03, 28.38, 28.16, 27.85, 27.46, 26.68]


""" Main compilying funcitons. Here we will control:

- Light curve & parameter space input (i.e distance, coordinates, luminosity)
- Light curve model type
- Light curve model parameters (per band)
- Fetch and query data
- Store/secure light curve data (as .json files)
- In the same .json as above, store all variable information
"""


# query = {"coordinates":[120, -12], "distance":0.055, "luminosity": -13, "enable_MW_ebv": True,
#        "model": "fourier_sine", "model_theta": all_models, #list of model parameters
#         "sample_bands": 'all', "survey_duration":365*3}


def calc_lc_band(table, query, band='r', type='mag'):
    """Will take the astropy table and return important properties about the light curve"""

    phot_bands = ['u', 'g', 'r', 'i', 'z', 'y']
    phot_index = np.arange(0, 6, step=1)

    phase = table['mjd']-table['mjd'][0] # phase from first detection

    model_mag = models.photometricmodel(phase[table['filter']==band]).fourier_cos(*query['model_theta'][phot_index[phot_bands==band]]) # calculate the flux at each epoch
    lim_sigma = table[table['filter']==band]['5sig'] #5-sigma depth
    snr = m5snr(model_flux, lim_sigma) # SNR at that band
    sigma_err = 2.5*np.log10(1+1./snr) # magnitude error

    if type=='mag':
        return phase[table['filter']==band], model_mag, sigma_err
    elif type=='flux':
        return phase[table['filter']==band], model_mag, astro_tool.mag_err_2_flux_err(model_mag, sigma_err) # convert magnitude error into flux error (assume that model_mag is now in flux)


def run_sampler(query, nc_lcs=1000, store_path='../data/lc_output/'):
    """ Read the query dictionary and begin simulating light curves.

        Input
        -----
        query (dict): dictionary containing information.

        Currently supports the following query criteria:
        coordinates [list]: ra, dec in degrees
        distance [float]: distance in Mpc
        luminosity [float]: assumed luminosity of the transient

    """

    # Unpack & load cadences
    sim_table = generator.simulation.load_table() #complete simulation table
    sim_table_pointing = generator.fetch_cadence_info(sim_table, query['coordinates'][0], query['coordinates'][1], filter_band=query['sample_bands'])

    for i in tqdm(range(nc_lcs)):
        sim_table_gen = generator.simple_mjd_sampler(sim_table_pointing, time_separation=query['survey_duration'], mode='random') # should have all the epochs of observation withiin that time frame
        #phase = sim_table_gen['mjd'] - sim_table_gen['mjd'][0] # starting from the first detection

    return calc_lc_band(sim_table_gen)
    #Select each filter?




        # TODO: Clean
        #model_M = np.zeros(shape=(5000,6))
        #lim_5s = np.zeros(shape=(5000,6)) # no more than >1000 pointings per field (set to 5k to be safe); mark as nan the empty ones
        #snr = np.zeros(shape=(5000,6))
        #sigma_band = np.zeros(shape=(5000,6))
        #for index, filter in tqdm(enumerate(lsst_bands)):
    #        model_M[0:len(phase[sim_table_gen['filter']==filter]), index] = models.photometricmodel(phase[sim_table_gen['filter']==filter]).fourier_cos(*query['model_theta'][index])
    #        model_M[:,index][model_M[:,index]==0] = np.inf # mask as infinity TODO -- numpy mask instead
#            lim_5s[0:len(phase[sim_table_gen['filter']==filter],index] = sim_table_gen[sim_table_gen['filter']==filter]['5sig']#
    #        lim_5s[:,index][lim_5s[:,index]==0] = np.in
    #        snr[0:len(phase[sim_table_gen['filter']==filter],index] = m5snr(model_M[:, index], lim_5s[:,index]) # should be in linear order of filters
    #        snr[:,index][snr[:,index]==0] = np.inf
#            sigma_band[0:len(phase[sim_table_gen['filter']==filter],index] = 2.5*np.log10(1+1./snr[:,index])


        #return model_M, sigma_band

        #model_u = calc_lc_band()


        #SNR = m52snr(model_u, lim_sigma_u)

        # Issue:
        # Each photometric detection per band is not the same.

        # Calculate the theoretical flux of each model; model has already been scaled by distance and luminosity
        #model_u = models.photometricmodel(phase[sim_table_gen['filter']=='u']).fourier_cos(*query['model_theta'][0])
        #model_g = models.photometricmodel(phase[sim_table_gen['filter']=='g']).fourier_cos(*query['model_theta'][1])
        #model_r = models.photometricmodel(phase[sim_table_gen['filter']=='r']).fourier_cos(*query['model_theta'][2])
        #model_i = models.photometricmodel(phase[sim_table_gen['filter']=='i']).fourier_cos(*query['model_theta'][3])
        #model_z = models.photometricmodel(phase[sim_table_gen['filter']=='z']).fourier_cos(*query['model_theta'][4])
        #model_y = models.photometricmodel(phase[sim_table_gen['filter']=='y']).fourier_cos(*query['model_theta'][5])

        # what about the errors?
        #im_sigma_u = sim_table_gen[sim_table_gen['filter']=='u']['5sig']
        #lim_sigma_g = sim_table_gen[sim_table_gen['filter']=='g']['5sig']
        #lim_sigma_r = sim_table_gen[sim_table_gen['filter']=='r']['5sig']
        #lim_sigma_i = sim_table_gen[sim_table_gen['filter']=='i']['5sig']
        #lim_sigma_z = sim_table_gen[sim_table_gen['filter']=='z']['5sig']
        #lim_sigma_y = sim_table_gen[sim_table_gen['filter']=='y']['5sig']

        #snr_u = m52snr(model_u, lim_sigma_u) # median limiting 5-sigma depth?
        #snr_g = m52snr(model_g, lim_sigma_g) # median limiting 5-sigma depth?
        #snr_r = m52snr(model_r, lim_sigma_r) # median limiting 5-sigma depth?
        #snr_i = m52snr(model_i, lim_sigma_i) # median limiting 5-sigma depth?
        #snr_z = m52snr(model_z, lim_sigma_z) # median limiting 5-sigma depth?
        #snr_y = m52snr(model_y, lim_sigma_y) # median limiting 5-sigma depth?

        # The magnitude uncertainties that go with amplified mags
        #sigma_u = 2.5*np.log10(1+1./snr_u) # error in magnitude
        #sigma_g = 2.5*np.log10(1+1./snr_g) # error in magnitude
        #sigma_r = 2.5*np.log10(1+1./snr_r) # error in magnitude
        ##sigma_i = 2.5*np.log10(1+1./snr_i) # error in magnitude
        #sigma_z = 2.5*np.log10(1+1./snr_z) # error in magnitude
        #sigma_y = 2.5*np.log10(1+1./snr_y) # error in magnitude


        #X.append(phase[sim_table_gen['filter']=='r'])
        #Y.append(model_r)
        #YERR.append(sigma_r)
    #return X, Y, YERR




    # query = {"coordinates":[120, -12], "distance":0.055, "luminosity": -13, "enable_MW_ebv": True,
    #        "model": "fourier_sine", "model_theta": all_models, #list of model parameters
    #         "sample_bands": 'all', "survey_duration":365*3}
