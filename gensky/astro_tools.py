""" A series of astronomical tools """

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import numpy as np
from astropy import constants as const

def flux_2_mag(flux, zp):
    """Convert flux into magnitude.
    Input
    -----
    flux (float or numpy.ndarray): array of your flux estimates in mJy
    zp (float): zero-point of your optical system
    Output
    ------
    magnitude (float or numpy.ndarray): array of the AB magnitude
     """
    return (-2.5*np.log10(flux) + zp)


def mag_2_flux(mag, zp):
    """Convert AB magnitudes into flux (mJy) assuming standard relationship.
    Input
    -----
    mag (float or numpy.ndarray): array of AB magnitudes of your source
    zp (float): zero-point of your optical system
    Output
    -----
    flux (float or numpy.ndarray): array flux measurements in mJy
    """
    return (10**((zp-mag)/2.5))

def mag_err_2_flux_err(flux, mag_err):
    """
    Convert AB magnitude uncertainity to flux uncertainty (mJy).
    Input
    -----
    flux (float or numpy.ndarray): array of flux (mJy) estimates
    mag_err (float or numpy.ndarray): array of your magnitude uncertainty
    Output
    -----
    flux_err (float or numpy.ndarray): array of flux uncertainty
    """
    return ((flux*mag_err)/1.0857)

def flux_err_2_mag_err(flux, f_err):
    """
    Convert flux uncertainity to AB magnitude uncertainity.
    Input
    -----
    flux (float or numpy.ndarray): array of flux (mJy) estimates
    f_err (float or numpy.ndarray): array of your flux uncertainty
    Output
    -----
    mag_err (float or numpy.ndarray): array of AB magnitude uncertainty
    """
    return (1.0857 * f_err)/flux

def dist_mod_mag(app_mag, distance):
    """Calculate the absolute magnitude of a cosmological source using the distance modulus.
        Input
        -----
        app_mag (float or numpy.ndarray): array of apparent AB magnitudes observed
        distance (float or numpy.ndarray): distance of source in Mpc
        Output
        ------
        abs_mag (float or numpy.ndarray): array of absolute magnitude of source
    """
    return (app_mag - 5*np.log10(distance)-25)

def app_mag(abs_mag, distance):
    """
    Calculate the apparent magnitude (AB) of a source assuming a luminosity and distance using the distance modulus.
    Input
    -----
    abs_mag (float or numpy.ndarray): absolute magnitude of source
    distance (float or numpy.ndarray): distance of source in Mpc
    Output
    ------
    app_mag (float or numpy.ndarray): apparent AB magnitude
    """
    return abs_mag + 5*np.log10(distance) + 25

def redshift_to_distance(z):
    """
    Convert the observed redshift of a transient into luminosity distance in Mpc.
    Input
    -----
    z (float or numpy.ndarray): redshift of transient/source
    Output
    ------
    distance (float): distance of transient in Mpc
    """
    cosmo = FlatLambdaCDM(H0=67.7 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.30966, Neff=3.05, m_nu=[0., 0., 0.06]*u.eV, Ob0=0.049)
    return (cosmo.luminosity_distance(z).value)

def m52snr(m, m5):
    """
    Calculate the SNR for a star of magnitude m in an
    observation with 5-sigma limiting magnitude depth m5.
    Assumes gaussian distribution of photons and might not be
    strictly due in bluer filters. See table 2 and equation 5
    in astroph/0805.2366. Function taken from: https://sims-maf.lsst.io/modules.html
    Parameters
    ----------
    m : float or numpy.ndarray
        The magnitude of the star
    m5 : float or numpy.ndarray
        The m5 limiting magnitude of the observation
    Returns
    -------
    float or numpy.ndarray
        The SNR
    """
    snr = 5.*10.**(-0.4*(m-m5))
    return snr
