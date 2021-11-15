"""
- Generate light curves for a specified model
- Using the sampling algorithm at hand
- Specify the epochs range
- Specify the number of light curves per pointing
"""


# LSST constants from https://smtn-002.lsst.io/
lsst_bands = ['u', 'g', 'r', 'i', 'z', 'y']
m5_limiting = [23.87, 24.82, 24.36, 23.93, 23.36, 22.47]
filter_zeropoint = [27.03, 28.38, 28.16, 27.85, 27.46, 26.68]

# Mock initialize:

# Support the following modes:
# deep drill light curve (i.e a fixed ra,dec)
# sky scan drill light curve (i.e scan with ra and dec)


# Goal 1: Using a single ra, dec do light curve deep drilling:

# >> ./generator.py drilling_1.json
# >> output all light curves in  a seperate file

# What will the .json file contain?
# {coordinates: {ra, dec}, photometric_model: {'v19'}, bands:{all}, survey_dur:{2 yrs}}
#
#
#
