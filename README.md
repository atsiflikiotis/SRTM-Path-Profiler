# SRTM-Path-Profiler
A Python program that reads SRTM1 (1-arc/30m accuracy) data and returns the path profile between 2 specified points.

# Usage
1. Uncompress all needed .hgt files to /srtmdata folder. You can download SRTM1 data files from here: 
https://search.earthdata.nasa.gov/search
(an earthdata account required: https://urs.earthdata.nasa.gov//users/new)

2. Import srtm_path_profiler and run 'get_path_profile' function:
```python
import srtm_path_profiler as srtm
dist, atrdeg, di, hi = srtm.get_path_profile(phi_t, psi_t, phi_r, psi_r, coord_samples=700, fill_missing=False)
```
Input:
phi_t, psi_t, phi_r, psi_r are latitude and longitude coordinates of two points (for example transmitter and receiver), 
coord_samples is the number of samples to construct the path, 
fill_missing (optional, default: False): set to True only if you want to fill with 0's all the missing data

Returns:
dist: Distance in km
atrdeg: Bearing angle from point 1 to point 2 (transmitter to receiver)
di: numpy array with distance of the i-th terrain point
hi: Elevation of the i-th terrain point


3. If you want to plot the path after :
```python
srtm.plotpathprofile(di, hi, earthcurvature=True)
di, hi: see above
earthcurvature (optional, default: True): set to False to pot path profile on flat Earth
```

Packages that you'll need:
numpy
sxipy
matplotlib (for plot only)

For more info on .hgt SRTM data:
https://librenepal.com/article/reading-srtm-data-with-python/

An alternative tool for downloading SRTM1 data easily (still requires nasa earthdata account):
https://dwtkns.com/srtm30m/

