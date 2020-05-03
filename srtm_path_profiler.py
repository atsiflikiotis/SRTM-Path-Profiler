import numpy as np
from pathlib import Path
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt

# folder name where srtmdata files are located:
hgtfolder = 'srtmdata'
hgtpath = Path.joinpath(Path.cwd(), hgtfolder)
x = np.linspace(0, 1, 3601, dtype=np.dtype('float32'))


class MissingHgtError(Exception):
    pass


def plotpathprofile(di, hi, earthcurvature=True):
    ec = np.full_like(hi, 0)
    if earthcurvature:
        Ro = 3971   # earth radius (km)
        # define dh_i as the change (reduction) of the elevation of the i-th terrain point, due to earth curvature:
        # At distance di from the point of reference x0, the elevation angle to the i-th terrain point, theta_i,
        # is reduced (the more we walk away from x0, the more the longest objects seem lower) by a factor:
        # arctan(di/(2Ro)), and tan(theta_i) = di/2Ro = dh_i/di.
        # Setting the reference point x0 as the middle of the great-arc distance from tx to rx (x0 = dist/2):
        ec = -1000 * (di - di[-1] / 2) ** 2 / (2 * Ro)
        ec -= min(ec)   # set tx and rx to 0, and the middle of path is getting the maximum height.
    h = hi + ec
    plt.style.use('seaborn')
    plt.plot(di, h, color='grey')
    plt.fill_between(di, h, ec, color='darkgoldenrod', alpha=0.6)
    plt.fill_between(di, ec, 0, color='silver')
    plt.xlim(0, di[-1])
    plt.ylim(0, max(2.5 * max(h), 1.7 * h[-1]))
    plt.xlabel('Distance (km)')
    plt.ylabel('Elevation (m)')


def get_path_profile(phi_t, psi_t, phi_r, psi_r, coord_samples=700, fill_missing=False):
    [dist, atrdeg, di, phi_values, psi_values] = get_path_geometry(phi_t, psi_t, phi_r, psi_r, coord_samples)
    subpaths, ufnames, flagsmissing = get_path_sections(phi_values, psi_values, fill_missing)
    hi = get_elevations(subpaths, ufnames, flagsmissing)

    #     plt.plot(di, elevationscurved, di, yc)
    if not fill_missing:
        if flagsmissing.any():
            missing = np.array2string(ufnames[np.nonzero(flagsmissing)])
            print('Missing HGT files', f"Missing srtmdata files and filled with 0 values:\n{missing}")
    print(f'Distance: {dist:.5f}km\nBearing angle: {atrdeg:.1f}°\nMax elevation:{np.max(hi):.0f}m')
    return dist, atrdeg, di, hi


def get_elevations(subpaths, ufnames, flagsmissing):
    interpvalues = np.empty(0, dtype=int)
    for i in range(len(subpaths)):
        # open srtmdata filename according to subpath
        path = Path.joinpath(hgtpath, ufnames[i])
        if not flagsmissing[i]:
            with open(path, 'rb') as data:
                # data is stored in srtmdata files in big-endian format:
                # SRTM1 (1-arc sampling) is used, so each file contains 3601x3601 individual elevations
                # after reading all elevatin ins file, we reshape it into a rectangular numpy array
                elevationdata = np.fromfile(data, np.dtype('>i2'), 3601 ** 2).reshape(3601, 3601)

            # interpolate elevation values with scipy.interpolate.griddata
            f1 = RectBivariateSpline(x, x, elevationdata, kx=1, ky=1)
            latcorr = int(ufnames[i][1:3]) + 1
            loncorr = int(ufnames[i][4:7])
            interpvalues = np.hstack((interpvalues, f1.ev(latcorr - subpaths[i][0], subpaths[i][1] - loncorr)))
        else:
            interpvalues = np.hstack((interpvalues, np.zeros_like(subpaths[i][0])))

    return interpvalues


def get_path_sections(phi_values, psi_values, fill_missing):
    fnames = get_hgt_names(phi_values, psi_values)

    # check if all srtmdata files exist in srtmdata directory
    _, indices = np.unique(fnames, return_index=True)
    indices = np.sort(indices)
    ufnames = fnames[indices]

    flagsmissing = np.full_like(indices, 0)
    for i in range(len(indices)):
        file = fnames[indices[i]]
        path = Path.joinpath(hgtpath, file)
        flagsmissing[i] = not path.exists()
    if not fill_missing:
        if flagsmissing.any():
            missing = np.array2string(ufnames[np.nonzero(flagsmissing)])
            print('Missing HGT files', f"Missing srtmdata files from srtmdata directory: {missing}")
            raise MissingHgtError(f"Missing srtmdata files from srtmdata directory:\n{missing}")

    # create sub-paths for each difference srtmdata file:
    subpaths = [None] * len(indices)
    for i in range(len(indices) - 1):
        startidx = indices[i]
        endidx = indices[i+1]
        subpaths[i] = (phi_values[startidx:endidx], psi_values[startidx:endidx])
        # subpaths[i] = list(zip(phi_values[startidx:endidx], psi_values[startidx:endidx]))
    subpaths[-1] = (phi_values[indices[-1]:], psi_values[indices[-1]:])
    # subpaths[-1] = list(zip(phi_values[indices[-1]:], psi_values[indices[-1]:]))
    return subpaths, ufnames, flagsmissing


def get_hgt_names(phi_values, psi_values):
    n_s = np.where(phi_values >= 0, 'N', 'S')
    e_w = np.where(psi_values >= 0, 'E', 'W')
    lat = abs(phi_values).astype(int).astype(str)
    lat = np.char.zfill(lat, 2)
    lon = abs(psi_values).astype(int).astype(str)
    lon = np.char.zfill(lon, 3)
    fnames = np.char.add(np.char.add(n_s, lat), np.char.add(e_w, lon))
    fnames = np.char.add(fnames, '.hgt')
    return fnames


def get_path_geometry(phi_t, psi_t, phi_r, psi_r, coord_samples):
    phi_t = np.deg2rad(phi_t)
    psi_t = np.deg2rad(psi_t)
    phi_r = np.deg2rad(phi_r)
    psi_r = np.deg2rad(psi_r)

    # The angle subtended by the path at the centre of the Earth, δ, from the stations’ geographic
    # coordinates:
    delta = np.arccos(np.sin(phi_t) * np.sin(phi_r) +
                      np.cos(phi_t) * np.cos(phi_r) *
                      np.cos(psi_t - psi_r))

    # great circle distance:
    R = 6371
    dist = 6371 * delta  # km

    # bearing (azimuthal direction clockwise from true North) from from transmitter to receiver
    atr = np.arccos((np.sin(phi_r) - np.sin(phi_t) * np.cos(delta)) /
                    (np.sin(delta) * np.cos(phi_t)))

    if psi_t > psi_r:
        atr = 2 * np.pi - atr
    atrdeg = np.rad2deg(atr)

    di = np.linspace(0, dist, coord_samples)

    phi_values = np.arcsin(np.sin(phi_t) * np.cos(di / R) + np.cos(phi_t) * np.sin(di / R) * np.cos(atr))
    psi_values = psi_t + np.arctan2(np.sin(atr) * np.sin(di / R) * np.cos(phi_t),
                                    np.cos(di / R) - np.sin(phi_t) * np.sin(phi_values))
    psi_values = np.remainder(psi_values + 3 * np.pi, (2 * np.pi) - np.pi)

    phi_values = np.rad2deg(phi_values)
    psi_values = np.rad2deg(psi_values)

    return [dist, atrdeg, di, phi_values, psi_values]
