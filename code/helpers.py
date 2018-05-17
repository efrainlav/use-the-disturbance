"""
Helper functions for loading and plotting data.
"""
import tarfile

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def minmax(data, fields):
    """
    Get the minimum and maximum data values out of all fields given.
    Returns them in a dictionary with the 'vmin' and 'vmax' keys.
    """
    vmin = min(data[field].values.min() for field in fields)
    vmax = max(data[field].values.max() for field in fields)
    return dict(vmin=vmin, vmax=vmax)


def plot_field(ax, data, field, cmap=None, gridline_spacing=30, cb_pad=0.03,
               cb_aspect=50, cb_shrink=0.7, title=True, **kwargs):
    """
    Make a pcolormesh plot of the given data field.
    Set's the plot extent and includes ticks in longitude and latitude.
    """
    if title:
        ax.set_title(field)
    if 'add_colorbar' not in kwargs:
        kwargs['cbar_kwargs'] = dict(orientation='horizontal',
                                     aspect=cb_aspect, pad=cb_pad,
                                     shrink=cb_shrink)
    data[field].plot.pcolormesh(ax=ax, add_labels=False, cmap=cmap,
                                transform=ccrs.PlateCarree(), **kwargs)
    ax.coastlines()
    w, e, s, n = [data.longitude.min(), data.longitude.max(),
                  data.latitude.min(), data.latitude.max()]
    xlocs = np.arange(w, e + 0.01, gridline_spacing)
    ylocs = np.arange(s, n + 0.01, gridline_spacing)
    ax.gridlines(color="#cccccc55", xlocs=xlocs, ylocs=ylocs)


def load_icgem_gdf(fname, dtype='float64'):
    """
    Load data from an ICGEM .gdf file into an xarray.Dataset.

    Reads metdata from the header, like the grid area, number of points,
    height over the ellipsoid (if it's constant), etc.

    Stores the file header in the ``attrs`` attribute of the Dataset.
    
    Assumes that files are in tar files with gz compression.

    Parameters
    ----------
    fname : str
        The name of the .gdf.tar.gz file.
    dtype : str or numpy dtype object
        The data type used when loading the data from the file.

    Returns
    -------
    data : xarray.Dataset

    """
    with tarfile.open(fname, 'r:gz') as archive:
        gdf_file = archive.extractfile(archive.getmembers()[0])
        # Read the header and extract metadata
        header = []
        shape = [None, None]
        size = None
        height_over_ell = None
        fields = None
        is_field_names = False
        west, east, south, north = [None]*4
        for line in gdf_file:
            if line.strip()[:11] == 'end_of_head':
                break
            if not line.strip():
                # The field names will come after a blank line
                is_field_names = True
                continue
            header.append(line)
            parts = line.strip().split()
            if parts[0] == 'height_over_ell':
                height_over_ell = float(parts[1])
            elif parts[0] == 'latitude_parallels':
                nlat = shape[0] = int(parts[1])
            elif parts[0] == 'longitude_parallels':
                nlon = shape[1] = int(parts[1])
            elif parts[0] == 'number_of_gridpoints':
                size = int(parts[1])
            elif parts[0] == 'latlimit_south':
                south = float(parts[1])
            elif parts[0] == 'latlimit_north':
                north = float(parts[1])
            elif parts[0] == 'longlimit_west':
                west = float(parts[1])
            elif parts[0] == 'longlimit_east':
                east = float(parts[1])
            if is_field_names:
                # Skip the first two because they are the coordinate
                # names.
                fields = line.strip().split()[2:]
                is_field_names = False
        # Read the numerical values
        rawdata = np.loadtxt(gdf_file, ndmin=2, unpack=True, dtype=dtype)

    # Sanity checks
    assert all(n is not None for n in shape), "Couldn't read shape of grid."
    assert size is not None, "Couldn't read size of grid."
    assert shape[0]*shape[1] == size, \
        "Grid shape '{}' and size '{}' mismatch.".format(shape, size)
    assert fields is not None, "Couldn't read column names."
    assert len(fields) == rawdata.shape[0] - 2, \
        "Number of attributes ({}) and data columns ({}) mismatch".format(
            len(fields), rawdata.shape[0] - 2)
    assert all(i is not None for i in [west, east, south, north]), \
        "Couldn't read the grid area."

    if height_over_ell is not None:
        fields.append('height_over_ell')
        rawdata.append(height_over_ell*np.ones(size, dtype=dtype))

    # Build the xarray container
    dims = ['latitude', 'longitude']
    latitude = np.linspace(south, north, nlat, dtype=dtype)
    longitude = np.linspace(west, east, nlon, dtype=dtype)
    # Cartopy doesn't like 0-360 longitude
    longitude[longitude > 180] -= 360
    coords = {'latitude': latitude, 'longitude': longitude}
    attrs = {'file_header': ''.join(header)}
    data_vars = {}
    for name, value in zip(fields, rawdata[2:]):
        # Need to invert the data matrices in latitude "[::-1]"
        # because the ICGEM grid varies latitude from N to S
        # instead of S to N.
        data_vars[name] = xr.DataArray(
            value.reshape(shape)[::-1], coords=coords, dims=dims, name=name,
            attrs=attrs)

    return xr.Dataset(data_vars, attrs=attrs)
