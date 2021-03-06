# Data sets used in the paper

This folder contains the raw and processed data files used in the paper.
All gravity and topography were generated from spherical harmonic models using
the [ICGEM Calculation Service](http://icgem.gfz-potsdam.de).

* `eigen-6c4-geoid.tar.gz`: Global grid of geoidal heights (in meters) generated
  from the [EIGEN-6c4](https://doi.org/10.5880/icgem.2015.1) model based on the
  WGS84 reference ellipsoid.
* `eigen-6c4-gravity.tar.gz`: Global grid of the magnitude of the gravity vector
  (in mGal) generated from the [EIGEN-6c4](https://doi.org/10.5880/icgem.2015.1)
  model. Values are calculated on the Earth's surface and observation heights
  given are referenced to the geoid.
* `etopo1.tar.gz`: Global grid of topography and bathymetry (in meters) referenced
  to the geoid and generated from a spherical harmonic model of
  [ETOPO1](https://doi.org/10.7289/V5C8276M).
* `disturbance-anomaly-difference.nc`: Global grid in a netCDF file the
  absolute difference of the gravity disturbance and the free-air anomaly.
  Generated by [`code/difference.ipynb`](code/difference.ipynb).
