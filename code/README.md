# Source code for producing the results and figures

Python code (mostly) to generate the results for the paper.
Most of the code is in Jupyter notebooks (the `.ipynb` files). Some helper
functions are defined in `helpers.py` and imported by the notebooks.

* `difference.ipynb`: Loads the gravity data from ICGEM, calculates the
  free-air anomaly, the disturbance, and their difference. Generates the
  difference netCDF grid `disturbance-anomaly-difference.nc`.
* `plot-difference.ipynb`: Plots the difference grid using GMT.
