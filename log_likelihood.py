## Imports necessary libraries 
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck18
from astropy import constants as const
import astropy.cosmology.units as cu
from dark_emulator import model_hod
from numpy.linalg import inv
from scipy.stats import norm
import astropy.units as u
from glob import glob 
import numpy as np

# Instantiates the base class that computes lensing + clustering signals 
hod = model_hod.darkemu_x_hod({"fft_num":8})
