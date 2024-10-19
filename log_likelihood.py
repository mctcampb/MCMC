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

# Instantiates the base class that computes g-g lensing + clustering signals 
hod = model_hod.darkemu_x_hod({"fft_num":8})

# Sets the Hubble constant to the Planck 2018 value 
H0 = Planck18.H0.value 
# Defines little h 
h = H0 / 100
# Sets the fraction of baryons to the Planck 2018 value 
Omega_b = Planck18.Ob0
# Defines little "Omega_b" 
omega_b = Omega_b * (h ** 2)
# Sets little "Omega_nu" to Dark Emulator's default value 
omega_nu = 6.4e-4
# Defines the fraction of neutrinos 
Omega_nu = omega_nu / (h ** 2)
# Defines the Boltzmann constant in eV/K 
kB_evK = const.k_B.to(u.eV / u.K).value
# Computes the mass (in eV) of a neutrino species using Astropy's formulae. Note: we have 3 neutrino species; two are massless. This ensures consistency between 
# Dark Emulator and Astropy
m_nu_massive = kB_evK * Planck18.Tnu0.value * (((((3 * Omega_nu / (Planck18.Ogamma0 * Planck18.Neff * 0.22710731766)) - 2) ** 1.83) - 1) ** 0.54644808743) / 0.3173
# Sets the wCDM cosmology to a flat lambda-CDM cosmology 
omega = -1.0
# Sets the primordial curvature power spectrum's tilt to the Planck 2018 value 
ns =  0.9665
# Sets the upper bound of the integral for wp to 100 Mpc/h 
pimax = 100.0 
# Defines the number of jackknife patches used to estimate measurement covariance matrices 
njk = 150 
# Defines the mean galaxy abundance in (Mpc/h)^-3 
mu = 3.85e-4
# Defines the spread in galaxy abundances   
sigma = 0.85e-4 
# Defines a normalization factor for the prior on galaxy abundance  
norm_fact = norm.pdf(mu, mu, sigma)
# Defines the median redshift of lenses 
zl = 0.54 * cu.redshift 
# Arbitrarily sets the median redshift of sources to 0.7 
zs = 0.7 * cu.redshift 

## Defines flat prior lower and upper bounds for the parameters that are allowed to vary  
lnAs_range = [2.4752, 3.7128]
Omega_m_range = [0.286, 0.338]
logMmin_range = [11.5, 14.5]
sigma_sq_range = [0.001, 1.0]
logM1_range = [12.0, 16.5]
alpha_range = [0.01, 3.0]
kappa_range = [0.001, 3.0]

# Defines the path to the folder that contains measurement data vectors + covariance matrices  
folder_path = "/scratch/campb951/mcmc/data_files/"  
# Defines the file extension
file_pattern = "*.npy"
# Creates a list of measurement files 
file_list = glob(folder_path + file_pattern)
# Initializes a dictionary to organize measurement files 
ds_wp_vects = {}

