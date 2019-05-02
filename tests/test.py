import sys
sys.path.append('../src/')
from CUBElib import *


# Numerical derivative
new_cube_plus = Cube("o4_ngwf01_hplus_600")
new_cube_minus = Cube("o4_ngwf01_hminus_600")
sym_deriv = new_cube_plus.symmetric_derivative(new_cube_minus,0.0001)
#new_cube_minus._write_cube_file()
exit()
# Plot Isosurfac
#new_cube_plus.plot_isosurface(1e-4, alpha = 0.5)
#new_cube_minus.plot_isosurface(1e-4, alpha = 0.5)

# Analytical derivative
analytical_deriv = Cube("ngwf_derivative_atom00004_ngwf01.cube")
analytical_deriv.plot_density_map(normalized=False, title=r"$d\phi/dR_y$ (O4, NGWF 0001) h = 0.0001")


# Create new cube file for the difference between derivatives and compute it
# Compute the symmetric derivatives
sym_deriv = new_cube_plus.symmetric_derivative(new_cube_minus,0.0001)
deriv2 = Cube("o4_ngwf01_hminus_600")
deriv2._cube_grid = sym_deriv - analytical_deriv._cube_grid
print(deriv2._cube_grid.shape)
deriv2.plot_density_map(normalized=False, title=r"$d\phi/dR_y$ (O4, NGWF 0001) h = 0.0001")
