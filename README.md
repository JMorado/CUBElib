# CUBElib

### Introduction
Python library to read, manipulate and visualize files in the CUBE format.



```python
import sys
sys.path.append('../src/')
from CUBElib import *

new_cube_plus = Cube("o4_ngwf01_hplus_600")
new_cube_minus = Cube("o4_ngwf01_hminus_600")

# Plot Isosurface
new_cube_plus.plot_isosurface(1e-4, alpha = 0.5)
new_cube_minus.plot_isosurface(1e-4, alpha = 0.5)

# Create new cube file for the derivate and compute it
# Compute the symmetric derivative
sym_deriv = new_cube_plus.symmetric_derivative(new_cube_minus,0.0001)
deriv2 = Cube("o4_ngwf01_hminus_600")
deriv2._cube_grid = sym_deriv
deriv2.plot_density_map(normed=False)
```

![GitHub Logo](/figures/orbital_isosurface_example.png)
Format: ![Alt Text](url)

### Symmetric derivatives

Let us consider the electronic density of given orbital  <img src="https://latex.codecogs.com/gif.latex? \phi_\alpha " />  where <img src="https://latex.codecogs.com/gif.latex? \alpha " /> is the atom to which this orbital belongs. 

