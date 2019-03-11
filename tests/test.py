import sys
sys.path.append('../src/')
from CUBElib import *

new_cube = Cube("final_atom00001_ngwf02.cube")
new_cube.plot_isosurface(1e-3)
new_cube.plot_density_map()
