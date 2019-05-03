# TODO: 2D density map (see https://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html)
# TODO: make derivative hypothesis more intuitive (e.g., derivative(self, otherCube, h=0.0001))
import numpy as np

class Cube:
    def __init__(self, cube_file):
        """
        :param cube_file (str): name of the cube file to be read
        """

        self._origin = [None, None, None]
        self._n = [None, None, None]
        self._dx = [None, None, None]
        self._natoms = None
        self._cube_grid = None
        self._coordinates = None
        self._atomic_numbers = None
        self._nuclear_charge = None

        # Variables that control the units
        self._bohr_to_angstrom = 0.529177249
        self._units = None

        if self._cube_file_exists(cube_file):
            self._cube_file = cube_file
            self._read_cube_file()
        else:
            raise Exception("File {} wat not found.".format(str(cube_file)))

    def __str__(self):
        return "This object contains the CUBE file with name: " + str(self._file_name) + "."

    def __repr__(self):
        return "This object contains the CUBE file with name: " + str(self._file_name) + "."

    @classmethod
    def _cube_file_exists(self, cube_file):
        """
        :param cube_file: Name of cube file.
        :return: Boolean that indicates whether cube_file exists (True) or not (False).
        """
        import os.path

        return os.path.isfile(cube_file)

    def _read_cube_file(self):
        """
        TODO: write description
        """

        f = open(self._cube_file, 'r')
        file_data = f.readlines()
        self._natoms = int(file_data[2].split()[0])

        for i in range(3):
            self._origin[i] = float(file_data[2].split()[i+1])

        # Line 3, 4 and 5
        for i in range(3):
            self._n[i] = int(file_data[i+3].split()[0])
            self._dx[i] = float(file_data[i+3].split()[i+1])

        # If the sign of the number of voxels in a dimension is positive then the units are Bohr,
        # If negative then Angstroms.
        if np.sign(self._n[0]) > 0:
            self._units = "BOHR"
        else:
            self._units = "ANGSTROM"

        self._cube_grid = np.ones((self._n[0],self._n[1],self._n[2]))

        # Read the coordinates, atomic numbers and nuclear charge
        self._coordinates = [[None, None, None] for x in range(self._natoms)]
        self._atomic_numbers = [None for x in range(self._natoms)]
        self._nuclear_charge = [None for x in range(self._natoms)]

        for i in range(self._natoms):
            self._atomic_numbers[i] = int(file_data[i+6].split()[0])
            self._nuclear_charge[i] = float(file_data[i+6].split()[1])
            for j in range(3):
                self._coordinates[i][j] = float(file_data[i+6].split()[2+j])

        if self._n[2] % 6 == 0:
            nlines = self._n[2] // 6
        else:
            nlines = self._n[2] // 6 + 1

        # Read the CUBE file grid data
        line = self._natoms + 6
        for x in range(self._n[0]):
            for y in range(self._n[1]):
                z = 0
                for j in range(nlines):
                    line_data = file_data[line].split()
                    for i in range(len(line_data)):
                        self._cube_grid[x,y,z] = float(line_data[i])
                        z += 1

                    line += 1

        f.close()
        return self._cube_grid

    def _write_cube_file(self, file_name = "output.cube"):
        """
        """

        f = open(file_name, 'w')

        f.write(" CUBElib generated cube file \n")
        f.write(" by Joao Morado \n")
        f.write("{:4d} {:12.5f} {:12.5f} {:12.5f} \n".format(self._natoms, *self._origin))

        for i in range(3):
            dxTmp = [0.0, 0.0, 0.0]
            dxTmp[i] = self._dx[i]
            f.write("{:4d} {:12.5f} {:12.5f} {:12.5f} \n".format(self._n[i], *dxTmp))

        for i in range(self._natoms):
            f.write("{:4d} {:12.5f} {:12.5f} {:12.5f} {:12.5f} \n".format(self._atomic_numbers[i], self._nuclear_charge[i], *self._coordinates[i]))  # CHANGE TO self._dx

        # Write the CUBE file grid data
        if self._n[2] % 6 == 0:
            nlines = self._n[2] // 6
        else:
            nlines = self._n[2] // 6 + 1

        line = self._natoms + 6
        for x in range(self._n[0]):
            for y in range(self._n[1]):
                for j in range(nlines):
                    count = 0
                    line_data = self._cube_grid[x,y,j*6:(j+1)*6]
                    string_to_write = "{:13.5E}" * len(line_data)
                    string_to_write = string_to_write.format(*line_data)
                    string_to_write = string_to_write + "\n"
                    f.write(string_to_write)
        f.close()

        return self._cube_grid

    def dot_product(self, otherCube):
        """
        This method performs the dot product between self._cube_grid and otherCube._cube_grid
        :param otherCube:
        :return:
        """
        # Assert that dimensions coincide
        for i in range(3):
            assert self._n[i] == otherCube._n[i], "Size of dimension {} does not match.".format(i)

        product = self._cube_grid * otherCube._cube_grid
        dot = 0
        for x in range(self._n[0]):
            for y in range(self._n[1]):
                for z in range(self._n[2]):
                    dot += product[x,y,z] # self._cube_grid[x,y,z] * otherCube._cube_grid[x,y,z]

        if self._units == "BOHR":
            return dot * self._dx[0] * self._dx[1] * self._dx[2] #* self._bohr_to_angstrom ** 3
        elif self._units == "ANGSTROM":
            return dot * self._dx[0] * self._dx[1] * self._dx[2]
        else:
            return "Something is wrong... unit {} is not implemented.".format(self._units)

    def symmetric_derivative(self, otherCube, h=0.0001):
        """
        Computes the numerical derivative of the quantity represented in the CUBE files.
        :param otherCube:
        :param h:
        :return: 3D grid containing the derivatives
        """
        import copy

        derivative = (self._cube_grid - otherCube._cube_grid) / (2.0 * h / self._bohr_to_angstrom)

        newCube = copy.deepcopy(self)
        newCube._coordinates = (np.asarray(self._coordinates) + np.asarray(otherCube._coordinates)) / 2.0
        newCube._cube_grid = derivative

        return newCube

    def grid_difference(self, otherCube):
        """
        Computes the difference between
        :param otherCube
        :return: 3D grid containing the difference between self and otherCube grids.
        """

        return self._cube_grid - otherCube._cube_grid

    def plot_isosurface(self, iso_value, alpha = 0.5):
        """

        :return:
        """
        def isosurface(grid,iso_value,dx):
            """
            Returns the isosurface of value iso_value of grid.

            """
            from skimage import measure

            verts, faces, normals, values = measure.marching_cubes_lewiner(grid, iso_value, spacing=dx)

            return verts, faces

        from matplotlib import pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        # Compute the isosurface contour
        verts, faces = isosurface(self._cube_grid,iso_value,self._dx)

        # Create the plot instance and plot the data
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Fancy indexing: 'verts[faces]' to generate a collection of triangles
        mesh = Poly3DCollection(verts[faces]+self._origin,alpha=alpha)

        # Customize the plot
        face_color = [0.5, 0.5, 1]
        mesh.set_facecolor(face_color)
        mesh.set_edgecolor('k')
        ax.add_collection3d(mesh)
        ax.set_title("Orbital isosurface (iso=" + str(iso_value) + ")")
        ax.set_xlim(40, 60)
        ax.set_ylim(40, 60)
        ax.set_zlim(40, 60)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.tight_layout()

        return plt.show()

    def plot_density_map(self, normalized = True, alpha=0.3, threshold=0.05, title = "Orbital Electronic Density"):
        """
        Plots the electronic density map of the orbital.
        """
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import pyplot as plt
        import copy

        X, Y, Z = np.mgrid[0:self._n[0], 0:self._n[1], 0:self._n[2]] * self._dx[0]
        X = X + self._origin[0]
        Y = Y + self._origin[1]
        Z = Z + self._origin[2]

        # Deepcopy in order to avoid modifying the original grid
        if normalized:
            cube_grid = copy.deepcopy(self._cube_grid) / self._cube_grid.max()
        else:
            cube_grid = copy.deepcopy(self._cube_grid)

        cube_grid[abs(cube_grid) <= threshold] = np.NaN

        # Create the plot instance and plot the data
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        density_scatter = ax.scatter(X, Y, Z, c=cube_grid.flatten(), alpha=alpha, cmap="seismic")

        # Customize the plot
        ax.set_title(title)
        ax.set_xlim(X.min()-5,X.max()+5)
        ax.set_ylim(Y.min()-5,Y.max()+5)
        ax.set_zlim(Z.min()-5,Z.max()+5)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        cbar = fig.colorbar(density_scatter,fraction=0.046, pad=0.04)
        cbar.set_label(r'$\AA^{-1}$', rotation=270)

        return plt.show()
