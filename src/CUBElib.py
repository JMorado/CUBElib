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

        if self._cube_file_exists(cube_file):
            self._cube_file = cube_file
            self._read_cube_file()
        else:
            raise Exception("File {} wat not found.".format(str(cube_file)))


    def __str__(self):
        return "This object contains the CUBE file with name: " + str(self._file_name) + "."

    def __repr__(self):
        return "This object contains the CUBE file with name: " + str(self._file_name) + "."

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
        import numpy as np

        f = open(self._cube_file, 'r')
        file_data = f.readlines()

        self._natoms = int(file_data[2].split()[0])

        for i in range(3):
            self._origin[i] = float(file_data[2].split()[i+1])

        # Line 3, 4 and 5
        for i in range(3):
            self._n[i] = int(file_data[i+3].split()[0])
            self._dx[i] = float(file_data[i+3].split()[i+1])


        self._cube_grid = np.ones((self._n[0],self._n[1],self._n[2]))

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

        return self._cube_grid



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


    def plot_density_map(self, alpha=0.3, threshold=0.0):
        """
        Plots the electronic density map of the orbital.
        """
        from matplotlib import pyplot as plt
        import numpy as np

        X, Y, Z = np.mgrid[0:self._n[0], 0:self._n[1], 0:self._n[2]] * self._dx[0]
        X = X + self._origin[0]
        Y = Y + self._origin[1]
        Z = Z + self._origin[2]

        cube_grid = self._cube_grid / self._cube_grid.max()
        cube_grid[cube_grid <= threshold] = np.NaN

        # Create the plot instance and plot the data
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        density_scatter = ax.scatter(X, Y, Z, c=cube_grid.flatten(), alpha=alpha, cmap="hot")

        # Customize the plot
        ax.set_title("Orbital Electronic Density")
        ax.set_xlim(X.min()-5,X.max()+5)
        ax.set_ylim(Y.min()-5,Y.max()+5)
        ax.set_zlim(Z.min()-5,Z.max()+5)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        fig.colorbar(density_scatter)
        plt.tight_layout()

        return plt.show()
