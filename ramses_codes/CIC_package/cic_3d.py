import numpy as np

def cloud_in_cell_3D(posx, posy, posz, masses, grid_size, box_size):
    """
    Optimized Cloud-In-Cell (CIC) mass assignment in 3D using Numba.
    
    Parameters:
    posx, posy, posz : ndarray
        (N,) arrays of particle positions in x, y, and z directions.
    masses : ndarray
        (N,) array of particle masses.
    grid_size : int
        Number of grid points along each axis.
    box_size : float
        Physical size of the cubic box.
    
    Returns:
    density_grid : ndarray
        (grid_size, grid_size, grid_size) array of mass density values.
    """
    dx = box_size / grid_size  # Grid cell size
    density_grid = np.zeros((grid_size, grid_size, grid_size), dtype=np.float64)

    num_particles = len(posx)

    # Loop over each particle in parallel
    for p in range(num_particles):
        # Compute the index of the lower-left grid point
        i = int(np.floor(posx[p] / dx))
        j = int(np.floor(posy[p] / dx))
        k = int(np.floor(posz[p] / dx))


        # Compute distances to the lower-left grid point
        dx_p = (posx[p] / dx) - i
        dy_p = (posy[p] / dx) - j
        dz_p = (posz[p] / dx) - k

        # Compute weights for the eight surrounding grid points
        weights = np.array([
            (1 - dx_p) * (1 - dy_p) * (1 - dz_p),
            dx_p * (1 - dy_p) * (1 - dz_p),
            (1 - dx_p) * dy_p * (1 - dz_p),
            dx_p * dy_p * (1 - dz_p),
            (1 - dx_p) * (1 - dy_p) * dz_p,
            dx_p * (1 - dy_p) * dz_p,
            (1 - dx_p) * dy_p * dz_p,
            dx_p * dy_p * dz_p
        ], dtype=np.float64)

        weights /= np.sum(weights)  # Normalize weights to sum exactly to 1
        
        # Get neighboring grid points (with periodic boundary handling)
        neighbors = np.array([
            (i % grid_size, j % grid_size, k % grid_size),
            ((i+1) % grid_size, j % grid_size, k % grid_size),
            (i % grid_size, (j+1) % grid_size, k % grid_size),
            ((i+1) % grid_size, (j+1) % grid_size, k % grid_size),
            (i % grid_size, j % grid_size, (k+1) % grid_size),
            ((i+1) % grid_size, j % grid_size, (k+1) % grid_size),
            (i % grid_size, (j+1) % grid_size, (k+1) % grid_size),
            ((i+1) % grid_size, (j+1) % grid_size, (k+1) % grid_size)
        ])

        # Distribute mass to the grid points
        for n in range(8):
            ii, jj, kk = neighbors[n]
            density_grid[ii, jj, kk] += masses[p] * weights[n]

    # Convert mass to density by dividing by the grid cell volume
    density_grid /= dx**3

    return density_grid
