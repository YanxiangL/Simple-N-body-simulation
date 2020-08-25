# Simple-N-body-simulation
This is a simple N-body simulation which uses Newton's law of graivty to calculate the forces between particles in the simulation. In addition, this code also allows users to use a simplified version of the tree code to approximate the motion of particles in the simulation with less time. This code will also calculate the two-point correlation function using the Peebles-Hauser estimator. The matter power spectrum of the simulation is obtained by Fourier transforming the two-point correlation function. 
## Scale of the simulation
The default scale in the simulation is:
One unit length in the simulation corresponds to 1 kpc (proper distance) of physical scale. 
One unit mass in the simulation corresponds to 10^5 solar masses (about the same mass of a star cluster) in reality.
One unit time in the simulation corresponds to 1 Myr (million year) in reality. 

## System requirement
This is a python code. The modules require to run this code are numpy (https://numpy.org/), pylab (https://scipy.github.io/old-wiki/pages/PyLab), matplotlib (https://matplotlib.org/) and scipy (https://www.scipy.org/). I haven't check the dependency on these modules yet. If the code doesn't run properly, try updating the module. To run the simulation, you also have to download the ffmpeg (https://ffmpeg.org/) to the same diretatory as the code. 

## Input parameter
1. mass_min: This is the minimum mass resolution in the simulation. It should be a nonzero float in the unit of 10^5 solar masses. 
2. mass_max: This is the maximum mass of the particle in the simulation, it should be bigger than mass_min. 
3. bound_xyz: This is half length of the simulation box in each dimension, the possible positions are between -bound_xyz and bound_xyz
4. method: This determines the method used to calculate the simulation. Enters 0 for grid approximation method (highly recommand if the number of particles in the simulation is greater than 500) and 1 for exact method from Newton's law of graivty. 
5. grid_xyz: The number of grid cells in each direction. To achieve maximum efficiency, grid_xyz = 4-6 is enough when number of particles is between 100 and a few thousand. 
6. Np: This is the total number of particles in the simulation. It must be a nonzero integer.
7. Nt: Total number of time steps. It must be a nonzero integer. 
8. dt: The time step of the simulation in the unit of Myr. A general rule of thumb is dt should be on the same order of the ratio between the soften length and typical velocity of particles in the simulation. 
9. v_max: The maximum drift velocity in the unit of kpc/Myr. This parameter is only important when you initialise a random velocity distribution. 
10. H_0: Hubble constant in km/s/Mpc. For static universe, enter 0.0. In this simulation, the Hubble constant will not change with respect to time. Therefore, the simulation is most appropriate at small scale where recession velocity is neglectable or in the dark energy dominated era of the universe. Otherwise, you can set the Hubble parameter to zero and interpret the distance in the simulation as comoving distance. 
11. smooth: This is the soften length of the simulation. The soften length is the minimum distance between two particles in the simulation when the Newtonian force is calculated. This prevents the force to blow up to infinity when two particles are two close together. The main reason to use the soften length is the particle in the simulation is treated as point particles, while all physical objects will have finite size. The soften length should help the particles behave like they have a finite size. The soften length is also the minimum length resolution of the simulation. Any results obtain from the scale below the soften length is no longer accurate. The power spectrum will be also ploted up to the scale of the soften length in Fourier space. The soften length should be chosen as the size of the most massive object in the simulation. For example, if the maximum mass is set to 10^5, this corresponds to 10^10 solar masses which is typically the mass of the dwarf galaxy. The typical length of the dwarf galaxy is on the order of kpc. Therefore, the soften length in this simulation should be set to 1.0 since any particle interaction within 1 kpc is not able to be resolved in the simulation. 

## Default setting
The initial condition for the position is set to random distribution of particles within a 10 kpc cube. The velocity of particles follow a 3 dimentional Gaussian distribution with the the mean velocity set to 0 and standard deviation set to 0.1 kpc/Myr in each dimension. The mass of the particle is assigned randomly between the mass_min and mass_max provides by the user. 

## Advanced settings.
Initial conditions:
The initial condition can be modified inside the code. To modify the initial position, change the variable called "position" to the desired initial condition. To modify the inital velocity, change the variable called "velocity" to the desired initial condition. To change the mass distribution of the particle, find a variable called "mass" and change it to the desired distribution.

For legends, axis label, title and size of the plot. You can also change these in the code directly. 

## Exact method
The exact method will use Newton's law of graivty to calculate the force between particles. When we calculate the force of a particle, the code will loop over all other particles to calculate the total force of all other particles acting on the desired particle. This process will repeat for all particles until the force of all particles are found. The position and velcity of the particle is advanced in time by using the Leapfrog method. The leapfrog method is similar to the Euler's method except we advanced the position half a time step in front of the velocity. 

## Grid approximation
The grid approximation is a simplified version of tree code, it takes less time to complete the simulation by sacrificing the accuracy of the simulation. The user will determine the number of grid cells in each direction the code will automatically find the index of particles inside each grid and the center of mass of all particles in each grid. If the distance between the particle and the center of mass of a grid cell is more than the length of the grid. The length of the grid is the main diagonal of the cube. The acceleration of the particle is then calculated between the particle and the center of mass of the grid. The soften length is the length of the grid since this calculation won't resolve the below the length of the grid. 
