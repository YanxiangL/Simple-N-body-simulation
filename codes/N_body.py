# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 19:56:46 2020

@author: Yansiang
"""

import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import CubicSpline
from scipy import integrate
import time
from scipy.interpolate import InterpolatedUnivariateSpline

# This is a simple skeleton code just to give some ideas.
# It plots collisionles particles moving at random in a cubic box in the main panel
# and shows the distribution of their separations in one of two other panels.

def find_particle(position, x, y, z):
    """
    This function finds the index of particles in each grid cell. 
    
    Parameters
    ----------
    position : array_like 
    This is the position of all particles in the simulation.
    x : array_like
    The array contains the grid point in the x direction.
    y : array_like
    The array contains the grid point in the y direction.
    z : array_like
    The array contains the grid point in the z direction.

    Returns
    -------
    num : List
        Number of particles in each grid.
    particle_in_cell : List
        The index of particles in each grid cell.

    """
    particle_in_cell =[]
    num = []
    #Initialise the output list
    limit_x = len(x)-1
    limit_y = len(y)-1
    limit_z = len(z)-1
    #Number of computation requires for each dimension. 
    counter_x = 0
    counter_y = 0
    counter_z = 0
    #This counts the number of cells that have been computed in each direction. 
    for i in range(limit_x*limit_y*limit_z):
        position_x = position[0,:]
        position_y = position[1,:]
        position_z = position[2,:]
        #The poistion of the particle in each direction. 
        xx = np.where(position_x < x[counter_x+1], position_x, 0)
        yy = np.where(position_y < y[counter_y+1], position_y, 0)
        zz = np.where(position_z < z[counter_z+1], position_z, 0)
        #Find the particles in the position array which are to the right of the desired grid. For such particle 
        #replace the position with zero. 
        xxx = np.where(xx > x[counter_x], xx, 0)
        yyy = np.where(yy > y[counter_y], yy, 0)
        zzz = np.where(zz > z[counter_z], zz, 0)
        #Find the particles in the position array which are to the left of the desired grid. For such particle 
        #replace the position with zero.
        
        index_x = np.nonzero(xxx)[0]
        index_y = np.nonzero(yyy)[0]
        index_z = np.nonzero(zzz)[0]
        #Find the index of the particle which are nonzero. These particles are located in the desired grid. 
        #print(index_x, index_y, index_z)
            
        xy = np.intersect1d(index_x, index_y, assume_unique=True)
        xyz = np.intersect1d(xy, index_z, assume_unique=True)
        #The codes above finds the index of particle in the desired grid. The index in each array is unique. 
        if (len(xyz != 0)):
            num.append(len(xyz))
            particle_in_cell.append(xyz)
            #Append the particle in the grid and the number of particle in the grid if the number of particles in
            #the grid is nonzero. 
            
        counter_x += 1
        #Move to the grid at the right
        if (counter_x == limit_x):
            #This means it completes calculate the particles in the grid in a row. Advance to the next row. 
            counter_x = 0
            counter_y += 1
        if (counter_y == limit_y):
            #This moves to the next layer of the xy-plane. 
            counter_y = 0
            counter_z += 1
    return num, particle_in_cell
                
def apply_boundary(p, Nd, Np):
    """
    This function applies the periodic boundary condition to the position of particles. 

    Parameters
    ----------
    p : array_like
        Position of all particle in the array.
    Nd : int
        Number of dimensions of the simulation.
    Np : int
        Number of particles in the simulation.

    Returns
    -------
    p : array_like
        The position of particles after applying the boundary condition.

    """
    # Modify to apply your chosen boundary conditions
    
    position_x = p[0,:]
    position_y = p[1,:]
    position_z = p[2,:]
    #The position of particles in the x, y and z position. 
    
    #The following lines will find the particle outside the simulation and move it back to the simulation
    #based on the periodic boundary condition.
    xx = np.where(position_x > x_bound, (position_x-x_bound)-x_bound, position_x)
    xxx = np.where(xx < -x_bound, x_bound - (-x_bound - xx), xx)
    yy = np.where(position_y > y_bound, (position_y-y_bound)-y_bound, position_y)
    yyy = np.where(yy < -y_bound, y_bound - (-y_bound - yy),yy)
    zz = np.where(position_z > z_bound, (position_z-z_bound)-z_bound, position_z)
    zzz = np.where(zz < -z_bound, z_bound - (-z_bound - zz), zz)
    
    p = np.concatenate((xxx, yyy, zzz), axis = 0)
    p = np.reshape(p, (Nd, Np))
    #Reconstruct the array for position of particles. 
    return p

def center_of_mass(particle_in_cell, num, mass, position):
    """
    This function calculates the center of mass of all particles in the same grid cell. 

    Parameters
    ----------
    particle_in_cell : List
        The list that contain the index of particles in each grid cell.
    num : List
        The list contains the number of particles in each grid cell. 
    mass : array_like
        The mass of all the particles.
    position : array_like
        The position of all the partiles.

    Returns
    -------
    result : List
        The center of mass position in each grid cell.
    total_mass : List
        The total mass of all particles in each grid cell.

    """
    result = []
    total_mass = []
    #Initialise the output lists
    position_x = position[0,:]
    position_y = position[1,:]
    position_z = position[2,:]
    
    for i in range(len(particle_in_cell)):
        COM_x = 0.0
        COM_y = 0.0
        COM_z = 0.0
        M_total = 0.0
        #Initialise the center of mass position and the total mass of particles in the grid
        for j in range(num[i]):
            COM_x += mass[particle_in_cell[i][j]]*position_x[particle_in_cell[i][j]]
            COM_y += mass[particle_in_cell[i][j]]*position_y[particle_in_cell[i][j]]
            COM_z += mass[particle_in_cell[i][j]]*position_z[particle_in_cell[i][j]]
            M_total += mass[particle_in_cell[i][j]]
            #Calculate the center off mass
        result.append(np.array([COM_x/M_total, COM_y/M_total, COM_z/M_total]))
        total_mass.append(M_total)
    return result, total_mass

def position2grid(particle_index, particle_in_cell):
    """
    This function matches the index of the particle to the index of grid cell it is in. 

    Parameters
    ----------
    particle_index : array_like
        Index of all particles in the simulation.
    particle_in_cell : List
        The index of particle in each grid cell.

    Returns
    -------
    result : List
        It matches the index of particle to the index of the grid it is in.

    """
    result = []
    for i in range(len(particle_index)):
        for j in range(len(particle_in_cell)):
            size = (np.intersect1d(np.array([i]), particle_in_cell[j]))
            #Determine whether the particle is in the grid[j]
            if (size.size > 0):#If the particle is in grid[j], the size of the array will be nonzero. 
            #Since the index of particle is also unique, we are certain that when the size of array is not zero. 
            #we find the cell of which the particle is in
                break
        result.append(np.array([i,j]))
    return result

def accel_grid(COM, total_mass, mass, index, p, particle, particle_in_cell, num, grid_length, smooth):
    """
    This uses the center of mass to calculate the acceleration of a particle. 

    Parameters
    ----------
    COM : List
        Center of mass of all grid cells.
    total_mass : List
        Mass of all particles in a single grid cell.
    mass : array_like
        The mass of each individual particle.
    index : int
        The index of the grid cell.
    p : array_like
        Position of all the particle.
    particle : int
        Index of the particle.
    particle_in_cell : List
        The list contains the index of particles in each grid cell.
    num : List
        Number of particles in each grid cell.
    grid_length : The length of the gridcell
        A reference length. If the distance between particle and center of mass of any grid is below the 
        reference length. We will calulate the interaction in particle-particle basis (similar to P3M).

    Returns
    -------
    float
        The acceleration in the x-direction.
    float
        The acceleration in the y-direction.
    float
        The acceleration in the z-direction.

    """
    G = 4.452*10**(-7)  #in unit of kpc^3/10^5 solar masses/Myr^2
    smooth_grid = grid_length #The smoothen scale which is set to length of the grid. 
    
    accel_x = 0.0
    accel_y = 0.0
    accel_z = 0.0
    
    for i in range(len(COM)):
        r_dash = np.sqrt((COM[i][0]-p[0, particle])**2 + (COM[i][1]-p[1, particle])**2 + (COM[i][2]-p[2,particle])**2)
        #The distance between the particle and the center of mass of particles.               
        if (r_dash <= grid_length):
            #If less than the grid size, calculate the force using individual particle.
            accel_x += accel_particle(p, particle, mass, i, particle_in_cell, smooth)[0]
            accel_y += accel_particle(p, particle, mass, i, particle_in_cell, smooth)[1]
            accel_z += accel_particle(p, particle, mass, i, particle_in_cell, smooth)[2]
        else:
            #Larger than the gridsize, calculate the force with center of mass.
            r = np.sqrt((COM[i][0]-p[0, particle])**2 + (COM[i][1]-p[1, particle])**2 + (COM[i][2]-p[2,particle])**2 + smooth_grid**2)
            accel_x += G*total_mass[i]*(COM[i][0]-p[0, particle])/r**3
            accel_y += G*total_mass[i]*(COM[i][1]-p[1, particle])/r**3
            accel_z += G*total_mass[i]*(COM[i][2]-p[2, particle])/r**3
    return accel_x, accel_y, accel_z

#Acceleration of the particles
def accel_particle(p, particle, mass, index, particle_in_cell, smooth):
    """
    This calculates the acceleration of particle on a particle-particle basis.

    Parameters
    ----------
    p : array_like
        The psition of all particles.
    particle : int
        The index of the particle.
    mass : array_like
        The mass of all particles.
    index : int
        The index of the grid of which the particle is in.
    particle_in_cell : List
        The index of particles in each grid cell.

    Returns
    -------
    accel_x: float
        The acceleration in the x-direction.
    accel_y: float
        The acceleration in the y-direction.
    accel_z: float
        The acceleration in the z-direction.

    """
    G = 4.452*10**(-7)  #in unit of kpc^3/10^5 solar masses/Myr^2
    #smooth = 1.0 #The smoothen scale is 100 pc which is bigger than the size of globular cluster (around 0.01 kpc, smallest possible
    #mass) and the size of a dwarf galaxy (around 1 kpc, largest possible mass)
    
    accel_x = 0.0
    accel_y = 0.0
    accel_z = 0.0
    
    total = particle_in_cell[index]
    #This is the collection of all particles in a specific grid. 
    
    for i in range(len(total)):
        if (total[i] != particle):
            #Calculate the force on the particle individually. 
            r = np.sqrt((p[0,total[i]]-p[0, particle])**2 + (p[1,total[i]]-p[1, particle])**2 + (p[2,total[i]]-p[2,particle])**2 + smooth**2)
            accel_x += G*mass[total[i]]*(p[0,total[i]]-p[0, particle])/r**3
            accel_y += G*mass[total[i]]*(p[1,total[i]]-p[1, particle])/r**3
            accel_z += G*mass[total[i]]*(p[2,total[i]]-p[2, particle])/r**3
            
    return accel_x, accel_y, accel_z

def recession_vel(position, H_0):
    """
    This calculates the recession velocity of the particle at a given position

    Parameters
    ----------
    position : array_like
        The position of every particle in the simulation
    H_0 : float
        The value of Hubble constant in km/s/Mpc

    Returns
    -------
    v_rec : array_like
        The recession velocity of the particle at each position.

    """
    v_rec = position*Hubble_convert(H_0) #return the recession velocity in kpc/Myr
    return v_rec

def Hubble_convert(H_0):
    """
    Converts the Hubble parameter from km/s/Mpc to Myr^-1

    Parameters
    ----------
    H_0 : float
        The Hubble parameter in km/s/Mpc.

    Returns
    -------
    result : float
        The Hubble parameter in Myr^-1.

    """
    result = H_0*1000.0*3.1536*10**13/(3.09*10**16)/10**6 #This formula convert the Hubble parameter from
    #km/s/Mpc to Myr^-1 in order to match the unit convention in this program
    return result

#Acceleration of the particles
def acceleration(p, num, Np, mass, smooth):
    """
    This uses the exact method to calculate the force on a particle. 

    Parameters
    ----------
    p : array_like
        Position of all particles.
    num : List
        The number of particle in each grid.
    Np : int
        Total number of particles.
    mass : array_like
        The mass of each particle.

    Returns
    -------
    accel_x: float
        The acceleration in the x-direction.
    accel_y: float
        The acceleration in the y-direction.
    accel_z: float
        The acceleration in the z-direction.

    """
    G = 4.452*10**(-7)  #in unit of kpc^3/10^5 solar masses/Myr^2
    accel_x = 0.0
    accel_y = 0.0
    accel_z = 0.0
    for i in range(Np):
        if (i != num):
            r = np.sqrt((p[0,i]-p[0, num])**2 + (p[1,i]-p[1, num])**2 + (p[2,i]-p[2,num])**2 + smooth**2)
            accel_x += G*mass[i]*(p[0,i]-p[0, num])/r**3
            accel_y += G*mass[i]*(p[1,i]-p[1, num])/r**3
            accel_z += G*mass[i]*(p[2,i]-p[2, num])/r**3
    return accel_x, accel_y, accel_z

def separation(p): 
    """
    This code clulates the seperation between two particles

    Parameters
    ----------
    p : array_like
        The position of all particles in the simulation.

    Returns
    -------
    float
        The separation between all particles in the simulation.

    """
    # Function to find separations from position vectors
    s = (p[:,None,:] - p[:,:,None]) # find N x N x Nd matrix of particle separations
    return np.sum(s**2,axis=0)**0.5 # return N x N matrix of scalar separations

#Basic user interface. Ask the user to input the following parameter
print('To modify the initial condition, please modify the code directly.')
mass_min = eval(input("The minimum nonzero resolution mass (in 10^5 solar masses): "))
mass_max = eval(input("The maximum mass of the dark matter halo (in 10^5 solar masses): "))
bound_xyz = eval(input("The length of simulation in the xy plane (in kpc): "))
method = eval(input("Evaluation method, enter 0 for grid approximation and 1 for exact solution: "))
grid_xyz = eval(input("Total number of grids in the x or y position: "))
Np = eval(input("Total number of particles: "))
Nt = eval(input("Total number of time steps: "))
dt = eval(input("Time step (in Myr): "))
v_max =  eval(input("The maximum drift velocity (in kpc/Myr): "))
H_0 = eval(input("The Hubble parameter (in km/s/Mpc). For static universe, enter 0.0: "))
smooth = eval(input("The soften length of the simulation (in kpc): "))

t_0 = time.time()
# For reproducibility, set a seed for randomly generated inputs. Change to your favourite integer.
np.random.seed(4080)

# Set the number of spatial dimensions (at least 2)
Nd = 3

# Set how long the animation should dispay each timestep (in milliseconds).
frame_duration = 100

#boundary for x, y and z
x_bound = bound_xyz
y_bound = bound_xyz
z_bound = bound_xyz

# Set initial positions at random within box
# position_xy = (np.random.normal(loc = 0.0, scale = 4, size = (2, Np)))
# position_z = np.random.random(size = (1, Np))*z_bound

# position = np.concatenate((position_xy, position_z), axis = 0)

#position = np.random.normal(loc = 0.0, scale = 4, size = (Nd, Np))

position_1 = np.random.normal(loc = -15.0, scale = 5, size = (3, 100))
position_2 = np.random.normal(loc = 15.0, scale = 5, size = (3, 100))

position = np.concatenate((position_1, position_2), axis = 1)


# Set initial velocities to be random fractions of the maximum
#velocity = v_max*(1-2*np.random.random((Nd,Np)))

velocity_1 = np.full((Nd, 100), v_max)
velocity_2 = np.full((Nd, 100), -v_max)
velocity = np.concatenate((velocity_1, velocity_2), axis = 1)

mass = 10**(np.random.random(size=(Np))*(np.log10(mass_max)-np.log10(mass_min)) + np.log10(mass_min)) 

position += velocity/2.0*dt #first step of Leapfrog method. 

position = apply_boundary(position, Nd, Np)#Apply the periodic boundary condition

position_new  = np.reshape(np.concatenate((position[0,:], position[1,:])), (2, Np)).T
#Position_new is a 2xNp matrix. The first column is the x position of particles and the second column
#is the y position. 

separation_max = np.sqrt(3.0)*x_bound*2.0
#The maximum separation possible in the simulation. 

# Set the axes on which the points will be shown
plt.ion() # Set interactive mode on
fig = figure(figsize=(12,6)) # Create frame and set size
subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95,wspace=0.15,hspace=0.2)

ax1 = subplot(121) # For normal 2D projection

ax1.set_xlabel('kpc')
ax1.set_ylabel('kpc')
ax1.set_title('Collision of two galaxies')

# # Create command which will plot the positions of the particles
scat = plt.scatter(position_new[:,0], position_new[:,1], s= (np.log10(mass))**2)


ax2 = subplot(222) # Create second set of axes as the top right panel in a 2x2 grid
xlabel('Separation (kpc)')
ylabel('Correlation function')
dx= 0.5 # Set width of x-axis bins
xb = np.arange(0, separation_max+dx,dx)
xb[0] = 1e-6 # Shift first bin edge by a fraction to avoid showing all the zeros (a cheat, but saves so much time!)
line, = ax2.plot([],[],drawstyle='steps-post', label='data') # Define a command that plots a line in this panel
smooth_line, =  ax2.plot([],[],drawstyle='steps-post',label='Spline') #This plots the spline interpolation of the 
#correlation function
ax2.legend(loc='upper right')


ax4 = plt.subplot(224) # Create last set of axes as the bottom right panel in a 2x2 grid
ax4.set_xscale('log')
ax4.set_yscale('log')
#Set both the x and y axis on a log scale. 
plane, = ax4.plot([], [], drawstyle = 'steps-post', label='data') #This plots the power spectrum
smooth_plane, = ax4.plot([], [], drawstyle = 'steps-post', label='spline') #This plots the spline interpolation 
#of the power spectrum. 
xlabel('Wavenumber (kpc^-1)')
ylabel('Power spectrum (kpc^3)')
ax4.legend(loc='best')


# Define procedure to update positions at each timestep
def update(i):
    global position,velocity, dx, mass, smooth # Get positions and velocities and bin width
    N = position.shape[1]
    year = i*dt #Shows how many year has passed since the initial condition. 
    scat.set_label('%lf Myrs'%year)
    ax1.legend(loc='upper right')#Display the time in the lower right corner.
    accel = np.empty(shape = position.shape)
    for i in range(N):
        accel[0, i], accel[1, i], accel[2, i] = acceleration(position, i, Np, mass, smooth)
        
    velocity += accel
    
    position += (velocity+recession_vel(position,H_0))*dt # Increment positions according to their velocites
    #The total velocity is the sum of the peculiar velocity and the recession velocity.
    position = apply_boundary(position, Nd, Np) # Apply boundary conditions
    ax1.set_xlim(-x_bound-x_bound*Hubble_convert(H_0)*year,x_bound+x_bound*Hubble_convert(H_0)*year)  # Set x-axis limits
    ax1.set_ylim(-y_bound-y_bound*Hubble_convert(H_0)*year,y_bound+y_bound*Hubble_convert(H_0)*year)  # Set y-axis limits
    #points.set_data(position[0,:], position[1,:]) # Show 2D projection of first 2 position coordinates
    scat.set_offsets(np.reshape(np.concatenate((position[0,:], position[1,:])), (2, Np)).T)#This line of code basically


    
    DD = np.ravel(tril(separation(position)))#The separation of all particles in the dataset. 
    
    factor = Np**2/((2*x_bound)**3) #This is the number density of pair of particles in the simulation. Since we use
    #periodic boundary condition, we can also consider our simulation in a sphere. 
    
    
    h_DD, x_DD = histogram(DD,bins=xb) #The number of pairs of galaxies in each bin. 
    
    h = np.zeros(len(h_DD)) #Correlation function
    x_max = 0.0#The mximum value in the x-axis
    for i in range(len(h_DD)):
        h[i] = h_DD[i]/((4.0/3.0*np.pi*(xb[i+1]**3-xb[i]**3)*factor))-1.0 # calculate the correlation function
        #using the estimator
        if (h[i] > 0):
            x_max = x_DD[i] #Find the largest separation where the correlation function is greater than 0. 
    
    line.set_data(x_DD[:-1],h) # Set the new data for the line in the 2nd panel
    ax2.set_xlim(0, x_max)
    ax2.set_ylim(-1, np.amax(h)+5)
    
    variable_x = x_DD[:-1]
    cs = CubicSpline(variable_x, h) #Find a spline interpolation between the bin and the correlation function
    x = np.linspace(xb[0], np.sqrt(3.0)*2.0*x_bound, num=10000)
    smooth_plot = cs(x) #Plot the correlation function with the spline interpolation. 
    smooth_line.set_data(x, smooth_plot)
    
    k = 2.0*np.pi/(x_DD[:-1]) #The wavenumber. 
    
    k_min = np.amin(k) #The minimum wavenumber
    k_max = 2.0*np.pi/smooth #the maximum wavenumber is fixed to the wavenumber of the soften length. This is 
    #because any scale below the soften length is in accurate. 
    
    k_order_max = int(np.floor(np.log10(k_max))) + 1
    k_order_min = int(np.ceil(np.log10(np.amin(k)))) - 1
    
    #all possible k are between 10^k_order_min  and 10^k_order_max
    
    N_op = k_order_max - k_order_min
    #Number of steps to go from k_order_min to k_order_max
    
    x_axis = []
    order = k_order_min
    
    #The following for loop will put each interval from 10^n to 10^(n+1) into 100 smaller intervals. This will 
    #help to smooth the power spectrum and make the spline interpolation easier later. 
    for i in range(N_op):
    
        segment = np.linspace(10**order, 10**(order+1), 100, endpoint=False)
        order += 1
    
        x_axis.append(segment)
    
      
    x_axis = np.array(x_axis)
    size_x, size_y = x_axis.shape

    x_axis = np.reshape(x_axis, (size_x*size_y))
    
    
    
    PS = [] #The power spectrum
    k_eff = [] #Wavenumbers between k_min and k_max
    for i in range(len(x_axis)):
        if ((x_axis[i] < k_min) or (x_axis[i] > k_max)):
            continue
        
        PS.append(np.trapz(variable_x**2*h*np.sin(x_axis[i]*variable_x)/(x_axis[i]*variable_x)*2.0*np.pi))
        k_eff.append(x_axis[i])
        
    PS = np.array(PS)
    k_eff = np.array(k_eff)
    
    ax4.set_xlim(k_min, k_max)
    ax4.set_ylim(1, np.amax(PS)) #Since the y axis is set to log scale. The minimum value of y cannot be less than zero.
    
    cs_ps = CubicSpline(k_eff, PS) #Spline interpolate the power spectrum. 
    k_final = np.linspace(k_min, k_max, 10000)
    PS_spline =  cs_ps(k_final) #The spline interpolation of the power spectrum. 
    
    smooth_plane.set_data(k_final, PS_spline)
    plane.set_data(k_eff,PS)
    return scat, plane, smooth_line, line, smooth_plane # Plot the points and the line
    
# Create animation
# https://matplotlib.org/api/_as_gen/matplotlib.animation.FuncAnimation.html
ani = animation.FuncAnimation(fig, update, frames=Nt,interval = frame_duration)
#plt.show()

ani.save("panels.mp4")

t_1 = time.time()
print(t_1-t_0)


