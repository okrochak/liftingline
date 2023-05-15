# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import functions as fcn

# Switch to TeX style for plots 
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})

airfoil = 'polarDU95W180.txt'
data1 = pd.read_csv(airfoil, header=0,
                    names=["alpha", "cl", "cd", "cm"],  sep='\s+')
polar_alpha = data1['alpha'][:]
polar_CL = data1['cl'][:]
polar_CD = data1['cd'][:]
delta_mu = 0.01
mu = np.arange(0.2, 1 + delta_mu / 2, delta_mu)
delta_Phi = np.pi/10
Phi = np.arange(0.2, 2*np.pi + delta_Phi / 2, delta_Phi)
TSR = 8

# Plot figures?
doPlot = True

# blade shape
R = 50
pitch = 2 # degrees
chord_distribution = 3 * (1 - mu) + 1 # meters
twist_distribution = 14 * (1 - mu) - pitch  # degrees

# define flow conditions
U_inf = 10  # unperturbed wind speed in m/s
N_B = 3 # number of blades

# Lifting line model inputs
Ncp = 15 # number of segments per blade = number of control points
psi = np.pi/2 # azimuthal position of the first rotor blade in radians
Loutlet = 1 # the distance from the rotor  to the domain boundary, in rotor diameters
dx = 0.1  # discretization distance for the vortex ring, in meters.


mu_tip = 1
mu_root = 0.2
'''Run the BEM model for the inputs '''
# TSR_distribution = np.linspace(4,10,24)
CT_distribution = []
CP_distribution = []
results = np.zeros([len(mu)-1, 11])

Omega = U_inf * TSR / R
# solve BEM model
for i in range(len(mu)-1):
        chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
        twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

        results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)

mu = results[:,0]
a_axial = results[:,1]
a_tan = results[:,2]
df_axial = results[:,5]
df_tan = results[:,6] 
alpha = results[:,10]
delta_r = (mu[1:] - mu[:-1]) * R
CT = np.sum(delta_r * df_axial[1:] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
CP = np.sum(delta_r * df_tan[1:] * mu[1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))
if doPlot == 1:
        fig1 = plt.figure(figsize=(8, 4))
        plt.title('Flow Field downstream of the rotor')
        plt.plot(mu, U_inf*(1-2*a_axial), 'k-', label=r'$U_{ax}$')
        plt.plot(mu, U_inf*2*a_tan, 'k--', label=r'$U_{tan}$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.savefig("figures/flowfield",bbox_inches='tight')

class VortexRing:
        def __init__(self, Rvec, theta, alpha, chord, Uax, Utan):   
                print("Initialized a new VortexRing object")
                self.r1 = Rvec[0] # position closest to the hub
                self.r2 = Rvec[1] # position closest to the tip
                self.theta = theta # angle in radians from the x-axis in the rotor plane
                self.alpha = alpha # angle of attack
                self.chord = chord # chord length
                self.Uax = Uax # downstream axial velocity
                self.Utan = Utan # tangential downstream velocity
                self.Gamma = 1 # circulation acting on this segment
                # define the downstream circulation lines
                self.coords_cyl = np.hstack((np.flip(fcn.downstreamLine(self.r1, theta, chord, alpha, Uax, Utan, R, Loutlet, dx),1), \
                                             fcn.downstreamLine(self.r2, theta, chord, alpha, Uax, Utan, R, Loutlet, dx)))
                #self.coords_cyl = fcn.downstreamLine(self.r1, theta, chord, alpha, Uax, Utan, R, Loutlet, dx)
                self.coords_cart = fcn.cyl2cart(self.coords_cyl)

### Discretize the rotor blades
print("Initialize the vortexRing system")
vortexSystem ={} 
vortexSystem["r_coord"] = np.linspace(mu_root*R, mu_tip*R,Ncp + 1)
vortexSystem["mu"] = (np.linspace(mu_root, mu_tip,Ncp + 1)[1:] + np.linspace(mu_root, mu_tip,Ncp + 1)[:-1]) / 2
vortexSystem["chord"] = 3 * (1 - vortexSystem["mu"]) + 1 # meters
vortexSystem["alpha"] = np.deg2rad(np.interp(vortexSystem["mu"],mu,alpha))
vortexSystem["U_ax"] = np.interp(vortexSystem["mu"],mu,U_inf*(1-2*a_axial))
vortexSystem["U_tan"] = np.interp(vortexSystem["mu"],mu,Omega*R*(1+2*a_tan))

for j in range(N_B):
        theta = psi + j*(2/3)*np.pi
        for i in range(Ncp):
                # Create a dictionary storing all vortex ring variables
                Rvec = np.array([vortexSystem["r_coord"][i],vortexSystem["r_coord"][i+1]])
                vortexSystem[f"inst{j}_{i}"] = VortexRing(Rvec, theta, vortexSystem["alpha"][i], vortexSystem["chord"][i], vortexSystem["U_ax"][i], vortexSystem["U_tan"][i])

if doPlot == 1:
        plt.figure(2)
        ax = plt.axes(projection='3d')
        for j in range(N_B):
                # for i in range(Ncp): 
                for i in range(Ncp-1,Ncp):
                        coords_plot = vortexSystem[f"inst{j}_{i}"].coords_cart
                        # Data for a three-dimensional line
                        ax.plot3D(coords_plot[0,:]/R, coords_plot[1,:]/R, coords_plot[2,:]/R, 'gray')
        plt.show()
# Data for three-dimensional scattered points

print("Solution converged")