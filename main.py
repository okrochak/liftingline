# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import functions as fcn
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

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
Ncp = 9 # number of segments per blade = number of control points
psi = np.pi/10 # np.pi/3 # azimuthal position of the first rotor blade in radians
Loutlet = 1# the distance from the rotor  to the domain boundary, in rotor diameters
dx = 0.5  # discretization distance for the vortex ring, in meters.
spacing = 1 # 0 - regular, 1 -  cosine
averageFactor = 1 # average the induction factors in radial direction?

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
if averageFactor == 1:
        a_axial[:] = np.mean(a_axial); a_tan[:] = np.mean(a_tan)
chord_distribution = 3 * (1 - mu) + 1 # meters
delta_r = (mu[1:] - mu[:-1]) * R
CT = np.sum(delta_r * df_axial[1:] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
CP = np.sum(delta_r * df_tan[1:] * mu[1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))


class VortexRing:
        def __init__(self, Rvec, theta, alphaVec, chordVec, UaxVec, UtanVec):   
                print("Initialized a new VortexRing object")
                self.r = Rvec # position closest to the hub / closest to the tip
                self.theta = theta
                self.alpha = alphaVec # angle of attack
                self.chord = chordVec # chord length
                self.Uax = UaxVec # downstream axial velocity
                self.Utan = UtanVec # tangential downstream velocity
                self.Gamma = 1 # circulation acting on this segment
                # define the downstream circulation lines
                self.coords_cyl = np.hstack((np.flip(fcn.downstreamLine(Rvec[0], theta, chordVec[0], alphaVec[0], UaxVec[0], UtanVec[0], R, Loutlet, dx),1), \
                                             fcn.downstreamLine(Rvec[1], theta, chordVec[1], alphaVec[1], UaxVec[1], UtanVec[1], R, Loutlet, dx)))
                self.coords_cart = fcn.cyl2cart(self.coords_cyl)
                self.control_point = fcn.cyl2cart(np.array([[0.5*0.5*(self.chord[0] + self.chord[1])], [0.5*(self.r[0] + self.r[1])], [theta] ]))
                self.blade = fcn.cyl2cart(np.array([[-0.25 * self.chord[0],-0.25 * self.chord[1], 0.75 * self.chord[1], 0.75 * self.chord[0]], \
                                      [self.r[0], self.r[1], self.r[1], self.r[0]], [theta, theta, theta, theta]])) 

### Discretize the rotor blades
print("Initialize the vortexRing system")
vortexSystem ={} 
if spacing == 0:
        # regular spacing
        vortexSystem["mu_coord"] = np.linspace(mu_root*R, mu_tip*R,Ncp + 1)
elif spacing == 1:
       # cosine spacing
        vortexSystem["mu_coord"] = mu_root + ((1 + np.flip(np.cos(np.linspace(0,np.pi,Ncp+1)))) / 2) * (mu_tip - mu_root)
vortexSystem["mu"] = (vortexSystem["mu_coord"][1:] + vortexSystem["mu_coord"][:-1]) / 2

vortexSystem["chord"] = np.interp(vortexSystem["mu_coord"],mu,chord_distribution)
vortexSystem["alpha"] = np.deg2rad(np.interp(vortexSystem["mu_coord"],mu,alpha))
vortexSystem["U_ax"] = np.interp(vortexSystem["mu_coord"],mu,U_inf*(1-2*a_axial))
vortexSystem["U_tan"] = np.interp(vortexSystem["mu_coord"],mu,Omega*mu*R*(1+2*a_tan))



# Create the system of vortex filaments with unit circulation
for j in range(N_B):
        theta = psi + j*(2/3)*np.pi
        for i in range(Ncp):
                # Create a dictionary storing all vortex ring variables
                Rvec = np.array([vortexSystem["mu_coord"][i],vortexSystem["mu_coord"][i+1]])*R
                alphaVec = np.array([vortexSystem["alpha"][i],vortexSystem["alpha"][i+1]])
                chordVec = np.array([vortexSystem["chord"][i],vortexSystem["chord"][i+1]])
                UaxVec = np.array([vortexSystem["U_ax"][i],vortexSystem["U_ax"][i+1]])
                UtanVec = np.array([vortexSystem["U_tan"][i],vortexSystem["U_tan"][i+1]])
                vortexSystem[f"inst{j}_{i}"] = VortexRing(Rvec, theta, alphaVec, chordVec, UaxVec, UtanVec)


# Assemble the induction matrix
controlPoints = {}
controlPoints["coords"] = np.empty((3*Ncp,3))
controlPoints["gamma"] = np.ones((3*Ncp))
controlPoints["vel_ind"] = np.zeros((3*Ncp,3*Ncp,3))
for j1 in range(N_B):
        for i1 in range(Ncp):
                c1 = j1*Ncp+i1 # row index in induction matrix
                controlPoints["coords"][c1,:] = vortexSystem[f"inst{j1}_{i1}"].control_point.flatten() # assemble global control points matrix.
                # now calculate induced velocity from each element
                for j2 in range(N_B):
                        for i2 in range(Ncp):
                                c2 = j2*Ncp+i2 # column index in induction matrix
                                coords_cp = controlPoints["coords"][j1*Ncp+i1,:]
                                vort_fil = vortexSystem[f"inst{j1}_{i1}"].coords_cart
                                r = vortexSystem[f"inst{j1}_{i1}"].coords_cyl[1,0]
                                for n in range(vort_fil.shape[1]-1):
                                        x1 = vort_fil[:,n]; x2 = vort_fil[:,n+1]
                                        # finally assemble the induction matrix
                                        controlPoints["vel_ind"][c1,c2,:] += fcn.velocity3d_vortex_filament(controlPoints["gamma"][c2], x1, x2, coords_cp, r)


if doPlot == 1:

        fig1 = plt.figure(figsize=(8, 4))
        plt.title('Flow Field downstream of the rotor')
        plt.plot(vortexSystem["mu_coord"], vortexSystem["U_ax"] , 'k-', label=r'$U_{ax}$')
        plt.plot(vortexSystem["mu_coord"], vortexSystem["U_tan"] , 'k--', label=r'$U_{tan}$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.savefig("figures/flowfield",bbox_inches='tight')
        
        plt.figure(2,constrained_layout=False)
        ax = plt.axes(projection='3d')
        for j in range(N_B):
                # setup blade geometry
                for i in range(Ncp): 
                #for i in range(Ncp-1,Ncp):
                        coords_plot = vortexSystem[f"inst{j}_{i}"].coords_cart
                        blade_plot = vortexSystem[f"inst{j}_{i}"].blade
                        # Data for a three-dimensional line
                        ax.plot3D(coords_plot[0,:]/R, coords_plot[1,:]/R, coords_plot[2,:]/R, 'gray')
                        print(i,j)
                        print(blade_plot)

                        ax.plot_trisurf(blade_plot[0,:]/R, blade_plot[1,:]/R, blade_plot[2,:]/R, linewidths=0,edgecolor='Gray', color='gray')
                        # ax.plot_trisurf(Poly3DCollection(blade_plot, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.20))

                        # compute the coordinates of the blades
        ax.scatter(controlPoints["coords"][:,0]/R,controlPoints["coords"][:,1]/R,controlPoints["coords"][:,2]/R)
        # the beam
        ax.plot3D([0,0], [0, 0], [0, -1], 'gray', linewidth=10)
        ax.set_box_aspect((5, 1, 1))
        plt.show()

print("The Induction matrix is assembled")
print("Solution converged")