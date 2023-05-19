# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import functions as fcn
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

# Switch to TeX style for plots 
# plt.rcParams['text.usetex'] = True
# plt.rcParams.update({'font.size': 18})

delta_mu = 0.01
mu = np.arange(0.2, 1 + delta_mu / 2, delta_mu)
TSR = 8

# Plot figures?
doPlot = 1

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
psi =  np.pi/11 # np.pi/3 # azimuthal position of the first rotor blade in radians
Loutlet = 1 # the distance from the rotor  to the domain boundary, in rotor diameters
dx = 0.5 # discretization time step for the vortex ring.
spacing = 0 # 0 - regular, 1 -  cosine
averageFactor = 1 # average the induction factors in radial direction?
rho = 1.225
mu_tip = 1
mu_root = 0.2
r_vortex = 1 #radius of solid body rotation
# Read the polar
airfoil = 'polarDU95W180.txt'
data1 = pd.read_csv(airfoil, header=0,
                    names=["alpha", "cl", "cd", "cm"],  sep='\s+')
polar_alpha = data1['alpha'][:]
polar_CL = data1['cl'][:]
polar_CD = data1['cd'][:]

'''1. Run the BEM model for the inputs '''
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
                twist = np.deg2rad(14 * (1 - (Rvec)/R) - pitch)  # twist in rad
                self.twist = np.mean(twist)
                # define the downstream circulation lines
                self.coords_cyl = np.hstack((np.flip(fcn.downstreamLine(Rvec[0], theta, chordVec[0], twist[0], UaxVec[0], UtanVec[0], R, Loutlet, dx),1), \
                                             fcn.downstreamLine(Rvec[1], theta, chordVec[1], twist[1], UaxVec[1], UtanVec[1], R, Loutlet, dx)))
                self.coords_cart = fcn.cyl2cart(self.coords_cyl)
                self.control_point = fcn.cyl2cart(np.array([[0*0.5*(self.chord[0] + self.chord[1])], [0.5*(self.r[0] + self.r[1])], [theta] ]))
                self.blade = fcn.cyl2cart(np.array([[-0.25 * self.chord[0],-0.25 * self.chord[1], 0.75 * self.chord[1], 0.75 * self.chord[0]], \
                                      [self.r[0], self.r[1], self.r[1], self.r[0]], [theta, theta, theta, theta]])) 
                radVec = fcn.cyl2cart(np.array([0, self.r[1], theta])) - fcn.cyl2cart(np.array([0, self.r[0], theta])) #  vector in radial direction x r theta
                
                theta1c = theta - np.sin(twist)*np.mean(chordVec) / np.mean(Rvec) # theta coordinate of the TE.
                
                # chordVec = fcn.cyl2cart(np.array([np.cos(twist)*np.mean(chordVec), np.mean(Rvec), theta1c])) - fcn.cyl2cart(np.array([0, np.mean(Rvec), theta])) 
                self.radVec = radVec/np.linalg.norm(radVec) # unit vector in radial direction in system (turbine) coordinates
                # self.chordVec = chordVec/np.linalg.norm(chordVec) #  unit vector aligned with the local chord
                self.xVec = np.array([1, 0, 0]) # unit vector in x-direction
                self.thetaVec = np.cross(self.radVec.flatten(),self.xVec.flatten())  # unit vector aligned with the local angle theta
### Discretize the rotor blades
print("Initialize the vortexRing system")
vortexSystem ={} 
if spacing == 0:
        # regular spacing
        vortexSystem["mu_coord"] = np.linspace(mu_root, mu_tip,Ncp + 1)
elif spacing == 1:
       # cosine spacing
        vortexSystem["mu_coord"] = mu_root + ((1 + np.flip(np.cos(np.linspace(0,np.pi,Ncp+1)))) / 2) * (mu_tip - mu_root)
vortexSystem["mu"] = (vortexSystem["mu_coord"][1:] + vortexSystem["mu_coord"][:-1]) / 2

vortexSystem["chord"] = np.interp(vortexSystem["mu_coord"],mu,chord_distribution)
vortexSystem["alpha"] = np.deg2rad(np.interp(vortexSystem["mu_coord"],mu,alpha))
vortexSystem["U_ax"] = np.interp(vortexSystem["mu_coord"],mu,U_inf*(1-a_axial))
vortexSystem["U_tan"] = np.interp(vortexSystem["mu_coord"],mu,Omega*mu*R*(1+a_tan))


'''2. Assemble the vortex ring system for the wake geometry '''
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

# this must be close to zero, orthogonal vectors
print("Dot product of r_hat and x_hat:" , np.dot(vortexSystem[f"inst{j}_{i}"].radVec.flatten(), vortexSystem[f"inst{j}_{i}"].xVec.flatten()))

'''3. Assemble the induction matrix'''
controlPoints = {}
controlPoints["coords"] = np.empty((N_B*Ncp,3))
controlPoints["gamma"] = np.ones((N_B*Ncp))
controlPoints["Uin"] = np.empty((3,N_B*Ncp))
controlPoints["Uin_blade"] = np.empty((3,N_B*Ncp))
controlPoints["matrix"] = np.zeros((N_B*Ncp,N_B*Ncp,3))
controlPoints["r_hat"] = np.empty((N_B*Ncp,3))
controlPoints["theta_hat"] = np.empty((N_B*Ncp,3))
controlPoints["chord"] = np.empty((N_B*Ncp))
controlPoints["twist"] = np.empty((N_B*Ncp))
controlPoints["r"] = np.empty((N_B*Ncp))

for j1 in range(N_B):
        for i1 in range(Ncp):
                c1 = j1*Ncp+i1 # row index in induction matrix, i.e control point index
                controlPoints["coords"][c1,:] = vortexSystem[f"inst{j1}_{i1}"].control_point.flatten() # assemble global control points matrix.
                controlPoints["r_hat"][c1,:] = vortexSystem[f"inst{j1}_{i1}"].radVec.flatten()
                # controlPoints["x_hat"][c1,:] = vortexSystem[f"inst{j1}_{i1}"].xVec.flatten()
                controlPoints["theta_hat"][c1,:] = vortexSystem[f"inst{j1}_{i1}"].thetaVec.flatten()
                controlPoints["chord"][c1] = np.mean(vortexSystem[f"inst{j1}_{i1}"].chord)
                controlPoints["twist"][c1] = np.mean(vortexSystem[f"inst{j1}_{i1}"].twist)
                controlPoints["r"][c1] = np.mean(vortexSystem[f"inst{j1}_{i1}"].r)
                # now calculate induced velocity from each element
                for j2 in range(N_B):
                        for i2 in range(Ncp):
                                c2 = j2*Ncp+i2 # column index in induction matrix , i.e., vortexRing object 
                                coords_cp = controlPoints["coords"][c1,:]
                                vort_fil = vortexSystem[f"inst{j2}_{i2}"].coords_cart
                                # r = vortexSystem[f"inst{j1}_{i1}"].coords_cyl[1,0]
                                for n in range(vort_fil.shape[1]-1):
                                        x1 = vort_fil[:,n]; x2 = vort_fil[:,n+1]
                                        # finally assemble the induction matrix
                                        controlPoints["matrix"][c1,c2,:] += fcn.velocity3d_vortex_filament(1, x1, x2, coords_cp, r_vortex)
print("The Induction matrix is assembled")

'''4. Begin the iteration loop '''
diff = 1000
tol = 0.01
iter = 0
controlPoints["gamma"] = np.ones((N_B*Ncp)) #-3 works well?
circulation_history = np.empty((N_B*Ncp,0))
while diff > tol:
        # Calculate the new induction velocities
        print("Iteration number:", iter)
        controlPoints["Uin"][0,:] = controlPoints["matrix"][:,:,0]@controlPoints["gamma"] # x-velocity
        controlPoints["Uin"][1,:] = controlPoints["matrix"][:,:,2]@controlPoints["gamma"] # y-velocity
        controlPoints["Uin"][2,:] = controlPoints["matrix"][:,:,1]@controlPoints["gamma"] # z-velocity, all in turbine (global) frame of reference

        '''Calculate the new circulation (lift) from the induced velocity'''
        # First decompose the induced velocities vector in the components in the blade local coordinate system
        controlPoints["Uin_blade"][0,:] = U_inf + controlPoints["Uin"][0,:]   #  x-direction velocity acting on the blade
        controlPoints["Uin_blade"][1,:] = Omega * controlPoints["r"] + np.einsum('ij,ji->i',controlPoints["theta_hat"],controlPoints["Uin"])# azimuthal velocity

        # Now we can easily calculate the magnitude of the in-plane velocity
        controlPoints["|V|bl"] = np.sqrt(((controlPoints["Uin_blade"][0,:]**2)+(controlPoints["Uin_blade"][1,:]**2)))
        # And the angle of attack
        controlPoints["alpha"] = np.arctan(controlPoints["Uin_blade"][0,:]/controlPoints["Uin_blade"][1,:]) - controlPoints["twist"] # phi - twist
        # Update the Gamma 

        CL = np.interp(np.rad2deg(controlPoints["alpha"]), polar_alpha, polar_CL)
        controlPoints["gamma_upd"] = 0.5*CL*(controlPoints["|V|bl"])*controlPoints["chord"]
        diff = np.mean(np.abs(controlPoints["gamma_upd"] - controlPoints["gamma"])) # record the difference
        controlPoints["gamma"] = 0.5*(controlPoints["gamma_upd"] + controlPoints["gamma"])
        circulation_history = np.append(circulation_history,controlPoints["gamma"][:,np.newaxis],axis=1)
        iter += 1

# Validation plots
print("Solution converged")
fig1 = plt.figure(figsize=(8, 4))
plt.title('Converged circulation solution')
fac = np.pi * U_inf**2 / (Omega*N_B)
plt.plot(controlPoints['r'][0:Ncp]/R,controlPoints['gamma'][0:Ncp]/fac,label=r'$\Gamma_1 /(\pi U^2 / \Omega N_B))$')
plt.plot(controlPoints['r'][Ncp:2*Ncp]/R,controlPoints['gamma'][0:Ncp]/fac,label=r'$\Gamma_2 /(\pi U^2 / \Omega N_B))$')
plt.plot(controlPoints['r'][2*Ncp:3*Ncp]/R,controlPoints['gamma'][0:Ncp]/fac,label=r'$\Gamma_3 /(\pi U^2 / \Omega N_B))$')

plt.xlabel(r'$r/R$')
plt.legend()
plt.grid()
plt.savefig("figures/circulation_converged",bbox_inches='tight')
# should be periodic and the same, as the norm of a vector is the first invariant :)
# plt.plot(np.linalg.norm(controlPoints["Uin"],axis=0))
# plt.plot(np.linalg.norm(controlPoints["Uin_blade"],axis=0))

fig2 = plt.figure(figsize=(8, 4))
plt.title('Converged circulation solution')
plt.plot(controlPoints['r'][0:Ncp]/R,np.rad2deg(controlPoints['alpha'][0:Ncp]),label=r'$\alpha [deg]$')
plt.xlabel(r'$r/R$')
plt.legend()
plt.grid()
plt.show()

#fig3 = plt.figure(2,constrained_layout=True)
#ax = plt.axes(projection='3d')
#ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
#          controlPoints["r_hat"][:,0], controlPoints["r_hat"][:,1], controlPoints["r_hat"][:,2])

if doPlot == 1:
        fig4 = plt.figure(figsize=(8, 4))
        plt.title('Flow Field downstream of the rotor')
        plt.plot(vortexSystem["mu_coord"], vortexSystem["U_ax"] , 'k-', label=r'$U_{ax}$')
        plt.plot(vortexSystem["mu_coord"], vortexSystem["U_tan"] , 'k--', label=r'$U_{tan}$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.ylabel(r'$U$')
        plt.legend()
        plt.savefig("figures/flowfield",bbox_inches='tight')
        
        fig5 = plt.figure(2,constrained_layout=False,figsize=(8, 4))
        plt.title('Wake Geometry')
        ax = plt.axes(projection='3d')
        for j in range(N_B):
                # setup blade geometry
                for i in range(Ncp): 
                        coords_plot = vortexSystem[f"inst{j}_{i}"].coords_cart
                        blade_plot = vortexSystem[f"inst{j}_{i}"].blade
                        # Data for a three-dimensional line
                        ax.plot3D(coords_plot[0,:], coords_plot[1,:], coords_plot[2,:], 'gray')

                        # print(i,j)
                        # print(blade_plot)

                        ax.plot_trisurf(blade_plot[0,:], blade_plot[1,:], blade_plot[2,:], linewidths=0,edgecolor='Gray', color='gray')
                        # ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
                        #          controlPoints["r_hat"][:,0], controlPoints["r_hat"][:,1], controlPoints["r_hat"][:,2],color = 'red')
                        # ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
                        #          controlPoints["c_hat"][:,0], controlPoints["c_hat"][:,1], controlPoints["c_hat"][:,2],color = 'black')
                        ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
                                 controlPoints["theta_hat"][:,0]*5, controlPoints["theta_hat"][:,1]*5, controlPoints["theta_hat"][:,2]*5,color = 'blue')
                        ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
                                controlPoints["Uin"][0,:]*0, controlPoints["Uin"][1,:], controlPoints["Uin"][2,:],color = 'black')
                        ax.quiver(np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]), \
                                 np.array([0,0,3]),np.array([0,3,0]),np.array([3,0,0]),color = 'green')
                        # compute the coordinates of the blades
                        fig5.subplots_adjust(left=0, right=1, bottom=0, top=1)

        ax.scatter(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2])
        # the beam
        ax.plot3D([0,0], [0, 0], [0, -50], 'gray', linewidth=10)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')
        # ax.set_ylim([-R,R])
        # ax.set_zlim([-R,R])
        # ax.set_ylim(0, R/3)
        # ax.set_zlim(0, R/3)
        # ax.set_xlim(-R/6, R/6)
        # ax.set_box_aspect((5, 1, 1))
        plt.savefig("figures/wake",bbox_inches='tight')
        plt.show()

