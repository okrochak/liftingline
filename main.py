# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import functions as fcn
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.cm import ScalarMappable

# Switch to TeX style for plots 
##plt.rcParams['text.usetex'] = True
#plt.rcParams.update({'font.size': 18})

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
Ncp = 10 # number of segments per blade = number of control points, don't go above 25
psi =  np.pi/11 # np.pi/3 # azimuthal position of the first rotor blade in radians
Loutlet = 2 # the distance from the rotor to the domain boundary, in rotor diameters
dx = 0.15 # discretization time step for the vortex ring.
spacing = 1 # 0 - regular, 1 -  cosine
averageFactor = 0 # average the induction factors in radial direction? 0 - no, 1 - yes
convFac = 0.95 # what part of original guess to take - higher is more stable ,  
rho = 1.225
mu_tip = 1
mu_root = 0.2
r_vortex = 0.4 #radius of solid body rotation - too small = unstable? I think must be above dx
# Read the polar
airfoil = 'polarDU95W180.txt'
data1 = pd.read_csv(airfoil, header=0,
                    names=["alpha", "cl", "cd", "cm"],  sep='\s+')
polar_alpha = data1['alpha'][:]
polar_CL = data1['cl'][:]
polar_CD = data1['cd'][:]

'''1. Run the BEM model for the inputs '''
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
phi_BEM = results[:,4]
df_axial_BEM = results[:,5]
df_tan_BEM = results[:,6] 
gamma_BEM = results[:,7] 
alpha_BEM = results[:,10]
chord_distribution = 3 * (1 - mu) + 1 # meters
delta_r = (mu[1:] - mu[:-1]) * R
CT = np.sum(delta_r * df_axial_BEM[1:] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
print('CT BEM = ' + str(CT))
CP = np.sum(delta_r * df_tan_BEM[1:] * mu[1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))
print('CP BEM = ' + str(CP))


# Present options to user
print('1 - Lifting Line Solution')
print('2 - Variable Convection Speed of Frozen Wake')
print('3 - Blade Discretisation - Regular/Cosine')
print('4 - Azimuthal Discretisation')
print('5 - Wake Length')

choice = input('Select an option: ')






if choice == '1':
        # Compute Lifting Line Solution
        [vortexSystem,controlPoints,df_axial,df_tan] = fcn.lifting_line(averageFactor,spacing,mu,mu_root,mu_tip,chord_distribution,R,Ncp,alpha_BEM,U_inf,a_axial,Omega,a_tan,psi,polar_alpha,polar_CL,polar_CD,N_B,convFac,r_vortex,pitch,Loutlet,dx)


        # Validation plots
        print("Solution converged")
        #fig1 = plt.figure(figsize=(8, 4))
        #plt.title('Converged circulation solution')
        #fac = np.pi * U_inf**2 / (Omega*N_B)
        #plt.plot(controlPoints['r'][0:Ncp]/R,controlPoints['gamma'][0:Ncp]/fac,label=r'$\Gamma_1 /(\pi U^2 / \Omega N_B))$')
        #plt.plot(controlPoints['r'][Ncp:2*Ncp]/R,controlPoints['gamma'][0:Ncp]/fac,label=r'$\Gamma_2 /(\pi U^2 / \Omega N_B))$')
        #plt.plot(controlPoints['r'][2*Ncp:3*Ncp]/R,controlPoints['gamma'][0:Ncp]/fac,label=r'$\Gamma_3 /(\pi U^2 / \Omega N_B))$')
        #plt.plot(mu, gamma_BEM/fac, label='BEM Solution')
        #plt.xlabel(r'$r/R$')
        #plt.legend()
        #plt.grid()
        #plt.savefig("figures/circulation_converged",bbox_inches='tight')
        # should be periodic and the same, as the norm of a vector is the first invariant :)
        # plt.plot(np.linalg.norm(controlPoints["Uin"],axis=0))
        # plt.plot(np.linalg.norm(controlPoints["Uin_blade"],axis=0))

        mu_LLT = controlPoints['r'][0:Ncp]/R

        fig1 = plt.figure(figsize=(8, 4))
        plt.title(r'Angle of Attack and Inflow Angle for $\lambda=$'+str(TSR))
        plt.plot(mu_LLT, np.rad2deg(controlPoints['alpha'][0:Ncp]), '-k',label=r'$\alpha_{LLT}$')
        plt.plot(mu, alpha_BEM, '--k', label=r'$\alpha_{BEM}$')
        plt.plot(mu_LLT, np.rad2deg(controlPoints['phi'][0:Ncp]), '-r', label=r'$\phi_{LLT}$')
        plt.plot(mu, np.rad2deg(phi_BEM), '--r', label=r'$\phi_{BEM}$')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        fig2 = plt.figure(figsize=(8, 4))
        plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega N_B}$ for $\lambda=$'+str(TSR))
        fac = np.pi * U_inf**2 / (Omega*N_B)
        plt.plot(mu_LLT,controlPoints['gamma'][0:Ncp]/fac, '-k', label='LLT Solution')
        plt.plot(mu, gamma_BEM/fac, '--k', label='BEM Solution')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        fig3 = plt.figure(figsize=(8, 4))
        plt.title(r'Thrust and Azimuthal Loading, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$ for $\lambda=$'+str(TSR))
        fac = 0.5 * rho * U_inf**2 * R
        plt.plot(mu_LLT, df_axial[0:Ncp] / fac,'-k', label=r'$dT_{LLT}$')
        plt.plot(mu, df_axial_BEM / fac, '--k', label=r'$dT_{BEM}$')
        plt.plot(mu_LLT, df_tan[0:Ncp] / fac, '-r', label=r'$dQ_{LLT}$')
        plt.plot(mu, df_tan_BEM / fac, '--r', label=r'$dQ_{BEM}$')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        #plt.show()

        print('CT BEM = ' + str(CT))
        print('CP BEM = ' + str(CP))

        delta_r = (mu_LLT[1:] - mu_LLT[:-1]) * R
        CT = np.sum(delta_r * df_axial[1:Ncp] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP = np.sum(delta_r * df_tan[1:Ncp] * mu_LLT[1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))

        print('CT LLT = ' + str(CT))
        print('CP LLT = ' + str(CP))

        #fig3 = plt.figure(2,constrained_layout=True)
        #ax = plt.axes(projection='3d')
        #ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
        #          controlPoints["r_hat"][:,0], controlPoints["r_hat"][:,1], controlPoints["r_hat"][:,2])

        #fig3 = plt.figure(2,constrained_layout=True)
        #ax = plt.axes(projection='3d')
        #ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
        #          controlPoints["r_hat"][:,0], controlPoints["r_hat"][:,1], controlPoints["r_hat"][:,2])
        gamma = controlPoints["gamma"]/ fac
        gamma_norm = (gamma - gamma.min()) / (gamma.max() - gamma.min())
        # gamma_norm = gamma / fac
        if doPlot == 1:
                # fig4 = plt.figure(figsize=(8, 4))
                # plt.title('Flow Field downstream of the rotor')
                # plt.plot(vortexSystem["mu_coord"], vortexSystem["U_ax"] , 'k-', label=r'$U_{ax}$')
                # plt.plot(vortexSystem["mu_coord"], vortexSystem["U_tan"] , 'k--', label=r'$U_{tan}$')
                # plt.grid()
                # plt.xlabel(r'$r/R$')
                # plt.ylabel(r'$U$')
                # plt.legend()
                # plt.savefig("figures/flowfield",bbox_inches='tight')
                
                fig5 = plt.figure(5,constrained_layout=False,figsize=(8, 4))
                plt.title('Wake Geometry')
                colormap = plt.cm.get_cmap('viridis')
                colors = colormap(gamma_norm)

                ax = plt.axes(projection='3d')
                for j in range(N_B):
                        # setup blade geometry
                        for i in range(Ncp): 
                                coords_plot = vortexSystem[f"inst{j}_{i}"].coords_cart
                                blade_plot = vortexSystem[f"inst{j}_{i}"].blade
                                # Data for a three-dimensional line
                                line = ax.plot3D(coords_plot[0,:], coords_plot[1,:], coords_plot[2,:], c=colors[Ncp*j+i])
                        
                                # print(i,j)
                                # print(blade_plot)

                                ax.plot_trisurf(blade_plot[0,:], blade_plot[1,:], blade_plot[2,:], linewidths=0,edgecolor='Gray', color='gray')
                                # ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
                                #          controlPoints["r_hat"][:,0], controlPoints["r_hat"][:,1], controlPoints["r_hat"][:,2],color = 'red')
                                # ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
                                #          controlPoints["c_hat"][:,0], controlPoints["c_hat"][:,1], controlPoints["c_hat"][:,2],color = 'black')
                                # ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
                                #          controlPoints["theta_hat"][:,0]*5, controlPoints["theta_hat"][:,1]*5, controlPoints["theta_hat"][:,2]*5,color = 'blue')
                                # ax.quiver(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2], \
                                #         controlPoints["Uin"][0,:]*0, controlPoints["Uin"][1,:], controlPoints["Uin"][2,:],color = 'black')
                                # ax.quiver(np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]), \
                                #          np.array([0,0,3]),np.array([0,3,0]),np.array([3,0,0]),color = 'green')
                                # compute the coordinates of the blades
                                fig5.subplots_adjust(left=0, right=1, bottom=0, top=1)
                scalar_map = ScalarMappable(cmap=colormap)
                scalar_map.set_array(colors)
                cbar = plt.colorbar(scalar_map, shrink=0.6)
                custom_labels = [round(np.min(gamma),2), round(0.5*(np.max(gamma)+np.min(gamma)),2), round(np.max(gamma),2)]  # Custom labels for colorbar ticks
                cbar.set_ticks([0, 0.5, 1])  # Set the tick positions
                cbar.set_ticklabels(custom_labels)
                cbar.set_label(r'$\Gamma / (\pi U_{\infty}^2 / \Omega N_B  )$')

                ax.scatter(controlPoints["coords"][:,0],controlPoints["coords"][:,1],controlPoints["coords"][:,2],color='red')
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
                plt.savefig("figures/wake",bbox_inches='tight',dpi = 300)
                plt.show()






if choice == '2':
        averageFactor = [0,1] # induction factors average - 0: non-averaged, 1: averaged

        vortex = []
        control = []
        axial = []
        tan = []
        mu_LL = []

        for i in range(len(averageFactor)):
                
                # Compute Lifting Line Solution
                [vortexSystem,controlPoints,df_axial,df_tan] = fcn.lifting_line(averageFactor[i],spacing,mu,mu_root,mu_tip,chord_distribution,R,Ncp,alpha_BEM,U_inf,a_axial,Omega,a_tan,psi,polar_alpha,polar_CL,polar_CD,N_B,convFac,r_vortex,pitch,Loutlet,dx)
                vortex.append(vortexSystem)
                control.append(controlPoints)
                axial.append(df_axial)
                tan.append(df_tan)

                # Validation plots
                print('Solution ' + str(i+1) + ' converged')

                mu_LLT = control[i]['r'][0:Ncp]/R
                mu_LL.append(mu_LLT)

        fig1 = plt.figure(figsize=(8, 4))
        plt.title(r'Angle of Attack and Inflow Angle for $\lambda=$'+str(TSR))
        plt.plot(mu_LL[0], np.rad2deg(control[0]['alpha'][0:Ncp]), '-r',label=r'$\alpha_{LLT} - non-averaged$')
        plt.plot(mu_LL[1], np.rad2deg(control[1]['alpha'][0:Ncp]), '-g',label=r'$\alpha_{LLT} - averaged$')
        plt.plot(mu, alpha_BEM, '-k', label=r'$\alpha_{BEM}$')
        plt.plot(mu_LL[0], np.rad2deg(control[0]['phi'][0:Ncp]), '--r', label=r'$\phi_{LLT} - non-averaged$')
        plt.plot(mu_LL[1], np.rad2deg(control[1]['phi'][0:Ncp]), '--g', label=r'$\phi_{LLT} - averaged$')
        plt.plot(mu, np.rad2deg(phi_BEM), '--k', label=r'$\phi_{BEM}$')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        fig2 = plt.figure(figsize=(8, 4))
        plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega N_B}$ for $\lambda=$'+str(TSR))
        fac = np.pi * U_inf**2 / (Omega*N_B)
        plt.plot(mu_LL[0],control[0]['gamma'][0:Ncp]/fac, '-r', label='LLT - non-averaged')
        plt.plot(mu_LL[1],control[1]['gamma'][0:Ncp]/fac, '-g', label='LLT - averaged')
        plt.plot(mu, gamma_BEM/fac, '-k', label='BEM Solution')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        fig3 = plt.figure(figsize=(8, 4))
        plt.title(r'Thrust and Azimuthal Loading, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$ for $\lambda=$'+str(TSR))
        fac = 0.5 * rho * U_inf**2 * R
        plt.plot(mu_LL[0], axial[0][0:Ncp] / fac,'-r', label=r'$dT_{LLT}$ - non-averaged')
        plt.plot(mu_LL[1], axial[1][0:Ncp] / fac,'-g', label=r'$dT_{LLT}$ - averaged')
        plt.plot(mu, df_axial_BEM / fac, '-k', label=r'$dT_{BEM}$')
        plt.plot(mu_LL[0], tan[0][0:Ncp] / fac, '--r', label=r'$dQ_{LLT}$ - non-averaged')
        plt.plot(mu_LL[1], tan[1][0:Ncp] / fac, '--g', label=r'$dQ_{LLT}$ - averaged')
        plt.plot(mu, df_tan_BEM / fac, '--k', label=r'$dQ_{BEM}$')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        plt.show()

        print('CT BEM = ' + str(CT))
        print('CP BEM = ' + str(CP))

        delta_r_nav = (mu_LL[0][1:] - mu_LL[0][:-1]) * R
        delta_r_av = (mu_LL[1][1:] - mu_LL[1][:-1]) * R
        CT_nav = np.sum(delta_r_nav * axial[0][1:Ncp] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP_nav = np.sum(delta_r_nav * tan[0][1:Ncp] * mu_LL[0][1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))
        CT_av = np.sum(delta_r_av * axial[1][1:Ncp] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP_av = np.sum(delta_r_av * tan[1][1:Ncp] * mu_LL[1][1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))

        print('CT LLT non_averaged = ' + str(CT_nav))
        print('CT LLT averaged = ' + str(CT_av))
        print('CP LLT non_averaged = ' + str(CP_nav))
        print('CP LLT averaged = ' + str(CP_av))






if choice == '3':
        spacing = [0,1] # induction factors average - 0: non-averaged, 1: averaged

        vortex = []
        control = []
        axial = []
        tan = []
        mu_LL = []

        for i in range(len(spacing)):
                
                # Compute Lifting Line Solution
                [vortexSystem,controlPoints,df_axial,df_tan] = fcn.lifting_line(averageFactor,spacing[i],mu,mu_root,mu_tip,chord_distribution,R,Ncp,alpha_BEM,U_inf,a_axial,Omega,a_tan,psi,polar_alpha,polar_CL,polar_CD,N_B,convFac,r_vortex,pitch,Loutlet,dx)
                vortex.append(vortexSystem)
                control.append(controlPoints)
                axial.append(df_axial)
                tan.append(df_tan)

                # Validation plots
                print('Solution ' + str(i+1) + ' converged')

                mu_LLT = control[i]['r'][0:Ncp]/R
                mu_LL.append(mu_LLT)

        fig1 = plt.figure(figsize=(8, 4))
        plt.title(r'Angle of Attack and Inflow Angle for $\lambda=$'+str(TSR))
        plt.plot(mu_LL[0], np.rad2deg(control[0]['alpha'][0:Ncp]), '-r',label=r'$\alpha_{LLT} - regular$')
        plt.plot(mu_LL[1], np.rad2deg(control[1]['alpha'][0:Ncp]), '-g',label=r'$\alpha_{LLT} - cosine$')
        plt.plot(mu, alpha_BEM, '-k', label=r'$\alpha_{BEM}$')
        plt.plot(mu_LL[0], np.rad2deg(control[0]['phi'][0:Ncp]), '--r', label=r'$\phi_{LLT} - regular$')
        plt.plot(mu_LL[1], np.rad2deg(control[1]['phi'][0:Ncp]), '--g', label=r'$\phi_{LLT} - cosine$')
        plt.plot(mu, np.rad2deg(phi_BEM), '--k', label=r'$\phi_{BEM}$')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        fig2 = plt.figure(figsize=(8, 4))
        plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega N_B}$ for $\lambda=$'+str(TSR))
        fac = np.pi * U_inf**2 / (Omega*N_B)
        plt.plot(mu_LL[0],control[0]['gamma'][0:Ncp]/fac, '-r', label='LLT - regular')
        plt.plot(mu_LL[1],control[1]['gamma'][0:Ncp]/fac, '-g', label='LLT - cosine')
        plt.plot(mu, gamma_BEM/fac, '-k', label='BEM Solution')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        fig3 = plt.figure(figsize=(8, 4))
        plt.title(r'Thrust and Azimuthal Loading, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$ for $\lambda=$'+str(TSR))
        fac = 0.5 * rho * U_inf**2 * R
        plt.plot(mu_LL[0], axial[0][0:Ncp] / fac,'-r', label=r'$dT_{LLT}$ - regular')
        plt.plot(mu_LL[1], axial[1][0:Ncp] / fac,'-g', label=r'$dT_{LLT}$ - cosine')
        plt.plot(mu, df_axial_BEM / fac, '-k', label=r'$dT_{BEM}$')
        plt.plot(mu_LL[0], tan[0][0:Ncp] / fac, '--r', label=r'$dQ_{LLT}$ - regular')
        plt.plot(mu_LL[1], tan[1][0:Ncp] / fac, '--g', label=r'$dQ_{LLT}$ - cosine')
        plt.plot(mu, df_tan_BEM / fac, '--k', label=r'$dQ_{BEM}$')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        plt.show()

        print('CT BEM = ' + str(CT))
        print('CP BEM = ' + str(CP))

        delta_r_reg = (mu_LL[0][1:] - mu_LL[0][:-1]) * R
        delta_r_cos = (mu_LL[1][1:] - mu_LL[1][:-1]) * R
        CT_reg = np.sum(delta_r_reg * axial[0][1:Ncp] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP_reg = np.sum(delta_r_reg * tan[0][1:Ncp] * mu_LL[0][1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))
        CT_cos = np.sum(delta_r_cos * axial[1][1:Ncp] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP_cos = np.sum(delta_r_cos * tan[1][1:Ncp] * mu_LL[1][1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))

        print('CT LLT regular = ' + str(CT_reg))
        print('CT LLT cosine = ' + str(CT_cos))
        print('CP LLT regular = ' + str(CP_reg))
        print('CP LLT cosine = ' + str(CP_cos))






if choice == '5':
        Loutlet = [0.5,2,5] # induction factors average - 0: non-averaged, 1: averaged

        vortex = []
        control = []
        axial = []
        tan = []
        mu_LL = []

        for i in range(len(Loutlet)):
                
                # Compute Lifting Line Solution
                [vortexSystem,controlPoints,df_axial,df_tan] = fcn.lifting_line(averageFactor,spacing,mu,mu_root,mu_tip,chord_distribution,R,Ncp,alpha_BEM,U_inf,a_axial,Omega,a_tan,psi,polar_alpha,polar_CL,polar_CD,N_B,convFac,r_vortex,pitch,Loutlet[i],dx)
                vortex.append(vortexSystem)
                control.append(controlPoints)
                axial.append(df_axial)
                tan.append(df_tan)

                # Validation plots
                print('Solution ' + str(i+1) + ' converged')

                mu_LLT = control[i]['r'][0:Ncp]/R
                mu_LL.append(mu_LLT)

        fig1 = plt.figure(figsize=(8, 4))
        plt.title(r'Angle of Attack and Inflow Angle for $\lambda=$'+str(TSR))
        plt.plot(mu_LL[0], np.rad2deg(control[0]['alpha'][0:Ncp]), '-r',label=r'$\alpha_{LLT}: L_{wake} =$ %.1fD' %Loutlet[0])
        plt.plot(mu_LL[1], np.rad2deg(control[1]['alpha'][0:Ncp]), '-g',label=r'$\alpha_{LLT}: L_{wake} =$ %.1fD' %Loutlet[1])
        plt.plot(mu_LL[2], np.rad2deg(control[2]['alpha'][0:Ncp]), '-b',label=r'$\alpha_{LLT}: L_{wake} =$ %.1fD' %Loutlet[2])
        plt.plot(mu, alpha_BEM, '-k', label=r'$\alpha_{BEM}$')
        plt.plot(mu_LL[0], np.rad2deg(control[0]['phi'][0:Ncp]), '--r', label=r'$\phi_{LLT}: L_{wake} =$ %.1fD' %Loutlet[0])
        plt.plot(mu_LL[1], np.rad2deg(control[1]['phi'][0:Ncp]), '--g', label=r'$\phi_{LLT}: L_{wake} =$ %.1fD' %Loutlet[1])
        plt.plot(mu_LL[2], np.rad2deg(control[2]['phi'][0:Ncp]), '--b', label=r'$\phi_{LLT}: L_{wake} =$ %.1fD' %Loutlet[2])
        plt.plot(mu, np.rad2deg(phi_BEM), '--k', label=r'$\phi_{BEM}$')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        fig2 = plt.figure(figsize=(8, 4))
        plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega N_B}$ for $\lambda=$'+str(TSR))
        fac = np.pi * U_inf**2 / (Omega*N_B)
        plt.plot(mu_LL[0],control[0]['gamma'][0:Ncp]/fac, '-r', label=r'$L_{wake} =$ %.1fD' %Loutlet[0])
        plt.plot(mu_LL[1],control[1]['gamma'][0:Ncp]/fac, '-g', label=r'$L_{wake} =$ %.1fD' %Loutlet[1])
        plt.plot(mu_LL[2],control[2]['gamma'][0:Ncp]/fac, '-b', label=r'$L_{wake} =$ %.1fD' %Loutlet[2])
        plt.plot(mu, gamma_BEM/fac, '-k', label='BEM Solution')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        fig3 = plt.figure(figsize=(8, 4))
        plt.title(r'Thrust and Azimuthal Loading, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$ for $\lambda=$'+str(TSR))
        fac = 0.5 * rho * U_inf**2 * R
        plt.plot(mu_LL[0], axial[0][0:Ncp] / fac,'-r', label=r'$dT_{LLT}: L_{wake} =$ %.1fD' %Loutlet[0])
        plt.plot(mu_LL[1], axial[1][0:Ncp] / fac,'-g', label=r'$dT_{LLT}: L_{wake} =$ %.1fD' %Loutlet[1])
        plt.plot(mu_LL[2], axial[2][0:Ncp] / fac,'-b', label=r'$dT_{LLT}: L_{wake} =$ %.1fD' %Loutlet[2])
        plt.plot(mu, df_axial_BEM / fac, '-k', label=r'$dT_{BEM}$')
        plt.plot(mu_LL[0], tan[0][0:Ncp] / fac, '--r', label=r'$dQ_{LLT}: L_{wake} =$ %.1fD' %Loutlet[0])
        plt.plot(mu_LL[1], tan[1][0:Ncp] / fac, '--g', label=r'$dQ_{LLT}: L_{wake} =$ %.1fD' %Loutlet[1])
        plt.plot(mu_LL[2], tan[2][0:Ncp] / fac, '--b', label=r'$dQ_{LLT}: L_{wake} =$ %.1fD' %Loutlet[2])
        plt.plot(mu, df_tan_BEM / fac, '--k', label=r'$dQ_{BEM}$')
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.grid()

        plt.show()

        print('CT BEM = ' + str(CT))
        print('CP BEM = ' + str(CP))

        delta_r_1 = (mu_LL[0][1:] - mu_LL[0][:-1]) * R
        delta_r_2 = (mu_LL[1][1:] - mu_LL[1][:-1]) * R
        delta_r_3 = (mu_LL[2][1:] - mu_LL[2][:-1]) * R
        CT_1 = np.sum(delta_r_1 * axial[0][1:Ncp] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP_2 = np.sum(delta_r_1 * tan[0][1:Ncp] * mu_LL[0][1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))
        CT_3 = np.sum(delta_r_2 * axial[1][1:Ncp] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP_1 = np.sum(delta_r_2 * tan[1][1:Ncp] * mu_LL[1][1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))
        CT_2 = np.sum(delta_r_3 * axial[2][1:Ncp] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP_3 = np.sum(delta_r_3 * tan[2][1:Ncp] * mu_LL[2][1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))

        print('CT LLT 0.5D = ' + str(CT_1))
        print('CT LLT 2D = ' + str(CT_2))
        print('CT LLT 5D = ' + str(CT_3))
        print('CP LLT 0.5D = ' + str(CP_1))
        print('CP LLT 2D = ' + str(CP_2))
        print('CP LLT 5D = ' + str(CP_3))
       
