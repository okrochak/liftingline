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

# define the blade geometry
delta_mu = 0.01
mu = np.arange(0.2, 1 + delta_mu / 2, delta_mu)
delta_Phi = np.pi/10
Phi = np.arange(0.2, 2*np.pi + delta_Phi / 2, delta_Phi)

# blade shape
R = 50
pitch = 2 # degrees
chord_distribution = 3 * (1 - mu) + 1 # meters
twist_distribution = 14 * (1 - mu) - pitch  # degrees

# define flow conditions
U_inf = 10  # unperturbed wind speed in m/s
N_B = 3 # number of blades

mu_tip = 1
mu_root = 0.2

print('1 - Tip Speed Ratio Variation')
print('2 - Yaw Angle Variation')
print('3 - Twist Angle Optimization')
print('4 - Influence of the Tip Correction')
print('5 - Influence of Number of Annuli with Constant Spacing')
print('6 - Stagnation Pressure Distribution')
print('7 - Airfoil Operational Point')

choice = input('Select an option: ')

'''Run the lifting Line Model '''


TSR_distribution = [8] # tip speed ratio
# TSR_distribution = np.linspace(4,10,24)
CT_distribution = []
CP_distribution = []

for j in range(len(TSR_distribution)):
        TSR = TSR_distribution[j]
        Omega = U_inf * TSR / R

        # solve BEM model
        results = np.zeros([len(mu)-1, 10])
        results_unc = np.zeros([len(mu)-1, 10])

        for i in range(len(mu)-1):
                chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                        R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
                
                results_unc[i,:] = fcn.solveStreamtube_unc(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                        R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
        
        delta_r = (mu[1:] - mu[:-1]) * R
        CT = np.sum(delta_r * results[:, 5] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
        CP = np.sum(delta_r * results[:, 6] * results[:, 0] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))

        CT_distribution.append(CT)
        CP_distribution.append(CP)

        print('CT = ', CT)
        print('CP = ', CP)

        fig1 = plt.figure(figsize=(12, 6))
        plt.title(r'Axial and Azimuthal Inductions for $\lambda=$' + str(TSR))
        plt.plot(results[:, 0], results[:, 1], 'k-', label=r'$a$ corr.')
        plt.plot(results[:, 0], results[:, 2], 'k--', label=r'$a^,$ corr.')
        plt.plot(results_unc[:, 0], results_unc[:, 1], 'r-', label=r'$a$ unc.')
        plt.plot(results_unc[:, 0], results_unc[:, 2], 'r--', label=r'$a^,$ unc.')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.savefig("figures/TSR/inductionTSR" + str(TSR),bbox_inches='tight')

#         fig1 = plt.figure(figsize=(8, 4))
#         plt.title(r'Angle of Attack and Inflow Angle for $\lambda=$' + str(TSR))
#         plt.plot(results[:, 0], results[:, 3], 'k-', label=r'$\alpha$')
#         plt.plot(results[:, 0], np.rad2deg(results[:, 4]), 'k--', label=r'$\phi$')
#         plt.grid()
#         plt.xlabel(r'$r/R$')
#         plt.legend()
#         plt.savefig("figures/TSR/anglesTSR" + str(TSR),bbox_inches='tight')

#         fig1 = plt.figure(figsize=(8, 4))
#         plt.title(r'Thrust and Azimuthal Loading, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$ for $\lambda=$' + str(TSR))
#         plt.plot(results[:, 0], results[:, 5] / (0.5 * U_inf ** 2 * R), 'k-', label=r'dT')
#         plt.plot(results[:, 0], results[:, 6] / (0.5 * U_inf ** 2 * R), 'k--', label=r'dQ')
#         plt.grid()
#         plt.xlabel(r'$r/R$')
#         plt.legend()
#         plt.savefig("figures/TSR/loadingTSR" + str(TSR),bbox_inches='tight')
        print(j)
        if j == 0:
                gamma1 = results[:,7]
        elif j == 1:
                gamma2 = results[:,7]
        else:
                gamma3 = results[:,7]




# fig1 = plt.figure(figsize=(8, 4))
# plt.title('Thrust and Power Coefficients')
# plt.plot(TSR_distribution, CT_distribution, 'k-', label=r'$C_T$')
# plt.plot(TSR_distribution, CP_distribution, 'k--', label=r'$C_P$')
# plt.grid()
# plt.xlabel('TSR')
# plt.legend()
# plt.savefig("figures/TSR/coeffs",bbox_inches='tight')
# # plt.show()


fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega N_B}$ for varying $\lambda$')
plt.plot(results[:, 0], gamma1/(np.pi * U_inf ** 2 / (N_B * Omega)), 'k-', label=r'$\Gamma $ for $\lambda = 6$')
plt.plot(results[:, 0], gamma2/(np.pi * U_inf ** 2 / (N_B * Omega)), 'r-', label=r'$\Gamma $ for $\lambda = 8$')
plt.plot(results[:, 0], gamma3/(np.pi * U_inf ** 2 / (N_B * Omega)), 'b-', label=r'$\Gamma $ for $\lambda = 10$')
plt.grid()
plt.xlabel(r'$r/R$')
plt.legend()
plt.savefig("figures/TSR/circulation_TSR",bbox_inches='tight')

fig1 = plt.figure(figsize=(12, 6))
plt.title(r'$C_T$ and $C_P$ vs. $\lambda$' )
plt.plot(TSR_distribution, CT_distribution, 'k-', label=r'$C_T$')
plt.plot(TSR_distribution, CP_distribution, 'r-', label=r'$C_P$')
plt.grid()
plt.xlabel(r'$\lambda$')
plt.legend()
plt.savefig("figures/TSR/loadingLambda",bbox_inches='tight')







if choice == '1':
        TSR_distribution = [8] # tip speed ratio
        # TSR_distribution = np.linspace(4,10,24)
        CT_distribution = []
        CP_distribution = []

        for j in range(len(TSR_distribution)):
                TSR = TSR_distribution[j]
                Omega = U_inf * TSR / R

                # solve BEM model
                results = np.zeros([len(mu)-1, 10])
                results_unc = np.zeros([len(mu)-1, 10])

                for i in range(len(mu)-1):
                    chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                    twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                    results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                            R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
                    
                    results_unc[i,:] = fcn.solveStreamtube_unc(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                            R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
                
                delta_r = (mu[1:] - mu[:-1]) * R
                CT = np.sum(delta_r * results[:, 5] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
                CP = np.sum(delta_r * results[:, 6] * results[:, 0] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))

                CT_distribution.append(CT)
                CP_distribution.append(CP)

                print('CT = ', CT)
                print('CP = ', CP)

                fig1 = plt.figure(figsize=(12, 6))
                plt.title(r'Axial and Azimuthal Inductions for $\lambda=$' + str(TSR))
                plt.plot(results[:, 0], results[:, 1], 'k-', label=r'$a$ corr.')
                plt.plot(results[:, 0], results[:, 2], 'k--', label=r'$a^,$ corr.')
                plt.plot(results_unc[:, 0], results_unc[:, 1], 'r-', label=r'$a$ unc.')
                plt.plot(results_unc[:, 0], results_unc[:, 2], 'r--', label=r'$a^,$ unc.')
                plt.grid()
                plt.xlabel(r'$r/R$')
                plt.legend()
                plt.savefig("figures/TSR/inductionTSR" + str(TSR),bbox_inches='tight')

        #         fig1 = plt.figure(figsize=(8, 4))
        #         plt.title(r'Angle of Attack and Inflow Angle for $\lambda=$' + str(TSR))
        #         plt.plot(results[:, 0], results[:, 3], 'k-', label=r'$\alpha$')
        #         plt.plot(results[:, 0], np.rad2deg(results[:, 4]), 'k--', label=r'$\phi$')
        #         plt.grid()
        #         plt.xlabel(r'$r/R$')
        #         plt.legend()
        #         plt.savefig("figures/TSR/anglesTSR" + str(TSR),bbox_inches='tight')

        #         fig1 = plt.figure(figsize=(8, 4))
        #         plt.title(r'Thrust and Azimuthal Loading, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$ for $\lambda=$' + str(TSR))
        #         plt.plot(results[:, 0], results[:, 5] / (0.5 * U_inf ** 2 * R), 'k-', label=r'dT')
        #         plt.plot(results[:, 0], results[:, 6] / (0.5 * U_inf ** 2 * R), 'k--', label=r'dQ')
        #         plt.grid()
        #         plt.xlabel(r'$r/R$')
        #         plt.legend()
        #         plt.savefig("figures/TSR/loadingTSR" + str(TSR),bbox_inches='tight')
                print(j)
                if j == 0:
                        gamma1 = results[:,7]
                elif j == 1:
                        gamma2 = results[:,7]
                else:
                        gamma3 = results[:,7]




        # fig1 = plt.figure(figsize=(8, 4))
        # plt.title('Thrust and Power Coefficients')
        # plt.plot(TSR_distribution, CT_distribution, 'k-', label=r'$C_T$')
        # plt.plot(TSR_distribution, CP_distribution, 'k--', label=r'$C_P$')
        # plt.grid()
        # plt.xlabel('TSR')
        # plt.legend()
        # plt.savefig("figures/TSR/coeffs",bbox_inches='tight')
        # # plt.show()


        fig1 = plt.figure(figsize=(12, 6))
        plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega N_B}$ for varying $\lambda$')
        plt.plot(results[:, 0], gamma1/(np.pi * U_inf ** 2 / (N_B * Omega)), 'k-', label=r'$\Gamma $ for $\lambda = 6$')
        plt.plot(results[:, 0], gamma2/(np.pi * U_inf ** 2 / (N_B * Omega)), 'r-', label=r'$\Gamma $ for $\lambda = 8$')
        plt.plot(results[:, 0], gamma3/(np.pi * U_inf ** 2 / (N_B * Omega)), 'b-', label=r'$\Gamma $ for $\lambda = 10$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.savefig("figures/TSR/circulation_TSR",bbox_inches='tight')

        fig1 = plt.figure(figsize=(12, 6))
        plt.title(r'$C_T$ and $C_P$ vs. $\lambda$' )
        plt.plot(TSR_distribution, CT_distribution, 'k-', label=r'$C_T$')
        plt.plot(TSR_distribution, CP_distribution, 'r-', label=r'$C_P$')
        plt.grid()
        plt.xlabel(r'$\lambda$')
        plt.legend()
        plt.savefig("figures/TSR/loadingLambda",bbox_inches='tight')


elif choice == '2':
        TSR = 8
        Omega = U_inf * TSR / R

        yaw_distribution = [0, 15, 30]

        CT_distribution = []
        CP_distribution = []

        #for i in range(len(yaw_distribution)):
                #yaw = yaw_distribution[i]*np.pi/180
        yaw = 15*np.pi/180
                # solve BEM model
        results = np.zeros([len(mu)-1, len(Phi)-1, 10])

        for i in range(len(mu)-1):
                chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                for j in range(len(Phi)-1):
                        results[i, j, :] = fcn.solveStreamtube_yaw(U_inf, mu[i], mu[i+1], mu_root, mu_tip, Omega, R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD, yaw, Phi[i], Phi[i+1], TSR)

                a_bar = np.sum(results[i,:,1]) 
                aline_bar = np.sum(results[i,:,2]) 
                mu_bar = np.mean(results[i,:,0])       
                #delta_r = (mu[1:] - mu[:-1]) * R
                #CT = np.sum(delta_r * results[:, 5] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
                #CP = np.sum(results[:, 7] / (0.5 * U_inf ** 3 * np.pi * R ** 2))

                #CT_distribution.append(CT)
                #CP_distribution.append(CP)

                #print('CT is', CT)
                #print('CP is ', CP)

        fig1 = plt.figure(figsize=(12, 6))
        plt.title(r'Axial and Azimuthal Inductions for $\theta=$' + str(yaw) + ' degrees')
        plt.plot(mu_bar, a_bar, 'k-', label=r'$a$')
        plt.plot(mu_bar, aline_bar, 'k--', label=r'$a^,$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
                #plt.savefig("figures/Yaw/inductionYaw" + str(yaw))

                #fig1 = plt.figure(figsize=(12, 6))
                #plt.title(r'Angle of Attack and Inflow Angle for $\theta=$' + str(yaw) + ' degrees')
                #plt.plot(results[:, 0], results[:, 3], 'k-', label=r'$\alpha$')
                #plt.plot(results[:, 0], results[:, 4], 'k--', label=r'$\phi$')
                #plt.grid()
                #plt.xlabel(r'$r/R$')
                #plt.legend()
                #plt.savefig("figures/Yaw/anglesYaw" + str(yaw))

                #fig1 = plt.figure(figsize=(12, 6))
                #plt.title(r'Thrust and Azimuthal Loading, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$, for $\theta=$' + str(yaw) + ' degrees')
                #plt.plot(results[:, 0], results[:, 5] / (0.5 * U_inf ** 2 * R), 'k-', label=r'dT')
                #plt.plot(results[:, 0], results[:, 6] / (0.5 * U_inf ** 2 * R), 'k--', label=r'dQ')
                #plt.grid()
                #plt.xlabel(r'$r/R$')
                #plt.legend()
                #plt.savefig("figures/Yaw/loadingYaw" + str(yaw))

        #fig1 = plt.figure(figsize=(12, 6))
        #plt.title('Thrust and Power Coefficients')
        #plt.plot(yaw_distribution, CT_distribution, 'k-', label=r'$C_T$')
        #plt.plot(yaw_distribution, CP_distribution, 'k--', label=r'$C_P$')
        #plt.grid()
        #plt.xlabel(r'$\theta$') 
        #plt.legend()
        # plt.savefig("figures/Yaw/coeffs" + str(yaw))
        #plt.show()

elif choice == '3':
        TSR = 8 
        tol = 0.005 # within which range of CTs to look [0.75-tol : 0.75+tol]
        a = np.linspace(9, 11, 5)     # start with bigger range, restrict based on the output
        b = np.linspace(0, 5, 5) #Optimal a =  12.0  b =  0.75 c =  -15.25 d =  3.0
        c = np.linspace(-18, -15, 5)
        d = np.linspace(2, 6, 5) 
        # phi_new = phi_original + a mu^2 + b mu + c
        N = len(a)*len(b)*len(c)*len(d)
        CT_distribution = np.empty((N,1))
        CP_distribution = np.empty((N,1))
        index = np.empty((N,4))

        for i1 in range(len(a)):
                for i2 in range(len(b)):
                        for i3 in range(len(c)):
                                for i4 in range(len(d)):
                                        # assume a twist distribution
                                        twist_distribution = 14*(1-mu) + a[i1]*mu**3 + b[i2]*mu**2  + c[i3]*mu + d[i4] - pitch
                                        Omega = U_inf * TSR / R
                                        # solve BEM model
                                        results = np.zeros([len(mu)-1, 10])

                                        for i in range(len(mu)-1):
                                                chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                                                twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                                                results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                                                R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
                                        delta_r = (mu[1:] - mu[:-1]) * R
                                        T = np.sum(delta_r * N_B * results[:, 5])
                                        Q = np.sum(delta_r * N_B * Omega * R * results[:, 6] * results[:, 0])
                                        CT = T / (0.5 * U_inf ** 2 * np.pi * R ** 2)
                                        CP = Q  / (0.5 * U_inf ** 3 * np.pi * R ** 2)
                                        
                                        counter = len(a)*len(b)*len(c)*i1 + len(b)*len(c)*i2 + len(c)*i3 + i4
                                        CT_distribution[counter] = CT
                                        CP_distribution[counter] = CP
                                        index[counter] = np.array([i1, i2 ,i3, i4])

        # Iteration loop finished,  pick the best result
        ind = np.where((CT_distribution < 0.75+tol) & (CT_distribution > 0.75-tol))[0]
        maxId = np.argmax(CP_distribution[ind])
        CP_max = CP_distribution[ind[maxId]]
        CT_max = CT_distribution[ind[maxId]]
        amaxId, bmaxId, cmaxId, dmaxid = index[ind[maxId]]
        amax = a[int(amaxId)]; bmax = b[int(bmaxId)]; cmax = c[int(cmaxId)]; dmax = d[int(dmaxid)]
        twist_initial = 14 * (1 - mu) - pitch
        twist_opt =  14 * (1 - mu) - pitch + amax*mu**3 + bmax*mu**2 + cmax*mu + dmax 

        print("Optimal a = ", str(amax), " b = ", str(bmax), "c = ", str(cmax), "d = ", str(dmax))
        print("Thrust coefficient = ", str(CT_max))
        print("Power coefficient = ", str(CP_max))


        twist_distribution = twist_opt
        Omega = U_inf * TSR / R
        # solve BEM model
        results = np.zeros([len(mu)-1, 10])
        results_opt = np.zeros([len(mu)-1, 10])
        results_def = np.zeros([len(mu)-1, 10])

        for i in range(len(mu)-1):
                chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                results_opt[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)


        twist_distribution = twist_initial
        Omega = U_inf * TSR / R
        # solve BEM model
        results = np.zeros([len(mu)-1, 10])

        for i in range(len(mu)-1):
                chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                results_def[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)

        plt.figure(figsize=(12, 6))
        plt.title('Optimized Twist Distribution')
        plt.plot(mu, twist_initial, 'k--', label=r'$\phi_{old}$')
        plt.plot(mu, twist_opt, 'k-', label=r'$\phi_{opt}$')
        plt.grid()
        plt.xlabel(r'$r/R$') 
        plt.legend()
        #plt.savefig("figures/Twist/opt_distribution.png")

        ## adapt to optimal case
        fig1 = plt.figure(figsize=(12, 6))
        plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega N_B}$')
        plt.plot(results_def[:, 0], results_def[:, 7]/(np.pi * U_inf ** 2 / (N_B * Omega)), 'k-', label=r'$\Gamma_{1}$')
        plt.plot(results_opt[:, 0], results_opt[:, 7]/(np.pi * U_inf ** 2 / (N_B * Omega)), 'r-', label=r'$\Gamma_{2}$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.show()                   
        
if choice == '4':
        # plot Prandtl tip, root and combined correction for a number of blades and induction 'a', over the non-dimensioned radius
        mu = np.arange(0.2, 1, .01)
        a = np.zeros(np.shape(mu)) + 0.3
        Prandtl, Prandtltip, Prandtlroot = fcn.PrandtlTipRootCorrection(mu, mu_root, mu_tip, 8, N_B, a)

        fig1 = plt.figure(figsize=(12, 6))
        plt.plot(mu, Prandtl, 'r-', label='Prandtl')
        plt.plot(mu, Prandtltip, 'g.', label='Prandtl tip')
        plt.plot(mu, Prandtlroot, 'b.', label='Prandtl root')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        #plt.show()
        plt.savefig("figures/prandtl.png",bbox_inches='tight')

if choice == '5':
        delta_mu_distribution = [0.1, 0.01, 0.001, 0.0001]

        TSR = 8
        Omega = U_inf * TSR / R

        N_distribution = []
        CT_distribution = []
        CP_distribution = []
        
        results_total = [] 

        for j in range(len(delta_mu_distribution)):
                delta_mu = delta_mu_distribution[j]
                mu = np.arange(0.2, 1 + delta_mu / 2, delta_mu)

                chord_distribution = 3 * (1 - mu) + 1 # meters
                twist_distribution = 14 * (1 - mu) - pitch  # degrees

                N = len(mu) - 1
                N_distribution.append(N)

                # solve BEM model
                results = np.zeros([len(mu)-1, 10])

                for i in range(len(mu)-1):
                    chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                    twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                    results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                            R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
                        
                delta_r = (mu[1:] - mu[:-1]) * R
                CT = np.sum(delta_r * results[:, 5] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
                CP = np.sum(delta_r * results[:, 6] * results[:, 0] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))

                CT_distribution.append(CT)
                CP_distribution.append(CP)

                print('CT = ', CT)
                print('CP = ', CP)

                results_total.append(results)

                # fig1 = plt.figure(figsize=(12, 6))
                # plt.title(r'Axial and Azimuthal Inductions for $TSR=$' + str(TSR))
                # plt.plot(results[:, 0], results[:, 1], 'k-', label=r'$a$')
                # plt.plot(results[:, 0], results[:, 2], 'k--', label=r'$a^,$')
                # plt.grid()
                # plt.xlabel(r'$r/R$')
                # plt.legend()
                # #plt.savefig("figures/TSR/inductionTSR" + str(TSR))

                # fig1 = plt.figure(figsize=(12, 6))
                # plt.title(r'Angle of Attack and Inflow Angle for $N=$' + str(N))
                # plt.plot(results[:, 0], results[:, 3], 'k-', label=r'$\alpha$')
                # plt.plot(results[:, 0], results[:, 4], 'k--', label=r'$\phi$')
                # plt.grid()
                # plt.xlabel(r'$r/R$')
                # plt.legend()
                # #plt.savefig("figures/TSR/anglesTSR" + str(TSR))

                # fig1 = plt.figure(figsize=(12, 6))
                # plt.title(r'Thrust and Azimuthal Loading, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$ for $N=$' + str(N))
                # plt.plot(results[:, 0], results[:, 5] / (0.5 * U_inf ** 2 * R), 'k-', label=r'dT')
                # plt.plot(results[:, 0], results[:, 6] / (0.5 * U_inf ** 2 * R), 'k--', label=r'dQ')
                # plt.grid()
                # plt.xlabel(r'$r/R$')
                # plt.legend()
                # #plt.savefig("figures/TSR/loadingTSR" + str(TSR))

                # fig1 = plt.figure(figsize=(12, 6))
                # plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega N_B}$ for $N=$' + str(N))
                # plt.plot(results[:, 0], results[:, 7]/(np.pi * U_inf ** 2 / (N_B * Omega)), 'k-', label=r'$\Gamma$')
                # plt.grid()
                # plt.xlabel(r'$r/R$')
                # plt.legend()
                # #plt.show()                                      




        fig1 = plt.figure(figsize=(12, 6))
        plt.title(r'Angle of attack for varying N')
        plt.plot(results_total[0][:, 0], results_total[0][:, 3], 'k:', label=r'$\alpha, N = 10$')
        plt.plot(results_total[1][:, 0], results_total[1][:, 3], 'r-', label=r'$\alpha, N = 100$')
        plt.plot(results_total[2][:, 0], results_total[2][:, 3], 'g--', label=r'$\alpha, N = 1000$')
        plt.plot(results_total[3][:, 0], results_total[3][:, 3], 'b-.', label=r'$\alpha, N = 10000$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.show()
        plt.savefig("figures/Discretization/alphaN")
 
        fig1 = plt.figure(figsize=(12, 6))
        plt.title(r'Axial Inductions for varying N')
        plt.plot(results_total[0][:, 0], results_total[0][:, 1], 'k:', label=r'$a, N = 10$')
        plt.plot(results_total[1][:, 0], results_total[1][:, 1], 'r-', label=r'$a, N = 100$')
        plt.plot(results_total[2][:, 0], results_total[2][:, 1], 'g--', label=r'$a, N = 1000$')
        plt.plot(results_total[3][:, 0], results_total[3][:, 1], 'b-.', label=r'$a, N = 10000$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.show()
        plt.savefig("figures/Discretization/a_inductionN")
                
        fig1 = plt.figure(figsize=(12, 6))
        plt.title(r'Azimuthal forces for varying N')
        plt.plot(results_total[0][:, 0], results_total[0][:, 6], 'k:', label=r'$F_{azim}, N = 10$')
        plt.plot(results_total[1][:, 0], results_total[1][:, 6], 'r-', label=r'$F_{azim}, N = 100$')
        plt.plot(results_total[2][:, 0], results_total[2][:, 6], 'g--', label=r'$F_{azim}, N = 1000$')
        plt.plot(results_total[3][:, 0], results_total[3][:, 6], 'b-.', label=r'$F_{azim}, N = 10000$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.show()
        plt.savefig("figures/Discretization/F_azimN")

  


if choice == '6':
        TSR = 8
        Omega = U_inf * TSR / R

        # solve BEM model
        results = np.zeros([len(mu)-1, 10])

        A_rotor = []

        for i in range(len(mu)-1):
                chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                            R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
            
                A = np.pi * ((mu[i+1] * R) ** 2 - (mu[i] * R) ** 2)
                A_rotor.append(A)


        A_rotor = np.array(A_rotor, dtype=np.float32)

        mu_rotor = results[:, 0]

        df_axial = results[:, 5]
        dF_axial = df_axial * R * delta_mu * N_B
        CdT = dF_axial / (0.5 * A_rotor * U_inf ** 2)
        a = 0.5 - 0.5 * np.sqrt(1 - CdT)

        U_rotor_axial = U_inf * (1 - a)

        p_inf = 0 #101325 # Pa
        p_stag_upwind = p_inf + 0.5 * U_inf ** 2
        A_upwind = U_rotor_axial * A_rotor / U_inf
        mu_upwind = np.zeros(len(A_upwind))
        mu_upwind[0] = 0.2 # assumption
        for i in range(len(mu_upwind) - 1):
                mu_upwind[i + 1] = np.sqrt(A_upwind[i] / (np.pi * R ** 2) + mu_upwind[i] ** 2) 

        p_stag_rotor_upwind = p_stag_upwind * np.ones(len(mu_rotor))

        delta_r = (mu[1:] - mu[:-1]) * R
        T = delta_r * df_axial * N_B
        delta_p_rotor = T / A_rotor
        p_stag_rotor_downwind = p_stag_rotor_upwind - delta_p_rotor

        U_downwind = np.sqrt(2 * (p_stag_rotor_downwind - p_inf * np.ones(len(mu_rotor))))
        A_downwind = U_rotor_axial * A_rotor / U_downwind

        mu_downwind = np.zeros(len(A_downwind))
        mu_downwind[0] = 0.2 # assumption
        for i in range(len(mu_downwind) - 1):
                mu_downwind[i + 1] = np.sqrt(A_downwind[i] / (np.pi * R ** 2) + mu_downwind[i] ** 2) 

        # fig1 = plt.figure(figsize=(12, 6))
        # plt.title('Stagnation Pressure at Upwind Infinity')
        # plt.plot(mu_upwind, p_stag_upwind * np.ones(len(mu_upwind)), 'k-', label=r'Stagnation Pressure')
        # #plt.plot(mu_upwind, p_inf * np.ones(len(mu_upwind)), 'k--', label=r'Atmospheric Pressure')
        # plt.grid()
        # plt.ylabel(r'$p_{\mathrm{stag}} - p_{\mathrm{atm}} \quad \mathrm{[Pa]}$')
        # plt.xlabel(r'$r/R$')
        # #plt.legend()

        # fig1 = plt.figure(figsize=(12, 6))
        # plt.title('Stagnation Pressure at Rotor (Upwind Side)')
        # plt.plot(mu_rotor, p_stag_rotor_upwind, 'k-', label=r'Stagnation Pressure')
        # #plt.plot(mu_rotor, p_inf * np.ones(len(mu_rotor)), 'k--', label=r'Atmospheric Pressure')
        # plt.grid()
        # plt.ylabel(r'$p_{\mathrm{stag}} - p_{\mathrm{atm}} \quad \mathrm{[Pa]}$')
        # plt.xlabel(r'$r/R$')
        # #plt.legend()

        # fig1 = plt.figure(figsize=(12, 6))
        # plt.title('Stagnation Pressure at Rotor (Downwind Side)')
        # plt.plot(mu_rotor, p_stag_rotor_downwind, 'k-', label=r'Stagnation Pressure')
        # #plt.plot(mu_rotor, p_inf * np.ones(len(mu_rotor)), 'k--', label=r'Atmospheric Pressure')
        # plt.grid()
        # plt.ylabel(r'$p_{\mathrm{stag}} - p_{\mathrm{atm}} \quad \mathrm{[Pa]}$')
        # plt.xlabel(r'$r/R$')
        # #plt.legend()

        # fig1 = plt.figure(figsize=(12, 6))
        # plt.title('Stagnation Pressure at Downwind Infinity')
        # plt.plot(mu_downwind, p_stag_rotor_downwind, 'k-', label=r'Stagnation Pressure')
        # #plt.plot(mu_downwind, p_inf * np.ones(len(mu_downwind)), 'k--', label=r'Atmospheric Pressure')
        # plt.grid()
        # plt.ylabel(r'$p_{\mathrm{stag}} - p_{\mathrm{atm}} \quad \mathrm{[Pa]}$')
        # plt.xlabel(r'$r/R$')
        # plt.legend()

        fig1 = plt.figure(figsize=(12, 6))
        plt.title('Stagnation Pressure')
        plt.plot(mu_upwind, p_stag_rotor_upwind, 'r-', label=r'far upstream')
        plt.plot(mu_rotor, p_stag_rotor_upwind, 'r--', label=r'upstream of rotor')
        plt.plot(mu_rotor, p_stag_rotor_downwind, 'k-', label=r'downstream of rotor')
        plt.plot(mu_downwind, p_stag_rotor_downwind, 'k--', label=r'far downstream')
        #plt.plot(mu_downwind, p_inf * np.ones(len(mu_downwind)), 'k--', label=r'Atmospheric Pressure')
        plt.grid()
        plt.legend(loc='center left')
        plt.ylabel(r'$p_{\mathrm{stag}} - p_{\mathrm{atm}} \quad \mathrm{[Pa]}$')
        plt.savefig('figures/stagnation_pressure.png')
        plt.xlabel(r'$r/R$')

        plt.show()  

if choice == '7':
        polar_CL_CD = polar_CL / polar_CD

        fig1 = plt.figure()
        plt.plot(polar_alpha, polar_CL, 'k-')
        plt.title('Lift Coefficient Polar')
        plt.grid()
        plt.ylabel(r'$C_L$')
        plt.xlabel(r'$\alpha$')

        fig1 = plt.figure()
        plt.plot(polar_alpha, polar_CD, 'k-')
        plt.title('Drag Coefficient Polar')
        plt.grid()
        plt.ylabel(r'$C_D$')
        plt.xlabel(r'$\alpha$')


        fig1 = plt.figure()
        plt.plot(polar_alpha, polar_CL_CD, 'k-')
        plt.title('Airfoil Efficiency Polar')
        plt.grid()
        plt.ylabel(r'$C_L / C_D$')
        plt.xlabel(r'$\alpha$')

        TSR = 8
        Omega = U_inf * TSR / R

        # solve BEM model
        results = np.zeros([len(mu)-1, 10])
        chord_interpolated = np.zeros([len(mu)-1, 1])

        for i in range(len(mu)-1):
                chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                chord_interpolated[i] = chord
                twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)

        
        alpha_opt = 8.73
        CL_CD_opt = 117.50

        plt.figure(figsize=(12, 6))
        plt.title('Lift Coefficient Distribution')
        plt.plot(results[:, 0], results[:, 8], 'k-')
        plt.grid()
        plt.xlabel(r'$r/R$') 
        plt.ylabel(r'$C_L$')

        plt.figure(figsize=(12, 6))
        plt.title('Chord Distribution')
        plt.plot(results[:, 0], chord_interpolated, 'k-')
        plt.grid()
        plt.xlabel(r'$r/R$') 
        plt.ylabel(r'$c$')

        plt.figure(figsize=(12, 6))
        plt.title('Operational Point Efficiency Distribution')
        plt.plot(results[:, 0], results[:, 8] / results[:, 9], 'k-')
        plt.grid()
        plt.xlabel(r'$r/R$') 
        plt.ylabel(r'$C_L/C_D$')

        plt.show()   

       
            



