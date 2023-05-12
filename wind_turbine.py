# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import functions as fcn
#from IPython.display import IFrame

airfoil = 'polarDU95W180.txt'
data1 = pd.read_csv(airfoil, header=0,
                    names=["alpha", "cl", "cd", "cm"],  sep='\s+')
polar_alpha = data1['alpha'][:]
polar_CL = data1['cl'][:]
polar_CD = data1['cd'][:]

# define the blade geometry
delta_mu = 0.01
mu = np.arange(0.2, 1 + delta_mu / 2, delta_mu)

# blade shape
R = 50
chord_distribution = 3 * (1 - mu) + 1 # meters
twist_distribution = 14 * (1 - mu)  # degrees

# define flow conditions
U_inf = 10  # unperturbed wind speed in m/s
N_B = 3 # number of blades

mu_tip = 1
mu_root = 0.2

print('1 - Tip Speed Ratio Variation')
print('2 - Yaw Angle Variation')
print('3 - Twist Angle Optimization')
choice = input('Select an option: ')

if choice == '1':
        TSR_distribution = [6, 8, 10] # tip speed ratio

        CT_distribution = []
        CP_distribution = []

        for i in range(len(TSR_distribution)):
                TSR = TSR_distribution[i]
                Omega = U_inf * TSR / R

                # solve BEM model
                results = np.zeros([len(mu)-1, 7])

                for i in range(len(mu)-1):
                        chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                        twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                        results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                        R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
                        
                T = np.sum(results[:, 5])
                Q = np.sum(results[:, 6])
                CT = T / (0.5 * U_inf ** 2 * np.pi * R ** 2)
                CP = Q * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2)

                CT_distribution.append(CT)
                CP_distribution.append(CP)

                print('CT = ', CT)
                print('CP = ', CP)

                fig1 = plt.figure(figsize=(12, 6))
                plt.title(r'Axial and Azimuthal Inductions for $TSR=$' + str(TSR) + ' degrees')
                plt.plot(results[:, 0], results[:, 1], 'k-', label=r'$a$')
                plt.plot(results[:, 0], results[:, 2], 'k--', label=r'$a^,$')
                plt.grid()
                plt.xlabel(r'$r/R$')
                plt.legend()
                plt.savefig("figures/TSR/inductionTSR" + str(TSR))

                fig1 = plt.figure(figsize=(12, 6))
                plt.title(r'Angle of Attack and Inflow Angle for $TSR=$' + str(TSR) + ' degrees')
                plt.plot(results[:, 0], results[:, 3], 'k-', label=r'$\alpha$')
                plt.plot(results[:, 0], results[:, 4], 'k--', label=r'$\phi$')
                plt.grid()
                plt.xlabel(r'$r/R$')
                plt.legend()
                plt.savefig("figures/TSR/anglesTSR" + str(TSR))

                fig1 = plt.figure(figsize=(12, 6))
                plt.title(r'Thrust and Azimuthal Loading for $TSR=$' + str(TSR) + ' degrees')
                plt.plot(results[:, 0], results[:, 5], 'k-', label=r'dT')
                plt.plot(results[:, 0], results[:, 6], 'k--', label=r'dQ')
                plt.grid()
                plt.xlabel(r'$r/R$')
                plt.legend()
                plt.savefig("figures/TSR/loadingTSR" + str(TSR))

        fig1 = plt.figure(figsize=(12, 6))
        plt.title('Thrust and Power Coefficients')
        plt.plot(TSR_distribution, CT_distribution, 'k-', label=r'$C_T$')
        plt.plot(TSR_distribution, CP_distribution, 'k--', label=r'$C_P$')
        plt.grid()
        plt.xlabel('TSR')
        plt.legend()
        plt.savefig("figures/TSR/coeffs")

        # plt.show()

elif choice == '2':
        TSR = 8
        Omega = U_inf * TSR / R

        yaw_distribution = [0, 15, 30]

        CT_distribution = []
        CP_distribution = []
        CP_yaw_distribution = []

        for i in range(len(yaw_distribution)):
                yaw = yaw_distribution[i]

                # solve BEM model
                results = np.zeros([len(mu)-1, 8])

                for i in range(len(mu)-1):
                        chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                        twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                        results[i, :] = fcn.solveStreamtube_yaw(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                                R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD, yaw)
                        
                T = np.sum(results[:, 5])
                Q = np.sum(results[:, 6])
                P_yaw = np.sum(results[:, 7])
                CT = T / (0.5 * U_inf ** 2 * np.pi * R ** 2)
                CP = Q * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2)
                CP_yaw = P_yaw / (0.5 * U_inf ** 3 * np.pi * R ** 2)

                CT_distribution.append(CT)
                CP_distribution.append(CP)
                CP_yaw_distribution.append(CP_yaw)

                print('CT is', CT)
                print('CP is ', CP)
                print('CP_yaw is ', CP_yaw)

                fig1 = plt.figure(figsize=(12, 6))
                plt.title(r'Axial and Azimuthal Inductions for $\theta=$' + str(yaw) + ' degrees')
                plt.plot(results[:, 0], results[:, 1], 'k-', label=r'$a$')
                plt.plot(results[:, 0], results[:, 2], 'k--', label=r'$a^,$')
                plt.grid()
                plt.xlabel(r'$r/R$')
                plt.legend()
                plt.savefig("figures/Yaw/inductionYaw" + str(yaw))

                fig1 = plt.figure(figsize=(12, 6))
                plt.title(r'Angle of Attack and Inflow Angle for $\theta=$' + str(yaw) + ' degrees')
                plt.plot(results[:, 0], results[:, 3], 'k-', label=r'$\alpha$')
                plt.plot(results[:, 0], results[:, 4], 'k--', label=r'$\phi$')
                plt.grid()
                plt.xlabel(r'$r/R$')
                plt.legend()
                plt.savefig("figures/Yaw/anglesYaw" + str(yaw))

                fig1 = plt.figure(figsize=(12, 6))
                plt.title(r'Thrust and Azimuthal Loading for $\theta=$' + str(yaw) + ' degrees')
                plt.plot(results[:, 0], results[:, 5], 'k-', label=r'dT')
                plt.plot(results[:, 0], results[:, 6], 'k--', label=r'dQ')
                plt.grid()
                plt.xlabel(r'$r/R$')
                plt.legend()
                plt.savefig("figures/Yaw/loadingYaw" + str(yaw))

        fig1 = plt.figure(figsize=(12, 6))
        plt.title('Thrust and Power Coefficients')
        plt.plot(yaw_distribution, CT_distribution, 'k-', label=r'$C_T$')
        plt.plot(yaw_distribution, CP_distribution, 'k--', label=r'$C_P$')
        plt.plot(yaw_distribution, CP_yaw_distribution, 'k-.', label=r'$C_P$')
        plt.grid()
        plt.xlabel(r'$\theta$') 
        plt.legend()
        plt.savefig("figures/Yaw/coeffs" + str(yaw))
        # plt.show()

elif choice == '3':
        TSR = 8 
        told = 0.02 # within which range of CTs to look [0.75 - tol : 0.75:tol]
        # maybe write a short prompt to ask for bounds?
        a = np.linspace(-2, -1, 10)     # start with bigger range, restrict based on the output
        b = np.linspace(-1, 0, 10)
        c = np.linspace(-3, -2, 10)
        # phi_new = phi_original + a mu^2 + b mu + c
        N = len(a)*len(b)*len(c)
        CT_distribution = np.empty((N,1))
        CP_distribution = np.empty((N,1))
        index = np.empty((N,3))

        for i1 in range(len(a)):
                for i2 in range(len(b)):
                        for i3 in range(len(c)):
                                # assume a twist distribution
                                twist_distribution = 14 * (1 - mu) + a[i1]*mu**2 + b[i2]*mu + c[i3]
                                Omega = U_inf * TSR / R
                                # solve BEM model
                                results = np.zeros([len(mu)-1, 7])

                                for i in range(len(mu)-1):
                                        chord = np.interp((mu[i] + mu[i+1]) / 2, mu, chord_distribution)
                                        twist = np.interp((mu[i] + mu[i+1]) / 2, mu, twist_distribution)

                                        results[i, :] = fcn.solveStreamtube(U_inf, mu[i], mu[i+1], delta_mu, mu_root, mu_tip, Omega,
                                                                        R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD)
                                        
                                T = np.sum(results[:, 5])
                                Q = np.sum(results[:, 6])
                                CT = T / (0.5 * U_inf ** 2 * np.pi * R ** 2)
                                CP = Q * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2)
                                
                                counter = len(a)*len(b)*i1 + len(b)*i2 + i3
                                CT_distribution[counter] = CT
                                CP_distribution[counter] = CP
                                index[counter] = np.array([i1, i2 ,i3])

        # Iteration loop finished,  pick the best result
        ind = np.where((CT_distribution < 0.75+told) & (CT_distribution > 0.75-told))[0]
        maxId = np.argmax(CP_distribution[ind])
        CP_max = CP_distribution[ind[maxId]]
        CT_max = CT_distribution[ind[maxId]]
        amaxId, bmaxId, cmaxId = index[ind[maxId]]
        amax = a[int(amaxId)]; bmax = b[int(bmaxId)]; cmax = c[int(cmaxId)]; 
        twist_initial = 14 * (1 - mu) 
        twist_opt = 14 * (1 - mu) + amax*mu**2 + bmax*mu + cmax

        print("Optimal a = ", str(amax), " b = ", str(bmax), "c = ", str(cmax))
        print("Thrust coefficient = ", str(CT_max))
        print("Power coefficient = ", str(CP_max))

        plt.figure(figsize=(12, 6))
        plt.title('Optimized Twist Distribution')
        plt.plot(mu, twist_initial, 'k--', label=r'$\phi_{old}$')
        plt.plot(mu, twist_opt, 'k-', label=r'$\phi_{opt}$')
        plt.grid()
        plt.xlabel(r'$r/R$') 
        plt.legend()
        plt.savefig("figures/Twist/opt_distribution.png")
        plt.show()
        print("end")