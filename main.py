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
TSR = 8

doPlot = True

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
'''Run the BEM model for the inputs '''
# TSR_distribution = np.linspace(4,10,24)
CT_distribution = []
CP_distribution = []
results = np.zeros([len(mu)-1, 10])

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
delta_r = (mu[1:] - mu[:-1]) * R
CT = np.sum(delta_r * df_axial[1:] * N_B / (0.5 * U_inf ** 2 * np.pi * R ** 2))
CP = np.sum(delta_r * df_tan[1:] * mu[1:] * N_B * R * Omega / (0.5 * U_inf ** 3 * np.pi * R ** 2))
if doPlot == 1:
        fig1 = plt.figure(figsize=(8, 4))
        plt.title('Flow Field downstream of the rotor')
        plt.plot(mu, U_inf*a_axial, 'k-', label=r'$U_{ax}$')
        plt.plot(mu, U_inf*a_tan, 'k--', label=r'$U_{tan}$')
        plt.grid()
        plt.xlabel(r'$r/R$')
        plt.legend()
        plt.savefig("figures/flowfield",bbox_inches='tight')

