# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def a_Glauert(CT, yaw):
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    
    yaw_crit = 15

    if yaw == 0:
        if CT < CT2:
            a = 0.5 - 0.5 * np.sqrt(1 - CT)
        elif CT >= CT2:
            a = 1 + (CT - CT1) / (4 * (np.sqrt(CT1) - 1))
    elif yaw > 0:
        #if yaw < yaw_crit:
            a = 0.5 * np.cos(np.deg2rad(yaw)) - 0.5 * np.sqrt(np.cos(np.deg2rad(yaw)) ** 2 - CT)
        #elif yaw >= yaw_crit:
        #    a = 0

    return a


def PrandtlTipRootCorrection(mu, mu_root, mu_tip, TSR, N_B, a):
    temp1 = - N_B / 2 * (mu_tip - mu) / mu * np.sqrt(1 + ((TSR * mu) ** 2)/((1 - a) ** 2))
    Ftip = np.array(2 / np.pi * np.arccos(np.exp(temp1)))
    Ftip[np.isnan(Ftip)] = 0

    temp1 = N_B / 2 * (mu_root - mu) / mu * np.sqrt(1 + ((TSR * mu) ** 2) / ((1 - a) ** 2))
    Froot = np.array(2 / np.pi * np.arccos(np.exp(temp1)))
    Froot[np.isnan(Froot)] = 0
    
    return Froot * Ftip, Ftip, Froot


def loadBladeElement(U_axial, U_tan, chord, twist, polar_alpha, polar_CL, polar_CD):
    W2 = U_axial ** 2 + U_tan ** 2
    phi = np.arctan2(U_axial, U_tan)
    alpha = np.rad2deg(phi) - twist

    CL = np.interp(alpha, polar_alpha, polar_CL)
    CD = np.interp(alpha, polar_alpha, polar_CD)
    L = 0.5 * W2 * chord * CL
    D = 0.5 * W2 * chord * CD

    df_axial = L * np.cos(phi) + D * np.sin(phi)
    df_tan = L * np.sin(phi) - D * np.cos(phi)

    gamma = 0.5 * np.sqrt(W2) * CL * chord

    return df_axial, df_tan, alpha, phi, gamma, CL, CD


def solveStreamtube(U_inf, mu_1, mu_2, delta_mu, mu_root, mu_tip, Omega, R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD):
    A = np.pi * ((mu_2 * R) ** 2 - (mu_1 * R) ** 2)  # area streamtube
    mu = (mu_1 + mu_2) / 2  # centroide

    # initialize variables
    a = 0.3  # axial induction
    aline = 0.0  # tangential induction factor

    N_iter = 100
    # error limit for iteration process, in absolute value of induction
    error_iter = 0.00001

    for i in range(N_iter):
        U_axial = U_inf * (1 - a)  # axial velocity at rotor
        U_tan = Omega * mu * R * (1 + aline)  # tangential velocity at rotor

        # calculate loads in blade segment in 2D (N/m)
        df_axial, df_tan, alpha, phi, gamma, CL, CD = loadBladeElement(U_axial, U_tan, chord, twist, polar_alpha, polar_CL, polar_CD)      
        
        # total force for each annuli
        dF_axial = df_axial * R * delta_mu * N_B

        # calculate thrust coefficient at the streamtube
        CdT = dF_axial / (0.5 * A * U_inf ** 2)

        # calculate new axial induction, accounting for Glauert's correction
        a_new = a_Glauert(CdT, 0)

        # calculate aximuthal induction
        aline_new = df_tan * N_B / (2 * np.pi * U_inf * (1 - a) * Omega * 2 * (mu * R) ** 2)

        # correct new axial induction with Prandtl's correction
        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(mu, mu_root, mu_tip, Omega * R / U_inf, N_B, a_new)
        if (Prandtl < 0.0001):
            Prandtl = 0.0001  # avoid divide by zero
        a_new = a_new / Prandtl  # correct estimate of axial induction
        aline_new = aline_new / Prandtl  # correct estimate of azimuthal induction with Prandtl's correction

        # for improving convergence, weigh current and previous iteration of axial induction
        a = 0.75 * a + 0.25 * a_new
        aline = 0.75 * aline + 0.25 * aline_new

        # test convergence of solution, by checking convergence of axial induction
        if (np.abs(a - a_new) < error_iter):
            #print("iterations")
            #print(i)
            break

    return [mu, a_new, aline, alpha, phi, df_axial, df_tan, gamma, CL, CD]


def solveStreamtube_unc(U_inf, mu_1, mu_2, delta_mu, mu_root, mu_tip, Omega, R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD):
    A = np.pi * ((mu_2 * R) ** 2 - (mu_1 * R) ** 2)  # area streamtube
    mu = (mu_1 + mu_2) / 2  # centroide

    # initialize variables
    a = 0.3  # axial induction
    aline = 0.0  # tangential induction factor

    N_iter = 100
    # error limit for iteration process, in absolute value of induction
    error_iter = 0.00001

    for i in range(N_iter):
        U_axial = U_inf * (1 - a)  # axial velocity at rotor
        U_tan = Omega * mu * R * (1 + aline)  # tangential velocity at rotor

        # calculate loads in blade segment in 2D (N/m)
        df_axial, df_tan, alpha, phi, gamma, CL, CD = loadBladeElement(U_axial, U_tan, chord, twist, polar_alpha, polar_CL, polar_CD)      
        
        # total force for each annuli
        dF_axial = df_axial * R * delta_mu * N_B

        # calculate thrust coefficient at the streamtube
        CdT = dF_axial / (0.5 * A * U_inf ** 2)

        # calculate new axial induction, accounting for Glauert's correction
        a_new = a_Glauert(CdT, 0)

        # calculate aximuthal induction
        aline_new = df_tan * N_B / (2 * np.pi * U_inf * (1 - a) * Omega * 2 * (mu * R) ** 2)

    
        # for improving convergence, weigh current and previous iteration of axial induction
        a = 0.75 * a + 0.25 * a_new
        aline = 0.75 * aline + 0.25 * aline_new

        # test convergence of solution, by checking convergence of axial induction
        if (np.abs(a - a_new) < error_iter):
            #print("iterations")
            #print(i)
            break

    return [mu, a_new, aline, alpha, phi, df_axial, df_tan, gamma, CL, CD]


def solveStreamtube_yaw(U_inf, mu_1, mu_2, mu_root, mu_tip, Omega, R, N_B, chord, twist, polar_alpha, polar_CL, polar_CD, yaw, Phi_1, Phi_2, TSR):
    mu = (mu_1 + mu_2) / 2  # radial centroide
    Phi = (Phi_1 + Phi_2)/2 # angular centroide
    Area = Phi * ((mu_2 * R) ** 2 - (mu_1 * R) ** 2)  # area streamtube segment
    sigma = (N_B*chord)/(2*np.pi*mu*R) # blade solidity

    # initialize variables
    a = 0.0  # axial induction
    aline = 0.0  # tangential induction factor

    F = 0.5*(mu+0.4*mu**3+0.4*mu**5)
    xi = (0.6*a+1)*yaw
    K = 2*np.tan(xi/2)

    N_iter = 100
    # error limit for iteration process, in absolute value of induction
    error_iter = 0.00001

    for i in range(N_iter):
        U_axial = U_inf *(np.cos(yaw)-a*(1+F*K*np.sin(Phi)))+Omega*mu*R*aline*np.cos(Phi)*np.sin(xi)*(1+np.sin(Phi)*np.sin(xi))   # axial velocity at rotor
        U_tan = Omega*mu*R*(1+aline*np.cos(xi)*(1+np.sin(Phi)*np.sin(xi)))+U_inf*np.cos(Phi)*(a*np.tan(xi/2)*(1+F*K*np.sin(Phi))-np.sin(yaw))  # tangential velocity at rotor

        # calculate loads in blade segment in 2D (N/m)
        df_axial, df_tan, alpha, phi = loadBladeElement(U_axial, U_tan, chord, twist,
                                                        polar_alpha, polar_CL, polar_CD)      
        
        # total force for each annuli
        dF_axial = df_axial * R * (mu_2 - mu_1) * (Phi_2-Phi_1) * N_B

        # calculate thrust coefficient at the streamtube
        CdT = dF_axial / (0.5 * Area * U_inf ** 2)

        # calculate new axial induction, accounting for Glauert's correction
        a = a_Glauert(CdT, 0)

        # correct new axial induction with Prandtl's correction
        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(mu, mu_root, mu_tip, Omega * R / U_inf, N_B, a)
        if (Prandtl < 0.0001):
            Prandtl = 0.0001  # avoid divide by zero
        a = a / Prandtl  # correct estimate of axial induction

        # compute a_new
        W2 = U_axial ** 2 + U_tan ** 2
        A = (-8*np.pi)/((np.cos(xi/2))**2)
        B = 8*np.pi*(np.cos(Phi)+np.tan(xi/2)*np.sin(yaw))
        C = sigma*W2*df_axial*(Phi_2-Phi_1)/(U_inf**2*0.5 * W2 * chord)

        a_new = (-B+np.sqrt(B**2-4*A*C))/(2*A)

        # for improving convergence, weigh current and previous iteration of axial induction
        a = 0.75 * a + 0.25 * a_new

        # calculate azimuthal induction
        C = sigma*W2*df_tan*(Phi_2-Phi_1)/(U_inf**2*0.5 * W2 * chord)
        aline = C/(4*Prandtl*(np.cos(yaw)-a)*TSR*mu*np.pi*(1+(np.cos(xi))**2))
        aline = aline / Prandtl  # correct estimate of azimuthal induction with Prandtl's correction

        #dP = df_axial * U_inf * (np.cos(np.deg2rad(yaw)) - a)

        # test convergence of solution, by checking convergence of axial induction
        if (np.abs(a - a_new) < error_iter):
            #print("iterations")
            #print(i)
            break

    return [mu, a_new, aline, alpha, phi, df_axial, df_tan]
