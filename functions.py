# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def a_Glauert(CT):
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    
    if CT < CT2:
        a = 0.5 - 0.5 * np.sqrt(1 - CT)
    elif CT >= CT2:
        a = 1 + (CT - CT1) / (4 * (np.sqrt(CT1) - 1))

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
        a_new = a_Glauert(CdT)

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

    return [mu, a_new, aline, alpha, phi, df_axial, df_tan, gamma, CL, CD, alpha]


def velocity3d_vortex_filament(GAMMA, XV1, XV2, XVP1, RV):
    # USES BIOS-SAVART LAW TO COMPUTE INDUCED VELOCITIES U, V AND W
    # AT A TARGET POINT XVP1 FROM A VORTEX FILAMENT OF RADIUS RV AND
    # PASSING THROUGH POINTS XV1 AND XV2
    X1 = XV1[0]; Y1 = XV1[1]; Z1 = XV1[2]
    X2 = XV2[0]; Y2 = XV2[1]; Z2 = XV2[2]

    XP = XVP1[0]; YP = XVP1[1]; ZP = XVP1[2]

    R1 = np.sqrt(np.power((XP - X1), 2) + np.power((YP - Y1), 2) + np.power((ZP - Z1), 2))
    R2 = np.sqrt(np.power((XP - X2), 2) + np.power((YP - Y2), 2) + np.power((ZP - Z2), 2))
    R1XR2_X = (YP - Y1) * (ZP - Z2) - (ZP - Z1) * (YP - Y2)
    R1XR2_Y = -(XP - X1) * (ZP - Z2) + (ZP - Z1) * (XP - X2)
    R1XR2_Z = (XP - X1) * (YP - Y2) - (YP - Y1) * (XP - X2)
    R1XR_SQR = np.power(R1XR2_X, 2) + np.power(R1XR2_Y, 2) + np.power(R1XR2_Z, 2)
    R0R1 = (X2 - X1) * (XP - X1) + (Y2 - Y1) * (YP - Y1) + (Z2 - Z1) * (ZP - Z1)
    R0R2 = (X2 - X1) * (XP - X2) + (Y2 - Y1) * (YP - Y2) + (Z2 - Z1) * (ZP - Z2)

    if R1XR_SQR < np.power(RV, 2):
        R1XR_SQR=np.power(RV, 2)
    if R1 < RV:
        R1 = RV
    if R2 < RV:
        R2 = RV

    K = GAMMA / 4 / np.pi / R1XR_SQR * (R0R1 / R1 - R0R2 / R2)

    U = K * R1XR2_X
    V = K * R1XR2_Y
    W = K * R1XR2_Z

    return np.array([U, V, W])


def downstreamLine(r, theta, chord, twist, Uax, Utan, R, Loutlet, dx):
    # this function discretizes a vorticity line
    N = int((Loutlet*R - chord)/ dx) # number of points along each line
    thetas = np.zeros((N)); rs = np.ones((N)) * r; xs = np.zeros((N))
    # compute the first two points along the chord
    thetas[0] = theta
    thetas[1] = theta - (np.cos(twist)*chord / r)
    xs[0] = 0; xs[1] = chord*np.sin(twist)
    # compute the rest of the downstream vortex lines
    dtheta = dx * (Utan/Uax) / r 
    xs[2:] = xs[1] + (np.arange(1,N-1) * dx)

    thetas[2:] = thetas[1] + (np.arange(1,N-1) * dtheta)
    coords = np.vstack([xs, rs, thetas])
    return coords #[x,r,theta]


def cyl2cart(arr): # array must be shaped as [ndim, npoints], with ndim = [x r theta]
    if arr.ndim == 1:
        arr = arr[:,np.newaxis]
    newarr = arr * 0
    newarr[0,:] = arr[0,:]
    newarr[1,:] = np.sin(arr[2,:])*arr[1,:] #y - coordinate
    newarr[2,:] = np.cos(arr[2,:])*arr[1,:] #z - cordinate
    return newarr #now x-y-z