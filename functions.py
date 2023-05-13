# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def velocity3d_vortex_filament(GAMMA, XV1, XV2, XVP1, RV):
    # USES BIOT-SAVART LAW TO COMPUTE INDUCED VELOCITIES U, V AND W
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

    return U, V, W
