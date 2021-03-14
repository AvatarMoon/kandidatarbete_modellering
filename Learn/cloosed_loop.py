import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
    
   
    # Colour-blind friendly palette (use nice colors)
cb_palette = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

def closed_loop(t,x):
    # Moste of the regular parameters
    b1 = 0.0059
    b2 = 0.1262
    b3 = 0.00005
    b4 = 0.4543
    b5 = 0.185
    b6 = 0.0102
    b7 = 0.03
    b8 = 0.022
    b9 = 0.022
    b10 = 0.022
    b11 = 0.02
    b12 = 28.66
    b13 = 0.0000095
    b14 = 0.0278
    b17 = 0.02 # random siffra
    b18 = 0.35
    b19 = 0.004
    b21 = 0.00876
    b22 = 0.0021
    b23 = 0.08
    b25 = 0.00026
    b27 = 0.014
    f = 0.9
    v = 15
    c = 0.1060
    #s = 0.03
    c0 = 1.8854
    c1 = 198
    c2 = 94
    e = 1
    Ge = 5
    c3 = 0.0554
    k8 = 0.5275
    # GLUT4 = 50.8472
    GtA = 135
    GLUT1 = 0.0283
    KmG4 = 146.851
    KmG1 = 1.082 
    Kgluc = 0.25
    Fat = 22
    #l = 0.006
    #m = 0.04
    #r = 0.04
    #p2u = 0.033
    # Ib =
    q1 = 0.0031
    q2 = 0.4054
    # Gb =
    H = 10
    """global b
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14,b17, b18, b19, b21, b22, b23, b25, b27 = b """
    S, L, G, C, I, W, E, M, A, Y, Q, H, INSA, GtA = x 

    #Stomach glucose [1]
    dS = b9*H-b8*S
    
    # Intestine glucose [2]
    dL = b8*S-b10*L
    
    # plasma glucose [3]
    dG = f*b10*L/v + f*b5*C/v - b1*G-b3*I*G
    
    # plasma insulin [4]
    dI = b4*G + c*W*G-b2*I
    
    # plasma incretin [5]
    dW = b6*L - b7*W + S
    
    # plasma glucagon [6] # 
    # u = np.heaviside(Ge, G)  # lägg till u enligt artikel
    dE = c0 + (c1/5)*(Ge -G)*(Ge - G) - c3*E  
    
    # liver glucose [7]
    dC = b23 - b25*I*e - b22*G + b21*E - b5*C
    
    # dynamics of glucose mass in muscle tissue [8]
    dM = 0.1*(v/f)*b3*G*I*e - b27*M

    # Adipose tissue glucose mass (A) [9]
    dA = k8*(GtA)/(KmG4 + GtA) + GLUT1*(GtA)/(KmG1 + GtA) - Kgluc*A

    # Dynamics of plasma leptin (Y) [10]
    dY = b13*A*Fat - b14*Y

    # Dynamics of ghrelin concentration n plasma (Q) [11]
    dQ = b12 - b11*Q

    # Glucose intake [12]
    dH = (b17*Q)/(b18*Y + 1) - b19*G*H - b9*H

    # Linking the whole body model with the cellular one [13]
    dINSA = -INSA + (I-2)

    # Linking the whole body model with the cellular one [14]
    dGtA = - q1*GtA + q2*(G - 2)

    return [dS, dL, dG, dI, dW, dE, dC, dM, dA, dY, dQ, dH, dINSA, dGtA]

    # Ranges 
    """ range_G = [4.5, 11] # mM
    range_I = [38, 400] #pM
    range_W = [5, 50]
    range_E = [28.68, 47.04]
    range_C = [0, 8]
    range_M = [2, 13]
    range_A = [30, 120]
    range_Y = [0, 0.6]
    range_Q = [8, 1146]
    # range H, S, L = none """

t_vec = np.linspace(0.1, 8, num=50)
time_span = [t_vec[0], t_vec[-1]] 
x0 = [10, 20, 15, 10, 20, 15, 10, 20, 15, 10, 20, 15, 10, 20]
sol = integrate.solve_ivp(closed_loop, time_span, x0, method="LSODA", t_eval=t_vec)

ymodel = sol.y[1]

print(ymodel)

plt.plot(t_vec, ymodel)
plt.title("Simulated model")
plt.show()
"""def solve(b, ):
    # initial condition
    
    # time points
    t_vec = np.linspace(0.1, 2.0, num=50)
    time_span = [t_vec[0], t_vec[-1]]
    # solve ODE-system
    

    y_model = sol.y[1]

    return y_model """



