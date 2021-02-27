import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
    
    #Stomach glucose
    dS = b9*H-b8*S
    # Colour-blind friendly palette (use nice colors)
cb_palette = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

def ode_closed_loop(t,x,b)
    # Moste of the regular parameters
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b18, b19 b21, b22, b23, b25, b27 = b 
    S, H, L, G, C, I, W, E, M, A, Y, Q = X

    
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
    u = np.heaviside(Ge, G)  
    dE = c0 + (c1/(c2 + I*e))*(Ge -G)*u(Ge - G) - c3*E  
    
    # liver glucose [7]
    dC = b23 - b25*I*e - b22*G + b21*E - b5*
    
    # dynamics of glucose mass in muscle tissue [8]
    dM = 0.1*(v/f)*b3*G*I*e – b27*M

    # Adipose tissue glucose mass (A) [9]
    dA = k8* GLUT4m*(GtA)/(KmG4 + Gt!) + GLUT1*(GtA)/(KmG1 + GtA) – Kgluc*A

    # Dynamics of plasma leptin (Y) [10]
    dY = b13*A*Fat – b14*Y

    # Dynamics of ghrelin concentration n plasma (Q) [11]
    dQ = b12*exp(IS)*exp(ml) – b11*Q

    # Glucose intake [12]
    dH = (b17*Q)/(b18*Y + 1) * exp (rl) - B19*G*H - b9*H

    # Linking the whole body model with the cellular one [13]
    dINSA = -p2U*INSA + P2U*(I-Ib)

    # Linking the whole body model with the cellular one [14]
    dGtA = - q1*GtA + q2*(G – Gb)

    #Ranges 
    range_G = [4.5, 11] # mM
    range_I = [38, 400] #pM
    range_W = [5, 50]
    range_E = [28.68, 47.04]
    range_C = [0, 8]
    range_M = [2, 13]
    range_A = [30, 120 ]
    range_Y = [0, 0.6]
    range_Q = [8, 1146]
    # range H, S, L = none



# initial condition
# time points
# solve ODE
sol = integrate.solve_ivp(model,t,C_start, method = 'RK45')

