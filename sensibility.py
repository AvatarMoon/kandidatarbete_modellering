import numpy as np 
import scipy.integrate as integrate 


def open_loop(t,x,b): 
    c0, c1, c2, c3, b1, b2, b3, b4, b5, b10, b21, b22, b23, b25, b27, b100, f, v = b

    Ge = 5
     
    # Says that the concentrations can't be lower than zero 
    x[x < 0] = 0.0 

    # Concentrations in the model as input 
    G, I, E, C, M, H = x 

    # Glucose plasma [1]  G0 closed loop = 5mM 
    dG = f*b10*H/v + f*b5*C/v - b1*G - b3*I*G 

    # Insulin plasma [2]  60 pM 
    dI = b4*G - b2*I    

    # Glucacon plasma [3] E0 closed loop = 34 pM 
    dE = c0 + (c1/(c2 + I))*(Ge - G)*np.heaviside(Ge-G,1) - c3*E  

    # GLucose liver [4] C0 closed loop = 3 mmol 
    dC = b23 - b25*I - b22*G + b21*E - b5*C 

    # Glucose musle [5] M0 = 2.5 mmol 
    dM = 0.1*(v/f)*b3*G*I - b27*M 

    # Glucos intake [6]  H0 = 200 mmol 
    dH = - b100*H*G 

    return [dG, dI, dE, dC, dM, dH] 

x0 = [5, 60, 34, 3, 2.5, 200]
t_vec = np.linspace(0, 240, num=20)
#start_values taken from openloop_optimering
# start_values = np.array(res.x)
start_values = np.array([1.885, 198, 94, 0.0554, 0.0059, 0.1262, 0.00005, 0.4543, 0.185, 0.022, 0.00876, 0.0021, 0.08, 0.00026, 0.014, 0.3, 0.9, 15])

def sensitivity(b,t,x):
    h = np.sqrt(np.finfo(np.float).eps) # Maskintoleransen, vår steglängd för finita differen 
    b_par = len(b)
    t_len = len(t)
    # Sensitivity analysis for each time step
    S = np.zeros([b_par, t_len * len(x)])
    time_span = [t[0], t[-1]]

    for n in range(len(b)):
        b1 = b.copy() 
        b2 = b.copy()  
        b1[n] += h 
        b2[n] -= h

        Sol_high = integrate.solve_ivp(open_loop, time_span, x, method='LSODA', args=(b1, ), t_eval = t)
        Sol_low = integrate.solve_ivp(open_loop, time_span, x, method='LSODA', args=(b2, ), t_eval= t)
        
        Sol_diff = (Sol_high.y-Sol_low.y)/(2*h)

        S[n,:] = Sol_diff.reshape(t_len*len(x))

    return S

S = sensitivity(start_values, t_vec, x0)

# Fisher matrix to make the covariance matrix
Fisher = 2 * S @ S.transpose()

print('Shape of S- and fisher matrix')
print(S.shape, Fisher.shape)

cov_mat = Fisher.transpose()

# Identification of the parameters
d_cov = np.diag(cov_mat)

var_coeff = np.square(d_cov)/start_values

print('Identification for each parameters')
print(var_coeff)


