import numpy as np
from scipy.optimize import newton
import pandas as pd
from numba import jit
import os
import ast
import matplotlib.pyplot as plt

Lr = 0.013
gr = 1.

eps = 0.01 # perturbation rate for local sensitivity analysis (e.g., eps=0.01 means 1% perturbation)
rank = 1 # specify the rank of the parameters to use (rank=1 is the best)

foldername = "foldername" # specify folder name in "param_estimation_result"


# Get a SEED_num
path_base = os.path.join("param_search", "param_estimation_result", foldername)
filename = os.path.join(path_base, "results.csv")
df = pd.read_csv(filename)
params = ast.literal_eval(df.iloc[rank-1]["params"])

def cross_time(param):
    # para: 
    # para[0]=gf, para[1]=gfs, para[2]=ade,
    # para[3]=ades, para[4]=Dif, para[5]=gam, para[6]=rv
    # para2:
    # para2[0]=Pf_u, para2[1]=Pr_u,
    # para2[2]=Pf_s, para2[3]=Pr_s,
    # para2[4]=lam1, para2[5]=lam2,
    # para2[6]=v11, para2[7]=v21,
    # para2[8]=v12, para2[9]=v22,
    para = np.array(param)
    para2 = [((para[2]+param[6]*para[4])*para[0]+para[4]*gr)/((para[2]+para[4]+param[6]*para[4])*para[2]), (param[6]*para[4]*para[0] + (para[2]+para[4])*gr)/((para[2]+para[4]+param[6]*para[4])*para[2]), 
               ((para[2]+param[6]*para[4])*para[1]+para[4]*gr)/((para[2]+param[6]*para[4])*para[3]+para[2]*para[4]), (param[6]*para[4]*para[1] + (para[3]+para[4])*gr)/((para[2]+param[6]*para[4])*para[3]+para[2]*para[4]), 
               -0.5*(para[3]+para[2]+para[4]+param[6]*para[4])-0.5*np.sqrt((para[3]-para[2]+para[4]-param[6]*para[4])**2+4.*param[6]*para[4]**2), -0.5*(para[3]+para[2]+para[4]+param[6]*para[4])+0.5*np.sqrt((para[3]-para[2]+para[4]-param[6]*para[4])**2+4.*param[6]*para[4]**2), 
                -0.5*(para[3]+para[2]+para[4]+param[6]*para[4])-0.5*np.sqrt((para[3]-para[2]+para[4]-param[6]*para[4])**2+4.*param[6]*para[4]**2)+para[2]+param[6]*para[4], param[6]*para[4], 
                para[4], -0.5*(para[3]+para[2]+para[4]+param[6]*para[4])+0.5*np.sqrt((para[3]-para[2]+para[4]-param[6]*para[4])**2+4.*param[6]*para[4]**2)+para[3]+para[4]]
    
    T_c_lis = []
    h = 0.0001
    t_learn = np.arange(0., 8.+h, h)
    t_rest = np.arange(8., 24.+h, h)

    p_ini =np.array([para2[0], para2[1], 0.])

    model_cross = Model_cross(para=para, para2=para2)

    ############ vs data_Q ############
    day = (int)(0)
    # day1~6 (b-A)
    for _ in range(6):
        T_c = solver_A(t_learn, t_rest, p_ini, para)
        T_c_lis.append(T_c)
        p_ini = model_cross.cross_b_Q(p_ini)
        day += (int)(1)
    
    # day7~8 (b-N)
    for _ in range(2):
        T_c = solver_N(t_learn, t_rest, p_ini, para)
        T_c_lis.append(T_c)
        p_ini = model_cross.cross_b_N(p_ini)
        day += (int)(1)
    
    # day9 (b-A)
    T_c = solver_A(t_learn, t_rest, p_ini, para)
    T_c_lis.append(T_c)

    ############ vs data_C1 ############
    day = (int)(0)
    # day1~5 (b-N)
    for _ in range(5):
        T_c = Lr/(para[5]*(para2[1]-para2[0]))
        T_c_lis.append(T_c)
        day += (int)(1)
    
    # day 6 (b-A)
    T_c = solver_A(t_learn, t_rest, p_ini, para)
    T_c_lis.append(T_c)
    p_ini = model_cross.cross_b_Q(p_ini)
    day += (int)(1)

    # day7~8 (b-N)
    for _ in range(2):
        T_c = solver_N(t_learn, t_rest, p_ini, para)
        T_c_lis.append(T_c)
        p_ini = model_cross.cross_b_N(p_ini)
        day += (int)(1)
    
    # day9 (b-A)
    T_c = solver_A(t_learn, t_rest, p_ini, para)
    T_c_lis.append(T_c)

    ############ vs data_C2 ############
    day = (int)(0)
    # day1~9 (b-N)
    for _ in range(9):
        T_c = Lr/(para[5]*(para2[1]-para2[0]))
        T_c_lis.append(T_c)
        day += (int)(1)
    
    return T_c_lis

# Definition of ODE system
@jit(nopython=True)
def system_sti(t, y, params): 
    pf, pr, l= y
    gf, gfs, ade, ades, Dif, gam, rv= params
    dpf = Dif*(pr-pf) -ades*pf + gfs
    dpr = rv*Dif*(pf-pr) -ade*pr + gr
    dl  = gam*(pr-pf)
    return np.array([dpf, dpr, dl])

@jit(nopython=True)
def system_nsti(t, y, params): 
    pf, pr, l= y
    gf, gfs, ade, ades, Dif, gam, rv= params
    dpf = Dif*(pr-pf) -ade*pf + gf
    dpr = rv*Dif*(pf-pr) -ade*pr + gr
    dl  = gam*(pr-pf)
    return np.array([dpf, dpr, dl])

# Definition of update function in RK4
def rk4_step_sti(t, y, h, params):
    k1 = system_sti(t, y, params)
    k2 = system_sti(t+0.5*h, y+0.5*h*k1, params)
    k3 = system_sti(t+0.5*h, y+0.5*h*k2, params)
    k4 = system_sti(t+h, y+h*k3, params)
    return y + (h/6.)*(k1+2.*k2+2.*k3+k4)

def rk4_step_nsti(t, y, h, params):
    k1 = system_nsti(t, y, params)
    k2 = system_nsti(t+0.5*h, y+0.5*h*k1, params)
    k3 = system_nsti(t+0.5*h, y+0.5*h*k2, params)
    k4 = system_nsti(t+h, y+h*k3, params)
    return y + (h/6.)*(k1+2.*k2+2.*k3+k4)

# ODE solver (RK4)
def solver_A(t_learn, t_rest, y_0, params):
    n_learn = len(t_learn)
    n_rest  = len(t_rest)
    y_1 = np.array([0.,0.,0.])
    for i in range(1, n_learn):
        h = t_learn[i] - t_learn[i-1]
        y_1 = rk4_step_sti(t_learn[i-1], y_0, h, params)
        y_0 = y_1
        if y_1[2] >= Lr:
            return t_learn[i]
    for i in range(1, n_rest):
        h = t_rest[i] - t_rest[i-1]
        y_1 = rk4_step_nsti(t_rest[i-1], y_0, h, params)
        y_0 = y_1
        if y_1[2] >= Lr:
            return t_rest[i]
    return None

def solver_N(t_learn, t_rest, y_0, params):
    n_learn = len(t_learn)
    n_rest  = len(t_rest)
    y_1 = np.array([0.,0.,0.])
    for i in range(1, n_learn):
        h = t_learn[i] - t_learn[i-1]
        y_1 = rk4_step_nsti(t_learn[i-1], y_0, h, params)
        y_0 = y_1
        if y_1[2] >= Lr:
            return t_learn[i]
    for i in range(1, n_rest):
        h = t_rest[i] - t_rest[i-1]
        y_1 = rk4_step_nsti(t_rest[i-1], y_0, h, params)
        y_0 = y_1
        if y_1[2] >= Lr:
            return t_rest[i]
    return None

@jit(nopython=True)
def cross_b_N_numba(T, p_ini, para, para2):
    pf0 = para2[0]+(p_ini[0]-p_ini[1]-para2[0]+para2[1])*np.exp(-(para[2]+para[4]+para[6]*para[4])*T)/(1+para[6]) + (para[6]*p_ini[0]+p_ini[1]-para[6]*para2[0]-para2[1])*np.exp(-para[2]*T)/(1+para[6])
    pr0 = para2[1]-para[6]*(p_ini[0]-p_ini[1]-para2[0]+para2[1])*np.exp(-(para[2]+para[4]+para[6]*para[4])*T)/(1+para[6]) + (para[6]*p_ini[0]+p_ini[1]-para[6]*para2[0]-para2[1])*np.exp(-para[2]*T)/(1+para[6])
    return [pf0, pr0, 0.]

@jit(nopython=True)
def cross_b_Q_numba(t_on, t_off, p_ini, para, para2):
    pf0 = para2[2] + para2[6]*(para2[9]*(p_ini[0]-para2[2])-para2[8]*(p_ini[1]-para2[3]))*np.exp(para2[4]*t_on)/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[8]*(-para2[7]*(p_ini[0]-para2[2])+para2[6]*(p_ini[1]-para2[3]))*np.exp(para2[5]*t_on)/(para2[6]*para2[9]-para2[8]*para2[7])
    pr0 = para2[3] + para2[7]*(para2[9]*(p_ini[0]-para2[2])-para2[8]*(p_ini[1]-para2[3]))*np.exp(para2[4]*t_on)/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[9]*(-para2[7]*(p_ini[0]-para2[2])+para2[6]*(p_ini[1]-para2[3]))*np.exp(para2[5]*t_on)/(para2[6]*para2[9]-para2[8]*para2[7])
    return cross_b_N_numba(t_off, [pf0, pr0], para, para2)

class Model_cross:
    def __init__(self, para, para2):
        self.para = para
        self.para2 = para2
    
    def cross_b_N(self, p_ini):
        return cross_b_N_numba(24., p_ini, self.para, self.para2)
    
    def cross_b_Q(self, p_ini):
        return cross_b_Q_numba(8., 16., p_ini, self.para, self.para2)

params_peps = [x for x in params]
params_meps = [x for x in params]
tmp = None
result = []


for i in range(len(params)):
    params_peps[i] = (1.+eps)*params[i]
    params_meps[i] = (1.-eps)*params[i]
    tmp = (np.array(cross_time(params_peps))-np.array(cross_time(params_meps)))/(2.*eps*np.array(cross_time(params)))
    params_peps[i] = params[i]
    params_meps[i] = params[i]
    result.append({"S_num": f"S{i}", "Values": tmp})
    

# make folder
path_base = os.path.join("param_search", "sensitivity_result", foldername)
os.makedirs(path_base, exist_ok=True)


df = pd.DataFrame(result)
filename = os.path.join(path_base, "result.csv")
df.to_csv(filename, index=False)

base_params = {"base_params" :params}
filename = os.path.join(path_base, "base_params.csv")
pd.DataFrame([base_params]).to_csv(filename, index=False)

##### plot #####
path_base = os.path.join("param_search", "sensitivity_result", foldername)

# Non-normalized local sensitvity
filename = os.path.join(path_base, "result.csv")
df = pd.read_csv(filename)
data = df["Values"].apply(lambda s: np.fromstring(s.strip("[]"), sep=" "))

param_num = len(df)
labels = ["gf", "gf_s", "a_de", "a_de^s", "Dif", "gam", "rv"]
Q_means = []
C1_means = []
C2_means = []

for i in range(param_num):
    Q_means.append(np.abs(data[i])[0:9].sum()/9.)
    C1_means.append(np.abs(data[i])[9:18].sum()/9.)
    C2_means.append(np.abs(data[i])[18:27].sum()/9.)

x = np.arange(len(labels))
width = 0.23

fig, ax = plt.subplots()
rects1 = ax.bar(x - width, Q_means, width, label='Q')
rects2 = ax.bar(x , C1_means, width, label='C1')
rects3 = ax.bar(x + width, C2_means, width, label='C2')

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('local sensitvity')
ax.legend()

filename = os.path.join(path_base, "result.png")
plt.savefig(filename, dpi=300, bbox_inches="tight")
plt.show()
