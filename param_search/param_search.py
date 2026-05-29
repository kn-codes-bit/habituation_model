import random
import numpy as np
from deap import base, creator, tools, algorithms
from scipy.optimize import newton
from scipy.optimize import minimize
import pandas as pd
from numba import jit
import datetime
import time
import os
import multiprocessing as mp
from functools import partial

Lr = 0.013
gr = 1.

# =========================
# hyperparameters for GA
# =========================
POP_SIZE = 200 # Population size
NGEN = 10000  # Number of generations
CXPB = 0.5  # Cross-probability
MUTPB = 0.3 # Mutation Probability
SEED_num = 400
# mutate
mu = 0.0
sigma = 0.04 
indpb = 0.2  # Probability of mutation for each parameter
# select
tournsize = 2
N_local_ref = 20
# =========================
# parameter for Parallelization
# =========================
N_CPU = 1
N_PROC = min(SEED_num, N_CPU)

# =========================
# parameter range（θ > 0）
# =========================
LOW  = np.array([1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3])
HIGH = np.array([100., 100., 100., 100., 100., 100., 1.])

LOG_LOW  = np.log10(LOW)
LOG_HIGH = np.log10(HIGH)

# read data
df = pd.read_csv(os.path.join("param_search/data/Q_data.csv"))
data_Q=df.to_numpy()
df = pd.read_csv(os.path.join("param_search/data/C1_data.csv"))
data_C1=df.to_numpy()
df = pd.read_csv(os.path.join("param_search/data/C2_data.csv"))
data_C2=df.to_numpy()


# para: 
# para[0]=gf, para[1]=gfs, para[2]=ade,
# para[3]=ades, para[4]=Dif, para[5]=gam, para[6]=rv
# para2:
# para2[0]=Pf_u, para2[1]=Pr_u,
# para2[2]=Pf_s, para2[3]=Pr_s,
# para2[4]=lam1, para2[5]=lam2,
# para2[6]=v11, para2[7]=v21,
# para2[8]=v12, para2[9]=v22,

def loss_function(individual):
    para = np.array(individual)
    para2 = [((para[2]+para[6]*para[4])*para[0]+para[4]*gr)/((para[2]+para[4]+para[6]*para[4])*para[2]), (para[6]*para[4]*para[0] + (para[2]+para[4])*gr)/((para[2]+para[4]+para[6]*para[4])*para[2]), 
               ((para[2]+para[6]*para[4])*para[1]+para[4]*gr)/((para[2]+para[6]*para[4])*para[3]+para[2]*para[4]), (para[6]*para[4]*para[1] + (para[3]+para[4])*gr)/((para[2]+para[6]*para[4])*para[3]+para[2]*para[4]), 
               -0.5*(para[3]+para[2]+para[4]+para[6]*para[4])-0.5*np.sqrt((para[3]-para[2]+para[4]-para[6]*para[4])**2+4.*para[6]*para[4]**2), -0.5*(para[3]+para[2]+para[4]+para[6]*para[4])+0.5*np.sqrt((para[3]-para[2]+para[4]-para[6]*para[4])**2+4.*para[6]*para[4]**2), 
                -0.5*(para[3]+para[2]+para[4]+para[6]*para[4])-0.5*np.sqrt((para[3]-para[2]+para[4]-para[6]*para[4])**2+4.*para[6]*para[4]**2)+para[2]+para[6]*para[4], para[6]*para[4], 
                para[4], -0.5*(para[3]+para[2]+para[4]+para[6]*para[4])+0.5*np.sqrt((para[3]-para[2]+para[4]-para[6]*para[4])**2+4.*para[6]*para[4]**2)+para[3]+para[4]]
    
    loss = 0.

    p_ini =[para2[0], para2[1]]

    model1 = Model_F1(p_ini=p_ini, para=para, para2=para2)
    model0 = Model_F0(p_ini=p_ini, para=para, para2=para2)
    model_cross = Model_cross(para=para, para2=para2)

    ############ vs data_Q ############
    data_size = len(data_Q[:,0])
    day = (int)(0)
    # day1~6 (b-Q)
    Tc_ini = Lr/(para[5]*(para2[3]-para2[2]))
    if Tc_ini<0:
        Tc_ini = Lr/(para[5]*(para2[1]-para2[0]))
    for _ in range(6):
        model1.p_ini=p_ini
        T_c = newton(model1.F1, Tc_ini, fprime=model1.dF1)
        loss += np.sum((data_Q[:,day]/60.-np.full(data_size, T_c))**2)/data_size
        p_ini = model_cross.cross_b_A(p_ini)
        Tc_ini = T_c
        day += (int)(1)
    
    # day7~8 (b-N)
    Tc_ini = Lr/(para[5]*(para2[1]-para2[0]))
    for _ in range(2):
        model0.p_ini=p_ini
        T_c = newton(model0.F0, Tc_ini, fprime=model0.dF0)
        loss += np.sum((data_Q[:,day]/60.-np.full(data_size, T_c))**2)/data_size
        p_ini = model_cross.cross_b_N(p_ini)
        Tc_ini = T_c
        day += (int)(1)
    
    # day9 (b-Q)
    Tc_ini = Lr/(para[5]*(para2[3]-para2[2]))
    if Tc_ini<0:
        Tc_ini = Lr/(para[5]*(para2[1]-para2[0]))
    model1.p_ini=p_ini
    T_c = newton(model1.F1, Tc_ini, fprime=model1.dF1)
    loss += np.sum((data_Q[:,day]/60.-np.full(data_size, T_c))**2)/data_size

    ############ vs data_C1 ############
    data_size = len(data_C1[:,0])
    day = (int)(0)
    # day1~5 (b-N)
    for _ in range(5):
        T_c = Lr/(para[5]*(para2[1]-para2[0]))
        loss += np.sum((data_C1[:,day]/60.-np.full(data_size, T_c))**2)/data_size
        day += (int)(1)
    
    # day 6 (b-Q)
    p_ini =[para2[0], para2[1]]
    Tc_ini = Lr/(para[5]*(para2[3]-para2[2]))
    if Tc_ini<0:
        Tc_ini = Lr/(para[5]*(para2[1]-para2[0]))
    model1.p_ini=p_ini
    T_c = newton(model1.F1, Tc_ini, fprime=model1.dF1)
    loss += np.sum((data_C1[:,day]/60.-np.full(data_size, T_c))**2)/data_size
    p_ini = model_cross.cross_b_A(p_ini)
    day += (int)(1)

    # day7~8 (b-N)
    Tc_ini = Lr/(para[5]*(para2[1]-para2[0]))
    for _ in range(2):
        model0.p_ini=p_ini
        T_c = newton(model0.F0, Tc_ini, fprime=model0.dF0)
        loss += np.sum((data_C1[:,day]/60.-np.full(data_size, T_c))**2)/data_size
        p_ini = model_cross.cross_b_N(p_ini)
        Tc_ini = T_c
        day += (int)(1)
    
    # day9 (b-Q)
    Tc_ini = Lr/(para[5]*(para2[3]-para2[2]))
    if Tc_ini<0:
        Tc_ini = Lr/(para[5]*(para2[1]-para2[0]))
    model1.p_ini=p_ini
    T_c = newton(model1.F1, Tc_ini, fprime=model1.dF1)
    loss += np.sum((data_C1[:,day]/60.-np.full(data_size, T_c))**2)/data_size

    ############ vs data_C2 ############
    data_size = len(data_C2[:,0])
    day = (int)(0)
    # day1~9 (b-N)
    for _ in range(9):
        T_c = Lr/(para[5]*(para2[1]-para2[0]))
        loss += np.sum((data_C2[:,day]/60.-np.full(data_size, T_c))**2)/data_size
        day += (int)(1)
    
    return loss

@jit(nopython=True)
def F0_numba(Tc, p_ini, para, para2):
    return Lr - (para[5]/((1.+para[6])*para[4] + para[2]))*((gr-para[0])*Tc + ((p_ini[1]-p_ini[0])-(para2[1]-para2[0]))*(1.-np.exp(-((1.+para[6])*para[4]+para[2])*Tc)))

@jit(nopython=True)
def dF0_numba(Tc, p_ini, para, para2):
    return -(para[5]/((1.+para[6])*para[4] + para[2]))*(gr-para[0] + ((1.+para[6])*para[4]+para[2])*((p_ini[1]-p_ini[0])-(para2[1]-para2[0]))*np.exp(-((1.+para[6])*para[4]+para[2])*Tc))

@jit(nopython=True)
def F1_numba(Tc, p_ini, para, para2):
    return Lr -para[5]*(para2[3]-para2[2])*Tc+(para[5]/(para2[6]*para2[9]-para2[8]*para2[7]))*((para2[9]*(p_ini[0]-para2[2])-para2[8]*(p_ini[1]-para2[3]))*(para2[6]-para2[7])*(1.-np.exp(para2[4]*Tc))/(-para2[4]) + (-para2[7]*(p_ini[0]-para2[2])+para2[6]*(p_ini[1]-para2[3]))*(para2[8]-para2[9])*(1.-np.exp(para2[5]*Tc))/(-para2[5]))

@jit(nopython=True)
def dF1_numba(Tc, p_ini, para, para2):
    return -para[5]*(para2[3]-para2[2])+(para[5]/(para2[6]*para2[9]-para2[8]*para2[7]))*((para2[9]*(p_ini[0]-para2[2])-para2[8]*(p_ini[1]-para2[3]))*(para2[6]-para2[7])*np.exp(para2[4]*Tc) + (-para2[7]*(p_ini[0]-para2[2])+para2[6]*(p_ini[1]-para2[3]))*(para2[8]-para2[9])*np.exp(para2[5]*Tc))

@jit(nopython=True)
def cross_b_N_numba(T, p_ini, para, para2):
    pf0 = para2[0]+(p_ini[0]-p_ini[1]-para2[0]+para2[1])*np.exp(-(para[2]+para[4]+para[6]*para[4])*T)/(1+para[6]) + (para[6]*p_ini[0]+p_ini[1]-para[6]*para2[0]-para2[1])*np.exp(-para[2]*T)/(1+para[6])
    pr0 = para2[1]-para[6]*(p_ini[0]-p_ini[1]-para2[0]+para2[1])*np.exp(-(para[2]+para[4]+para[6]*para[4])*T)/(1+para[6]) + (para[6]*p_ini[0]+p_ini[1]-para[6]*para2[0]-para2[1])*np.exp(-para[2]*T)/(1+para[6])
    return [pf0, pr0]

@jit(nopython=True)
def cross_b_A_numba(t_on, t_off, p_ini, para, para2):
    pf0 = para2[2] + para2[6]*(para2[9]*(p_ini[0]-para2[2])-para2[8]*(p_ini[1]-para2[3]))*np.exp(para2[4]*t_on)/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[8]*(-para2[7]*(p_ini[0]-para2[2])+para2[6]*(p_ini[1]-para2[3]))*np.exp(para2[5]*t_on)/(para2[6]*para2[9]-para2[8]*para2[7])
    pr0 = para2[3] + para2[7]*(para2[9]*(p_ini[0]-para2[2])-para2[8]*(p_ini[1]-para2[3]))*np.exp(para2[4]*t_on)/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[9]*(-para2[7]*(p_ini[0]-para2[2])+para2[6]*(p_ini[1]-para2[3]))*np.exp(para2[5]*t_on)/(para2[6]*para2[9]-para2[8]*para2[7])
    return cross_b_N_numba(t_off, [pf0, pr0], para, para2)

class Model_F0:
    def __init__(self, p_ini, para, para2):
        self.p_ini = p_ini
        self.para = para
        self.para2 = para2
    
    def F0(self, Tc):
        return F0_numba(Tc, self.p_ini, self.para, self.para2)
    
    def dF0(self, Tc):
        return dF0_numba(Tc, self.p_ini, self.para, self.para2)

class Model_F1:
    def __init__(self, p_ini, para, para2):
        self.p_ini = p_ini
        self.para = para
        self.para2 = para2
    
    def F1(self, Tc):
        return F1_numba(Tc, self.p_ini, self.para, self.para2)
    
    def dF1(self, Tc):
        return dF1_numba(Tc, self.p_ini, self.para, self.para2)

class Model_cross:
    def __init__(self, para, para2):
        self.para = para
        self.para2 = para2
    
    def cross_b_N(self, p_ini):
        return cross_b_N_numba(24., p_ini, self.para, self.para2)
    
    def cross_b_A(self, p_ini):
        return cross_b_A_numba(8., 16., p_ini, self.para, self.para2)

# =========================
# z in [0,1] → parameter 
# =========================
def decode_individual(ind):

    z = np.asarray(ind)

    log_theta = LOG_LOW + z * (LOG_HIGH - LOG_LOW)

    theta = 10.**log_theta

    return theta

def evaluate_log_individual(ind):
    try:
        theta = decode_individual(ind)
        if np.any(LOW>theta) or np.any(HIGH<theta):
            return (1e20,)
        
        loss = loss_function(theta)

        if not np.isfinite(loss):
            return (1e20,)

        return (loss,)

    except Exception as e:
        return (1e20,)

def evaluate_log_individual_local(ind):
    try:
        theta = decode_individual(ind)
        if np.any(LOW>theta) or np.any(HIGH<theta):
            return (1e20,)
        
        loss = loss_function(theta)

        if not np.isfinite(loss):
            return (1e20,)

        return (loss,)

    except Exception as e:
        return (1e20,)

def clip_log(operator):
    def wrapper(*args, **kwargs):
        offspring = operator(*args, **kwargs)

        # "mutate" returns (ind,)
        # "mate" returns (ind1, ind2)
        for ind in offspring:
            for i in range(len(ind)):
                ind[i] = min(max(ind[i], 0.0), 1.0)

        return offspring
    return wrapper

def save_logbook(logbook, seed, path_GA_his):
    filename = os.path.join(path_GA_his, f"seed_{seed:03d}.csv")
    df = pd.DataFrame(logbook)
    df.to_csv(filename,index=False)


# =========================
# Fitness / Individual
# =========================
if not hasattr(creator, "FitnessMin"):
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))    # Setting minimization problem
if not hasattr(creator, "Individual"):
    creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

toolbox.register(
    "individual",
    tools.initCycle,
    creator.Individual,
    [
        lambda: random.uniform(0.0, 1.0)
        for _ in range(len(LOG_LOW))
    ],
    n=1
)

toolbox.register(
    "population",
    tools.initRepeat,
    list,
    toolbox.individual
)

toolbox.register("evaluate", evaluate_log_individual)

toolbox.register("mate", tools.cxBlend, alpha=0.5)
toolbox.decorate("mate", clip_log)

toolbox.register(
    "mutate",
    tools.mutGaussian,
    mu=mu,
    sigma=sigma,  
    indpb=indpb    # Probability of mutation for each parameter
)
toolbox.decorate("mutate", clip_log)

toolbox.register(
    "select",
    tools.selTournament,
    tournsize=tournsize
)

def run_one_seed(RANDOM_SEED, path_GA_his):

    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    pop = toolbox.population(n=POP_SIZE)

    # =========================
    # Hall of Fame（log-space）
    # =========================
    hof = tools.HallOfFame(N_local_ref)

    # =========================
    # GAhistory（fitness is tuple）
    # =========================
    stats = tools.Statistics(lambda ind: ind.fitness.values[0])
    stats.register("avg", np.mean)
    stats.register("min", np.min)

    # =========================
    # GA run
    # =========================
    pop, logbook = algorithms.eaSimple(
        population=pop,
        toolbox=toolbox,
        cxpb=CXPB,
        mutpb=MUTPB,
        ngen=NGEN,
        stats=stats,
        halloffame=hof,
        verbose=False
    )
    save_logbook(logbook, RANDOM_SEED, path_GA_his)

    # =========================
    # best param（decode）
    # =========================
    def objective_local(ind):
        return evaluate_log_individual_local(ind)[0]

    best_res = None
    best_loss = np.inf
    method_name = "Nelder-Mead"
    for ind in hof:
        z0 = np.array(ind)

        res = minimize(
            objective_local,
            x0=z0,
            method=method_name,
            options={"maxiter": 5000}
        )

        if res.fun < best_loss:
            best_loss = res.fun
            best_res = res

    best_z = best_res.x
    best_params = decode_individual(best_z)
    return{
        "seed": RANDOM_SEED,
        "GA_loss": float(hof[0].fitness.values[0]),
        "local_loss": float(best_res.fun),
        "success": bool(best_res.success),
        "message": str(best_res.message),
        "params_z": best_z.tolist(),
        "params": best_params.tolist()
    }

if __name__ == "__main__":

    # warm up for numba
    loss_function([0.002164826579262351, 7.500441829972404, 0.027774541234562047, 0.22239744342485224, 0.40302120018610277, 0.004977511927778899, 0.1])
    
    # make folder for results
    now = datetime.datetime.now()
    path_base = os.path.join("param_search/param_estimation_result")
    path_base = os.path.join(path_base, now.strftime('%Y%m%d_%H%M%S'))
    os.makedirs(path_base, exist_ok=True)
    path_GA_his = os.path.join(path_base, "GA_loss_history")
    os.makedirs(path_GA_his, exist_ok=True)

    run_one_seed_with_path = partial(run_one_seed, path_GA_his=path_GA_his)
    
    hyperparams = {
        "parameter_bounds": {
            "LOW": LOW.tolist(),
            "HIGH": HIGH.tolist()
        },
        "GA":{
        "POP_SIZE": POP_SIZE,
        "NGEN": NGEN,
        "CXPB": CXPB,
        "MUTPB": MUTPB,
        "SEED_num": SEED_num,
        "mut_mu": mu,
        "mut_sigma": sigma,
        "mut_indpb": indpb, 
        "tournsize": tournsize, 
        "N_local_ref ": N_local_ref
        }
    }
    filename = os.path.join(path_base, "hyperparams.csv")
    pd.DataFrame([hyperparams]).to_csv(filename, index=False)
    
    result_file = os.path.join(path_base, "results.csv")
    t0 = time.perf_counter()
    rows = []
    with mp.Pool(processes=N_PROC) as pool:
        for i, res in enumerate(
            pool.imap_unordered(run_one_seed_with_path, range(SEED_num)),
            start=1
        ):
            rows.append(res)

            # Save the results so far
            df = pd.DataFrame(rows)
            df_sorted = df.sort_values(
                by="local_loss",
                ascending=True
            )
            df_sorted.to_csv(result_file, index=False)
            print(f"[{i}/{SEED_num}] seeds finished "
                  f"({100*i/SEED_num:.1f}%)")

    # Output to CSV
    df = pd.DataFrame(rows)
    df_sorted = df.sort_values(
        by="local_loss",
        ascending=True
    )
    print(df_sorted.head())
    df_sorted.to_csv(result_file, index=False)

    now = datetime.datetime.now()
    print(now.strftime('%Y%m%d_%H%M%S'))
    t1 = time.perf_counter()
    print(f"elapsed: {t1 - t0:.3f} sec")
