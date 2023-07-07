import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize
import warnings
import seaborn as sns
import scipy.stats as st
from scipy.stats import ranksums
import pickle
import theano
import time
import pymc3 as pm
from pymc3 import HalfNormal as pm_pool_initialize
from pymc3 import HalfCauchy as pm_initialize
warnings.filterwarnings("ignore", category=DeprecationWarning)
def load_data(small_sample_path,large_sample_path):
    df_small = pd.read_csv(small_sample_path)
    df_large = pd.read_csv(large_sample_path)

    data = pd.concat((df_large, df_small), axis=0)
    data['eco_cause'] = data['eco_id'].apply(str) + '_' + data['cause']
    return data

def initial_parameters():
    mu_a_fm = pm_pool_initialize("mu_a_fm", 1.0)
    mu_b_fm = pm.Normal("mu_b_fm", mu=0.0, sigma=100)

    mu_a_fa = pm_pool_initialize("mu_a_fa", 1.0)
    mu_b_fa = pm.Normal("mu_b_fa", mu=0.0, sigma=100)

    mu_a_popu = pm_pool_initialize("mu_a_popu", 1.0)
    mu_b_popu = pm.Normal("mu_b_popu", mu=0.0, sigma=100)

    a_fm = pm_initialize("a_fm", beta=mu_a_fm, dims="eco")
    b_fm = pm.Normal("b_fm", mu=mu_b_fm, sigma=100, dims="ecocause")
    a_fa = pm_initialize("a_fa", beta=mu_a_fa, dims="eco")
    b_fa = pm.Normal("b_fa", mu=mu_b_fa, sigma=100, dims="ecocause")
    a_popu = pm_initialize("a_popu", beta=mu_a_popu, dims="eco")
    b_popu = pm.Normal("b_popu", mu=mu_b_popu, sigma=100, dims="ecocause")

    return a_fm,b_fm,a_fa,b_fa,a_popu,b_popu

def Bayesian_inference(model_fpath,df_small_sample_path,df_large_sample_path):
    data=load_data(df_small_sample_path,df_large_sample_path)
    ecocause_idxs, ecocauses = pd.factorize(data.eco_cause)
    eco_idxs, ecos = pd.factorize(data.eco_id)
    coords = {
        "ecocause": ecocauses,
        "eco": ecos,
        "obs_id": np.arange(len(ecocause_idxs)),
    }

    print(ecocauses, len(ecocauses))
    data["fire_class"] = data["fire_class"].astype(theano.config.floatX)

    with pm.Model(coords=coords) as fire_model:
        ecocause_idx = pm.Data("ecocause_id", ecocause_idxs, dims="obs_id")
        eco_idx = pm.Data("eco_id", eco_idxs, dims="obs_id")
        a_fm,b_fm,a_fa,b_fa,a_popu,b_popu=initial_parameters()
        fm = pm.Data("vpd", data.vpd.values, dims="obs_id")
        fa = pm.Data("NPP", data.NPP.values, dims="obs_id")
        popu = pm.Data("Popu", data.Popu.values, dims="obs_id")
        likelihood_fuel_mositure = pm.invlogit(a_fm[eco_idx] * (fm - b_fm[ecocause_idx]))
        likelihood_fuel_availability = pm.invlogit(a_fa[eco_idx] * (fa - b_fa[ecocause_idx]))
        likelihood_popu_suppression = pm.invlogit(-a_popu[eco_idx] * popu + b_popu[ecocause_idx])
        likelihood = likelihood_fuel_mositure * likelihood_fuel_availability * likelihood_popu_suppression
        # Bernoulli random vector with probability of success
        # given by sigmoid function and actual data as observed
        pm.Bernoulli(name='y', p=likelihood, observed=data.fire_class, dims="obs_id")

    with fire_model:
        # do the sampling
        trace = pm.sample(tune=2000,
                          draws=1000, target_accept=0.9)

    with open(model_fpath, 'wb') as buff:
        pickle.dump({'model': fire_model, 'trace': trace}, buff)
        print('save trace successfully!')

if __name__ == '__main__':
    model_fpath = f'model_para/eco_thresholds.para'
    small_sample_path = f'data_for_small_fires.csv'
    large_sample_path = f'data_for_large_fires.csv'
    Bayesian_inference(model_fpath,small_sample_path,large_sample_path)