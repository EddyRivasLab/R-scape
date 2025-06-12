import math
import sys
import numpy as np
import jax
import jax.numpy as jnp
import jax.nn as jnn
import jax.scipy as jsp
import functools

import lib.probability as prob


def G5_param_tornado(verbose):
    # paired emission probabilities 16x1 matrix
    e_pair = np.array([-6.625427,-6.793430,-6.335155,-1.874087,
                       -6.478658,-7.308380,-1.316858,-6.830704, 
                       -6.429442,-1.265831,-6.722807,-2.817939, 
                       -1.794277,-6.818125,-2.804975,-6.307268]) # UA, UC, UG, UU
    e_pair = prob.logpNorm(e_pair);
    pe_pair = jnp.exp(e_pair)

    # unpaired emission probabilities 4x1 matrix
    e_single = np.array([-1.012312,-1.753221,-1.518464,-1.405849]) # A, C, G, U
    e_single = prob.logpNorm(e_single);
    pe_single = jnp.exp(e_single)

    # transition probabilities (t1, t2, t3)
    log_t = np.array([-0.725970,-1.365468,-1.341384]) # aS | aSa'S | e
    log_t = prob.logpNorm(log_t);    
    t = jnp.exp(log_t)

    if verbose:
        print("G5 param tornado")
        print("     transitions S", t, log_t)
        print("     emissions single", pe_single, e_single)
        print("     emissions paired", pe_pair, e_pair)

    return log_t, t, e_single, pe_single, e_pair, pe_pair

def G5_param_uniform(K, verbose):
    # paired emission probabilities 4x4 matrix
    log_val = -2.0 * np.log(K)
    e_pair = np.array([log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val]) # UA, UC, UG, UU
    e_pair = prob.logpNorm(e_pair);
    pe_pair = np.exp(e_pair)

    # unpaired emission probabilities 4x1 matrix
    log_val = -np.log(K)
    e_single = np.array([log_val,log_val,log_val,log_val]) # A, C, G, U
    e_single = prob.logpNorm(e_single);
    pe_single = np.exp(e_single)

    # transition probabilities (t1, t2, t3)
    log_val = -np.log(3.0)
    log_t = np.array([log_val,log_val,log_val]) # aS | aSa'S | e
    log_t = prob.logpNorm(log_t);

    t = np.exp(log_t)

    if verbose:
        print("G5 param uniform")
        print("     transitions S", t)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)

    return log_t, t, e_single, pe_single, e_pair, pe_pair

def G5_param_naive(verbose):
    # paired emission probabilities 6 WC cases only
    pe_pair = np.array([0.0,0.0,0.0,1.0,
                       0.0,0.0,1.0,0.0,
                       0.0,1.0,0.0,1.0,
                       1.0,0.0,1.0,0.0,]) # UA, UC, UG, UU
    pe_pair = prob.pNorm(pe_pair);
    e_pair = np.log(pe_pair)

    # unpaired emission probabilities 4x1 matrix
    pe_single = np.array([1.0,1.0,1.0,1.0,]) # A, C, G, U
    pe_single = prob.pNorm(pe_single);
    e_single = np.log(pe_single)

    # transition probabilities (t1, t2, t3)
    t = np.array([1.0,1.0,1.0]) # aS | aSa'S | e
    t = prob.pNorm(t);
    log_t = np.log(t)

    if verbose:
        print("G5 param naive")
        print("     transitions S", t)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)

    return log_t, t, e_single, pe_single, e_pair, pe_pair



def G5_normalize_params(uparams, scaled):
    if scaled:
        return {
            "t":         prob.pNorm(uparams['t']),
            "t2":        prob.pNorm(uparams['t2']),
            "pe_single": prob.pNorm(uparams['pe_single']),
            "pe_pair":   prob.pNorm(uparams['pe_pair'])
        }
    else:
        return {
            "log_t":    prob.logpNorm(uparams['log_t']),
            "e_single": prob.logpNorm(uparams['e_single']),
            "e_pair":   prob.logpNorm(uparams['e_pair'])
        }
        

def G5_write_paramfile(rundir, epoch, params):
    param_file = rundir / f"param_i{epoch}.param"
    
    with open(param_file, "a") as f:
        f.write(f"3\n")
        f.write(f"0 3 {params['log_t'][0]} {params['log_t'][1]} {params['log_t'][2]}\n")
        f.write(f"2\n")
        f.write(f"0 e1_2_0_0 2 0 1 (WW_C 0 1)\n")
        f.write(f"{params['e_pair'][0]}  {params['e_pair'][1]}  {params['e_pair'][2]}  {params['e_pair'][3]}\n")
        f.write(f"{params['e_pair'][4]}  {params['e_pair'][5]}  {params['e_pair'][6]}  {params['e_pair'][7]}\n")
        f.write(f"{params['e_pair'][8]}  {params['e_pair'][9]}  {params['e_pair'][10]} {params['e_pair'][11]}\n")
        f.write(f"{params['e_pair'][12]} {params['e_pair'][13]} {params['e_pair'][14]} {params['e_pair'][15]}\n")
        f.write(f"1 e1_1_0_0 1 0 0\n")
        f.write(f"{params['e_single'][0]} {params['e_single'][1]} {params['e_single'][2]} {params['e_single'][3]}\n")
        f.write(f"0\n")
    return param_file
