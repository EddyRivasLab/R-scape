import math
import sys
import numpy as np
import jax
import jax.numpy as jnp
import jax.nn as jnn
import jax.scipy as jsp
import functools

import lib.probability as prob

def G6X_param_tornado(verbose):

    # paired emission probabilities 4x4 matrix
    e_pair = np.array([-6.716494,-6.893637,-6.396199,-1.875715,  #AA AC AG AU
                       -6.556362,-7.468432,-1.316826,-6.944485,  #CA CC CG CU
                       -6.497039,-1.265024,-6.815633,-2.819071,  #GA GC GG GU
                       -1.796519,-6.921051,-2.806054,-6.366497]) #UA UC UG UU
    e_pair = prob.logpNorm(e_pair);
    pe_pair = jnp.exp(e_pair)

    # unpaired emission probabilities 4x1 matrix
    e_single = np.array([-1.013576,-1.753141,-1.518408,-1.405715]) # A, C, G, U
    e_single = prob.logpNorm(e_single);
    pe_single = jnp.exp(e_single)

    # transition probabilities (t0=tS[2],t1=tL[3],t2=tF[3])
    log_t0 = np.array([-0.293521,-1.368194]) # S -> LS | e
    log_t1 = np.array([-1.396876,-9.210340,-0.283914]) # L -> aFa' | aa' | a
    log_t2 = np.array([-0.984961,-9.210340,-0.467213]) # F -> aFa' | aa' | LS

    log_t0 = prob.logpNorm(log_t0);   
    log_t1 = prob.logpNorm(log_t1);   
    log_t2 = prob.logpNorm(log_t2);   
    t0 = np.exp(log_t0)
    t1 = np.exp(log_t1)
    t2 = np.exp(log_t2)

    if verbose:
        print("G6X param tornado")
        print("     transitions S t0", t0)
        print("     transitions L t1", t1)
        print("     transitions F t2", t2)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)

    return log_t0, t0, log_t1, t1, log_t2, t2, e_single, pe_single, e_pair, pe_pair

def G6X_param_uniform(K, verbose):
    # paired emission probabilities 4x4 matrix
    log_val = -2.0 * np.log(K)
    e_pair = np.array([log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val]) # UA, UC, UG, UU
    e_pair = prob.logpNorm(e_pair);
    pe_pair = jnp.exp(e_pair)

    # unpaired emission probabilities 4x1 matrix
    log_val = -np.log(K)
    e_single = np.array([log_val,log_val,log_val,log_val]) # A, C, G, U
    e_single = prob.logpNorm(e_single);
    pe_single = jnp.exp(e_single)

    # transition probabilities S
    log_val = -np.log(2.0)
    log_t0 = np.array([log_val,log_val])         # S -> LS | end
    log_t0 = prob.logpNorm(log_t0);   
    t0 = np.exp(log_t0)
    # transition probabilities L
    log_val = -np.log(3.0)
    log_t1 = np.array([log_val,log_val,log_val]) # L -> aFa' | aa' | a
    log_t1 = prob.logpNorm(log_t1);   
    t1 = np.exp(log_t1)
    # transition probabilities F
    log_val = -np.log(3.0)
    log_t2 = np.array([log_val,log_val,log_val]) # F -> aFa' | aa' | LS
    log_t2 = prob.logpNorm(log_t2);   
    t2 = np.exp(log_t2)

    if verbose:
        print("G6X param uniform")
        print("     transitions S t0", t0)
        print("     transitions L t1", t1)
        print("     transitions F t2", t2)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)

    return log_t0, t0, log_t1, t1, log_t2, t2, e_single, pe_single, e_pair, pe_pair



