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
    e_pair = np.array([-6.704189,-6.887298,-6.393479,-1.874734,
                       -6.546295,-7.470841,-1.317225,-6.928318,
                       -6.493723,-1.266180,-6.809987,-2.819615,
                       -1.794874,-6.914457,-2.806629,-6.363942]) # UA, UC, UG, UU
    e_pair = prob.logpNorm(e_pair);
    pe_pair = jnp.exp(e_pair)

    # unpaired emission probabilities 4x1 matrix
    e_single = np.array([-1.012587,-1.753798,-1.518921,-1.406257]) # A, C, G, U
    e_single = prob.logpNorm(e_single);
    pe_single = jnp.exp(e_single)

    # transition probabilities (t1, t2, t3)
    log_t = np.array([-0.726176,-1.365860,-1.341766]) # aS | aSa'S | e
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



