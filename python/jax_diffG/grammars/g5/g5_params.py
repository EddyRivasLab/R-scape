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



