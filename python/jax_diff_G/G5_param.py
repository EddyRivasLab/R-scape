import math
import sys
import numpy as np
import jax
import jax.numpy as jnp
import jax.nn as jnn
import jax.scipy as jsp
import functools


def G5_param_tornado(verbose):
    # paired emission probabilities 16x1 matrix
    e_pair = np.array([-6.704189,-6.887298,-6.393479,-1.874734,
                       -6.546295,-7.470841,-1.317225,-6.928318,
                       -6.493723,-1.266180,-6.809987,-2.819615,
                       -1.794874,-6.914457,-2.806629,-6.363942]) # UA, UC, UG, UU
    pe_pair = np.exp(e_pair)

    # unpaired emission probabilities 4x1 matrix
    e_single = np.array([-1.012587,-1.753798,-1.518921,-1.406257]) # A, C, G, U
    pe_single = np.exp(e_single)

    # transition probabilities (t1, t2, t3)
    log_t = np.array([-0.726176,-1.365860,-1.341766]) # aS | aSa'S | e
    t = np.exp(log_t)

    if verbose:
        print("G5 param tornado")
        print("     transitions S", t)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)

    return log_t, t, e_single, pe_single, e_pair, pe_pair

def G5_param_uniform(K, verbose):
    # paired emission probabilities 4x4 matrix
    log_val = -2.0 * np.log(K)
    e_pair = np.array([log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val,
                       log_val,log_val,log_val,log_val]) # UA, UC, UG, UU
    pe_pair = np.exp(e_pair)

    # unpaired emission probabilities 4x1 matrix
    log_val = -np.log(K)
    e_single = np.array([log_val,log_val,log_val,log_val]) # A, C, G, U
    pe_single = np.exp(e_single)

    # transition probabilities (t1, t2, t3)
    log_val = -np.log(3.0)
    log_t = np.array([log_val,log_val,log_val]) # aS | aSa'S | e
    t = np.exp(log_t)

    if verbose:
        print("G5 param uniform")
        print("     transitions S", t)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)

    return log_t, t, e_single, pe_single, e_pair, pe_pair



