import os
import math
import sys
import numpy as np
import random
from pathlib import Path

import jax
import jax.numpy as jnp
import jax.nn as jnn
import jax.scipy as jsp
import functools
from copy import deepcopy
import pprint
import re
import matplotlib.pyplot as plt

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

def G5_param_random(K, verbose):
    random.seed()
    
    # paired emission probabilities 4x4 matrix
    pe_pair = []
    for i in range(K*K):
        pe_pair.append(random.random())
    pe_pair = np.array(pe_pair)
    pe_pair = prob.pNorm(pe_pair);
    e_pair = jnp.log(pe_pair)

    # unpaired emission probabilities 4x1 matrix
    pe_single = []
    for i in range(K):
        pe_single.append(random.random())
    pe_single = np.array(pe_single)
    pe_single = prob.pNorm(pe_single);
    e_single = jnp.log(pe_single)

    # transition probabilities S
    t = np.array([random.random(),random.random(),random.random()]) # S -> a S | a S b S | e
    t = prob.pNorm(t);   
    log_t = np.log(t)

    if verbose:
        print("G5 param uniform")
        print("     transitions S t", t)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)

    return log_t, t, e_single, pe_single, e_pair, pe_pair



def G5_normalize_params(uparams, scaled):
    if scaled:
        return {
            "t":         prob.pNorm(uparams['t']),
            "pe_single": prob.pNorm(uparams['pe_single']),
            "pe_pair":   prob.pNorm(uparams['pe_pair'])
        }
    else:
        return {
            "log_t":    prob.logpNorm(uparams['log_t']),
            "e_single": prob.logpNorm(uparams['e_single']),
            "e_pair":   prob.logpNorm(uparams['e_pair'])
        }
        

def G5_plot_params(outdir, epoch, params, params_ref, param_ref_name):
    default_color = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    if (os.path.exists(outdir)):
        nts = ['A', 'C', 'G', 'U']

        pair = [['AA', 'AC', 'AG', 'AU'],
                ['CA', 'CC', 'CG', 'CU'],
                ['GA', 'GC', 'GG', 'GU'],
                ['UA', 'UC', 'UG', 'UU']]
        
        S_rules = ['a S', 'a S b', 'e']
        
        t     = np.exp(params['log_t'])
        t_ref = np.exp(params_ref['log_t'])
        t0max = np.max([t, t_ref])
              
        e_single     = np.exp(params['e_single'])
        e_single_ref = np.exp(params_ref['e_single'])
        singlemax    = np.max([e_single, e_single_ref])
        
        e_pair     = np.exp(params['e_pair'])
        e_pair_ref = np.exp(params_ref['e_pair'])
        pairmax    = np.max([e_pair, e_pair_ref])
        
        fig = plt.figure(constrained_layout=True)
        gs  = fig.add_gridspec(2, 4)
        
        ax_pairA = fig.add_subplot(gs[0, 0])
        ax_pairC = fig.add_subplot(gs[0, 1])
        ax_pairG = fig.add_subplot(gs[0, 2])
        ax_pairU = fig.add_subplot(gs[0, 3])
        ax_pairC.set_title('Pair Probabilities P(ab) [sum_{ab} P(ab) = 1] a,b = {A,C,G,U}')
        
        ax_single = fig.add_subplot(gs[1, 0])
        ax_single.set_title('Unpaired P(a)')
        
        ax_t = fig.add_subplot(gs[1, 1])
        ax_t.set_title('S')
 
        ax_pairA.plot(pair[0], e_pair_ref[0:4], color=default_color[1])
        ax_pairA.plot(pair[0], e_pair[0:4],     color=default_color[0])
        #ax_pairA.set_xlabel(pair[0])
        ax_pairA.set_ylim(0,pairmax)
        
        ax_pairC.plot(pair[1], e_pair_ref[4:8], color=default_color[1])
        ax_pairC.plot(pair[1], e_pair[4:8],     color=default_color[0])
        #ax_pairC.set_xlabel(pair[1])
        ax_pairC.set_ylim(0,pairmax)
        ax_pairC.set_yticklabels([]) 
        
        ax_pairG.plot(pair[2], e_pair_ref[8:12], color=default_color[1])
        ax_pairG.plot(pair[2], e_pair[8:12],     color=default_color[0])
        #ax_pairG.set_xlabel(pair[2])
        ax_pairG.set_ylim(0,pairmax)
        ax_pairG.set_yticklabels([]) 
      
        ax_pairU.plot(pair[3], e_pair_ref[12:16], color=default_color[1])
        ax_pairU.plot(pair[3], e_pair[12:16],     color=default_color[0])
        #ax_pairU.set_xlabel(pair[3])
        ax_pairU.set_ylim(0,pairmax)
        ax_pairU.set_yticklabels([])
        
        ax_single.plot(nts, e_single_ref, color=default_color[1])
        ax_single.plot(nts, e_single,     color=default_color[0])
        ax_single.set_ylim(0,singlemax)
        
        ax_t.plot(S_rules, t_ref, color=default_color[1])
        ax_t.plot(S_rules, t,     color=default_color[0])
        ax_t.set_ylim(0,1)
        
        fig.suptitle(f'G5 grammar {str(Path(param_ref_name).stem)}')
        plt.savefig(outdir / f"g5_params_i{epoch}.pdf",  format="pdf")
        plt.clf()
        plt.close()

def G5_read_paramfile(param_file, scaled):
    
    #1
    #0  3  -1.005410 -1.155520 -1.140023 
    #2
    #0 e1_2_0_0 2 0 1 (WW_C 0 1)
    #-5.208315372467041  -6.664770126342773  -5.874554634094238  -1.5640732049942017
    #-5.951663017272949  -5.347833633422852  -2.1773369312286377  -6.041226387023926
    #-5.891550064086914  -1.4333771467208862  -5.636483669281006 -3.331711530685425
    #-1.4301977157592773 -4.504973411560059 -3.259310483932495 -2.421502113342285
    #1 e1_1_0_0 1 0 0
    #-0.8405029773712158 -1.5922355651855469 -1.5958921909332275 -1.8182547092437744
    #0

    comment_pattern = r'\#'
    
    t_pattern = r'0\s+3\s+(\S+)\s+(\S+)\s+(\S+)'

    single_pattern = r'e1_1_0_0'
    pair_pattern   = r'e1_2_0_0'
    
    emit_pattern = r'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)'

    n = 0
    is_single = False
    is_pair   = False       

    e_pair = []
    
    fp = open(param_file, "r")
    lines = fp.readlines()
    for line in lines:
        comment_match = re.search(comment_pattern, line)
        if comment_match: continue
        
        t_match = re.search(t_pattern, line)
       
        single_match = re.search(single_pattern, line)
        pair_match   = re.search(pair_pattern, line)
        
        emit_match = re.search(emit_pattern, line)

        if single_match:
            is_single = True
            is_pair   = False
            continue
        if pair_match:
            is_single = False
            is_pair   = True
            continue
            
        if t_match: t = np.array([float(t_match.group(1)), float(t_match.group(2)), float(t_match.group(3))])
  
        if emit_match and is_single and not is_pair:
            e_single = np.array([float(emit_match.group(1)),float(emit_match.group(2)),float(emit_match.group(3)),float(emit_match.group(4))])
            
        if emit_match and is_pair and not is_single:
            e_pair.append([float(emit_match.group(1)),float(emit_match.group(2)),float(emit_match.group(3)),float(emit_match.group(4))])

    e_pair = np.array(e_pair).flatten()
    if scaled: params = {"t":     deepcopy(t), "pe_single": deepcopy(e_single), "pe_pair": deepcopy(e_pair)  }
    else:      params = {"log_t": deepcopy(t), "e_single":  deepcopy(e_single),  "e_pair": deepcopy(e_pair)  }

    print(f"{pprint.pformat(params)}")

    return G5_normalize_params(params, scaled)

def G5_write_paramfile(rundir, epoch, params):
    param_file = rundir / f"param_i{epoch}.param"

    # e1_1_0_0 appears first
    with open(param_file, "a") as f:
        f.write(f"1\n")
        f.write(f"0 3 {params['log_t'][0]} {params['log_t'][1]} {params['log_t'][2]}\n")
        f.write(f"2\n")
        f.write(f"0 e1_1_0_0 1 0 0\n")
        f.write(f"{params['e_single'][0]} {params['e_single'][1]} {params['e_single'][2]} {params['e_single'][3]}\n")
        f.write(f"1 e1_2_0_0 2 0 1 (WW_C 0 1)\n")
        f.write(f"{params['e_pair'][0]}  {params['e_pair'][1]}  {params['e_pair'][2]}  {params['e_pair'][3]}\n")
        f.write(f"{params['e_pair'][4]}  {params['e_pair'][5]}  {params['e_pair'][6]}  {params['e_pair'][7]}\n")
        f.write(f"{params['e_pair'][8]}  {params['e_pair'][9]}  {params['e_pair'][10]} {params['e_pair'][11]}\n")
        f.write(f"{params['e_pair'][12]} {params['e_pair'][13]} {params['e_pair'][14]} {params['e_pair'][15]}\n")
        f.write(f"0\n")
    return str(param_file)
