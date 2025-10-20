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

def G6S_param_tornado(verbose):

    # stacked emission probabilities 16x16 matrix
    e_stck = np.array([
        [-3.056579,-4.657551,-2.720715,-2.028327,-3.969687,-4.657551,-1.834305,-5.340216,-3.565989,-1.623116,-9.210340,-3.565989,-1.576619,-9.210340,-2.587373,-3.969687], # e1_2_2_0
        [-2.909155,-4.289962,-2.909155,-2.909155,-4.289962,-3.886932,-1.811766,-9.210340,-4.289962,-0.800555,-4.975840,-3.886932,-2.686378,-9.210340,-4.289962,-2.909155], # e1_2_2_1
        [-3.643576,-3.929985,-3.643576,-2.634411,-9.210340,-9.210340,-1.260587,-3.929985,-9.210340,-1.814211,-9.210340,-1.319754,-2.547516,-9.210340,-3.929985,-3.643576], # e1_2_2_2
        [-7.250949,-7.250949,-7.250949,-1.577124,-7.563546,-9.210340,-1.513192,-7.056399,-7.250949,-1.251385,-7.056399,-3.091547,-1.731993,-7.364546,-2.822618,-6.821176], # e1_2_2_3
        [-9.210340,-3.946057,-3.255500,-2.850900,-3.032875,-9.210340,-1.552876,-3.946057,-9.210340,-2.078642,-4.634045,-3.388661,-1.385894,-4.232017,-2.245551,-2.696996], # e1_2_2_4
        [-9.210340,-3.593669,-3.593669,-2.902342,-4.283186,-4.283186,-1.804945,-9.210340,-4.283186,-2.092424,-1.987145,-9.210340,-3.189417,-3.593669,-2.210106,-1.399682], # e1_2_2_5
        [-7.287139,-7.287139,-7.093659,-1.827732,-6.447138,-8.215549,-1.130815,-6.510543,-7.218401,-1.561791,-7.400009,-3.436954,-1.589205,-7.783657,-2.750302,-7.093659], # e1_2_2_6
        [-2.469738,-9.210340,-2.064667,-2.651823,-4.255605,-4.941727,-1.697185,-2.389786,-9.210340,-1.489648,-9.210340,-3.565989,-2.469738,-9.210340,-9.210340,-2.469738], # e1_2_2_7
        [-3.646816,-4.737787,-2.955588,-2.350326,-2.955588,-3.493213,-1.886410,-4.050365,-5.419580,-1.570736,-4.737787,-3.360093,-1.529933,-3.646816,-2.801712,-4.737787], # e1_2_2_8
        [-7.624153,-8.060281,-6.945091,-2.136234,-6.990936,-8.851494,-1.343546,-7.932172,-6.742927,-1.075761,-7.289362,-2.953436,-1.963297,-7.389003,-2.492520,-7.389003], # e1_2_2_9
        [-9.210340,-3.094498,-9.210340,-1.917373,-3.094498,-2.976960,-1.917373,-4.474192,-2.341896,-2.179533,-3.227714,-4.474192,-1.709861,-3.563178,-2.976960,-3.381445], # e1_2_2_10
        [-8.158023,-6.876121,-8.158023,-1.913489,-9.210340,-8.158023,-1.209059,-8.158023,-6.876121,-1.206563,-7.192112,-2.628947,-1.894716,-8.158023,-3.553528,-7.075332], # e1_2_2_11
        [-6.697439,-6.819461,-6.958468,-1.702906,-6.538466,-7.696541,-1.375779,-7.425010,-6.697439,-1.354681,-6.641593,-3.305696,-1.588497,-7.425010,-2.843509,-7.211684], # e1_2_2_12
        [-3.971242,-9.210340,-4.374058,-2.078642,-2.876173,-4.374058,-1.896453,-9.210340,-3.684887,-1.183027,-3.126981,-4.374058,-1.673443,-4.374058,-9.210340,-3.280751], # e1_2_2_13
        [-6.874182,-7.475025,-8.005760,-1.948015,-6.670568,-8.436684,-1.291512,-7.475025,-5.692991,-1.628859,-7.475025,-2.021943,-1.837985,-6.066796,-2.465219,-6.670568], # e1_2_2_14
        [-4.948646,-4.948646,-4.262573,-2.189363,-9.210340,-5.627793,-1.595056,-9.210340,-3.350551,-1.144818,-4.262573,-3.859477,-1.966398,-3.859477,-3.350551,-2.658831]  # e1_2_2_15
        ]}
    K = e_stck.shape[0]
    e_stck = jnp.array([prob.logpNorm(e_stck[ab])  for ab in range(K)])
    pe_stck = np.exp(e_stck)

    # paired emission probabilities 4x4 matrix
    e_pair = np.array([-6.344996,-6.476250,-5.859319,-1.953408,  #AA AC AG AU
                       -6.441786,-6.821846,-1.348678,-6.908021,  #CA CC CG CU
                       -6.219902,-1.172563,-6.574419,-2.397389,  #GA GC GG GU
                       -1.976776,-6.430557,-3.129559,-5.944866]) #UA UC UG UU
    e_pair = prob.logpNorm(e_pair)
    pe_pair = jnp.exp(e_pair)

    # unpaired emission probabilities 4x1 matrix
    e_single = np.array([-1.013294,-1.752666,-1.518000,-1.405201]) # A, C, G, U
    e_single = prob.logpNorm(e_single)
    pe_single = jnp.exp(e_single)

    # transition probabilities (t0=tS[2],t1=tL[3],t2=tF[3])
    log_t0 = np.array([-0.157572,-1.922886]) # S -> LS   | L
    log_t1 = np.array([-2.139511,-0.124784]) # L -> aFa' | a
    log_t2 = np.array([-0.292894,-1.369245]) # F -> aFa' | LS

    log_t0 = prob.logpNorm(log_t0)
    log_t1 = prob.logpNorm(log_t1)  
    log_t2 = prob.logpNorm(log_t2)   
    t0 = np.exp(log_t0)
    t1 = np.exp(log_t1)
    t2 = np.exp(log_t2)

    if verbose:
        print("G6S param tornado")
        print("     transitions S t0", t0)
        print("     transitions L t1", t1)
        print("     transitions F t2", t2)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)
        print("     emissions stacked", pe_stck)
       
    return log_t0, t0, log_t1, t1, log_t2, t2, e_single, pe_single, e_pair, pe_pair, e_stck, pe_stck

def G6S_param_uniform(K, verbose):
    # stacked emission probabilities 4x4x4x4 matrix
    log_val = -2.0 * np.log(K)
    e_stck = np.array([
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val], 
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val],
        [log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val,log_val]])

    K = e_stck.shape[0]
    e_stck = jnp.array([prob.logpNorm(e_stck[ab])  for ab in range(K)])
    pe_stck = np.exp(e_stck)

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
    log_t0 = np.array([log_val,log_val]) # S -> LS | L
    log_t0 = prob.logpNorm(log_t0);   
    t0 = np.exp(log_t0)
    # transition probabilities L
    log_val = -np.log(2.0)
    log_t1 = np.array([log_val,log_val]) # L -> aFa' | a
    log_t1 = prob.logpNorm(log_t1);   
    t1 = np.exp(log_t1)
    # transition probabilities F
    log_val = -np.log(2.0)
    log_t2 = np.array([log_val,log_val]) # F -> aFa' | LS
    log_t2 = prob.logpNorm(log_t2);   
    t2 = np.exp(log_t2)

    if verbose:
        print("G6S param uniform")
        print("     transitions S t0", t0)
        print("     transitions L t1", t1)
        print("     transitions F t2", t2)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)
        print("     emissions stacked", pe_stck)

    return log_t0, t0, log_t1, t1, log_t2, t2, e_single, pe_single, e_pair, pe_pair, e_stck, pe_stck

def G6S_param_random(K, verbose):
    random.seed()

    # stacked emission probabilities 4x4x4x4 matrix
    pe_stck = []
    for i in range(K*K):
        this_pair = []
        
        for j in range(K*K):
            this_pair.append(random.random())
        this_pair = np.array(this_pair)
        this_pair = prob.pNorm(this_pair);
        
        pe_stck.append(this_pair)
    e_stck = jnp.log(pe_stck)

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
    t0 = np.array([random.random(),random.random()]) # S -> LS | L
    t0 = prob.pNorm(t0);   
    log_t0 = np.log(t0)
    # transition probabilities L
    t1 = np.array([random.random(),random.random()]) # L -> aFa' | a
    t1 = prob.pNorm(t1);   
    log_t1 = np.log(t1)
    # transition probabilities F
    t2 = np.array([random.random(),random.random()]) # F -> aFa' | LS
    t2 = prob.pNorm(t2);   
    log_t2 = np.log(t2)

    if verbose:
        print("G6S param uniform")
        print("     transitions S t0", t0)
        print("     transitions L t1", t1)
        print("     transitions F t2", t2)
        print("     emissions single", pe_single)
        print("     emissions paired", pe_pair)
        print("     emissions stacked", pe_stck)

    return log_t0, t0, log_t1, t1, log_t2, t2, e_single, pe_single, e_pair, pe_pair, e_stck, pe_stck


   
    e_stck = jnp.array([prob.logpNorm(e_stck[ab])  for ab in range(K)])
 
def G6S_normalize_params(uparams, scaled):
    if scaled:
         K = uparams['pe_stck'].shape[0]
         
        return {
            "t0":         prob.pNorm(uparams['t0']),
            "t1":         prob.pNorm(uparams['t1']),
            "t2":         prob.pNorm(uparams['t2']),
            "pe_single":  prob.pNorm(uparams['pe_single']),
            "pe_pair":    prob.pNorm(uparams['pe_pair']),
            "pe_stck":    jnp.array([prob.pNorm(uparams['pe_stck'][ab])  for ab in range(K)])
        }
    else:
        K = uparams['e_stck'].shape[0]
        
        return {
            "log_t0":    prob.logpNorm(uparams['log_t0']),
            "log_t1":    prob.logpNorm(uparams['log_t1']),
            "log_t2":    prob.logpNorm(uparams['log_t2']),
            "e_single":  prob.logpNorm(uparams['e_single']),
            "e_pair":    prob.logpNorm(uparams['e_pair']),
            "e_stck":    jnp.array([prob.logpNorm(uparams['e_stck'][ab]) for ab in range(K)])
        }
        

def G6S_plot_params(outdir, epoch, params, params_ref, param_ref_name, pair_ymax):
    default_color = plt.rcParams['axes.prop_cycle'].by_key()['color']

    if (os.path.exists(outdir)):
        nts = ['A', 'C', 'G', 'U']

        pair = [['AA', 'AC', 'AG', 'AU'],
                ['CA', 'CC', 'CG', 'CU'],
                ['GA', 'GC', 'GG', 'GU'],
                ['UA', 'UC', 'UG', 'UU']]
        
        S_rules = [' LS ',    ' L ']
        L_rules = [' a F b ', ' a ']
        F_rules = [' a F b ', ' LS ']
        
        t0     = np.exp(params['log_t0'])
        t0_ref = np.exp(params_ref['log_t0'])
        t0max  = np.max([t0, t0_ref])
        
        t1     = np.exp(params['log_t1'])
        t1_ref = np.exp(params_ref['log_t1'])
        t1max  = np.max([t1, t1_ref])
        
        t2     = np.exp(params['log_t2'])
        t2_ref = np.exp(params_ref['log_t2'])
        t2max  = np.max([t2, t2_ref])
        
        e_single     = np.exp(params['e_single'])
        e_single_ref = np.exp(params_ref['e_single'])
        singlemax    = np.max([e_single, e_single_ref])
        
        e_pair     = np.exp(params['e_pair'])
        e_pair_ref = np.exp(params_ref['e_pair'])
        pairmax    = np.max([e_pair, e_pair_ref])
        if (pair_ymax > pairmax): pairmax = pair_ymax
        
        e_stck     = jnp.array(np.exp(params['e_stck'][ab]) for ab in range(K)])
        e_stck_ref = jnp.array(np.exp(params_ref['e_stck'][ab]) for ab in range(K)])
        pairmax    = np.max([e_stck[ab], e_stck_ref[ab]] for ab in range(K)]) 
        if (pair_ymax > pairmax): pairmax = pair_ymax
        
        fig = plt.figure(constrained_layout=True)
        gs  = fig.add_gridspec(2, 4)
        
        ax_pairA = fig.add_subplot(gs[0, 0])
        ax_pairC = fig.add_subplot(gs[0, 1])
        ax_pairG = fig.add_subplot(gs[0, 2])
        ax_pairU = fig.add_subplot(gs[0, 3])
        ax_pairC.set_title('Pair Probabilities P(ab) [sum_{ab} P(ab) = 1] a,b = {A,C,G,U}')
        
        ax_single = fig.add_subplot(gs[1, 0])
        ax_single.set_title('Unpaired P(a)')
        
        ax_t0 = fig.add_subplot(gs[1, 1])
        ax_t0.set_title('S')
        ax_t1 = fig.add_subplot(gs[1, 2])
        ax_t1.set_title('L')
        ax_t2 = fig.add_subplot(gs[1, 3])
        ax_t2.set_title('F')

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
        
        ax_t0.plot(S_rules, t0_ref, color=default_color[1])
        ax_t0.plot(S_rules, t0,     color=default_color[0])
        ax_t0.set_ylim(0,1)
        
        ax_t1.plot(L_rules, t1_ref, color=default_color[1])
        ax_t1.plot(L_rules, t1,     color=default_color[0])
        ax_t1.set_ylim(0,1)
        
        ax_t2.plot(F_rules, t2_ref, color=default_color[1])
        ax_t2.plot(F_rules, t2,     color=default_color[0])
        ax_t2.set_ylim(0,1)
        
        fig.suptitle(f'G6S grammar {str(Path(param_ref_name).stem)}')
        plt.savefig(outdir / f"g6s_params_i{epoch}.pdf",  format="pdf")
        plt.clf()
        plt.close()
        
def G6S_plot_params_alone(outdir, epoch, params):
    default_color = plt.rcParams['axes.prop_cycle'].by_key()['color']

    if (os.path.exists(outdir)):
        nts = ['A', 'C', 'G', 'U']

        pair = [['AA', 'AC', 'AG', 'AU'],
                ['CA', 'CC', 'CG', 'CU'],
                ['GA', 'GC', 'GG', 'GU'],
                ['UA', 'UC', 'UG', 'UU']]
        
        S_rules = [' LS ',    ' L ']
        L_rules = [' a F b ', ' a ']
        F_rules = [' a F b ', ' LS ']
        
        t0 = np.exp(params['log_t0'])        
        t1 = np.exp(params['log_t1'])       
        t2 = np.exp(params['log_t2'])
        
        e_single = np.exp(params['e_single'])
        e_pair   = np.exp(params['e_pair'])
        pairmax   = np.max([e_pair])
        
        fig = plt.figure(constrained_layout=True)
        gs  = fig.add_gridspec(2, 4)
        
        ax_pairA = fig.add_subplot(gs[0, 0])
        ax_pairC = fig.add_subplot(gs[0, 1])
        ax_pairG = fig.add_subplot(gs[0, 2])
        ax_pairU = fig.add_subplot(gs[0, 3])
        ax_pairC.set_title('Pair Probabilities P(ab) [sum_{ab} P(ab) = 1] a,b = {A,C,G,U}')
        
        ax_single = fig.add_subplot(gs[1, 0])
        ax_single.set_title('Unpaired P(a)')
        
        ax_t0 = fig.add_subplot(gs[1, 1])
        ax_t0.set_title('S')
        ax_t1 = fig.add_subplot(gs[1, 2])
        ax_t1.set_title('L')
        ax_t2 = fig.add_subplot(gs[1, 3])
        ax_t2.set_title('F')

        ax_pairA.plot(pair[0], e_pair[0:4],     color=default_color[0])
        #ax_pairA.set_xlabel(pair[0])
        ax_pairA.set_ylim(0,pairmax)
        
        ax_pairC.plot(pair[1], e_pair[4:8],     color=default_color[0])
        #ax_pairC.set_xlabel(pair[1])
        ax_pairC.set_ylim(0,pairmax)
        ax_pairC.set_yticklabels([]) 
        
        ax_pairG.plot(pair[2], e_pair[8:12],     color=default_color[0])
        #ax_pairG.set_xlabel(pair[2])
        ax_pairG.set_ylim(0,pairmax)
        ax_pairG.set_yticklabels([]) 
      
        ax_pairU.plot(pair[3], e_pair[12:16],     color=default_color[0])
        #ax_pairU.set_xlabel(pair[3])
        ax_pairU.set_ylim(0,pairmax)
        ax_pairU.set_yticklabels([])
        
        ax_single.plot(nts, e_single,     color=default_color[0])
        
        ax_t0.plot(S_rules, t0,     color=default_color[0])
        ax_t0.set_ylim(0,1)
        
        ax_t1.plot(L_rules, t1,     color=default_color[0])
        ax_t1.set_ylim(0,1)
        
        ax_t2.plot(F_rules, t2,     color=default_color[0])
        ax_t2.set_ylim(0,1)
        
        fig.suptitle(f'G6S grammar')
        plt.savefig(outdir / f"g6s_params_i{epoch}.pdf",  format="pdf")
        plt.clf()
        plt.close()

def G6S_read_paramfile(param_file, scaled):
    
    #3
    #0 2 -0.5707096457481384 -0.8326951265335083
    #1 2 -1.2575370073318481 -0.3345690071582794
    #2 2 -0.3449913263320923 -1.231777548789978
    #18
    #  0 e1_1_0_0 1 0 0 
    #-1.013294 -1.752666 -1.518000 -1.405201 
    #1 e1_2_0_0 2 0 1 (WW_C 0 1)
    # -6.288652 -6.412247 -5.824289 -1.952719 
    # -6.379887 -6.732561 -1.348315 -6.811076 
    # -6.170024 -1.172263 -6.504039 -2.396299 
    # -1.976071 -6.369329 -3.127263 -5.906763 
    #2 e1_2_2_0 2 2 1 (WW_C 0 1)
    # -3.052000 -4.606606 -2.719562 -2.031480 
    # -3.947375 -4.606606 -1.838221 -5.235189 
    # -3.553476 -1.627709 -7.314246 -3.553476 
    # -1.581343 -7.314246 -2.587294 -3.947375 
    #3 e1_2_2_1 2 2 1 (WW_C 0 1)
    # -2.905732 -4.245267 -2.905732 -2.905732 
    # -4.245267 -3.860315 -1.817814 -7.043635 
    # -4.245267 -0.809640 -4.879288 -3.860315 
    # -2.685784 -7.043635 -4.245267 -2.905732 
    #4 e1_2_2_2 2 2 1 (WW_C 0 1)
    # -3.632591 -3.913690 -3.632591 -2.633679 
    # -7.553935 -7.553935 -1.264281 -3.913690 
    # -7.553935 -1.816794 -7.553935 -1.323356 
    # -2.547276 -7.553935 -3.913690 -3.632591 
    #5 e1_2_2_3 2 2 1 (WW_C 0 1)
    # -7.114444 -7.114444 -7.114444 -1.576683 
    # -7.381260 -8.498326 -1.512782 -6.942721 
    # -7.114444 -1.251083 -6.942721 -3.089326 
    # -1.731468 -7.212801 -2.820934 -6.730289 
    #6 e1_2_2_4 2 2 1 (WW_C 0 1)
    # -7.297693 -3.924042 -3.248175 -2.848537 
    # -3.028524 -7.297693 -1.557791 -3.924042 
    # -7.297693 -2.081661 -4.583498 -3.379219 
    # -1.391232 -4.200367 -2.247728 -2.696058 
    #7 e1_2_2_5 2 2 1 (WW_C 0 1)
    # -6.475659 -3.563123 -3.563123 -2.897521 
    # -4.203358 -4.203358 -1.817703 -6.475659 
    # -4.203358 -2.102228 -1.998129 -6.475659 
    # -3.175937 -3.563123 -2.218436 -1.415405 
    #8 e1_2_2_6 2 2 1 (WW_C 0 1)
    # -7.147671 -7.147671 -6.977363 -1.827134 
    # -6.384550 -7.894351 -1.130536 -6.443988 
    # -7.087621 -1.561341 -7.245101 -3.433812 
    # -1.588742 -7.563783 -2.748739 -6.977363 
    #9 e1_2_2_7 2 2 1 (WW_C 0 1)
    # -2.471456 -7.017148 -2.069534 -2.651657 
    # -4.211511 -4.845947 -1.703997 -2.392230 
    # -7.017148 -1.497281 -7.017148 -3.549064 
    # -2.471456 -7.017148 -7.017148 -2.471456 
    #10 e1_2_2_8 2 2 1 (WW_C 0 1)
    # -3.633434 -4.685472 -2.952277 -2.351640 
    # -2.952277 -3.482695 -1.889795 -4.027051 
    # -5.312589 -1.575071 -4.685472 -3.351729 
    # -1.534371 -3.633434 -2.799852 -4.685472 
    #11 e1_2_2_9 2 2 1 (WW_C 0 1)
    # -7.434281 -7.780056 -6.844328 -2.135404 
    # -6.885685 -8.312929 -1.343189 -7.681776 
    # -6.659871 -1.075496 -7.149952 -2.951514 
    # -1.962604 -7.236048 -2.491320 -7.236048 
    #12 e1_2_2_10 2 2 1 (WW_C 0 1)
    # -7.182310 -3.088884 -7.182310 -1.921846 
    # -3.088884 -2.972957 -1.921846 -4.426345 
    # -2.343983 -2.182653 -3.220036 -4.426345 
    # -1.715181 -3.548929 -2.972957 -3.371021 
    #13 e1_2_2_11 2 2 1 (WW_C 0 1)
    # -7.834803 -6.775569 -7.834803 -1.912898 
    # -8.471646 -7.834803 -1.208842 -7.834803 
    # -6.775569 -1.206347 -7.056579 -2.627581 
    # -1.894138 -7.834803 -3.549864 -6.953882 
    #14 e1_2_2_12 2 2 1 (WW_C 0 1)
    # -6.616918 -6.728942 -6.855119 -1.702395 
    # -6.469388 -7.491295 -1.375426 -7.264901 
    # -6.616918 -1.354336 -6.565285 -3.302933 
    # -1.588047 -7.264901 -2.841789 -7.080419 
    #15 e1_2_2_13 2 2 1 (WW_C 0 1)
    # -3.943350 -7.107656 -4.327914 -2.082798 
    # -2.873304 -4.327914 -1.901572 -7.107656 
    # -3.666285 -1.190607 -3.120488 -4.327914 
    # -1.679526 -4.327914 -7.107656 -3.271551 
    #16 e1_2_2_14 2 2 1 (WW_C 0 1)
    # -6.771789 -7.295392 -7.717062 -1.947418 
    # -6.586298 -8.021215 -1.291293 -7.295392 
    # -5.660561 -1.628477 -7.295392 -2.021286 
    # -1.837470 -6.019926 -2.464092 -6.586298 
    #17 e1_2_2_15 2 2 1 (WW_C 0 1)
    # -4.892611 -4.892611 -4.236652 -2.190924 
    # -7.509773 -5.515291 -1.598415 -7.509773 
    # -3.343356 -1.148982 -4.236652 -3.843897 
    # -1.968761 -3.843897 -3.343356 -2.657993 
    #0
    comment_pattern = r'\#'
    
    t0_pattern = r'0\s+2\s+(\S+)\s+(\S+)'
    t1_pattern = r'1\s+2\s+(\S+)\s+(\S+)'
    t2_pattern = r'2\s+2\s+(\S+)\s+(\S+)'

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
        
        t0_match = re.search(t0_pattern, line)
        t1_match = re.search(t1_pattern, line)
        t2_match = re.search(t2_pattern, line)
        
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
            
        if t0_match: t0 = np.array([float(t0_match.group(1)), float(t0_match.group(2))])
        if t1_match: t1 = np.array([float(t1_match.group(1)), float(t1_match.group(2))])
        if t2_match: t2 = np.array([float(t2_match.group(1)), float(t2_match.group(2))])

        if emit_match and is_single and not is_pair:
            e_single = np.array([float(emit_match.group(1)),float(emit_match.group(2)),float(emit_match.group(3)),float(emit_match.group(4))])
            
        if emit_match and is_pair and not is_single:
            e_pair.append([float(emit_match.group(1)),float(emit_match.group(2)),float(emit_match.group(3)),float(emit_match.group(4))])

    e_pair = np.array(e_pair).flatten()
    if scaled: params = {"t0":     deepcopy(t0),     "t1": deepcopy(t1),     "t2": deepcopy(t2), "pe_single": deepcopy(e_single), "pe_pair": deepcopy(e_pair) }
    else:      params = {"log_t0": deepcopy(t0), "log_t1": deepcopy(t1), "log_t2": deepcopy(t2), "e_single":  deepcopy(e_single),  "e_pair": deepcopy(e_pair)  }

    print(f"{pprint.pformat(params)}")

    return G6S_normalize_params(params, scaled)

def G6S_write_paramfile(rundir, epoch, params):
    param_file = rundir / f"param_i{epoch}.param"

    # e1_2_0_0 appears first
    with open(param_file, "a") as f:
        f.write(f"3\n")
        f.write(f"0 2 {params['log_t0'][0]} {params['log_t0'][1]}\n")
        f.write(f"1 2 {params['log_t1'][0]} {params['log_t1'][1]}\n")
        f.write(f"2 2 {params['log_t2'][0]} {params['log_t2'][1]}\n")
        f.write(f"2\n")
        f.write(f"0 e1_2_0_0 2 0 1 (WW_C 0 1)\n")
        f.write(f"{params['e_pair'][0]}  {params['e_pair'][1]}  {params['e_pair'][2]}  {params['e_pair'][3]}\n")
        f.write(f"{params['e_pair'][4]}  {params['e_pair'][5]}  {params['e_pair'][6]}  {params['e_pair'][7]}\n")
        f.write(f"{params['e_pair'][8]}  {params['e_pair'][9]}  {params['e_pair'][10]} {params['e_pair'][11]}\n")
        f.write(f"{params['e_pair'][12]} {params['e_pair'][13]} {params['e_pair'][14]} {params['e_pair'][15]}\n")
        f.write(f"1 e1_1_0_0 1 0 0\n")
        f.write(f"{params['e_single'][0]} {params['e_single'][1]} {params['e_single'][2]} {params['e_single'][3]}\n")
        f.write(f"0\n")
    return str(param_file)
