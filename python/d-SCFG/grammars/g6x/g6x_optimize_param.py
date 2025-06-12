import os
import shutil

import math
import numpy as np
import functools
import pdb
import argparse
from pathlib import Path
from copy import deepcopy
import pprint
from tqdm import tqdm
import matplotlib.pyplot as plt

import jax
import jax.numpy as jnp
import jax.scipy as jsp
import jax.nn    as jnn
import jax.random as random
from jax import grad, value_and_grad
import optax

import lib.seqio       as seqio
import lib.probability as prob
from   lib.utils import bcolors, tree_stack, plot_losses

from grammars.g6x import g6x_inside
from grammars.g6x import g6x_params

from grammars.tornado import tornado_fold

def get_argparse():
    parser = argparse.ArgumentParser(description="Optimize SCFG parameters via JAX")

    #outdir arguments
    parser.add_argument('run_name', help='name of output directory')

    # optional arguments
    #output parameters
    parser.add_argument('--outdir', type=str, default="experiments/",
                        help='path to output directory')
    parser.add_argument('--verbose', action='store_true',
                        help="verbose")

    # grammar parameters
    parser.add_argument('--scaled', action='store_true',
                        help="uses scaled probability parameters in the inside algorithm")
    parser.add_argument('--force_WCF', action='store_true',
                        help="force Watson-Crick-Franklin base pairs")
    parser.add_argument('--min_hairpin', type=int, default=0,
                        help="minimum hairpin size")

    # data parameters
    parser.add_argument('--train_data', type=str, default="data/trna1415_annote_1of1.fa",
                        help='input fasta file')
    parser.add_argument('--bymin', action='store_true',
                        help="sequence padding my min length, default is by max length")

    # test data (optional)
    parser.add_argument('--test_data', type=str, default="",
                        help='input fasta file')
    parser.add_argument('--fold_method', type=str, default="mea",
                        help='folding method: mea, cyk')
    parser.add_argument('--grm_file', type=str, default="../../lib/tornado/grammars/g6.grm",
                        help='file with grammar definition')
    parser.add_argument('--postgrm_file', type=str, default="../../lib/tornado/grammars/g6.grm",
                        help='file with postgrammar definition')
   

    init_param = parser.add_mutually_exclusive_group()
    init_param.add_argument('--init_tornado', action='store_true')
    init_param.add_argument('--init_naive',   action='store_true')
    init_param.add_argument('--init_uniform', action='store_true')

    # optimization parameters
    parser.add_argument('--n_epoch', type=int, default=50,
                        help="# of gradient descent iterations")
    parser.add_argument('--batch_size', type=int, default=200,
                        help="batch size for SGD")
    parser.add_argument('--lr', type=float, default=0.1,
                        help="learning rate")
    parser.add_argument('--grad_accum', action='store_true',
                        help="uses gradient accumulation rather than SGD")

    return parser

def loss_fn_g6x_logs(g6x_inside_jax, mask, log_psq, log_psq2, log_t0, log_t1, log_t2, e_single, e_pair):

    def sq_loss(mask, log_psq, log_psq2):
        return g6x_inside_jax(mask, log_psq, log_psq2, log_t0, log_t1, log_t2, e_single, e_pair)

    # vectorize single_sequence_loss across all sequences
    ns = log_psq.shape[0]
    all_losses = jax.vmap(sq_loss)(mask, log_psq, log_psq2)
    return -jnp.sum(all_losses) / ns  # minimize negative log likelihood

def loss_fn_g6x_probs(g6x_inside_jax, mask, psq, psq2, t0, t1, t2, pe_single, pe_pair):

    def sq_loss(mask, psq, psq2):
        return g6x_inside_jax(mask, psq, psq2, t0, t1, t2, pe_single, pe_pair)

    # vectorize single_sequence_loss across all sequences
    ns = psq.shape[0]
    all_losses = jax.vmap(sq_loss)(mask, psq, psq2)
    return -jnp.sum(jnp.log(all_losses)) / ns # minimize negative log likelihood
  

def optimize_param_G6X(args, rundir, K, mask, psq, log_psq, verbose):

    # initialization parameters
    init_tornado = args['init_tornado']
    init_naive   = args['init_naive']
    init_uniform = args['init_uniform']

    # grammar parameters
    scaled      = args['scaled']
    min_hairpin = args['min_hairpin']
    force_WCF   = args['force_WCF']
    scale       = K
    if scaled: print("optimize G6X params (scaled)")
    else:      print("optimize G6X params (logsumexp)")

    # optimization parameters
    n_epoch    = args['n_epoch']
    batch_size = args['batch_size']
    l_rate     = args['lr']

    # test data (optional)
    test_file = args['test_data']
    if os.path.exists(test_file):
        print(f"Test sequences:", test_file)
    fold_method  = args['fold_method']
    grm_file     = args['grm_file']
    postgrm_file = args['postgrm_file']
  
               
    def loss_fn_g6x(unnorm_params, mask, psq, psq2):
        params = g6x_params.G6X_normalize_params(unnorm_params, scaled)
        
        if scaled:
            loss = loss_fn_g6x_probs(g6x_inside_jax, mask, psq, psq2, params['t0'], params['t1'], params['t2'], params['pe_single'], params['pe_pair'])
        else:
            log_psq  = jnp.log(psq)
            log_psq2 = jnp.log(psq2)
            
            loss = loss_fn_g6x_logs(g6x_inside_jax, mask, log_psq, log_psq2, params['log_t0'], params['log_t1'], params['log_t2'], params['e_single'], params['e_pair'])
            
            return loss
        
    grad_fn_g6x = value_and_grad(loss_fn_g6x)
    grad_fn_g6x = jax.jit(grad_fn_g6x)
 
    # initialized to pre-trained ML values
    if   init_tornado: log_t0, t0, log_t1, t1, log_t2, t2, e_single, pe_single, e_pair, pe_pair = g6x_params.G6X_param_tornado(verbose)
    elif init_naive:   log_t0, t0, log_t1, t1, log_t2, t2, e_single, pe_single, e_pair, pe_pair = g6x_params.G6X_param_naive(verbose)
    else:              log_t0, t0, log_t1, t1, log_t2, t2, e_single, pe_single, e_pair, pe_pair = g6x_params.G6X_param_uniform(K, verbose)

    if scaled: params = {"t0":     deepcopy(t0),     "t1":     deepcopy(t1),     "t2":     deepcopy(t2),     "pe_single": deepcopy(pe_single), "pe_pair": deepcopy(pe_pair) }
    else:      params = {"log_t0": deepcopy(log_t0), "log_t1": deepcopy(log_t1), "log_t2": deepcopy(log_t2), "e_single":  deepcopy(e_single),  "e_pair":  deepcopy(e_pair)  }

    if   init_tornado: print("tornado initialization:")
    elif init_naive:   print("naive initialization:")
    else:              print("uniform initialization:")
    print(f"input params: {pprint.pformat(params)}")

    if scaled: g6x_inside_jax = g6x_inside.G6X_Inside_JAX_scaled(scale, verbose, K, min_hairpin)
    else:      g6x_inside_jax = g6x_inside.G6X_Inside_JAX(              verbose, K, min_hairpin)

    # initialize optimizer
    optimizer = optax.adam(learning_rate=l_rate)
    opt_state = optimizer.init(params)

    # Training loop with mini-batch SGD
    loss_file = rundir / "loss.txt"

    key = random.PRNGKey(0)
    losses  = []
    acc_sen = [0 for _ in range(n_epoch)]
    acc_ppv = [0 for _ in range(n_epoch)]
    acc_f1  = [0 for _ in range(n_epoch)]
    
    plot_every = 1
    for epoch in range(n_epoch):
        key, subkey = random.split(key)
        perm = random.permutation(subkey, len(psq))
        
        psq_sh  = jnp.array([psq[i]  for i in perm])
        mask_sh = jnp.array([mask[i] for i in perm])

        epoch_loss = 0.0

        for i in range(0, len(psq), batch_size):
            psq_b  = psq_sh [i:i+batch_size]
            mask_b = mask_sh[i:i+batch_size]
            psq2_b = np.einsum('...ia,...jb->...ijab', psq_b, psq_b)
            
            # compute gradients
            loss, grads = grad_fn_g6x(params, mask_b, psq_b, psq2_b)
            epoch_loss += loss
 
            # update parameters and normalize
            updates, opt_state = optimizer.update(grads, opt_state, params)
            params = optax.apply_updates(params, updates)

        norm_params = g6x_params.G6X_normalize_params(params, scaled)
        norm_params_str = pprint.pformat(norm_params)

        print(f"\n----- EPOCH {epoch} -----")
        print(f"- Loss: {epoch_loss}")
        print(f"- Params: {norm_params_str}")

        with open(loss_file, "a") as f:
            f.write(f"{epoch_loss}\n")
        losses.append(epoch_loss)

        if epoch and epoch == 0 or epoch % plot_every == 0:
            param_file = g6x_params.G6X_write_paramfile(rundir, epoch, norm_params)

            # run TORNADO with those parameters
            if os.path.exists(test_file):
                sen, ppv, f1 = tornado_fold.tornado_fold(fold_method, param_file, grm_file, postgrm_file, test_file)
                print(f"Epoch {epoch} sen {sen} ppv {ppv} f1 {f1}")
                
                acc_sen[epoch] = sen
                acc_ppv[epoch] = ppv
                acc_f1[epoch]  = f1

            # print the losses
            plot_losses(rundir, epoch, losses, acc_sen, acc_ppv, acc_f1)

    # plot after all iterations
    plot_losses(rundir, epoch, losses, acc_sen, acc_ppv, acc_f1)

    
def main(args):

    # train data
    train_file = args['train_data']
    bymin      = args['bymin']
    print(f"Train sequences:", train_file)
    mask, psq, log_psq, psq_name = seqio.read_fasta(train_file, bymin)

    # data parameters
    K = 4

    # ouput directory
    outdir = Path(args['outdir'])
    assert(outdir.exists())
 
    run_name = args['run_name']
    rundir = outdir / run_name
    with os.scandir(rundir) as it:
        if any(it): shutil.rmtree(rundir)
    rundir.mkdir(parents=False, exist_ok=True)

    # other option
    verbose = args['verbose']

    # run paramater optimizer
    optimize_param_G6X(args, rundir, K, mask, psq, log_psq, verbose)


    
if __name__ == "__main__":
    parser = get_argparse()
    args = vars(parser.parse_args())

    main(args)
