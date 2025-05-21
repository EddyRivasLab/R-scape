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
from   lib.utils import bcolors, tree_stack

from grammars.g6x import g6x_inside
from grammars.g6x import g6x_params

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
    parser.add_argument('--data', type=str, default="data/trna1415_annote_1of1.fa",
                        help='input fasta file')
    parser.add_argument('--bymin', action='store_true',
                        help="sequence padding my min length, default is by max length")

    init_param = parser.add_mutually_exclusive_group()
    init_param.add_argument('--init_tornado', action='store_true')
    init_param.add_argument('--init_naive',   action='store_true')
    init_param.add_argument('--init_uniform', action='store_true')

    # optimization parameters
    parser.add_argument('--n_epoch', type=int, default=500,
                        help="# of gradient descent iterations")
    parser.add_argument('--batch_size', type=int, default=200,
                        help="batch size for SGD")
    parser.add_argument('--lr', type=float, default=0.1,
                        help="learning rate")
    parser.add_argument('--grad_accum', action='store_true',
                        help="uses gradient accumulation rather than SGD")

    return parser

def loss_fn_logs(g6x_inside_jax, mask, log_psq, log_psq2, log_t, e_single, e_pair):

    def sq_loss(mask, log_psq, log_psq2):
        return g6x_inside_jax(mask, log_psq, log_psq2, log_t, e_single, e_pair)

    # vectorize single_sequence_loss across all sequences
    ns = log_psq.shape[0]
    all_losses = jax.vmap(sq_loss)(mask, log_psq, log_psq2)
    return -jnp.sum(all_losses) / ns  # minimize negative log likelihood

def loss_fn_probs(g6x_inside_jax, mask, psq, psq2, t, pe_single, pe_pair):

    def sq_loss(mask, psq, psq2):
        return g6x_inside_jax(mask, psq, psq2, t, pe_single, pe_pair)

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

    def process_params_g6x(uparams):
        if scaled:
            return {
                "t0":         prob.pNorm(uparams['t0']),
                "t1":         prob.pNorm(uparams['t1']),
                "t2":         prob.pNorm(uparams['t2']),
                "pe_single": prob.pNorm(uparams['pe_single']),
                "pe_pair":   prob.pNorm(uparams['pe_pair'])
            }
        else:
            return {
                "log_t0":    prob.logpNorm(uparams['log_t0']),
                "log_t1":    prob.logpNorm(uparams['log_t1']),
                "log_t2":    prob.logpNorm(uparams['log_t2']),
                "e_single": prob.logpNorm(uparams['e_single']),
                "e_pair":   prob.logpNorm(uparams['e_pair'])
            }
        

    def loss_fn(unprocessed_params, mask, psq, psq2):
        params = process_params_g6x(unprocessed_params)
        
        if scaled:
            loss = loss_fn_probs(g6x_inside_jax, mask, psq, psq2, params['t0'], params['t1'], params['t2'], params['pe_single'], params['pe_pair'])
        else:
            log_psq  = jnp.log(psq)
            log_psq2 = jnp.log(psq2)
            
            loss = loss_fn_logs(g6x_inside_jax, mask, log_psq, log_psq2, params['log_t0'], params['log_t1'], params['log_t2'], params['e_single'], params['e_pair'])
            
        return loss
    grad_fn = value_and_grad(loss_fn)
    grad_fn = jax.jit(grad_fn)

        
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
    losses = []
    plot_every = 5
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
            loss, grads = grad_fn(params, mask_b, psq_b, psq2_b)
            epoch_loss += loss
 
            # update parameters and normalize
            updates, opt_state = optimizer.update(grads, opt_state, params)
            params = optax.apply_updates(params, updates)

        processed_params = process_params_g6x(params)
        processed_params_str = pprint.pformat(processed_params)

        print(f"\n----- EPOCH {epoch} -----")
        print(f"- Loss: {epoch_loss}")
        print(f"- Params: {processed_params_str}")

        with open(loss_file, "a") as f:
            f.write(f"{epoch_loss}\n")

        losses.append(epoch_loss)

        if epoch and epoch % plot_every == 0:
            plt.plot(losses)
            plt.xlabel("Iteration")
            plt.ylabel("Loss")
            plt.savefig(rundir / f"losses_i{epoch}.png")
            plt.clf()

    # plot after all iterations
    plt.plot(losses)
    plt.xlabel("Iteration")
    plt.ylabel("Loss")
    plt.savefig(rundir / f"losses_i{epoch}.png")
    plt.clf()


    
def main(args):

    # load data
    filedir = args['data']
    bymin   = args['bymin']
    
    filename = 'data/trna1415_annote_1of1.fa'
    print(f"loading sequences from:", filename)
    mask, psq, log_psq, psq_name = seqio.read_fasta(filename, bymin)

    # data parameters
    K = 4

    # ouput parameters
    outdir = Path(args['outdir'])
    assert(outdir.exists())
 
    run_name = args['run_name']
    rundir = outdir / run_name
    rundir.mkdir(parents=False, exist_ok=False)
    print(f"output dir:", rundir)

    # other option
    verbose = args['verbose']

    # run paramater optimizer
    optimize_param_G6X(args, rundir, K, mask, psq, log_psq, verbose)


    
if __name__ == "__main__":
    parser = get_argparse()
    args = vars(parser.parse_args())

    main(args)
