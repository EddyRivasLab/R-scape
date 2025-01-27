import math
import numpy as np
import functools
import jax
import jax.numpy as jnp
import jax.scipy as jsp
import jax.nn    as jnn
import jax.random as random
import optax

from jax import grad, value_and_grad

from G5 import (
    G5_Inside_JAX,
    G5_Inside_std,
    G5_Inside_JAX_scaled,
    G5_Inside_std_scaled,
)
from G5_param import (
    G5_param_tornado,
    G5_param_uniform,
)
from seq_io import (
    read_fasta,
    pad_seq,
)

def loss_fn(g5_inside_jax, mask, log_psq, log_psq2, log_t, e_single, e_pair):
    
    def sq_loss(mask_s, log_psq_s, log_psq2_s):
        return g5_inside_jax(mask_s, log_psq_s, log_psq2_s, log_t, e_single, e_pair)
    
    # vectorize single_sequence_loss across all sequences
    ns = log_psq.shape[0]
    all_losses = jax.vmap(sq_loss)(mask, log_psq, log_psq2)
    return -jnp.sum(all_losses) / ns  # minimize negative log likelihood

def loss_fn_scaled(g5_inside_jax, mask, psq, psq2, t, pe_single, pe_pair):

    def sq_loss(mask_s, psq_s, psq2_s):
        return g5_inside_jax(mask_s, psq_s, psq2_s, t, pe_single, pe_pair)
    
    # vectorize single_sequence_loss across all sequences
    ns = psq.shape[0]
    all_losses = jax.vmap(sq_loss)(mask, psq, psq2)
    return -jnp.sum(jnp.log(all_losses)) / ns # minimize negative log likelihood


def optimize_param_G5(mask, psq, log_psq, verbose: int, init_uniform: int, scaled: int, K: int=4, min_hairpin: int = 0):

    scale = K

    # initialized to pre-trained ML values
    if init_uniform: log_t, t, e_single, pe_single, e_pair, pe_pair = G5_param_uniform(K, verbose=1)
    else:            log_t, t, e_single, pe_single, e_pair, pe_pair = G5_param_tornado(verbose=1)

    if scaled: g5_inside_jax = G5_Inside_JAX_scaled(scale, verbose, K, min_hairpin)
    else:      g5_inside_jax = G5_Inside_JAX(              verbose, K, min_hairpin)

    # initialize optimizer
    l_rate = 1.0 
    optimizer = optax.adam(learning_rate=l_rate)
    if scaled: opt_state = optimizer.init((t,     pe_single, pe_pair))
    else:      opt_state = optimizer.init((log_t, e_single,  e_pair))
    
    # training loop with mini-batch
    n_epoch    = 50
    batch_size = 200

    key = random.PRNGKey(0)
    losses = []
    for epoch in range(n_epoch):
        key, subkey = random.split(key)
        perm = random.permutation(subkey, len(psq))
        
        psq_sh  = jnp.array([psq[i] for i in perm])
        mask_sh = jnp.array([mask[i] for i in perm])
  
        for i in range(0, len(psq), batch_size):
            psq_b  = psq_sh [i:i+batch_size]
            mask_b = mask_sh[i:i+batch_size]
            psq2_b = np.einsum('...ia,...jb->...ijab', psq_b, psq_b)
            
            log_psq_b  = jnp.log(psq_b)
            log_psq2_b = jnp.log(psq2_b)

            # compute gradients
            if scaled:
                loss, grads = value_and_grad(loss_fn_scaled, argnums=(4,5,6))(g5_inside_jax, mask_b, psq_b,     psq2_b,     t,     pe_single, pe_pair)
            else:
                loss, grads = value_and_grad(loss_fn,        argnums=(4,5,6))(g5_inside_jax, mask_b, log_psq_b, log_psq2_b, log_t, e_single,  e_pair)
                
            # update parameters
            if scaled:
                updates, opt_state = optimizer.update(grads, opt_state, (t, pe_single, pe_pair))
                t, pe_single, pe_pair = optax.apply_updates((t, pe_single, pe_pair), updates)

                # normalize the probabilities
                t         = jnn.softmax(t,         axis=-1)
                pe_single = jnn.softmax(pe_single, axis=-1)
                pe_pair   = jnn.softmax(pe_pair,   axis=-1)
                
            else:
                updates, opt_state = optimizer.update(grads, opt_state, (log_t, e_single, e_pair))
                log_t, e_single, e_pair = optax.apply_updates((log_t, e_single, e_pair), updates)

                # normalize the probabilities
                log_t    = jnp.log(jnn.softmax(np.exp(log_t),    axis=-1))
                e_single = jnp.log(jnn.softmax(np.exp(e_single), axis=-1))
                e_pair   = jnp.log(jnn.softmax(np.exp(e_pair),   axis=-1))
               
        if scaled:
            print("t:", t)
            print("pe_single:", pe_single)
            print("pe_pair:", pe_pair)
        else:
            print("log_t:", log_t)
            print("e_single:", e_single)
            print("e_pair:", e_pair)
                
            
        losses.append(loss)

        if epoch % 1 == 0:
            print(f"Epoch {epoch}, Loss: {loss:.4f}")
   
def main():

    verbose = 0
    
    filepath = 'data/trna1415_annote_1of1.fa'
    mask, psq, log_psq, psq_name = read_fasta(filepath)
    
    print("optimize param G5 (logsumexp)")
    #optimize_param_G5(mask, psq, log_psq, verbose, init_uniform=0, scaled=0, K=4, min_hairpin=0)
    optimize_param_G5(mask, psq, log_psq, verbose, init_uniform=1, scaled=0, K=4, min_hairpin=0)

    #print("optimize param G5 (scaled)")
    #optimize_param_G5(mask, psq, log_psq, verbose, init_uniform=0, scaled=1, K=4, min_hairpin=0)
    
        
if __name__ == "__main__":
    main()
