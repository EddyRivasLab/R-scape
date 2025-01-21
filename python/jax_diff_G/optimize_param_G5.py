import math
import numpy as np
import functools
import jax
import jax.numpy as jnp
import jax.scipy as jsp
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
    G5_param_random,
}

from seq_io import (
    read_fasta,
    pad_sq,
}

def loss_fn(log_t, e_single, e_pair, psqs, masks):
    
    def sq_loss(p_seq, mask):
        return g5_inside_jax(p_seq, log_t, e_single, e_pair, mask)
    
    # vectorize single_sequence_loss across all sequences
    all_losses = vmap(sq_loss)(psqs, masks)
    return -jnp.log(jnp.sum(all_losses)) # minimize negative log likelihood





def optimize_param_G5(psq, int_random: int, scaled: int, verbose: int, K: int=4, min_hairpin: int = 0):

    scale = K

    # initialized to pre-trained ML values
    if init_random: log_t, t, e_single, pe_single, e_pair, pe_pair = G5_param_random(verbose)
    else:           log_t, t, e_single, pe_single, e_pair, pe_pair = G5_param_tornado(verbose)

    if scaled: g5_inside_jax = G5_Inside_JAX_scaled(scale, verbose, K, min_hairpin)
    else:      g5_inside_jax = G5_Inside_JAX(              verbose, K, min_hairpin)


    # initialize optimizer
    l_rate = 0.05 
    optimizer = optax.adam(learning_rate=l_rate)
    opt_state = optimizer.init((e_single, e_pair, log_t))
    
    # training loop with mini-batch
    n_epoch = 50
    batch_size = 200 # 

    key = random.PRNGKey(0)
    losses = []
    for epoch in range(n_epoch):
        key, subkey = random.split(key)
        perm = random.permutation(subkey, len(psq))
        psq_sh = [psq[i] for i in perm]
        
        for i in range(0, len(psq), batch_size):
            batch_psq = sequences_shuffled[i:i+batch_size]
            padded_psq, masks = pad_sequences(batch_psq)

            # compute gradients
            loss, grads = value_and_grad(loss_fn, argnums=(0, 1, 2))(log_t, e_single, e_pair, padded_psq, masks)
            
            # update parameters
            updates, opt_state = optimizer.update(grads, opt_state, (e_single, e_pair, log_t))
            e_single, e_pair, log_t = optax.apply_updates((e_single, e_pair, log_t), updates)
    
    # if epoch % 10 == 0:
    print(f"Epoch {epoch}, Loss: {loss:.4f}")
    losses.append(loss)

        # gradients wrt parameters
        # jax_log
        grad_g5_inside = value_and_grad(g5_inside_jax, argnums=(2,3,4))
        logp_inside, logp_inside_grads = grad_g5_inside(log_psq, log_psq2, log_t, e_single, e_pair)
        print("logp_inside:", logp_inside)
        print("logp_inside_gradients:", logp_inside_grads)
    

def main():

    verbose = 0
    
    filepath = 'data/trna1415_annote_1of1.fa'
    psq, psq_name = read_fasta(filepath)
    
    
    print("optimize param G5 (scaled)")
    optimize_param_G5(psq, init_random=0, scaled=1, verbose, K=4, min_hairpin=0)

    print("optimize param G5 (logsumexp)")
    optimize_param_G5(psq, init_random=0, scaled=0, verbose, K=4, min_hairpin=0)
    optimize_param_G5(psq, init_random=1, scaled=0, verbose, K=4, min_hairpin=0)

        



if __name__ == "__main__":
    main()
