import math
import numpy as np
import jax
import jax.numpy as jnp

RNA_ALPHA = "ACGU"
K = len(RNA_ALPHA)

def dsq_to_psq(dsq, K):
    # Create a prob sequence
    psq = jnp.full((len(dsq), K), 0.0) 
    
    # Create a dictionary 
    nt_idx = {symbol: i for i, symbol in enumerate(RNA_ALPHA)}
 
    for i, a in enumerate(dsq):
       if a in nt_idx:
            idx = nt_idx[a]
            psq = psq.at[i,idx].set(1.0)
       else:
           if a == "N":
               psq = psq.at[i].set(1/K)
           elif a == "Y":
               psq = psq.at[i,nt_idx["C"]].set(2/K)
               psq = psq.at[i,nt_idx["U"]].set(2/K)
           elif a == "R":
               psq = psq.at[i,nt_idx["A"]].set(2/K)
               psq = psq.at[i,nt_idx["G"]].set(2/K)
    
    return psq

def pad_seq(sqs):
    L_max = max(len(sq) for sq in sqs)
        
    sq_padded = jnp.array([jnp.pad(sq,               ((0, L_max - len(sq)), (0,0)), mode='constant',constant_values=1/K) for sq in sqs])
    sq_mask   = jnp.array([jnp.pad(jnp.ones(len(sq)), (0, L_max - len(sq)),         mode='constant') for sq in sqs])
    
    return sq_padded, sq_mask

def read_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

        dseq = [] # discrete sequence
        name = []
        for line in lines:
            if line[0] == '>':
                dseq.append('')
                name.append(line[1:])
            else:
                dseq[-1] += line.strip()

    # the psqs
    psq = [dsq_to_psq(sq, K) for sq in dseq[:]]

    # pad the sequences
    psq, mask = pad_seq(psq)
    log_psq = jnp.log(psq)
 
    return mask, psq, log_psq, name
