import numpy as np
import jax
import jax.numpy as jnp

RNA_ALPHA = "ACGU"
K = len(RNA_ALPHA)

def dseq_to_pseq(dsq, K):

    # Create a prob sequence
    psq = np.zeros((len(dsq), K), dtype=float)
    
    # Create a dictionary 
    nt_idx = {symbol: i for i, symbol in enumerate(RNA_ALPHA)}
 
    for i, a in enumerate(dsq):
        if a in nt_index:
            idx = nt_to_index[nt]
            psq[i, idx] = 1.0
 
    return psq;

def pad_seq(sqs):
    L_max = max(len(sq) for sq in sqs)
    
    sq_padded = jnp.array([jnp.pad(sq,                ((0, L_max - len(sq)), (0, 0)), mode='constant') for sq in sqs])
    sq_mask   = jnp.array([jnp.pad(jnp.ones(len(sq)),  (0, L_max - len(sq)),          mode='constant') for sq in sqs])
    return sq_padded, sq_masks

def read_fasta(filepath):
    with open(path, 'r') as f:
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
    psq = [dsq_to_psq(sq) for sq in dsq[:]]

    return psq, name
