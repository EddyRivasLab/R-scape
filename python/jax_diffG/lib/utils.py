import math
import numpy as np

import jax
import jax.numpy as jnp

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def tree_stack(trees):
    return jax.tree.map(lambda *v: jnp.stack(v), *trees)
