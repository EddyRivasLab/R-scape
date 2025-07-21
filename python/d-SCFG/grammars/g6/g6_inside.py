import math
import sys
import numpy as np
import jax
import jax.numpy as jnp
import jax.nn    as jnn
import jax.scipy as jsp
import functools

from scipy.special           import logsumexp

# adapted from Ward et al. 2023
from lib.checkpoint import checkpoint_scan

checkpoint_every = 1
if checkpoint_every is None:
    scan = jax.lax.scan
else:
    scan = functools.partial(checkpoint_scan, checkpoint_every=checkpoint_every)


def G6_Inside_JAX_scaled(scale, verbose, K: int = 4, min_hairpin: int = 0):

    @jax.jit
    def g6_inside_jax_scaled(mask, psq, psq2, t0, t1, t2, pe_single, pe_pair):
        n   = psq.shape[0]
        ns  = jnp.arange(n)
        ss  = jnp.arange(K)
        ps  = jnp.arange(K*K)

        pe_single = pe_single * scale         # single emission probabilities
        pe_pair   = pe_pair   * scale * scale # pair   emission probabilities

        # G6/G6S 
        #
        #   S -> LS   | L
        #   L -> aFa' | a
        #   F -> aFa' | LS
        #
        def g6_S_i_js(i: int, S, L, F, t0, pe_single, pe_pair, K):

            """fills S[i,0..n-1] """
            g6_S_i_js_map = jax.vmap(lambda j: g6_S_ij_scaled(i, j, S, L, F, psq, psq2, t0, pe_single, pe_pair, K))(ns)
            S_i_js_ret = S.at[i].set(g6_S_i_js_map)
            return S_i_js_ret

        def g6_L_i_js(i: int, S, L, F, psq, psq2, t1, pe_single, pe_pair, K, min_hairpin):
            """fills L[i,0..n-1] """
            g6_L_i_js_map = jax.vmap(lambda j: g6_L_ij_scaled(i, j, S, L, F, psq, psq2, t1, pe_single, pe_pair, K, min_hairpin))(ns)
            L_i_js_ret = L.at[i].set(g6_L_i_js_map)
            return L_i_js_ret

        def g6_F_i_js(i: int, S, L, F, psq, psq2, t2, pe_single, pe_pair, K, min_hairpin):
            """fills F[i,0..n-1] """
            g6_F_i_js_map = jax.vmap(lambda j: g6_F_ij_scaled(i, j, S, L, F, psq, psq2, t2, pe_single, pe_pair, K, min_hairpin))(ns)
            F_i_js_ret = F.at[i].set(g6_F_i_js_map)
            return F_i_js_ret
                  
        def g6_fill(carry, i):
            S, L, F, psq, psq2, t0, t1, t2, pe_single, pe_pair = carry

            # order is: L F S
            L = g6_L_i_js(i, S, L, F, psq, psq2, t1, pe_single, pe_pair, K, min_hairpin)
            if verbose: jax.debug.print('^^FILL i={} L={}', i, L)
            
            F = g6_F_i_js(i, S, L, F, psq, psq2, t2, pe_single, pe_pair, K, min_hairpin)
            if verbose: jax.debug.print('^^FILL i={} F={}', i, F)
                
            S = g6_S_i_js(i, S, L, F, t0, pe_single, pe_pair, K)
            if verbose: jax.debug.print('^^FILL i={} S={}', i, S)
                
            return (S, L, F, psq, psq2, t0, t1, t2, pe_single, pe_pair), None

        # initialize
        S = jnp.zeros((n, n))        
        L = jnp.zeros((n, n))
        F = jnp.zeros((n, n))
        initial_carry = (S, L, F, psq, psq2, t0, t1, t2, pe_single, pe_pair)
        ies = jnp.arange(n-1, -1, -1)
        (S, L, F, _, _, _,  _,  _,  _, _), _ = scan(g6_fill, initial_carry, ies)
        
        len = jnp.count_nonzero(mask)
        return S[0,len-1]
                      
    return g6_inside_jax_scaled

def G6_Inside_JAX(verbose, K: int = 4, min_hairpin: int = 0):
    
    @jax.jit
    def g6_inside_jax(mask, log_psq, log_psq2, log_t0, log_t1, log_t2, e_single, e_pair):
        n   = log_psq.shape[0]
        ns  = jnp.arange(n)
        ss  = jnp.arange(K)
        ps  = jnp.arange(K*K)

        # G6/G6S Grammar 
        #
        #   S -> LS   | L
        #   L -> aFa' | a
        #   F -> aFa' | LS
        #
        def g6_logS_i_js(i: int, logS, logL, logF, log_psq, log_psq2, log_t0, e_single, e_pair, K):
            """fills logS[i,0..n-1] """
            g6_logS_i_js_map = jax.vmap(lambda j: g6_logS_ij(i, j, logS, logL, logF, log_psq, log_psq2, log_t0, e_single, e_pair, K))(ns)
            logS_i_js_ret = logS.at[i].set(g6_logS_i_js_map)
            return logS_i_js_ret

        def g6_logL_i_js(i: int, logS, logL, logF, log_psq, log_psq2, log_t1, e_single, e_pair, K, min_hairpin):
            """fills logL[i,0..n-1] """
            g6_logL_i_js_map = jax.vmap(lambda j: g6_logL_ij(i, j, logS, logL, logF, log_psq, log_psq2, log_t1, e_single, e_pair, K, min_hairpin))(ns)
            logL_i_js_ret = logL.at[i].set(g6_logL_i_js_map)
            return logL_i_js_ret

        def g6_logF_i_js(i: int, logS, logL, logF, log_psq, log_psq2, log_t2, e_single, e_pair, K, min_hairpin):
            """fills logF[i,0..n-1] """
            g6_logF_i_js_map = jax.vmap(lambda j: g6_logF_ij(i, j, logS, logL, logF, log_psq, log_psq2, log_t2, e_single, e_pair, K, min_hairpin))(ns)
            logF_i_js_ret = logF.at[i].set(g6_logF_i_js_map)
            return logF_i_js_ret
                  
        def g6_fill(carry, i):
            logS, logL, logF, log_psq, log_psq2, log_t0, log_t1, log_t2, e_single, e_pair = carry
            
            # order is: L F S
            logL = g6_logL_i_js(i, logS, logL, logF, log_psq, log_psq2, log_t1, e_single, e_pair, K, min_hairpin)
            if verbose: jax.debug.print('^^FILL i={} logL={}', i, logL)
            
            logF = g6_logF_i_js(i, logS, logL, logF, log_psq, log_psq2, log_t2, e_single, e_pair, K, min_hairpin)
            if verbose: jax.debug.print('^^FILL i={} logF={}', i, logF)
                
            logS = g6_logS_i_js(i, logS, logL, logF, log_psq, log_psq2, log_t0, e_single, e_pair, K)
            if verbose: jax.debug.print('^^FILL i={} logS={}', i, logS)
                
            return (logS, logL, logF, log_psq, log_psq2, log_t0, log_t1, log_t2, e_single, e_pair), None

        # initialize
        veryneg = np.log(sys.float_info.epsilon)
        logS = jnp.full((n, n),veryneg)     
        logL = jnp.full((n, n),veryneg)     
        logF = jnp.full((n, n),veryneg)     
        initial_carry = (logS, logL, logF, log_psq, log_psq2, log_t0, log_t1, log_t2, e_single, e_pair)
        ies = jnp.arange(n-1, -1, -1)
        (logS, logL, logF, _, _, _,  _,  _, _, _), _ = scan(g6_fill, initial_carry, ies)
        
        len = jnp.count_nonzero(mask)
        return logS[0,len-1]
                      
    return g6_inside_jax
    
def g6_logS_ij(i: int, j: int, logS, logL, logF, log_psq, log_psq2, log_t0, e_single, e_pair, K):
    """fills logS[i,j]"""
    n  = log_psq.shape[0]
    ns = jnp.arange(n)

    def logS_rule1_k(k: int):
        """ S-> L(i,k) S(k+1,j) """
        nonterm_1 = jax.lax.select(k>=j, np.log(1e-300), logL[i,k])
        nonterm_2 = jax.lax.select(k>=j, np.log(1e-300), logS[k+1,j])
            
        #  S -> LS 
        logS_rule1_k_val = log_t0[0] + nonterm_1 + nonterm_2
        return logS_rule1_k_val
    
    def logS_rule2_L():
        """ S -> L(i,j) """
        nonterm = jax.lax.select(i>j, np.log(1e-300), logL[i,j])
         
        # S -> L
        logS_rule2_L_val = log_t0[1] + nonterm
        return logS_rule2_L_val
    
    logS_rule1_val = jsp.special.logsumexp(jax.vmap(lambda k: logS_rule1_k(k))(ns))
    logS_rule2_val = logS_rule2_L()
    logS_ij_val    = jsp.special.logsumexp(jnp.array([logS_rule1_val, logS_rule2_val]))

    return logS_ij_val

def g6_S_ij_scaled(i: int, j: int, S, L, F, psq, psq2, t0, pe_single, pe_pair, K):
    """fills S[i,j]"""
    n  = psq.shape[0]
    ns = jnp.arange(n)

    def S_rule1_k(k: int):
        """ L(i,k) S(k+1,j) """
        nonterm_1 = jax.lax.select(k>=j, 0.0, L[i,k])
        nonterm_2 = jax.lax.select(k>=j, 0.0, S[k+1,j])
            
        #  S -> LS 
        S_rule1_k_val = t0[0] * nonterm_1 * nonterm_2
        return S_rule1_k_val
    
    def S_rule2_L():
        """ S -> L(i,j)  """
        nonterm = jax.lax.select(i>j, 0.0, L[i,j])
            
        #  S -> L
        S_rule2_L_val = t0[1] * nonterm
        return S_rule2_L_val

    S_rule1_val = jnp.sum(jax.vmap(lambda k: S_rule1_k(k))(ns))
    S_rule2_val = S_rule2_L()
    
    #  S -> LS + L
    S_ij_val    = S_rule1_val + S_rule2_val
    
    return S_ij_val

def g6_logL_ij(i: int, j: int, logS, logL, logF, log_psq, log_psq2, log_t1, e_single, e_pair, K, min_hairpin):
    """fills L[i,j]"""
    ss = jnp.arange(K)
    ps = jnp.arange(K*K)

    def logL_rule1_ab(ab: int):
        """ L(i,j) -> a F(i+1,j-1) b """
        nonterm = jax.lax.select(j>i+1, logF[i+1,j-1], np.log(1e-300))
        logL_rule1_ab_val = log_psq2[i,j,ab//K,ab%K] + log_t1[0] + e_pair[ab] + nonterm
        return logL_rule1_ab_val
    
    def logL_rule2_a(a: int):
        """ L(i,j) -> a """
        nonterm = jax.lax.select(j==i, 0.0, np.log(1e-300))
        logL_rule2_a_val = log_psq[i,a] + log_t1[1] + e_single[a] + nonterm
        return logL_rule2_a_val
  
    logL_rule1_val = jsp.special.logsumexp(jax.vmap(logL_rule1_ab)(ps))
    logL_rule2_val = jsp.special.logsumexp(jax.vmap(logL_rule2_a)(ss))
    
    #  L -> aFb + a
    logL_ij_val = jsp.special.logsumexp(jnp.array([logL_rule1_val, logL_rule2_val]))
    
    return logL_ij_val

def g6_L_ij_scaled(i: int, j: int, S, L, F, psq, psq2, t1, pe_single, pe_pair, K, min_hairpin):
    """fills L[i,j]"""
    ss = jnp.arange(K)
    ps = jnp.arange(K*K)

    def L_rule1_ab(ab: int):
        """ L(i,j) -> a F(i+1,j-1) b """
        nonterm = jax.lax.select(j>i+1, F[i+1,j-1], 0.0)
        L_rule1_ab_val = psq2[i,j,ab//K,ab%K] * t1[0] * pe_pair[ab] * nonterm
        return L_rule1_ab_val
    
    def L_rule2_a(a: int):
        """ L(i,j) -> a """
        nonterm = jax.lax.select(j==i, 1.0, 0.0)
        L_rule2_a_val = psq[i,a] * t1[1] * pe_single[a] *  nonterm
        return L_rule2_a_val
  
    L_rule1_val = jnp.sum(jax.vmap(L_rule1_ab)(ps))
    L_rule2_val = jnp.sum(jax.vmap(L_rule2_a)(ss))
    
    #  L -> aFb + a
    L_ij_val = L_rule1_val + L_rule2_val 
    
    return L_ij_val

def g6_logF_ij(i: int, j: int, logS, logL, logF, log_psq, log_psq2, log_t2, e_single, e_pair, K, min_hairpin):
    """fills F[i,j]"""
    n  = log_psq.shape[0]
    ns = jnp.arange(n)
    ps = jnp.arange(K*K)

    # rule1: F -> a F b
    def logF_rule1_ab(ab: int):
        """ F -> a F(i+1,j-1) b """
        nonterm = jax.lax.select(j>i+1, logF[i+1,j-1], np.log(1e-300))
        logF_rule1_ab_val = log_psq2[i,j,ab//K,ab%K] + log_t2[0] + e_pair[ab] + nonterm
        return logF_rule1_ab_val
    
    # rule2: F -> L S
    def logF_rule2_k(k: int):
        """ F -> LS """
        nonterm_1 = jax.lax.select(k>=j, np.log(1e-300), logL[i,k])
        nonterm_2 = jax.lax.select(k>=j, np.log(1e-300), logS[k+1,j])

        logF_rule2_k_val = log_t2[1] + nonterm_1 + nonterm_2
        return logF_rule2_k_val
  
    logF_rule1_val = jsp.special.logsumexp(jax.vmap(logF_rule1_ab)(ps))
    logF_rule2_val = jsp.special.logsumexp(jax.vmap(lambda k: logF_rule2_k(k))(ns))
    
    # F -> aFb + LS
    logF_ij_val = jsp.special.logsumexp(jnp.array([logF_rule1_val, logF_rule2_val]))
    
    return logF_ij_val

def g6_F_ij_scaled(i: int, j: int, S, L, F, psq, psq2, t2, pe_single, pe_pair, K, min_hairpin):
    """fills F[i,j]"""
    n  = psq.shape[0]
    ns = jnp.arange(n)
    ps = jnp.arange(K*K)

    # rule1: F -> a F b
    def F_rule1_ab(ab: int):
        """ F -> a F(i+1,j-1) b """
        nonterm = jax.lax.select(j>i+1, F[i+1,j-1], 0.0)
        F_rule1_ab_val = psq2[i,j,ab//K,ab%K] * t2[0] * pe_pair[ab] * nonterm
        return F_rule1_ab_val
    
    # rule2: F -> L S
    def F_rule2_k(k: int):
        """ F -> LS """
        nonterm_1 = jax.lax.select(k>=j, 0.0, L[i,k])
        nonterm_2 = jax.lax.select(k>=j, 0.0, S[k+1,j])

        F_rule2_k_val = t2[1] * nonterm_1 * nonterm_2
        return F_rule2_k_val
  
    F_rule1_val = jnp.sum(jax.vmap(F_rule1_ab)(ps))
    F_rule2_val = jnp.sum(jax.vmap(lambda k: F_rule2_k(k))(ns))
    
    #  F -> aFb + LS
    F_ij_val = F_rule1_val + F_rule2_val 
    
    return F_ij_val

def G6_Inside_std(mask, log_psq, log_psq2, log_t0, log_t1, log_t2, e_single, e_pair, verbose, K: int=4, min_hairpin: int=0):

    """Standard implementation of g6"""
    n = log_psq.shape[0]

    print("G6 param")
    print("     transistion S", log_t0)
    print("     transistion L", log_t1)
    print("     transistion F", log_t2)
    print("     emissions single", e_single)
    print("     emissions paired", e_pair)
    
    #   L -> aFa' | a
    #   F -> aFa' | LS
    #   S -> LS   | end
    logS = [[-math.inf for _ in range(n)] for _ in range(n)]
    logL = [[-math.inf for _ in range(n)] for _ in range(n)]
    logF = [[-math.inf for _ in range(n)] for _ in range(n)]

    for i in range(n-1, -1, -1):
        for j in range(i, n):
            
            log_emit_pair = -math.inf
            for ab in range(K*K):
                log_emit_pair = logsumexp([log_emit_pair,log_psq2[i][j][ab//K][ab%K] + e_pair[ab]])
 
            # L
            # ruleL1: L -> a F b
            term  = log_emit_pair + log_t1[0]  + (logF[i+1][j-1] if i < j-1  else -math.inf)
            logL[i][j] = logsumexp([logL[i][j], term])
 
            # ruleL2: L -> a
            for a in range(K):
                term = log_psq[i][a] + log_t1[1] + e_single[a] + (0.0 if i == j else -math.inf)
                logL[i][j] = logsumexp([logL[i][j], term])
                 
            # F
            # ruleF1: F -> a F b
            term  = log_emit_pair + log_t2[0] + (logF[i+1][j-1] if i < j-1  else -math.inf)
            logF[i][j] = logsumexp([logF[i][j], term])
            
            # ruleF2: F -> L S
            for k in range(i, j+1):
                nonterm_1 = logL[i][k]
                nonterm_2 = (logS[k+1][j] if k < j  else -math.inf)
                term = log_t2[1] + nonterm_1 + nonterm_2
                logF[i][j] = logsumexp([logF[i][j], term])

            # S 
            # ruleS1: S -> L S 
            for k in range(i, j+1):
                nonterm_1 = logL[i][k]
                nonterm_2 = (logS[k+1][j] if k < j  else -math.inf)
                term = log_t0[0] + nonterm_1 + nonterm_2
                logS[i][j] = logsumexp([logS[i][j], term])

            # ruleS2: S -> L
            nonterm = logL[i][j]
            term = log_t0[1] + nonterm
            logS[i][j] = logsumexp([logS[i][j], term])

            if verbose: print("log_standard", "i", i, "j", j, "logL", logL[i][j], "logF", logF[i][j], "logS", logS[i][j])
            
    len = jnp.count_nonzero(mask)
    return float(logS[0][len-1])

def G6_Inside_std_scaled(mask, psq, psq2, t0, t1, t2, pe_single, pe_pair, scale, verbose, K: int=4, min_hairpin: int=0):
    
    """Standard implementation of g6"""
    n = psq.shape[0]

    pe_single = pe_single * scale          # single emission probabilities
    pe_pair   = pe_pair   * scale * scale  # pair   emission probabilities

    print("G6 param")
    print("     transistion S", t0)
    print("     transistion L", t1)
    print("     transistion F", t2)
    print("     emissions single", pe_single)
    print("     emissions paired", pe_pair)
    
    #   L -> aFa' | a
    #   F -> aFa' | LS
    #   S -> LS   | L 
    #S = jnp.zeros((n, n))        
    #L = jnp.zeros((n, n))        
    #F = jnp.zeros((n, n))
    # Why does do not work here?
    S = [[0.0 for _ in range(n)] for _ in range(n)]
    L = [[0.0 for _ in range(n)] for _ in range(n)]
    F = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n-1, -1, -1):
        for j in range(i, n):

            emit_pair = 0
            for ab in range(K*K):
                emit_pair += psq2[i][j][ab//K][ab%K] * pe_pair[ab]

            # L
            # ruleL1: L -> a F b
            L[i][j] += emit_pair * t1[0] * (F[i+1][j-1] if i < j-1  else 0.0)
                
            # ruleL2: L -> a
            for a in range(K):
                L[i][j] += psq[i][a] * t1[1] * pe_single[a] * (1.0 if i == j else 0)
                
            # F
            # ruleF1: F -> a F b
            F[i][j] += emit_pair * t2[0] * (F[i+1][j-1] if i < j-1  else 0.0)
            
            # ruleF2: F -> L S 
            for k in range(i, j+1):
                nonterm_1 = L[i][k]
                nonterm_2 = (S[k+1][j] if k < j  else 0.0)
                F[i][j] += t2[1] * nonterm_1 * nonterm_2

            # S 
            # ruleS1: S -> L S 
            for k in range(i, j+1):
                nonterm_1 = L[i][k]
                nonterm_2 = (S[k+1][j] if k < j  else 0.0)
                S[i][j] += t0[0] * nonterm_1 * nonterm_2
                
            # ruleS2: S -> L
            nonterm = L[i][j]
            S[i][j] += t0[1] * nonterm
             
            if verbose: print("standard", "i", i, "j", j, "L", np.log(L[i][j]), "F", np.log(F[i][j]), "S", np.log(S[i][j]))
            
    #print("standard", "L", L, "F", F, "S", S)
            
    len = jnp.count_nonzero(mask)
    return (np.log(float(S[0][len-1])) - len * np.log(scale))

