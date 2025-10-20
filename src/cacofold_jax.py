import numpy as np
import jax
import jax.numpy as jnp
import jax.nn as jnn


def CACO_ProbNuss_Inside_JAX(K: int, min_hairpin: int = 0):
    
    @jax.jit
    def probnuss_inside_jax(log_psq, log_t, e_single, e_pair):
        # Nussinov Grammar
        #
        # S -> a S | a S b S | end
        #      t1       t2
        #
        L      = len(log_psq)
        psq    = jnn.softmax(log_psq, axis=1)
        L_arr  = jnp.arange(L)
        a_arr  = jnp.arange(K)
        ab_arr = jnp.arange(K*K)

        t         = jnn.softmax(log_t)     # transition probabilities
        pe_single = jnn.softmax(e_single)  # single emission probabilities
        pe_pair   = jnn.softmax(e_pair)    # pair   emission probabilities

        # initialize
        S_zero = jnp.zeros((L, L))
 
        def S_i_jarray(i_: int, S):
            """fills S[i,0..L-1]"""
            # Iterate backwards through i
            i = L - i_ - 1
            
            def S_ij(j: int):
                """fills S[i,j]"""

                def S_rule1_a(a):
                    """ a(i) S(i+1,j) """
                    nonterm_1 = jax.lax.select(i == j, 1.0, S[i+1,j])
                    S_rule1_a_val = psq[i,a] * t[0] * pe_single[a] * nonterm_1
                    return S_rule1_a_val
                
                def S_rule2_k(k: int):
                    def S_rule2_k_ab(ab):
                        """ a(i) S(i+1,k-1) b(k) S(k+1,j) """
                        nonterm_1 = jax.lax.select(k == i+1, 1.0, S[i+1,k-1])
                        nonterm_2 = jax.lax.select(k == j,   1.0, S[k+1,j])
                        
                        S_rule2_ab_val = psq[i,ab//4] * psq[k,ab%4] * t[1] * pe_pair[ab] * nonterm_1 * nonterm_2
                        return S_rule2_ab_val
                    
                    S_rule2_k_val = jnp.sum(jax.vmap(S_rule2_k_ab)(ab_arr))
                    S_rule2_k_val = jax.lax.select(jnp.logical_and(i + min_hairpin < k, k <= j), S_rule2_k_val, 0.0)
                    return S_rule2_k_val
                
                S_rule1_val = jnp.sum(jax.vmap(S_rule1_a)(a_arr))
                S_rule2_val = jnp.sum(jax.vmap(S_rule2_k)(L_arr))
                
                #  aS + a S b S 
                S_ij_val = jax.lax.select(j >= i, S_rule1_val + S_rule2_val, 0.0)
                return S_ij_val
                
            #S_i_jarray_ret = S.at[i].set(jax.vmap(S_ij)(L_arr))
            S_i_jarray_ret = S.at[i].set(jax.vmap(lambda j: nuss_S_ij(K,i,j,S,psq,t,pe_single,pe_pair,min_hairpin))(L_arr))
            return S_i_jarray_ret

        """S[0..L-1,0..L-1]"""
        S_iarray_jarray = jax.lax.fori_loop(0, L, S_i_jarray, S_zero)
        
        S_1_L_val = S_iarray_jarray[0, L-1]
        return S_1_L_val
                      
    return probnuss_inside_jax

def nuss_S_ij(K: int, i: int, j: int, S, psq, t, pe_single, pe_pair, min_hairpin):
    """fills S[i,j]"""
    L      = len(psq)
    L_arr  = jnp.arange(L)
    a_arr  = jnp.arange(K)
    ab_arr = jnp.arange(K*K)

    def S_rule1_a(a):
        """ a(i) S(i+1,j) """
        nonterm_1 = jax.lax.select(i == j, 1.0, S[i+1,j])
        S_rule1_a_val = psq[i,a] * t[0] * pe_single[a] * nonterm_1
        return S_rule1_a_val
    
    def S_rule2_k(k: int):
        def S_rule2_k_ab(ab):
            """ a(i) S(i+1,k-1) b(k) S(k+1,j) """
            nonterm_1 = jax.lax.select(k == i+1, 1.0, S[i+1,k-1])
            nonterm_2 = jax.lax.select(k == j,   1.0, S[k+1,j])
            
            S_rule2_ab_val = psq[i,ab//4] * psq[k,ab%4] * t[1] * pe_pair[ab] * nonterm_1 * nonterm_2
            return S_rule2_ab_val
        
        S_rule2_k_val = jnp.sum(jax.vmap(S_rule2_k_ab)(ab_arr))
        S_rule2_k_val = jax.lax.select(jnp.logical_and(i + min_hairpin < k, k <= j), S_rule2_k_val, 0.0)
        return S_rule2_k_val
    
    S_rule1_val = jnp.sum(jax.vmap(S_rule1_a)(a_arr))
    S_rule2_val = jnp.sum(jax.vmap(S_rule2_k)(L_arr))
    
    #  aS + a S b S 
    S_ij_val = jax.lax.select(j >= i, S_rule1_val + S_rule2_val, 0.0)
    return S_ij_val

def probnuss_inside_standard(log_psq, log_t, e_single, e_pair, K: int=4, min_hairpin: int = 0):
    
    """Standard implementation of probabilistic nussinov"""
    L   = len(log_psq)
    psq = jnn.softmax(log_psq, axis=1)

    t         = jnn.softmax(log_t)     # transition probabilities
    pe_single = jnn.softmax(e_single)  # single emission probabilities
    pe_pair   = jnn.softmax(e_pair)    # pair   emission probabilities

    S = [[0.0 for _ in range(L)] for _ in range(L)]
    for i in range(L-1, -1, -1):
        for j in range(i, L):

            for a in range(K):
                S[i][j] += psq[i][a] * t[0] * pe_single[a] * (1.0 if i == j else S[i+1][j])
                
            for k in range(i+min_hairpin+1, j+1):
                for p in range(K*K):
                    nonterm_1 = 1.0 if k == i+1 else S[i+1][k-1]
                    nonterm_2 = 1.0 if k == j   else S[k+1][j]
                    S[i][j] += psq[i][p//K] * psq[k][p%K] * t[1] * pe_pair[p] * nonterm_1 * nonterm_2

    return float(S[0][L-1])


def test_CACO_ProbNuss_Inside_JAX(L: int, min_hairpin: int = 0):
    K = 4
    epsilon = 1e-3
    
    # paired emission probabilities 4x4 matrix
    pe_pair = np.array([epsilon, epsilon, epsilon, 1.0/6.0,  # AA, AC, AG, AU
                        epsilon, epsilon, 1.0/6.0, epsilon,  # CA, CC, CG, CU
                        epsilon, 1.0/6.0, epsilon, 1.0/6.0,  # GA, GC, GG, GU
                        1.0/6.0, epsilon, 1.0/6.0, epsilon]) # UA, UC, UG, UU
    pe_pair = pe_pair / np.sum(pe_pair)
    e_pair = np.log(pe_pair)

    # unpaired emission probabilities 4x1 matrix
    pe_single = np.array([0.25, 0.25, 0.25, 0.25]) # A, C, G, U
    e_single = np.log(pe_single)

    # transition probabilities (t1, t2, t3)
    t = np.ones(3)# aS, aSa'S, e
    t[0] = 0.55 - epsilon
    t[1] = 0.45 - epsilon
    t[2] = 1.0 - t[0] - t[1]
    log_t = np.log(t)

    print("Prob Nussinov param")
    print("     transistion", t)
    print("     emissions signle", pe_single)
    print("     emissions paired", pe_pair)

    # a random probabilistic sequence"
    log_psq = np.random.normal(size=(L, K))
    psq     = np.exp(log_psq) / np.sum(np.exp(log_psq), axis=1, keepdims=True)
    print("The prob sequence of length", len(psq))
    print(psq)

    probnuss_inside_jax = CACO_ProbNuss_Inside_JAX(K, min_hairpin)
    inside_jax          = probnuss_inside_jax     (log_psq, log_t, e_single, e_pair)
    inside_standard     = probnuss_inside_standard(log_psq, log_t, e_single, e_pair, K, min_hairpin)
    
    assert np.allclose(
        inside_jax, inside_standard, atol=1e-5), f"jax={inside_jax} std={inside_standard}"
    print(f"Test passed", inside_standard, inside_jax)


def main():

    L = 200
    
    print("test Prob Nussinov inside")
    test_CACO_ProbNuss_Inside_JAX(L, min_hairpin=0)
    test_CACO_ProbNuss_Inside_JAX(L, min_hairpin=3)


if __name__ == "__main__":
    main()
