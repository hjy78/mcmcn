# mcmcn
Markov Chain Monte Carlo Sampler

**Description:** Simulates continuous distributions of random vectors using
    Markov chain Monte Carlo (MCMC). Users specify the distribution by
    an R function that evaluates the unnormalized density. Algorithms are random
    walk Metropolis algorithm (function mtrp), independence sampler algorithm
    (function mtrp_exp and mtrp_unif), Metropolis-Hastings algorithm 
    (function mtrp_expu) and Gibbs sampler algorithm (function gibbs_norm and
    gibbs_multinom).
