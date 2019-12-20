# mcmcn

Markov Chain Monte Carlo Sampler

Package: mcmcn

Type: Package

Title: Markov Chain Monte Carlo Sampler

Version: 1.0.0

Author@R: c(person("Jiayi", "Huang", , email = "hjy0708@mail.ustc.edu.cn", role = c("aut", "cre")),
            person("Yuecheng", "Zhang", , email = "an1999@mail.ustc.edu.cn", role = c("aut"))
            
Description: Simulates continuous distributions of random vectors using
    Markov chain Monte Carlo (MCMC). Users specify the distribution by
    an R function that evaluates the unnormalized density. Algorithms are random
    walk Metropolis algorithm (function mtrp), independence sampler algorithm
    (function mtrp_exp and mtrp_unif), Metropolis-Hastings algorithm 
    (function mtrp_expu) and Gibbs sampler algorithm (function gibbs_norm and
    gibbs_multinom).
    
Depends: R (>= 3.6.1)

License: GPL-2 | file LICENSE

Encoding: UTF-8

LazyData: true

RoxygenNote: 6.1.1

URL: http://github.com/hjy78/mcmcn

BugReports: http://github.com/hjy78/mcmcn/issues

Suggests: knitr, rmarkdown

VignetteBuilder: knitr
