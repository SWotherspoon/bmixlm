# Binary Mixtures of Linear Models

Fit a binary mixture of linear models by Gibbs sampling.  Assumes
responses are drawn two subpopulations, each described by a different
linear model, and it is unknown which responses are drawn from which
subpopulation.  The mixture model estimates the coefficient of the two
subpopulation models and the probability that an observation is drawn
from a particular subpopulation.

## Installing

The current version of bmixlm can be installed from GitHub using the 
remotes package. 
```r
# install.packages("remotes")
remotes::install_github("SWotherspoon/bmixlm")
```


## TODO

* **Multiple Chains**.  Ideally the sampler would draw multiple
  chains.  This would require rewriting the summary facilities to cope
  with multiple chains.  This is a low priority.

* **Parallelisation**.  At this point, the sampler is only capable of
  utilizing a single core on a multicore machine.  It would be
  relatively simple to introduce coarse grain parallelism by having
  the sampler draw multiple chains in parallel, using something like
  the multicore facility in the parallel package.  Unfortunately, at
  the time of writing there does not seem to be a good parallization
  solution that works equally well on all platforms.  This is a low
  priority.

