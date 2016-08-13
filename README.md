# Binary Mixtures of Linear Models

[![Travis-CI Build Status](https://travis-ci.org/.svg?branch=master)](https://travis-ci.org/)

Fit a binary mixture of linear models by Gibbs sampling.  Assumes
responses are drawn two subpopulations, each described by a different
linear model, and it is unknown which responses are drawn from which
subpopulation.  The mixture model estimates the coefficient of the two
subpopulation models and the probability that an observation is drawn
from a particular subpopulation.

## Installing

The package is easily installed from GitHub, using the devtools package.

```R
devtools::install_github("SWotherspoon/bmixlm")
```

If you don't have `devtools` installed already, install it first.

```R
install.packages("devtools")
```

(bmixlm otherwise does not need devtools for normal use.)


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

