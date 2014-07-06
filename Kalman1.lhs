% Fun with Filters
% Dominic Steinitz
% 3rd July 2014

Suppose we wish to estimate the mean of a sample drawn from a normal
distribution. In the Bayesian approach, we know the prior distribution
for the mean (it could be a non-informative prior) and then we update
this with our observations to create the posterior, the latter giving
us improved information about the distribution of the mean. In symbols

$$
p(\theta \,\vert\, x) \propto p(x \,\vert\, \theta)p(\theta)
$$

Typically, the samples are chosen to be independent, and all of the
data is used to perform the update but, given independence, there is
no particular reason to do that, updates can performed one at a time
and the result is the same; nor is the order of update
important. Being a bit imprecise, we have

$$
p(z \,\vert\, x, y) = p(z, x, y)p(x, y) = p(z, x, y)p(x)p(y) =
p((z \,\vert\, x) \,\vert\, y)
p((z \,\vert\, y) \,\vert\, x)
$$

The standard notation in Bayesian statistics is to denote the
parameters of interest as $\theta \in \mathbb{R}^p$ and the
observations as $x \in \mathbb{R}^n$. For reasons that will become
apparent, let us change notation and label the parameters as $x$ and
the observations as $y$.



