% Fun with (Kalman) Filters Part I
% Dominic Steinitz
% 3rd July 2014

---
bibliography: Kalman.bib
---

\newcommand{\condprob} [3] {#1 \left( #2 \,\vert\, #3 \right)}

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
p((z \,\vert\, x) \,\vert\, y) =
p((z \,\vert\, y) \,\vert\, x)
$$

The standard notation in Bayesian statistics is to denote the
parameters of interest as $\theta \in \mathbb{R}^p$ and the
observations as $x \in \mathbb{R}^n$. For reasons that will become
apparent, let us change notation and label the parameters as $x$ and
the observations as $y$.


Let us take a very simple example of a prior $X \sim {\cal{N}}(0,
\sigma^2)$ where $\sigma^2$ is known and then sample from a normal
distribution with mean $x$ and variance for the $i$-th sample $c_i^2$
where $c_i$ is known (normally we would not know the variance but
adding this generality would only clutter the exposition
unnecessarily).

$$
p(y_i \,\vert\, x) = \frac{1}{\sqrt{2\pi c_i^2}}\exp\bigg(\frac{(y_i - x)^2}{2c_i^2}\bigg)
$$

The likelihood is then

$$
p(\boldsymbol{y} \,\vert\, x) = \prod_{i=1}^n \frac{1}{\sqrt{2\pi c_i^2}}\exp\bigg(\frac{(y_i - x)^2}{2c_i^2}\bigg)
$$

As we have already noted, instead of using this with the prior to
calculate the posterior, we can update the prior with each observation
separately. Suppose that we have obtained the posterior given $i - 1$
samples (we do not know this is normally distributed yet but we soon
will):

$$
p(x \,\vert\, y_1,\ldots,y_{i-1}) = {\cal{N}}(\hat{x}_{i-1}, \hat{\sigma}^2_{i-1})
$$

Then we have

$$
\begin{aligned}
p(x \,\vert\, y_1,\ldots,y_{i}) &\propto p(y_i \,\vert\, x)p(x \,\vert\, y_1,\ldots,y_{i-1}) \\
&\propto \exp-\bigg(\frac{(y_i - x)^2}{2c_i^2}\bigg) \exp-\bigg(\frac{(x - \hat{x}_{i-1})^2}{2\hat{\sigma}_{i-1}^2}\bigg) \\
&\propto \exp-\Bigg(\frac{x^2}{c_i^2} - \frac{2xy_i}{c_i^2} + \frac{x^2}{\hat{\sigma}_{i-1}^2} - \frac{2x\hat{x}_{i-1}}{\hat{\sigma}_{i-1}^2}\Bigg) \\
&\propto \exp-\Bigg( x^2\Bigg(\frac{1}{c_i^2} + \frac{1}{\hat{\sigma}_{i-1}^2}\Bigg) - 2x\Bigg(\frac{y_i}{c_i^2} + \frac{\hat{x}_{i-1}}{\hat{\sigma}_{i-1}^2}\Bigg)\Bigg)
\end{aligned}
$$

Writing

$$
\frac{1}{\hat{\sigma}_{i}^2} \triangleq \frac{1}{c_i^2} + \frac{1}{\hat{\sigma}_{i-1}^2}
$$

and then completing the square we also obtain

$$
\frac{\hat{x}_{i}}{\hat{\sigma}_{i}^2} \triangleq \frac{y_i}{c_i^2} + \frac{\hat{x}_{i-1}}{\hat{\sigma}_{i-1}^2}
$$

More Formally
=============

Now let's be a bit more formal about conditional probability and use
the notation of $\sigma$-algebras to define ${\cal{F}}_i =
\sigma\{Y_1,\ldots, Y_i\}$ and $M_i \triangleq \mathbb{E}(X \,\vert\,
{\cal{F}}_i)$ where $Y_i = X + \epsilon_i$, $X$ is as before and
$\epsilon_i \sim {\cal{N}}(0, c_k^2)$. We have previously calculated
that $M_i = \hat{x}_i$ and that ${\cal{E}}((X - M_i)^2 \,\vert\, Y_1,
\ldots Y_i) = \hat{\sigma}_{i}^2$ and the tower law for conditional
probabilities then allows us to conclude ${\cal{E}}((X - M_i)^2) =
\hat{\sigma}_{i}^2$. By [Jensen's
inequality](http://en.wikipedia.org/wiki/Jensen%27s_inequality), we have

$$
{\cal{E}}(M_i^2) = {\cal{E}}({\cal{E}}(X \,\vert\, {\cal{F}}_i)^2)) \leq
{\cal{E}}({\cal{E}}(X^2 \,\vert\, {\cal{F}}_i))) =
{\cal{E}}(X^2) = \sigma^2
$$

Hence $M$ is bounded in $L^2$ and therefore converges in $L^2$ and
almost surely to $M_\infty \triangleq {\cal{E}}(X \,\vert\,
{\cal{F}}_\infty)$.

Bibliography
============