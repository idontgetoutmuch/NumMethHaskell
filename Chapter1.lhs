% Haskell & Gibbs Sampling
% Dominic Steinitz
% 15th June 2014

This is really intended as a draft chapter for our book. Given the
diverse natures of the intended intended audiences, it is probably a
bit light on explanation of the Haskell (use of monad transformers)
for those with a background in numerical methods. It is hoped that the
explanation of the mathematics is adequate for those with a background
in Haskell but not necessarily in numerical methods. As always, any
feedback is gratefully accepted.

Introduction
============

Imagine an insect, a grasshopper, trapped on the face of a clock which
wants to visit each hour an equal number of times. However, there is a
snag: it can only see the value of the hour it is on and the value of
the hours immediately anti-clockwise and immediately clockwise. For
example, if it is standing on 5 then it can see the 5, the 4, and the
6 but no others.

It can adopt the following strategy: toss a fair coin and move
anti-clockwise for a head and move clockwise for a tail. Intuition
tells us that over a large set of moves the grasshopper will visit each
hour (approximately) the same number of times.

```{.dia width='500'}
import Diagrams
dia = example 12
```

Can we confirm our intuition somehow? Suppose that the strategy has
worked and the grasshopper is now to be found with equal probability
on any hour. Then at the last jump, the grasshopper must either have
been at the hour before the one it is now on or it must have been at
the hour after the one it is now on. Let us denote the probability
that the grasshopper is on hour $n$ by $\pi(n)$ and the (conditional)
probability that the grasshopper jumps to state $n$ given it was in
state $m$ by $p(n \, |\, m)$. Then we have

$$
\pi'(n) = p(n \, |\, n - 1)\pi(n - 1) + p(n \, |\, n + 1)\pi(n + 1)
$$

Substituting in where $N$ is a normalising constant (12 in this case)
we obtain

$$
\pi'(n) = \frac{1}{2}\frac{1}{N} + \frac{1}{2}\frac{1}{N} = \frac{1}{N}
$$

This tells us that the required distribution is a fixed point of the
grasshopper's strategy. But does the strategy actually converge to the
fixed point? Let us perform an experiment.

First we import some modules from
[hmatrix](http://hackage.haskell.org/package/hmatrix).

> {-# LANGUAGE FlexibleContexts #-}

> module Chapter1 where

> import Data.Packed.Matrix
> import Numeric.LinearAlgebra.Algorithms
> import Numeric.Container

> import Data.Random
> import Control.Monad.State
> import qualified Control.Monad.Writer as W
> import qualified Control.Monad.Loops as ML
> import Data.Random.Source.PureMT

Let us use a clock with 5 hours to make the matrices sufficiently
small to fit on one page.

```{.dia width='500'}
import Diagrams
dia = example 5
```

Here is the strategy encoded as a matrix. For example the first row
says jump to position 1 with probablity 0.5 or jump to position 5 with
probability 0.5.

> eqProbsMat :: Matrix Double
> eqProbsMat = (5 >< 5)
>         [ 0.0, 0.5, 0.0, 0.0, 0.5
>         , 0.5, 0.0, 0.5, 0.0, 0.0
>         , 0.0, 0.5, 0.0, 0.5, 0.0
>         , 0.0, 0.0, 0.5, 0.0, 0.5
>         , 0.5, 0.0, 0.0, 0.5, 0.0
>         ]

We suppose the grasshopper starts at 1 o'clock.

> startOnOne :: Matrix Double
> startOnOne = ((1 >< 5) [1.0, 0.0, 0.0, 0.0, 0.0])

If we allow the grasshopper to hop 1000 times then we see that it is
equally likely to be found on any hour hand with a 20% probability.

    [ghci]
    startOnOne
    eqProbsMat
    take 1 $ drop 1000  $ iterate (<> eqProbsMat) startOnOne

In this particular case, the strategy does indeed converge.

Now suppose the grasshopper wants to visit each hour in proportion the
value of the number on the hour. Lacking pen and paper (and indeed
opposable thumbs), it decides to adopt the following strategy: toss a
fair coin as in the previous strategy but only move if the number is
larger than the one it is standing on; if, on the other hand, the
number is smaller then choose a number at random from between 0 and 1
and move if this value is smaller than the ratio of the proposed hour
and the hour on which it is standing otherwise stay put. For example,
if the grasshopper is standing on 5 and gets a tail then it will move
to 6 but if it gets a head then four fifths of the time it will move
to 4 but one fifth of the time it will stay where it is.

```{.dia width='500'}
import Diagrams
dia = clockSelf 12
```

Suppose that the strategy has worked (it is not clear that is
has) and the grasshopper is now to be found at 12 o'clock 12 times as
often as at 1 o'clock, at 11 o'clock 11 times as often as at 1
o'clock, etc. Then at the last jump, the grasshopper must either have
been at the hour before the one it is now on, the hour after the one
it is now on or the same hour it is now on. Let us denote the
probability that the grasshopper is on hour $n$ by $\pi(n)$.

$$
\pi'(n) = p(n \, |\, n - 1)\pi(n - 1) + p(n \, |\, n)\pi(n) + p(n \, |\, n + 1)\pi(n + 1)
$$

Substituting in at 4 say

$$
\begin{aligned}
\pi'(4) &= \frac{1}{2}\pi(3) +
           \frac{1}{2}\frac{1}{4}\pi(4) +
           \frac{1}{2}\frac{4}{5}\pi(5) \\
        &= \frac{1}{2}\bigg(\frac{3}{N} + \frac{1}{4}\frac{4}{N} + \frac{4}{5}\frac{5}{N}\bigg) \\
        &= \frac{1}{N}\frac{8}{2} \\
        &= \frac{4}{N} \\
        &= \pi(4)
\end{aligned}
$$

The reader can check that this relationship holds for all other
hours. This tells us that the required distribution is a fixed
point of the grasshopper's strategy. But does this strategy actually
converge to the fixed point?

Again, let us use a clock with 5 hours to make the matrices
sufficiently small to fit on one page.

```{.dia width='500'}
import Diagrams
dia = clockSelf 5
```

Here is the strategy encoded as a matrix. For example the first row
says jump to position 1 with probablity 0.5 or jump to position 5 with
probability 0.5.

> incProbsMat :: Matrix Double
> incProbsMat = scale 0.5 $
>   (5 >< 5)
>     [ 0.0,         1.0,     0.0,        0.0, 1.0
>     , 1.0/2.0, 1.0/2.0,     1.0,        0.0, 0.0
>     , 0.0,     2.0/3.0, 1.0/3.0,        1.0, 0.0
>     , 0.0,         0.0, 3.0/4.0,    1.0/4.0, 1.0
>     , 1.0/5.0,     0.0,     0.0,    4.0/5.0, 1.0/5.0 + 4.0/5.0
>     ]

We suppose the grasshopper starts at 1 o'clock.

If we allow the grasshopper to hop 1000 times then we see that it is
equally likely to be found on any hour hand $n$ with a probability of
$n$ times the probability of being found on 1.

    [ghci]
    incProbsMat
    take 1 $ drop 1000  $ iterate (<> incProbsMat) startOnOne

In this particular case, the strategy does indeed converge.

Surprisingly, this strategy produces the desired result and is known
as the Metropolis Algorithm. What the grasshopper has done is to
construct a (discrete) Markov Process which has a limiting
distribution (the stationary distribution) with the desired feature:
sampling from this process will result in each hour being sampled in
proportion to its value.

Markov Chain Theory
===================

Let us examine what is happening in a bit more detail.

The grasshopper has started with a very simple Markov Chain: one which jumps
clockwise or anti-clockwise with equal probability and then modified
it. But what is a Markov Chain?

A time homogeneous **Markov chain** is a countable sequence of random variables $X_0, X_1,
\ldots$ such that

$$
\mathbb{P} (X_{n+1} = j \,|\, X_0 = i_0, X_1 = i_1, \dots X_n = i) = \mathbb{P} (X_{n+1} = j \,|\, X_n = i)
$$

We sometimes say that a Markov Chain is discrete time stochastic
process with the above property.

So the very simple Markov Chain can be described by

$$
q(i, j) =
\begin{cases}
\mathbb{P} (X_{n+1} = j \,|\, X_n = i) = \frac{1}{2} & \text{if } j = i + 1 \mod N \\
\mathbb{P} (X_{n+1} = j \,|\, X_n = i) = \frac{1}{2} & \text{if } j = i - 1 \mod N \\
\mathbb{P} (X_{n+1} = j \,|\, X_n = i) = 0 & \text{otherwise }
\end{cases}
$$

The grasshopper knows that $\pi(i) = i/N$ so it can calculate $\pi(j)/\pi(i) = j/i$ without knowing $N$. This is important because now, without knowing $N$, the grasshopper can evaluate

$$
p(i, j) =
\begin{cases}
q(i,j)\bigg[\frac{\pi(j) q(j,i)}{\pi(i) q(i,j)} \land 1 \bigg] & \text{if } j \ne i \\
1 - \sum_{k : k \ne i} q(i,k) \bigg[\frac{\pi(k) q(k,i)}{\pi(i) q(i,k)} \land 1 \bigg] & \text{if } j = i
\end{cases}
$$

where $\land$ takes the maximum of its arguments. Simplifying the above by substituing in the grasshopper's probabilities and noting that $j = i \pm 1 \mod N$ is somewhat obscure way of saying jump clockwise or anti-clockwise we obtain

$$
q(i, j) =
\begin{cases}
\frac{1}{2} (\frac{j}{i} \land 1) & \text{if } j \text{ is 1 step clockwise} \\
\frac{1}{2} (\frac{j}{i} \land 1) & \text{if } j \text{ is 1 step anti-clockwise} \\
1 - \frac{1}{2}(\frac{j^c}{i} \land 1) - \frac{1}{2}(\frac{j^a}{i} \land 1) & \text{if } j = i \text{ and } j^c \text{ is one step clockwise and } j^a \text{ is one step anti-clockwise} \\
0 & \text{otherwise}
\end{cases}
$$

The Ergodic Theorem
-------------------

In most studies of Markov chains, one is interested in whether a chain
has a stationary distribution. What we wish to do is take a
distribution and create a chain with this distribution as its
stationary distribution. We will still need to show that our chain
does indeed have the correct stationary distribution and we state the
relevant theorem somewhat informally and with no proof.

Theorem (L3)
------------

An irreducible, aperiodic and positive recurrent Markov chain has a
unique stationary distribution.

Roughly speaking

* Irreducible means it is possible to get from any state to any other state.

* Aperiodic means that returning to a state having started at that
state occurs at irregular times.

* Positive recurrent means that the first time to hit a state is
finite (for every state and more pedantically except on
sets of null measure).

Note that the last condition is required when the state space is
infinite - see
[Skrikant](http://www.ifp.illinois.edu/~srikant/ECE534/Spring09/DTMC.pdf)'s
lecture notes for an example and also for a more formal definition of
the theorem and its proof.

Algorithm (L3)
--------------

Let $\pi$ be a probability distribution on the state space
$\Omega$ with $\pi(i) \gt 0$ for all $i$ and let $(Q, \pi_0)$ be an
ergodic Markov chain on $\Omega$ with transition probabilities $q(i,j)
\gt 0$ (the latter condition is slightly stronger than it need be but
we will not need fully general conditions).

Create a new (ergodic) Markov chain with transition probabilities

$$
p_{ij} =
\begin{cases}
q(i,j)\bigg[\frac{\pi(j) q(j,i)}{\pi(i) q(i,j)} \land 1 \bigg] & \text{if } j \ne i \\
1 - \sum_{k : k \ne i} q(i,k) \bigg[\frac{\pi(j) q(j,i)}{\pi(i) q(i,j)} \land 1 \bigg] & \text{if } j = i
\end{cases}
$$

where $\land$ takes the maximum of its arguments.

Calculate the value of interest on the state space e.g. the total
magnetization for each step produced by this new chain.

Repeat a sufficiently large number of times and take the average. This
gives the estimate of the value of interest.

Convergence (L3)
----------------

Let us first note that the Markov chain produced by this algorithm
almost trivially satisfies the detailed balance condition, for
example,

$$
\begin{aligned}
\pi(i) q(i,j)\bigg[\frac{\pi(j) q(j, i)}{\pi(i)q(i,j)} \land 1\bigg]
&= \pi(i)q(i,j) \land \pi(j)q(j,i) \\
&= \pi(j)q(j,i)\bigg[\frac{\pi(i) q(i, j)}{\pi(j)q(j,i)} \land 1\bigg]
\end{aligned}
$$

Secondly since we have specified that $(Q, \pi_0)$ is ergodic then
clearly $(P, \pi_0)$ is also ergodic (all the transition probabilities
are $\gt 0$).

So we know the algorithm will converge to the unique distribution we
specified to provide estimates of values of interest.

Gibbs Sampling
==============

Random Scan
-----------

For simplicity let us consider a model with two parameters and that we
sample from either parameter with equal probability. In this sampler,
We update the parameters in a single step.

$$
\begin{cases}
\text{Sample } \theta_1^{(i+1)} \sim \pi(\theta_1 \,\big|\, \theta_2^{(i)}) & \text{with probability } \frac{1}{2} \\
\text{Sample } \theta_2^{(i+1)} \sim \pi(\theta_2 \,\big|\, \theta_1^{(i)}) & \text{with probability } \frac{1}{2}
\end{cases}
$$

The transition density kernel is then given by

$$
q\big(\boldsymbol{\theta}^{(i+1)}, \boldsymbol{\theta}^{(i)}\big) =
\frac{1}{2}\pi(\theta_1^{(i+1)} \,\big|\, \theta_2^{(i)})\delta({\theta_2^{(i)},\theta_2^{(i+1)}}) +
\frac{1}{2}\pi(\theta_2^{(i+1)} \,\big|\, \theta_1^{(i)})\delta({\theta_1^{(i)},\theta_1^{(i+1)}})
$$

where $\delta$ is the Dirac delta function.

Detailed balance (L3)
---------------------

This sampling scheme satisifies the detailed balance condition. We have

$$
\begin{aligned}
\pi(\theta_1, \theta_2)
\bigg[
\frac{1}{2}\pi(\theta_1' \,\big|\, \theta_2)\delta({\theta_2,\theta_2'}) +
\frac{1}{2}\pi(\theta_2' \,\big|\, \theta_1)\delta({\theta_1,\theta_1'})\bigg] &= \\
\frac{1}{2}\bigg[\pi(\theta_1, \theta_2)
\pi(\theta_1' \,\big|\, \theta_2)\delta({\theta_2,\theta_2'}) +
\pi(\theta_1, \theta_2)
\pi(\theta_2' \,\big|\, \theta_1)\delta({\theta_1,\theta_1'})\bigg] &= \\
\frac{1}{2}\bigg[\pi(\theta_1, \theta_2')
\pi(\theta_1' \,\big|\, \theta_2)\delta({\theta_2,\theta_2'}) +
\pi(\theta_1', \theta_2)
\pi(\theta_2' \,\big|\, \theta_1)\delta({\theta_1,\theta_1'})\bigg] &= \\
\frac{1}{2}\bigg[
\pi(\theta_2')\pi(\theta_1 \,\big|\, \theta_2')
\frac{1}{2}\pi(\theta_1' \,\big|\, \theta_2)\delta({\theta_2,\theta_2'}) +
\pi(\theta_1')\pi(\theta_2 \,\big|\, \theta_1')
\pi(\theta_2' \,\big|\, \theta_1)\delta({\theta_1,\theta_1'})
\bigg] &= \\
\frac{1}{2}\bigg[
\pi(\theta_1', \theta_2')\pi(\theta_1 \,\big|\, \theta_2')
\delta({\theta_2',\theta_2}) +
\pi(\theta_1', \theta_2')\pi(\theta_2 \,\big|\, \theta_1')
\delta({\theta_1',\theta_1})
\bigg] &= \\
\pi(\theta_1', \theta_2')\bigg[
\frac{1}{2}\pi(\theta_1 \,\big|\, \theta_2')
\delta({\theta_2',\theta_2}) +
\frac{1}{2}\pi(\theta_2 \,\big|\, \theta_1')
\delta({\theta_1',\theta_1})
\bigg] &
\end{aligned}
$$

In other words

$$
\pi\big({\boldsymbol{\theta}}\big)q\big(\boldsymbol{\theta}', \boldsymbol{\theta}\big) =
\pi\big({\boldsymbol{\theta'}}\big)q\big(\boldsymbol{\theta}, \boldsymbol{\theta}'\big)
$$

Hand waving slightly, we can see that this scheme satisfies the
premises of the ergodic theorem and so we can conclude that there is a
unique stationary distribution and $\pi$ must be that distribution.

Systematic Scan
===============

Most references on Gibbs sampling do not describe the random scan but
instead something called a systematic scan.

Again for simplicity let us consider a model with two parameters. In this
sampler, we update the parameters in two steps.

$$
\begin{eqnarray}
\text{Sample } \theta_1^{(i+1)} & \sim & \pi(\theta_1 \,\big|\, \theta_2^{(i)}) \\
\text{Sample } \theta_2^{(i+1)} & \sim & \pi(\theta_2 \,\big|\, \theta_1^{(i+1)})
\end{eqnarray}
$$

We observe that this is **not** time-homegeneous; at each step the
transition matrix flips between the two transition matrices given by
the individual steps. Thus although, as we show below, each individual
transtion satisifies the detailed balance condition, we cannot apply
the ergodic theorem as it only applies to time-homogeneous processes.

The transition density kernel is then given by

$$
q\big(\boldsymbol{\theta}^{(i)}, \boldsymbol{\theta}^{(i+1)}\big) =
q_1\big(\boldsymbol{\theta}^{(i)}, \tilde{\boldsymbol{\theta}}\big)
q_2\big(\tilde{\boldsymbol{\theta}}, \boldsymbol{\theta}^{(i+1)}\big)
$$

where $\tilde{\boldsymbol{\theta}} = (\theta_1^{(i+1)}, \theta_2^{(i)})^\top$.

Thus

$$
q\big(\boldsymbol{\theta}, \boldsymbol{\theta}'\big) =
\pi(\theta_1' \,\big|\, \theta_2)
\pi(\theta_2' \,\big|\, \theta_1')
$$

Detailed balance (L3)
---------------------

Suppose that we have two states $\boldsymbol{\theta} = (\theta_1,
\theta_2)^\top$ and $\boldsymbol{\theta}' = (\theta_1',
\theta_2')^\top$ and that $\theta_2 \neq \theta_2'$. Then
$q_1\big(\boldsymbol{\theta}, \boldsymbol{\theta}'\big) =
0$. Trivially we have

$$
\pi\big({\boldsymbol{\theta}}\big)q_1\big(\boldsymbol{\theta}, \boldsymbol{\theta}'\big) =
\pi\big({\boldsymbol{\theta'}}\big)q_1\big(\boldsymbol{\theta}', \boldsymbol{\theta}\big)
$$

Now suppose that $\theta_2 = \theta_2'$

$$
\begin{aligned}
\pi(\theta_1, \theta_2)q_1((\theta_1, \theta_2), (\theta_1', \theta_2)) & =
\pi(\theta_1, \theta_2)\pi(\theta_1' \,\big|\, \theta_2) \\
& = \pi(\theta_1 \,\big|\, \theta_2)\pi(\theta_1', \theta_2) \\
& = \pi(\theta_1 \,\big|\, \theta_2')\pi(\theta_1', \theta_2') \\
& = \pi(\theta_1', \theta_2')q_1((\theta_1', \theta_2), (\theta_1, \theta_2))
\end{aligned}
$$

So again we have

$$
\pi\big({\boldsymbol{\theta}}\big)q_1\big(\boldsymbol{\theta}, \boldsymbol{\theta}'\big) =
\pi\big({\boldsymbol{\theta'}}\big)q_1\big(\boldsymbol{\theta}', \boldsymbol{\theta}\big)
$$

Similarly we can show

$$
\pi\big({\boldsymbol{\theta}}\big)q_2\big(\boldsymbol{\theta}, \boldsymbol{\theta}'\big) =
\pi\big({\boldsymbol{\theta'}}\big)q_2\big(\boldsymbol{\theta}', \boldsymbol{\theta}\big)
$$

But note that

$$
\begin{aligned}
\pi(\theta_1, \theta_2)
q_1((\theta_1, \theta_2), (\theta_1', \theta_2))
q_2((\theta_1', \theta_2), (\theta_1', \theta_2')) & =
\pi(\theta_1, \theta_2)
\pi(\theta_1' \,\big|\, \theta_2)
\pi(\theta_2' \,\big|\, \theta_1') \\
& = \pi(\theta_1', \theta_2) \pi(\theta_1 \,\big|\, \theta_2) \pi(\theta_2' \,\big|\, \theta_1') \\
& = \pi(\theta_1' \,\big|\, \theta_2) \pi(\theta_1 \,\big|\, \theta_2) \pi(\theta_2', \theta_1')
\end{aligned}
$$

whereas

$$
\begin{aligned}
\pi(\theta_1', \theta_2')
q_1((\theta_1', \theta_2'), (\theta_1, \theta_2'))
q_2((\theta_1, \theta_2'), (\theta_1, \theta_2)) & =
\pi(\theta_1', \theta_2')
\pi(\theta_1 \,\big|\, \theta_2')
\pi(\theta_2 \,\big|\, \theta_1) \\
& = \pi(\theta_1, \theta_2') \pi(\theta_1' \,\big|\, \theta_2') \pi(\theta_2 \,\big|\, \theta_1) \\
& = \pi(\theta_2' \,\big|\, \theta_1) \pi(\theta_1' \,\big|\, \theta_2') \pi(\theta_2, \theta_1)
\end{aligned}
$$

and these are not necessarily equal.

So the detailed balance equation is not satisfied, another sign that
we cannot appeal to the ergodic theorem.

An Example: The Bivariate Normal
================================

Let us demonstrate the Gibbs sampler with a distribution which we
actually know: the bivariate normal.

$$
\begin{bmatrix}
\theta_1 \\
\theta_2
\end{bmatrix}
\bigg| y \sim
N \begin{bmatrix}
\begin{bmatrix}
\theta_1 \\
\theta_2
\end{bmatrix} &
\begin{bmatrix}
1 & \rho \\
\rho & 1
\end{bmatrix}
\end{bmatrix}
$$

The conditional distributions are easily calculated to be
$$
\begin{aligned}
\theta_1 \,\vert\, \theta_2, y &\sim {\cal{N}}(y_1 + \rho(\theta_2 - y_2), 1 - \rho^2) \\
\theta_2 \,\vert\, \theta_1, y &\sim {\cal{N}}(y_2 + \rho(\theta_1 - y_1), 1 - \rho^2)
\end{aligned}
$$

Let's take a correlation of 0.8, a data point of (0.0, 0.0) and start
the chain at (2.5, 2.5).

> rho :: Double
> rho = 0.8
>
> y :: (Double, Double)
> y = (0.0, 0.0)
>
> y1, y2 :: Double
> y1 = fst y
> y2 = snd y
>
> initTheta :: (Double, Double)
> initTheta = (2.5, 2.5)

We pre-calculate the variance needed for the sampler.

> var :: Double
> var = 1.0 - rho^2

In Haskell and in the
[random-fu](http://hackage.haskell.org/package/random-fu) package,
sampling from probability distributions is implemented as a monad. We
sample from the relevant normal distributions and keep the trajectory
using a writer monad.

> gibbsSampler :: Double -> RVarT (W.Writer [(Double,Double)]) Double
> gibbsSampler oldTheta2 = do
>   newTheta1 <- rvarT (Normal (y1 + rho * (oldTheta2 - y2)) var)
>   lift $ W.tell [(newTheta1, oldTheta2)]
>   newTheta2 <- rvarT (Normal (y2 + rho * (newTheta1 - y1)) var)
>   lift $ W.tell [(newTheta1, newTheta2)]
>   return $ newTheta2

It is common to allow the chain to "burn in" so as to "forget" its
starting position. We arbitrarily burn in for 10,000 steps.

> burnIn :: Int
> burnIn = 10000

We sample repeatedly from the sampler using the monadic form of
iterate. Running the monadic stack is slightly noisy but nonetheless
straightforward. We use
[mersenne-random-pure64](http://hackage.haskell.org/package/mersenne-random-pure64)
(albeit indirectly via
[random-source](http://hackage.haskell.org/package/random-source))
as our source of entropy.

> runMCMC :: Int -> [(Double, Double)]
> runMCMC n =
>   take n $
>   drop burnIn $
>   snd $
>   W.runWriter (evalStateT (sample (ML.iterateM_ gibbsSampler (snd initTheta))) (pureMT 2))

We can look at the trajectory of our sampler for various run lengths.

```{.dia width='800'}
import Chapter1
import Chart
dia = (diagLine (runMCMC 16) 16 ||| diagLine (runMCMC 32) 32)
      ===
      (diagLine (runMCMC 64) 64 ||| diagLine (runMCMC 128) 128)
```

For bigger sample sizes, plotting the distribution sampled re-assures
us that we are indeed sampling from a bivariate normal distribution as
the theory predicted.

```{.dia width='800'}
import Chapter1
import Chart
dia = (diagPoint (runMCMC 128) 128 ||| diagPoint (runMCMC 512) 512)
      ===
      (diagPoint (runMCMC 2048) 2048 ||| diagPoint (runMCMC 8192) 8192)
```


Applications to Bayesian Statistics
===================================

Some of what is
[here](http://idontgetoutmuch.wordpress.com/2014/04/09/gibbs-sampling-in-r-haskell-jags-and-stan/)
and
[here](http://idontgetoutmuch.wordpress.com/2014/03/20/bayesian-analysis-a-conjugate-prior-and-markov-chain-monte-carlo/)
excluding JAGS and STAN (after all this is a book about Haskell).

Applications to Physics

Applications to Physics
=======================

Most of what is
[here](http://idontgetoutmuch.wordpress.com/2013/12/07/haskell-ising-markov-metropolis/).




