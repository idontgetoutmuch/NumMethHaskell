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
worked and the grasshopper is now to be found with equal probability on any
hour. Then at the last jump, the grasshopper must either have been at the
hour before the one it is now on or it must have been at the hour
after the one it is now on. Let us denote the probability that the
grasshopper is on hour $n$ by $\pi(n)$.

$$
\pi'(n) = p(n \, |\, n - 1)\pi(n - 1) + p(n \, |\, n + 1)\pi(n + 1)
$$

Substituting in

$$
\pi'(n) = \frac{1}{2}\frac{1}{N} + \frac{1}{2}\frac{1}{N} = \frac{1}{N}
$$

This tells us that the required distribution is a fixed point of the
grasshopper's strategy. But does the strategy actually converge to the
fixed point?


First we import some modules from
[hmatrix](http://hackage.haskell.org/package/hmatrix).

> import Data.Packed.Matrix
> import Numeric.LinearAlgebra.Algorithms
> import Numeric.Container

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

    [ghci]
    take 1 $ drop 1000  $ iterate (<> myMat) myMat
    ((1 >< 4) [0.1, 0.2, 0.3, 0.4]) <> myMat


Surprisingly, this strategy produces the desired result and is known
as the Metropolis Algorithm. What the grasshopper has done is to
construct a (discrete) Markov Process which has a limiting
distribution (the stationary distribution) with the desired feature:
sampling from this process will result in each hour being sampled in
proportion to its value.

Let us examine what is happening in a bit more detail.

Applications to Bayesian Statistics
===================================

Markov Chain Theory
===================

