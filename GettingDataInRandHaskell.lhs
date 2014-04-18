Imagine an insect trapped on the face of a clock which wants to visit
each hour an equal number of times. However, there is a snag: it can
only see the value of the hour it is on and the value of the hours
immediately anti-clockwise and immediately clockwise. For example, if
it is standing on 5 then it can see the 5, the 4, and the 6 but no
others.

It can adopt the following strategy: toss a fair coin and move
anti-clockwise for a head and move clockwise for a tail. Intuition
tells us that over a large set of moves the insect will visit each
hour (approximately) the same number of times.

```{.dia width='500'}
import Diagrams
dia = example
```

Now suppose it wants to visit each hour in proportion the value on the
hour. Lacking pen and paper (and indeed opposable thumbs), it decides
to adopt the following strategy: toss a fair coin as in the previous
strategy but only move if the number is larger than the one it is
standing on; if, on the other hand, the number is smaller then choose
a number at random from between 0 and 1 and move if this value is
smaller than the ratio of the proposed hour and the hour on which it
is standing otherwise stay put. For example, if the insect is standing
on 5 and gets a tail then it will move to 6 but if it gets a head then
four fifths of the time it will move to 4 but one fifth of the time it
stays where it is.

```{.dia width='500'}
import Diagrams
dia = example
```

Surprisingly, this strategy produces the desired result and is known
as the Metropolis Algorithm. What the insect has done is to construct
a (discrete) Markov Process which has a limiting distribution (the
stationary distribution) with the desired feature: sampling from this
process will result in each hour being sampled in proportion to its
value.

Let us examine what is happening in a bit more detail.

> import Data.Packed.Matrix
> import Numeric.LinearAlgebra.Algorithms
> import qualified Numeric.Container as M
>
> myMat :: Matrix Double
> myMat = (4 >< 4)
>         [ 0.0,                        0.5,             0.0,                    0.5
>         , 0.5* 1.0 / 2.0, 0.5 * 1.0 / 2.0,             0.5,                    0.0
>         , 0.0,            0.5 * 2.0 / 3.0, 0.5 * 1.0 / 3.0,                    0.5
>         , 0.5 * 1.0 / 4.0,            0.0, 0.5 * 3.0 / 4.0,  0.5 * 3.0 / 4.0 + 0.5 * 1.0 / 4.0
>         ]
> -- myMat = (4 >< 4)
> --         [ 0.5,                        0.5,             0.0,                    0.0
> --         , 0.5* 1.0 / 2.0, 0.5 * 1.0 / 2.0,             0.5,                    0.0
> --         , 0.0,            0.5 * 2.0 / 3.0, 0.5 * 1.0 / 3.0,                    0.5
> --         , 0.0,                        0.0, 0.5 * 3.0 / 4.0,  0.5 + 0.5 * 1.0 / 4.0
> --         ]

    [ghci]
    take 1 $ drop 1000  $ iterate (M.<> myMat) myMat
    ((1 >< 4) [0.1, 0.2, 0.3, 0.4])  M.<> myMat

