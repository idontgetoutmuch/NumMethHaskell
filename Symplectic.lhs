% Rao-Blackwellisation
% Dominic Steinitz
% 9th April 2017

---
bibliography: DynSys.bib
---

Introduction
============


Summary
-------

Back in January, a colleague pointed out to me that GHC did not
produce very efficient code for performing floating point *abs*. Below
is a write-up of my notes about hacking on GHC.

But maybe getting GHC to produce high performance numerical code is
"swimming uphill". Also below is a comparison of a "state of the art"
numerical language, [Julia](https://julialang.org), and an alternative
Haskell approach, using a Haskell domain-specific embedded language
(DSEL),
[accelerate](https://hackage.haskell.org/package/accelerate). The
problem is a real one albeit quite simple from a programming point of
view and should therefore be taken with some grains of salt.

A bit more background
---------------------

The reason for improving GHC's (numerical) performance is that I'd
like to have a fast way of performing statistical inference either via
Hamiltonian Monte Carlo or Sequential Monte
Carlo. [Stan](http://mc-stan.org) is my "go to" tool for offline
analysis but its inference engine is fixed and according to [Fast
Forward
Labs](http://blog.fastforwardlabs.com/2017/01/18/new-research-on-probabilistic-programming.html):
"it's very awkward to productize". What one would really like is the
ability to choose and compose inference methods like can be done in
[monad-bayes](https://github.com/adscib/monad-bayes). Although I
haven't compared monad-bayes and Stan in very much depth, the latter
seems faster in my very limited experience.

The view of one of the folks that's worked hard on making Haskell
a great tool for doing numerical applications is that trying to get
GHC to produce llvm on which tools like [polly](https://polly.llvm.org)
can work is "swimming up hill".

A lot of people who work on numerical problems use
[matlab](https://en.wikipedia.org/wiki/MATLAB) but it seems like if
you want to put the code into production then you have to re-write
it. [Julia](https://julialang.org) is an attempt to improve this
situation and also provide higher performance.

A sort of manifesto
-------------------


Symplectic Integrators
======================

> {-# OPTIONS_GHC -Wall                   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults #-}
>
> {-# LANGUAGE FlexibleContexts #-}
> {-# LANGUAGE TypeOperators    #-}

> module Symplectic (
>     runSteps
>   , runSteps'
>   , reallyRunSteps'
>   , inits
>   , h
>   , bigH1
>   , hs
>   , nabla1
>   , nabla2
>   , bigH2BodyH98
>   ) where
>
> import Data.Number.Symbolic
> import Numeric.AD
> import Prelude                            as P
> import Data.Array.Accelerate              as A   hiding ((^))
> import Data.Array.Accelerate.LLVM.Native  as CPU
> import Data.Array.Accelerate.Linear              hiding (trace)
> import Data.Array.Accelerate.Control.Lens
> import qualified Linear                   as L


The [Störmer-Verlet
scheme](http://www.unige.ch/~hairer/poly_geoint/week2.pdf) is an
implicit symplectic method of order 2.

$$
\begin{aligned}
p_{n+1 / 2} &= p_n         - \frac{h}{2}\nabla_q H(p_{n+1 / 2}, q_n) \\
q_{n+1}     &= q_n         + \frac{h}{2}(\nabla_p H(p_{n+1 / 2}, q_n) + \nabla_p H(p_{n+1 / 2}, q_{n+1}) \\
p_{n+1}     &= p_{n+1 / 2} - \frac{h}{2}\nabla_q H(p_{n+1 / 2}, q_{n+1})
\end{aligned}
$$

Let's assume that the Hamiltonian is of the form $H(p, q) = T(q) +
V(p)$ then this becomes an explicit scheme.

$$
\begin{aligned}
p_{n+1 / 2} &= p_n         - \frac{h}{2}\nabla_q T(q_n) \\
q_{n+1}     &= q_n         + \frac{h}{2}(\nabla_p V(p_{n+1 / 2}) + \nabla_q V(p_{n+1 / 2})) \\
p_{n+1}     &= p_{n+1 / 2} - \frac{h}{2}\nabla_q T(q_{n+1})
\end{aligned}
$$

Simple Pendulum
---------------

Consider the following Hamiltonian for the pendulum:

$$
H(q,p) = \frac{1}{2}p^2 - \cos q
$$

> bigH1 :: P.Floating a => a -> a -> a
> bigH1 q p = p^2 / 2 - cos q

> ham1 :: P.Floating a => [a] -> a
> ham1 [qq1, pp1] = 0.5 * pp1^2 - cos qq1
> ham1 _ = error "Hamiltonian defined for 2 generalized co-ordinates only"

Although it is trivial to find the derivatives in this case, let us
check using automatic symbolic differentiation

> q1, q2, p1, p2 :: Sym a
> q1 = var "q1"
> q2 = var "q2"
> p1 = var "p1"
> p2 = var "p2"

> nabla1 :: [[Sym Double]]
> nabla1 = jacobian ((\x -> [x]) . ham1) [q1, p1]

    [ghci]
    nabla1

which after a bit of simplification gives us $\nabla_q T(q) = \sin q$
and $\nabla_p V(p) = p$ or

> nablaQ, nablaP :: P.Floating a => a -> a
> nablaQ = sin
> nablaP = id

One step of the Störmer-Verlet

> oneStep :: P.Floating a => (a -> a) ->
>                                  (a -> a) ->
>                                  a ->
>                                  (a, a) ->
>                                  (a, a)
> oneStep nabalQQ nablaPP hh (qPrev, pPrev) = (qNew, pNew)
>   where
>     h2 = hh / 2
>     pp2 = pPrev - h2 * nabalQQ qPrev
>     qNew = qPrev + hh * nablaPP pp2
>     pNew = pp2 - h2 * nabalQQ qNew

> h :: Double
> h = 0.01

> hs :: (Double, Double) -> [(Double, Double)]
> hs = P.iterate (oneStep nablaQ nablaP h)

We can plot the result in phase space

FIXME: Plot here!

or we can check that the energy is conserved directly

    [ghci]
    P.map (P.uncurry bigH1) $ P.take 5 $               hs (pi/4, 0.0)
    P.map (P.uncurry bigH1) $ P.take 5 $ P.drop 1000 $ hs (pi/4, 0.0)


Two body problem
----------------

Newton's equations of motions for the two body problem are

$$
\ddot{q}_1 = -\frac{x}{(x^2 + y^2)^{3/2}}, \quad
\ddot{q}_2 = -\frac{y}{(x^2 + y^2)^{3/2}}
$$

And we can re-write those to use Hamilton's equations with this Hamiltonian

$$
H(p_1,p_2,q_1,q_2) = \frac{1}{2}(p_1^2 +p_2^2) - \frac{1}{\sqrt{q_1^2 + q_2^2}}
$$

The advantage of using this small example is that we know the exact
solution should we need it.

> ham2 :: P.Floating a => [a] -> a
> ham2 [qq1, qq2, pp1, pp2] = 0.5 * (pp1^2 + pp2^2) - recip (sqrt (qq1^2 + qq2^2))
> ham2 _ = error "Hamiltonian defined for 4 generalized co-ordinates only"

Again we can calculate the derivatives

> nabla2 :: [[Sym Double]]
> nabla2 = jacobian ((\x -> [x]) . ham2) [q1, q2, p1, p2]

    [ghci]
    (P.mapM_ . P.mapM_) putStrLn . P.map (P.map show) $ nabla2

which after some simplification becomes

$$
\begin{matrix}
\frac{q_1}{(q_1^2 + q_2^2)^{3/2}} \\
\frac{q_2}{(q_1^2 + q_2^2)^{3/2}} \\
p_1 \\
p_2
\end{matrix}
$$

Here's one step of Störmer-Verlet using Haskell 98.

> oneStepH98 :: Double -> V2 (V2 Double) -> V2 (V2 Double)
> oneStepH98 hh prev = V2 qNew pNew
>   where
>     h2 = hh / 2
>     hhs = V2 hh hh
>     hh2s = V2 h2 h2
>     pp2 = psPrev - hh2s * nablaQ' qsPrev
>     qNew = qsPrev + hhs * nablaP' pp2
>     pNew = pp2 - hh2s * nablaQ' qNew
>     qsPrev = prev ^. L._x
>     psPrev = prev ^. L._y
>     nablaQ' qs = V2 (qq1 / r) (qq2 / r)
>       where
>         qq1 = qs ^. L._x
>         qq2 = qs ^. L._y
>         r   = (qq1 ^ 2 + qq2 ^ 2) ** (3/2)
>     nablaP' ps = ps

And here is the same thing using accelerate.

> oneStep2 :: Double -> Exp (V2 (V2 Double)) -> Exp (V2 (V2 Double))
> oneStep2 hh prev = lift $ V2 qNew pNew
>   where
>     h2 = hh / 2
>     hhs = lift $ V2 hh hh
>     hh2s = lift $ V2 h2 h2
>     pp2 = psPrev - hh2s * nablaQ' qsPrev
>     qNew = qsPrev + hhs * nablaP' pp2
>     pNew = pp2 - hh2s * nablaQ' qNew
>     qsPrev :: Exp (V2 Double)
>     qsPrev = prev ^. _x
>     psPrev = prev ^. _y
>     nablaQ' :: Exp (V2 Double) -> Exp (V2 Double)
>     nablaQ' qs = lift (V2 (qq1 / r) (qq2 / r))
>       where
>         qq1 = qs ^. _x
>         qq2 = qs ^. _y
>         r   = (qq1 ^ 2 + qq2 ^ 2) ** (3/2)
>     nablaP' :: Exp (V2 Double) -> Exp (V2 Double)
>     nablaP' ps = ps

With initial values below the solution is an ellipse with eccentricity
$e$.

> e, q10, q20, p10, p20 :: Double
> e = 0.6
> q10 = 1 - e
> q20 = 0.0
> p10 = 0.0
> p20 = sqrt ((1 + e) / (1 - e))

> initsH98 :: V2 (V2 Double)
> initsH98 = V2 (V2 q10 q20) (V2 p10 p20)
>
> inits :: Exp (V2 (V2 Double))
> inits = lift initsH98

We can either keep all the steps of the simulation using accelerate
and the CPU

> nSteps :: Int
> nSteps = 10000
>
> dummyInputs :: Acc (Array DIM1 (V2 (V2 Double)))
> dummyInputs = A.use $ A.fromList (Z :. nSteps) $
>                P.replicate nSteps (V2 (V2 0.0 0.0) (V2 0.0 0.0))

> runSteps :: Array DIM1 (V2 (V2 Double))
> runSteps = CPU.run $ A.scanl (\s _x -> (oneStep2 h s)) inits dummyInputs

Or we can do the same in plain Haskell

> dummyInputsH98 :: [V2 (V2 Double)]
> dummyInputsH98 = P.replicate nSteps (V2 (V2 0.0 0.0) (V2 0.0 0.0))

> runStepsH98 :: [V2 (V2 Double)]
> runStepsH98= P.scanl (\s _x -> (oneStepH98 h s)) initsH98 dummyInputsH98

![](diagrams/symplectic.png)

We'd like to measure performance and running the above for many steps
might use up all available memory. Let's confine ourselves to looking
at the final result.

> runSteps' :: Int -> Exp (V2 (V2 Double)) -> Exp (V2 (V2 Double))
> runSteps' nSteps = A.iterate (lift nSteps) (oneStep2 h)

> reallyRunSteps' :: Int -> Array DIM1 (V2 (V2 Double))
> reallyRunSteps' nSteps = CPU.run $
>                          A.scanl (\s _x -> runSteps' nSteps s) inits
>                          (A.use $ A.fromList (Z :. 1) [V2 (V2 0.0 0.0) (V2 0.0 0.0)])

Performance
===========

Accelerate's LLVM
-----------------

Let's see what accelerate generates with

~~~~ {.haskell include="RunAccGPU.hs"}
~~~~

It's a bit verbose but we can look at the key "loop": while5.

~~~~ {.llvm .numberLines include="RunAccGPU.ll"}
~~~~

And then we can run it for $10^8$ steps and see how long it takes.

~~~~ {include="TimeAccGPU.txt"}
~~~~

Julia's LLVM
------------

Let's try the same problem in Julia.

~~~~ {.julia include="JuliaCPU.jl"}
~~~~

Again we can see how long it takes

~~~~ {include="TimeJulGPU.txt"}
~~~~

> bigH2BodyH98 :: (V2 Double, V2 Double) -> Double
> bigH2BodyH98 x = ke + pe
>   where
>     pe = let V2 q1' q2' = P.fst x in negate $ recip (sqrt (q1'^2 + q2'^2))
>     ke = let V2 p1' p2' = P.snd x in 0.5 * (p1'^2 + p2'^2)


> bigH2BodyH98' :: (V2 Double, V2 Double) -> Double
> bigH2BodyH98' x = ke + pe
>   where
>     pe = let V2 q1' q2' = P.fst x in negate $ recip (sqrt (q1'^2 + q2'^2))
>     ke = let V2 p1' p2' = P.snd x in 0.5 * (p1'^2 + p2'^2)

> main :: IO ()
> main = do
>   putStrLn $ show $ runSteps
>   putStrLn $ show $ runStepsH98
