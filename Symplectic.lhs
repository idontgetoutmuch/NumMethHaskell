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


Symplectic Integrator
=====================

> {-# OPTIONS_GHC -Wall                   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults #-}
>
> {-# LANGUAGE FlexibleContexts #-}
> {-# LANGUAGE TypeOperators    #-}

> import Data.Number.Symbolic
> import Numeric.AD
> import Prelude                            as P
> import Data.Array.Accelerate              as A   hiding ((^), iterate)
> import Data.Array.Accelerate.LLVM.Native  as CPU
> import Data.Array.Accelerate.Linear              hiding (trace)
> import Data.Array.Accelerate.Control.Lens
> import qualified Linear                   as L
>
> import Debug.Trace

Let us consider the following Hamiltonian

> bigH :: P.Floating a => a -> a -> a
> bigH p q = p^2 / 2 - cos q

Symplectic Euler

$$p_{n+1} = p_n - h\nabla_q H(p_{n+1}, q_n)$$

$$q_{n+1} = q_n + h\nabla_p H(p_{n+1}, q_n)$$

If e.g.

$$H(p, q) = \frac{1}{2}p^2 - \cos q$$

then the scheme is explicit.

Let's assume that the Hamiltonian is of the form $H(p, q) = T(q) +
V(p)$ then the midpoint rule

$$
y_{n+1} = y_n + hJ^{-1}\nabla H(\frac{y_{n+1} + y_n}{2})
$$

which we can re-write as

$$
\begin{aligned}
p_{n+1} &= p_n - h\nabla_q H(\frac{p_{n+1} + p_n}{2}, \frac{q_{n+1} + q_n}{2}) \\
q_{n+1} &= q_n + h\nabla_p H(\frac{p_{n+1} + p_n}{2}, \frac{q_{n+1} + q_n}{2})
\end{aligned}
$$

becomes

$$
\begin{aligned}
p_{n+1} &= p_n - h\nabla_q T(\frac{q_{n+1} + q_n}{2}) \\
q_{n+1} &= q_n + h\nabla_p V(\frac{p_{n+1} + p_n}{2})
\end{aligned}
$$

Störmer-Verlet

$$
\begin{aligned}
p_{n+1 / 2} &= p_n         - \frac{h}{2}\nabla_q H(p_{n+1 / 2}, q_n) \\
q_{n+1}     &= q_n         + \frac{h}{2}(\nabla_p H(p_{n+1 / 2}, q_n) + \nabla_p H(p_{n+1 / 2}, q_{n+1}) \\
p_{n+1}     &= p_{n+1 / 2} - \frac{h}{2}\nabla_q H(p_{n+1 / 2}, q_{n+1})
\end{aligned}
$$

$$
\begin{aligned}
p_{n+1 / 2} &= p_n         - \frac{h}{2}\nabla_q T(q_n) \\
q_{n+1}     &= q_n         + \frac{h}{2}(\nabla_p V(p_{n+1 / 2}) + \nabla_q V(p_{n+1 / 2}) \\
p_{n+1}     &= p_{n+1 / 2} - \frac{h}{2}\nabla_q T(q_{n+1})
\end{aligned}
$$

For our example we have $\nabla_q T(q) = \sin q$ and $\nabla_p V(p) = p$

> nablaQ, nablaP :: P.Floating a => a -> a
> nablaQ = sin
> nablaP = id

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
> hs = iterate (oneStep nablaQ nablaP h)

Two body problem

$$
\ddot{q}_1 = -\frac{q_1}{(q_1^2 + q_2^2)^{3/2}}, \quad
\ddot{q}_2 = -\frac{q_2}{(q_1^2 + q_2^2)^{3/2}}
$$

Hamiltonian

$$
H(p_1,p_2,q_1,q_2) = \frac{1}{2}(p_1^2 +p_2^2) - \frac{1}{\sqrt{q_1^2 + q_2^2}}
$$

> e, q10, q20, p10, p20 :: Double
> e = 0.6
> q10 = 1 - e
> q20 = 0.0
> p10 = 0.0
> p20 = sqrt ((1 + e) / (1 - e))
>
> ham :: P.Floating a => [a] -> a
> ham [pp1, pp2, qq1, qq2] = 0.5 * (pp1^2 + pp2^2) - recip (sqrt (qq1^2 + qq2^2))
> ham _ = error "Hamiltonian defined for 4 generalized co-ordinates only"

> q1, q2, p1, p2 :: Sym a
> q1 = var "q1"
> q2 = var "q2"
> p1 = var "p1"
> p2 = var "p2"

> testHam :: P.Floating a => [a] -> [[a]]
> testHam xs = jacobian ((\x -> [x]) . ham) xs

$$
\begin{matrix}
0.5*p1+0.5*p1 \\
0.5*p2+0.5*p2 \\
q1/(2.0*sqrt (q1*q1+q2*q2))/sqrt (q1*q1+q2*q2)/sqrt (q1*q1+q2*q2)+q1/(2.0*sqrt (q1*q1+q2*q2))/sqrt (q1*q1+q2*q2)/sqrt (q1*q1+q2*q2) \\
q2/(2.0*sqrt (q1*q1+q2*q2))/sqrt (q1*q1+q2*q2)/sqrt (q1*q1+q2*q2)+q2/(2.0*sqrt (q1*q1+q2*q2))/sqrt (q1*q1+q2*q2)/sqrt (q1*q1+q2*q2)
\end{matrix}
$$

$$
\begin{matrix}
p_1 \\
p_2 \\
\frac{q_1}{(q_1^2 + q_2^2)^{3/2}} \\
\frac{q_2}{(q_1^2 + q_2^2)^{3/2}}
\end{matrix}
$$

> oneStep2 :: Double -> Exp (V2 Double, V2 Double) -> Exp (V2 Double, V2 Double)
> oneStep2 h' prev = lift (qNew, pNew)
>   where
>     h2 = h' / 2
>     gs :: Exp (V2 Double)
>     gs = lift ((pure h') :: V2 Double)
>     hs' :: Exp (V2 Double)
>     hs' = (lift ((pure h2) :: V2 Double))
>     p2' = psPrev - hs' * nablaQ' qsPrev
>     qNew = qsPrev + gs * nablaP' p2'
>     pNew = p2' - hs' * nablaQ' qNew
>     qsPrev = A.fst prev
>     psPrev = A.snd prev
>     nablaQ' :: Exp (V2 Double) -> Exp (V2 Double)
>     nablaQ' qs = lift (V2 (q1' / r) (q2' / r))
>       where
>         q1' = qs ^. _x
>         q2' = qs ^. _y
>         r   = (q1' ^ 2 + q2' ^ 2) ** (3/2)
>     nablaP' :: Exp (V2 Double) -> Exp (V2 Double)
>     nablaP' ps = ps

> oneStepH98 :: Double -> (V2 Double, V2 Double) -> (V2 Double, V2 Double)
> oneStepH98 h' prev = (qNew, pNew)
>   where
>     h2 = h' / 2
>     gs = pure h'
>     hs' = pure h2
>     p2' = psPrev - hs' * nablaQ' qsPrev
>     qNew = qsPrev + gs * nablaP' p2'
>     pNew = p2' - hs' * nablaQ' qNew
>     qsPrev = P.fst prev
>     psPrev = P.snd prev
>     nablaQ' qs = V2 (q1' / r) (q2' / r)
>       where
>         q1' = qs ^. L._x
>         q2' = qs ^. L._y
>         r   = (q1' ^ 2 + q2' ^ 2) ** (3/2)
>     nablaP' ps = ps

> oneStep2' :: Double ->
>              Exp (V2 Double, V2 Double) ->
>              Exp (V2 Double, V2 Double) ->
>              Exp (V2 Double, V2 Double)
> oneStep2' h' prev dummies = lift (xs + qNew, pNew)
>   where
>     h2 = h' / 2
>     gs :: Exp (V2 Double)
>     gs = lift ((pure h') :: V2 Double)
>     hs' :: Exp (V2 Double)
>     hs' = (lift ((pure h2) :: V2 Double))
>     p2' = psPrev - hs' * nablaQ' qsPrev
>     qNew = qsPrev + gs * nablaP' p2'
>     pNew = p2' - hs' * nablaQ' qNew
>     qsPrev = A.fst prev
>     psPrev = A.snd prev
>     xs = A.fst dummies
>     nablaQ' :: Exp (V2 Double) -> Exp (V2 Double)
>     nablaQ' qs = lift (V2 (q1' / r) (q2' / r))
>       where
>         q1' = qs ^. _x
>         q2' = qs ^. _y
>         r   = (q1' ^ 2 + q2' ^ 2) ** (3/2)
>     nablaP' :: Exp (V2 Double) -> Exp (V2 Double)
>     nablaP' ps = ps

> symplecticEuler ::  Double -> Exp (V2 Double, V2 Double) -> Exp (V2 Double, V2 Double)
> symplecticEuler h' prev = lift(qNew, pNew)
>   where
>     qsPrev = A.fst prev
>     psPrev = A.snd prev
>     gs = lift ((pure h') :: V2 Double)
>     pNew = psPrev - gs * nablaQ' qsPrev
>     qNew = qsPrev + gs * nablaP' pNew
>     nablaQ' :: Exp (V2 Double) -> Exp (V2 Double)
>     nablaQ' qs = lift (V2 (q1' / r) (q2' / r))
>       where
>         q1' = qs ^. _x
>         q2' = qs ^. _y
>         r   = (q1' ^ 2 + q2' ^ 2) ** (3/2)
>     nablaP' :: Exp (V2 Double) -> Exp (V2 Double)
>     nablaP' ps = ps

> nSteps :: Int
> nSteps = 100

> dummyStart :: Exp (V2 Double, V2 Double)
> dummyStart = lift (V2 q10 q20, V2 p10 p20)

> dummyInputs :: Acc (Array DIM1 (V2 Double, V2 Double))
> dummyInputs = A.use $ A.fromList (Z :. nSteps) $
>               P.replicate nSteps (pure 0.0 :: V2 Double, pure 0.0 :: V2 Double)

> runSteps :: Acc (Array DIM1 (V2 Double, V2 Double))
> runSteps = A.scanl (\s _x -> (oneStep2 h s)) dummyStart dummyInputs

> runSteps' :: Acc (Array DIM1 (V2 Double, V2 Double))
> runSteps' = A.scanl (oneStep2' h) dummyStart dummyInputs

> reallyRunSteps :: (Array DIM1 (V2 Double, V2 Double))
> reallyRunSteps = run runSteps

> reallyRunSteps' :: (Array DIM1 (V2 Double, V2 Double))
> reallyRunSteps' = run runSteps'

> symplecticEuler98 ::  Double -> (V2 Double, V2 Double) -> (V2 Double, V2 Double)
> symplecticEuler98 h' prev = (qNew, pNew)
>   where
>     qsPrev = P.fst prev
>     psPrev = P.snd prev
>     gs = pure h'
>     pNew = psPrev - gs * nablaQ' qsPrev
>     qNew = qsPrev + gs * nablaP' pNew
>     nablaQ' qs = V2 (q1' / r) (q2' / r)
>       where
>         q1' = qs ^. L._x
>         q2' = qs ^. L._y
>         r   = (q1' ^ 2 + q2' ^ 2) ** (3/2)
>     nablaP' ps = ps

> dummyStartH98 :: (V2 Double, V2 Double)
> dummyStartH98 = (V2 q10 q20, V2 p10 p20)
>
> dummyInputsH98 :: [(V2 Double, V2 Double)]
> dummyInputsH98 = P.replicate nSteps (pure 0.0 :: V2 Double, pure 0.0 :: V2 Double)

> runStepsH98 :: [(V2 Double, V2 Double)]
> runStepsH98= P.scanl (\s _x -> (oneStepH98 h s)) dummyStartH98 dummyInputsH98

> runStepsSE :: [(V2 Double, V2 Double)]
> runStepsSE = P.scanl (\s _x -> (symplecticEuler98 0.01 s)) dummyStartH98 dummyInputsH98

> bigH2BodyH98 :: (V2 Double, V2 Double) -> Double
> bigH2BodyH98 x = ke + pe
>   where
>     pe = let V2 q1' q2' = P.fst x in negate $ recip (sqrt (q1'^2 + q2'^2))
>     ke = let V2 p1' p2' = P.snd x in 0.5 * (p1'^2 + p2'^2)

> main :: IO ()
> main = do
>   putStrLn $ show $ reallyRunSteps
