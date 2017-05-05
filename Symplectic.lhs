% Rao-Blackwellisation
% Dominic Steinitz
% 9th April 2017

---
bibliography: DynSys.bib
---

> {-# OPTIONS_GHC -Wall                   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults #-}
>
> {-# LANGUAGE FlexibleContexts #-}

> import Data.Number.Symbolic
> import Numeric.AD
> import Data.Array.Accelerate              as A   hiding ((^), iterate)
> import Data.Array.Accelerate.LLVM.Native  as CPU
> import Data.Array.Accelerate.Linear
> import Data.Array.Accelerate.Control.Lens

Let us consider the following Hamiltonian

> bigH :: Prelude.Floating a => a -> a -> a
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

> nablaQ, nablaP :: Prelude.Floating a => a -> a
> nablaQ = sin
> nablaP = id

> oneStep :: Prelude.Floating a => (a -> a) ->
>                                  (a -> a) ->
>                                  a ->
>                                  (a, a) ->
>                                  (a, a)
> oneStep nablaQ nablaP h (qPrev, pPrev) = (qNew, pNew)
>   where
>     h2 = h / 2
>     p2 = pPrev - h2 * nablaQ qPrev
>     qNew = qPrev + h * nablaP p2
>     pNew = p2 - h2 * nablaQ qNew

> h :: Double
> h = 0.1

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
> ham :: Prelude.Floating a => [a] -> a
> ham [p1, p2, q1, q2] = 0.5 * (p1^2 + p2^2) - recip (sqrt (q1^2 + q2^2))
> ham _ = error "Hamiltonian defined for 4 generalized co-ordinates only"

> q1, q2, p1, p2 :: Sym a
> q1 = var "q1"
> q2 = var "q2"
> p1 = var "p1"
> p2 = var "p2"

> testHam :: Prelude.Floating a => [a] -> [[a]]
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
>     gs = (lift ((pure h') :: V2 Double))
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



