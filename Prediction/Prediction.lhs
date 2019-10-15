% Haskell for Numerics?
% Dominic Steinitz
% 2nd June 2017

> {-# OPTIONS_GHC -Wall #-}

> {-# LANGUAGE DataKinds #-}
> {-# LANGUAGE QuasiQuotes #-}
> {-# LANGUAGE LambdaCase #-}
> {-# LANGUAGE GADTs #-}
> {-# LANGUAGE ScopedTypeVariables #-}
> {-# LANGUAGE OverloadedStrings #-}
> {-# LANGUAGE TypeOperators #-}

> {-# LANGUAGE DataKinds #-}
> {-# LANGUAGE KindSignatures #-}
> {-# LANGUAGE TypeFamilies #-}
> {-# LANGUAGE TypeOperators #-}
> {-# LANGUAGE FlexibleInstances #-}
> {-# LANGUAGE FlexibleContexts #-}
> {-# LANGUAGE ScopedTypeVariables #-}
> {-# LANGUAGE ConstraintKinds #-}
> {-# LANGUAGE ExistentialQuantification #-}
> {-# LANGUAGE RankNTypes #-}
> {-# LANGUAGE PolyKinds #-}

> {-# LANGUAGE DataKinds             #-}
> {-# LANGUAGE TypeOperators         #-}
> {-# LANGUAGE KindSignatures        #-}
> {-# LANGUAGE GADTs                 #-}
> {-# LANGUAGE Rank2Types            #-}
> {-# LANGUAGE ScopedTypeVariables   #-}
> {-# LANGUAGE MultiParamTypeClasses #-}
> {-# LANGUAGE FlexibleInstances     #-}
> {-# LANGUAGE TypeFamilies          #-}
> {-# LANGUAGE UndecidableInstances  #-}
> {-# LANGUAGE PolyKinds             #-}

> import           Numeric.Sundials.ARKode.ODE
> import           Numeric.LinearAlgebra
> import qualified Naperian as N
> import qualified Data.Foldable as F
> import           Control.Applicative ( liftA2 )
> import qualified GHC.TypeNats as M
>
> import Graphics.Rendering.Chart hiding (Matrix, Vector)
> import Graphics.Rendering.Chart.Backend.Diagrams
> import Diagrams.Backend.Cairo.CmdLine
> import Diagrams.Prelude hiding (render, Renderable, (*~), Time, Vector)
> import Diagrams.Backend.CmdLine

> import System.IO.Unsafe

> displayHeader :: FilePath -> Diagram B -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 700) fn
>              , DiagramLoopOpts False Nothing 0
>              )

> chart :: String ->
>          String ->
>          [[(Double, Double)]] ->
>          Renderable ()
> chart t l obss = toRenderable layout
>   where

>     actual x l c = plot_lines_values .~ [x]
>                    $ plot_lines_style  . line_color .~ opaque c
>                    -- $ plot_lines_title .~ l
>                    $ plot_lines_style  . line_width .~ 1.0
>                    $ def

>     ls = map (\n -> "Path " ++ show n) [1..]
>     cs = cycle [blue, green, red, brown, crimson]

>     actuals' :: [PlotLines Double Double]
>     actuals' = zipWith3 actual obss ls cs

>     layout = layout_title .~ t
>            $ layout_plots .~ (map toPlot actuals')
>            $ layout_y_axis . laxis_title .~ l
>            $ layout_y_axis . laxis_override .~ axisGridHide
>            $ layout_x_axis . laxis_title .~ "Time"
>            $ layout_x_axis . laxis_override .~ axisGridHide
>            $ def

> diagrmM :: String -> String -> [[(Double, Double)]] -> IO (Diagram Cairo)
> diagrmM t l xss = do
>   denv <- defaultEnv vectorAlignmentFns 600 500
>   return $ fst $ runBackendR denv (chart t l xss)

> main :: IO ()
> main = do
>   putStrLn "Hello"

We want to solve the heat equation

$$
\frac{\partial u}{\partial t} = \beta \bigg[\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial x^2}\bigg]
$$

Why? Because we want to solve e.g

$$
\frac{\mathrm{d}x}{\mathrm{d}t} = \theta x\bigg(1 - \frac{x}{k}\bigg)
$$

but with

$$
\mathrm{d}\theta = \mathrm{d}W_t
$$

Notice that there is nothing stochastic about the biology but we
express our uncertainty about the parameter by making it a
time-varying stochastic variable which says the further we go into the
future the less certain we are about it.

We are going to turn this into a Fokker-Planck equation which we can
then solve using e.g. the method of lines. But before turning to
Fokker-Planck, let's show that we can indeed solve a diffusion
equation using the method of lines.

Let us solve the heat equation over the unit square to some arbitrary
point in the future.


> x1, a, x2 :: Double
> x1 = 0
> a = 1.0
> x2 = a

> y1, y2 :: Double
> y1 = 0.0
> y2 = 1.0

> bigT :: Double
> bigT = 1000.0


> n :: Int
> n = 3

> dx :: Double
> dx = a / (fromIntegral n + 1)

> dy :: Double
> dy = a / (fromIntegral n + 1)

> beta, s :: Double
> beta = 1.0e-5
> s = beta / (dx * dx)

We discretize the space variables. I probably ought to show how this is done.

With one spatial dimension we have:

$$
u_{t}=k u_{x x} + f
$$

initial condition $u(0, x)=0$

stationary boundary conditions

$$
\frac{\partial u}{\partial t}(t, 0)=\frac{\partial u}{\partial t}(t, 1)=0
$$

$$
f(t, x)=\left\{\begin{array}{ll}{1} & {\text { if } x=1 / 2} \\ {0} & {\text { otherwise }}\end{array}\right.
$$

and we can discretize over this spatial dimension using:

$$
u_{x x}=\frac{u_{j+1}-2 u_{j}+u_{j-1}}{\Delta x^{2}}
$$

where

$$
u_{j}(t) \triangleq u\left(t, x_{j}\right), \quad x_{j} \triangleq j \Delta x, \quad 0 \leq j \leq n+1
$$

$$
\dot{u}_i = \sum_0^{n+1} A_{i\,j} u_j + B_i, \quad 0 \leq i \leq n+1
$$

where

$$
\begin{aligned}
A_{0\,j}     = 0, & \quad 0 \leq j \leq n+1, & \text{boundary condition} \\
A_{i\,i-1}   = 1  &                          &                           \\
A_{i\,i}     = 2  &                          &                           \\
A_{i\,i+1}   = 1  &                          &                           \\
A_{{n+1}\,j} = 0, & \quad 0 \leq j \leq n+1, & \text{boundary condition} \\
A_{i\,j}     = 0  & \quad \text{otherwise}   &                           \\
\end{aligned}
$$

Converting this to a system of ODEs is straightforward:

$$
\begin{bmatrix}
\dot{u_0} \\
\dot{u_1} \\
\dot{u_2}
\end{bmatrix}
=
\begin{bmatrix}
0 & 0 & 0 \\
1 & 2 & 1 \\
0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
u_0 \\
u_1 \\
u_2
\end{bmatrix}
+
\begin{bmatrix}
f_0 \\
f_1 \\
f_2
\end{bmatrix}
$$

where $f_j \triangleq f(t, x_j)$.

How about two spatial variables?

$$
\frac{\partial u}{\partial t}=k_{x} \frac{\partial^{2} u}{\partial x^{2}}+k_{y} \frac{\partial^{2} u}{\partial y^{2}}+h
$$

with initial condition $u(0, x, y) = 0$ and stationary boundary conditions

$$
\frac{\partial u}{\partial t}(t, 0, y)=\frac{\partial u}{\partial t}(t, 1, y)=\frac{\partial u}{\partial t}(t, x, 0)=\frac{\partial u}{\partial t}(t, x, 1)=0
$$

and a periodic heat source

$$
h(x, y)=\sin (\pi x) \sin (2 \pi y)
$$

This has analytic solution

$$
u(t, x, y)=\frac{1-e^{-\left(k_{x}+4 k_{y}\right) \pi^{2} t}}{\left(k_{x}+4 k_{y}\right) \pi^{2}} \sin (\pi x) \sin (2 \pi y)
$$

$$
u_{i\,j}(t) \triangleq u\left(t, x_{i}, y_{j}\right), \quad x_{i} \triangleq i \Delta x, \quad 0 \leq j \leq n+1, \quad  y_{j} \triangleq j \Delta y
$$


$$
\begin{align}
u_{x x} &= \frac{u_{i+1\,j}-2 u_{i\,j}+u_{i-1\,j}}{\Delta x^{2}} \\
u_{y y} &= \frac{u_{i\,j+1}-2 u_{i\,j}+u_{i\,j-1}}{\Delta y^{2}}
\end{align}
$$

$$
\dot{u}_{i\, j} = \frac{k_x}{(\Delta x)^2}({u_{i+1\,j}-2 u_{i\,j}+u_{i-1\,j}})
                + \frac{k_y}{(\Delta y)^2}({u_{i\,j+1}-2 u_{i\,j}+u_{i\,j-1}})
$$

$$
\dot{u}_{i\, j} = \sum_{k=0}^{n+1}\sum_{l=0}^{n+1}A_{i\,j\,k\,l} u_{k\,l}
$$

$$
\begin{align}
A_{0\,j\,k\,l} &= 0 \\
A_{i\,j\,i-1\,j} &= 1 \\
A_{i\,j\,i\,j} &= -2 \\
A_{i\,j\,i+1\,j} &= 1 \\
A_{n+1\,j\,k\,l} &= 0 \\
\end{align}
$$

$$
\begin{align}
\begin{bmatrix}
\dot{u}_{0,0} \\
\dot{u}_{0,1} \\
\dot{u}_{0,2} \\
\dot{u}_{1,0} \\
\dot{u}_{1,1} \\
\dot{u}_{1,2} \\
\dot{u}_{2,0} \\
\dot{u}_{2,1} \\
\dot{u}_{2,2} \\
\dot{u}_{3,0} \\
\dot{u}_{3,1} \\
\dot{u}_{3,2} \\
\end{bmatrix}
&=
\begin{bmatrix}
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & -2 &  0 & 0 &  1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  1 &  0 & 0 & -2 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
u_{0,0} \\
u_{0,1} \\
u_{0,2} \\
u_{1,0} \\
u_{1,1} \\
u_{1,2} \\
u_{2,0} \\
u_{2,1} \\
u_{2,2} \\
u_{3,0} \\
u_{3,1} \\
u_{3,2} \\
\end{bmatrix} \\
&+
\begin{bmatrix}
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & -2 &  0 & 0 &  1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  1 &  0 & 0 & -2 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 &  0 &  0 & 0 &  0 & 0 & 0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
u_{0,0} \\
u_{0,1} \\
u_{0,2} \\
u_{1,0} \\
u_{1,1} \\
u_{1,2} \\
u_{2,0} \\
u_{2,1} \\
u_{2,2} \\
u_{3,0} \\
u_{3,1} \\
u_{3,2} \\
\end{bmatrix}
\end{align}
$$

> tensor11 :: (Num a, N.Dimension f, N.Dimension g) => f a -> g a -> f (g a)
> tensor11 xs ys = liftA2 (liftA2 (*)) (fmap pure xs) (pure ys)

> tensor22 :: (Num a, N.Dimension e, N.Dimension f, N.Dimension g, N.Dimension h) =>
>             f (e a) -> g (h a) -> f (e (g (h a)))
> tensor22 xs ys = fmap N.transpose $ liftA2 (liftA2 (tensor11)) (fmap pure xs) (pure ys)


> foo :: (M.KnownNat n, M.KnownNat m, M.KnownNat p, (n M.* m) ~ p) =>
>        N.Hyper '[N.Vector n, N.Vector m] Double ->
>        N.Hyper '[N.Vector p] Double
> foo w = case z of
>   Nothing -> error "foo"
>   Just v  -> N.Prism $ N.Scalar v
>   where
>     z = N.fromList $
>         concat $ fmap F.toList $ F.toList $
>         N.point $ N.crystal $ N.crystal w

> bigAA1 :: Matrix Double
> bigAA1 = assoc (n * n, n * n) 0.0 [((i, j), f (i, j)) | i <- [0 .. n * n - 1]
>                                                       , j <- [i - n, i,  i + n]
>                                                       , j `elem` ([0 .. n * n -1] :: [Int])]
>  where
>    f (i, j) | i     == j = (-2.0)
>             | i - n == j = 1.0
>             | i + n == j = 1.0
>             | otherwise = error $ show (i, j)

> bigAA2 :: Matrix Double
> bigAA2 = diagBlock (replicate 4 bigB)
>  where
>    bigB :: Matrix Double
>    bigB = assoc (n, n) 0.0 [((i, j), f (i, j)) | i <- [0 .. n - 1]
>                                                , j <- [i-1..i+1]
>                                                , j `elem` ([0..n-1] :: [Int])]
>      where
>        f (i, j) | i     == j = (-2.0)
>                 | i - 1 == j = 1.0
>                 | i + 1 == j = 1.0
>                 | otherwise  = error "bigAA2"

> bigAA :: Matrix Double
> bigAA = cmap (*s) $ bigAA1 + bigAA2

Set up the initial conditions
-----------------------------

Since we know that a solution is the normal distribution the variance
of which gets larger as time increases, we can start with a normal
distribution.

> bb :: Vector Double
> bb = assoc (n * n) 0.0 []

> x, y :: Vector Double
> x = linspace n (x1 + dx, x2 - dx)
> y = linspace n (y1 + dy, y2 - dy)

> sigma :: Double
> sigma = 0.05

> bigUU0 :: Vector Double
> bigUU0 = vector $
>         map (\(u, v) -> exp (-((u - mux)^n2 + (v - muy)^n2) / (2 * sigma^n2)) / (2 * pi * sigma^n2))
>             [(u, v) | u <- toList x, v <- toList y]
>  where
>    mux = (x2 - x1) / 2
>    muy = (y2 - y1) / 2
>    n2 :: Int
>    n2 = 2

Now we can solve our system with a rather solver and solver parameters.

> sol :: Matrix Double
> sol = odeSolveV SDIRK_5_3_4' Nothing 1.0e-6 1.0e-10 (const bigUU') bigUU0 (vector [0, bigT])
>   where
>     bigUU' bigUU = bigAA #> bigUU + bb

For now you will have to open the html files yourself. I think I know how to fix this but one step at a time.

Fokker-Planck
=============

$$
\frac{\partial}{\partial t} p(t, \mathbf{x})+\sum_{k=1}^{n} \frac{\partial}{\partial x_{k}}\left(g_{k}(t, \mathbf{x}) p(t, \mathbf{x})\right)=\frac{1}{2} \sum_{j=1, k=1}^{n} \frac{\partial^{2}}{\partial x_{j} \partial x_{k}}\left[\left(\sigma(t, \mathbf{x}) \sigma^{T}(t, \mathbf{x})\right)_{j k} p(t, \mathbf{x})\right]
$$

$$
\frac{\partial}{\partial t} p(t, y, \theta)+\frac{\partial}{\partial y}(f(t, y, \theta ; k) p(t, y, \theta))=\frac{1}{2}\left[\sigma_{y}^{2} \frac{\partial^{2}}{\partial y^{2}} p(t, y, \theta)+\sigma_{\theta}^{2} \frac{\partial^{2}}{\partial \theta^{2}} p(t, y, \theta)\right]
$$

$$
\frac{\partial}{\partial t} p(t, y, \theta)+\frac{\partial}{\partial y}(f(t, y, \theta ; k) p(t, y, \theta))=0
$$

Since

$$
\dot{y}=\theta y\left(1-\frac{y}{k}\right)
$$

we have

$$
f(t, y, \theta ; k)=\theta y\left(1-\frac{y}{k}\right)
$$

$$
\frac{\partial p}{\partial t} +\frac{\partial}{\partial y}\bigg(\theta y\bigg(1 - \frac{y}{k}\bigg) p\bigg)=0
$$

$$
\frac{\partial p}{\partial t} + p\frac{\partial}{\partial y}\bigg(\theta y\bigg(1 - \frac{y}{k}\bigg) \bigg) + \theta y\bigg(1 - \frac{y}{k}\bigg)\frac{\partial p}{\partial y} = 0
$$

$$
\frac{\partial p}{\partial t} + p\theta\bigg(1 - \frac{y}{k} - \frac{y}{k}\bigg) + \theta y\bigg(1 - \frac{y}{k}\bigg)\frac{\partial p}{\partial y} = 0
$$

$$
\frac{\partial p}{\partial t} + \theta y\bigg(1 - \frac{y}{k}\bigg)\frac{\partial p}{\partial y} = - p\theta\bigg(1 - \frac{2y}{k}\bigg)
$$

We can solve the transport PDE with initial condition
$$
\left\{\begin{array}{l}{u_{t}+a u_{x}=0} \\ {u(x, 0)=\phi(x)}\end{array}\right.
$$

using the Method of Characteristics with $a = 2$ and $\phi (x) =
e^{-x^2}$ to give the solutuion illustrated

![](diagrams/transport.png)

$$
\begin{aligned}
\frac{d S}{d t} &=-\delta S(t) I(t) \\
\frac{d I}{d t} &=\delta S(t) I(t)-\gamma I(t) \\
\frac{d R}{d t} &=\quad \gamma I(t)
\end{aligned}
$$

$$
\begin{aligned}
s^{\prime}(t) &=-\beta s(t) i(t)-\mu s(t)+\mu \\
i^{\prime}(t)    &= \beta s(t) i(t)-\mu i(t)
\end{aligned}
$$

theta = c(0.0026, 0.5)

> nPupilsS = 762
> nPupilsI = 1

> delta, gamma :: Double
> delta = 0.0026 -- beta
> gamma = 0.5    -- mu

> s0, i0 :: Double
> s0 = nPupilsS
> i0 = nPupilsI

> bigC :: Double
> bigC = s0 + i0 - 1

> lambda :: Double
> lambda = delta - gamma + delta * bigC

> bigD :: Double
> bigD = d / n
>   where
>     d = lambda - i0 * delta
>     n = lambda * i0 * exp (delta * bigC / gamma)

> ii :: Double -> Double
> ii t = lambda / (delta + lambda * bigD * exp (-lambda * t) * exp (delta * bigC / gamma))
>
> ss :: Double -> Double
> ss t = 1 + bigC * exp (-gamma * t) - ii t

> sirSol = odeSolve f [i0, s0] (fromList [0.0, 0.1 .. 20.0])
>   where
>     f _t [i, s] = [ -delta * s * i - gamma * s + gamma
>                   ,  delta * s * i - gamma * i
>                   ]

$$
\begin{aligned}
0 &=-\beta s_{\infty} i_{\infty}-\mu s_{\infty}+\mu \\
0 &= \beta s_{\infty} i_{\infty}-\mu i_{\infty}
\end{aligned}
$$

$$
\begin{aligned}
0 &= -\mu i_{\infty} -\mu s_{\infty}+\mu
\end{aligned}
$$

$$
\begin{aligned}
s_{\infty} + i_{\infty} = 1
\end{aligned}
$$

$$
\begin{aligned}
\beta (1 - i_{\infty})i_{\infty} - \mu i_{\infty} &= 0 \\
(1 - i_{\infty}) i_{\infty} - \frac{\mu}{\beta} i_{\infty} &= 0
\end{aligned}
$$

$$
\begin{aligned}
i_{\infty} &= 1 - \frac{\mu}{\beta} \\
s_{\infty} &=     \frac{\mu}{\beta}
\end{aligned}
$$

$$
i(t)=\frac{\lambda}{\beta+\lambda\left(\frac{\lambda-i_{0} \beta}{\lambda i_{0} e} \frac{\beta-i_{0} \beta}{\mu}\right) e^{-\lambda t+\frac{\beta\left(s_{0}+i_{0}-1\right)}{\mu}}}
$$

If we let $t \rightarrow \infty$ then we obtain $i_{\infty} =
\frac{\lambda}{\beta}$.

$$
\left[ \begin{array}{c}{y_{i}} \\ {\theta_{i}}\end{array}\right]=\left[ \begin{array}{c}{\frac{k y_{i-1} \exp \theta_{i-1}\left(t_{i}-t_{i-1}\right)}{k+y_{i-1}\left(\exp \theta_{i-1}\left(t_{i}-t_{i-1}\right)-1\right)}} \\ {\theta_{i-1}}\end{array}\right]+\psi_{i-1}
$$

We just need to apply the MoC to the PDE for the probability distribution.

$$
\left[ \begin{array}{c}{y_{i}} \\ {\log \theta_{i}}\end{array}\right]=\left[ \begin{array}{c}{\frac{k y_{i-1} \exp \theta_{i-1}\left(t_{i}-t_{i-1}\right)}{k+y_{i-1}\left(\exp \theta_{i-1}\left(t_{i}-t_{i-1}\right)-1\right)}} \\ {\log \theta_{i-1}}\end{array}\right]+\psi_{i-1}
$$

$$
\left[ \begin{array}{c}{z_{i}}\end{array}\right]=\left[ \begin{array}{ll}{1} & {0}\end{array}\right] \left[ \begin{array}{l}{y_{i}} \\ {\theta_{i}}\end{array}\right]+\nu_{i}
$$

$$
\psi_{i} \sim \mathcal{N}(0, Q)
$$

$$
v_{i} \sim \mathcal{N}(0, R)
$$

$$
f(t, y, \theta ; k)=\frac{k y_{0} \exp \theta t}{k+y_{0}(\exp \theta t-1)}
$$

$$
\dot{y}=\theta y\left(1-\frac{y}{k}\right)
$$

$$
y=\frac{k y_{0} \exp \theta t}{k+y_{0}(\exp \theta t-1)}
$$
