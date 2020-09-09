\documentclass{article}
\usepackage{amsmath}
\usepackage{graphicx}
%include polycode.fmt
%options ghci

\begin{document}

\section{Introduction}

When I ran the London Marathon (never again) I was surprised to find I
 couldn't even start to run for about 20 minutes because of the
 density of all the participants.

\section{Some Incorrect Methods}

\subsection{Diffusion?}

$$
u_{t}+a u_{x}=0 \quad \text { for } x \in \mathbb{R}, t \geq 0
$$

$$
\frac{1}{h}(\psi(x-h)-\psi(x))=-\psi_{x}(x)+\mathcal{O}(h)
$$

$$
w_{i}^{\prime}(t)=\frac{a}{h}\left(w_{i-1}(t)-w_{i}(t)\right), \quad i=1,2, \ldots, m
$$

$$
A=\frac{a}{h}\left(\begin{array}{ccccc}
-1 & & & & 1 \\
1 & -1 & & \\
& \ddots & \ddots & \\
& & 1 & -1 \\
& & & 1 & -1
\end{array}\right)
$$



%if False
\begin{code}
{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE NumDecimals #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedStrings #-}

{-# OPTIONS_GHC -Wall            #-}

module WENO where

import Prelude as P

import           Numeric.LinearAlgebra
import           Numeric.Sundials
import           Numeric.Integration.TanhSinh

import           Control.Exception
import           Data.Coerce
import           Katip
import           Katip.Monadic
import           System.IO (stderr)
import           GHC.Int

import           Data.Csv
import           Data.Char
import qualified Data.ByteString.Lazy as BL

import           Control.Monad.Writer
import           Control.Monad.Loops

import           Data.Ratio
\end{code}
%endif

\begin{code}
bigN :: Int
bigN = 50

deltaX :: Double
deltaX = 1.0 / (fromIntegral bigN - 1)

deltaX' :: Double
deltaX' = 1.0 / (fromIntegral bigN)

bigNt :: Int
bigNt = 10

t0 :: Double
t0 = 0.0

tf :: Double
tf =1.0

deltaTT :: Double
deltaTT = (tf - t0) / (fromIntegral bigNt)

bigU0 :: Vector Double
bigU0 = assoc bigN 0.0 [(i, f i) | i <- [0 .. bigN - 1]]
  where
  f i = (sin (pi * fromIntegral i * deltaX)) ^ 100

bigA :: Matrix Double
bigA = assoc (bigN, bigN) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. bigN - 1]
                                                   , j <- [0 .. bigN - 1]
                        ]
 where
   f (i, j) | i       == j          = -1
            | i - 1 == j          =  1
            | i       == 0 &&
              j       == bigN - 1 =  1
            | otherwise               =  0

bigB :: Matrix Double
bigB = assoc (bigN, bigN) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. bigN - 1]
                                                   , j <- [0 .. bigN - 1]
                        ]
 where
   f (i, j) | i       == j          = -1
            | i - 1 == j          =  1
            | i       == 0 &&
              j       == bigN - 1 =  1
            | otherwise               =  0

bigC :: Matrix Double
bigC = assoc (bigN, bigN) 0.0 [ ((i, j), (1/(2*deltaX)) * f (i, j)) | i <- [0 .. bigN - 1]
                                                                    , j <- [0 .. bigN - 1]
                        ]
 where
   f (i, j)
            | i     == 0        =  0
            | i     == bigN - 1 =  0
            | i     == j        =  0
            | i - 1 == j        = -1
            | i + 1 == j        =  1
            | otherwise         =  0

bigD :: Matrix Double
bigD = assoc (bigN, bigN) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. bigN - 1]
                                                   , j <- [0 .. bigN - 1]
                        ]
 where
   f (i, j) | i     == j        =  1
            | otherwise         =  0

simpleAdvect :: OdeProblem
simpleAdvect = emptyOdeProblem
  { odeRhs = odeRhsPure $ \_t x -> coerce (bigA #> (coerce x))
  , odeJacobian = Just (\_t _x -> bigA)
  , odeEvents = mempty
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = bigU0
  , odeSolTimes = vector $ map (deltaTT *) [0] -- [0 .. 1]
  , odeTolerances = defaultTolerances
  }

sol :: IO (Matrix Double)
sol = do
  x <- runNoLoggingT $ solve (defaultOpts SDIRK_5_3_4) simpleAdvect
  case x of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

bigV :: Vector Double
bigV = assoc bigN 0.0 [(i, f i) | i <- [0 .. bigN - 1]]
  where
  f i = (sin (2 * pi * fromIntegral i * deltaX))

bigV' :: Vector Double
bigV' = assoc (bigN + 1) 0.0 [(i, f i) | i <- [0 .. bigN]]
  where
  f i = (sin (2 * pi * fromIntegral i * deltaX'))

sol' :: IO (Matrix Double)
sol' = do
  w <- runNoLoggingT $ solve (defaultOpts' HEUN_EULER_2_1_2) burgers
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

sol'' :: IO (Matrix Double)
sol'' = do
  -- w <- runNoLoggingT $ solve (defaultOpts' HEUN_EULER_2_1_2) burgersWeno
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  logEnv <- registerScribe "stdout" handleScribe defaultScribeSettings =<< initLogEnv "namespace" "devel"
  w <- runKatipT logEnv $ solve (defaultOpts' BOGACKI_SHAMPINE_4_2_3) burgersWeno
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

burgers :: OdeProblem
burgers = emptyOdeProblem
  { odeRhs = odeRhsPure $ \_t x -> coerce (cmap negate ((bigD #> (coerce x)) * (bigC #> (coerce x))))
  , odeJacobian = Nothing
  , odeEvents = mempty
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = bigV
  , odeSolTimes = vector $ map (0.025 *) [0 .. 10]
  , odeTolerances = defaultTolerances
  }

burgersWeno :: OdeProblem
burgersWeno = emptyOdeProblem
  { odeRhs = odeRhsPure $ \_t x -> coerce (rhs' bigN (coerce x))
  , odeJacobian = Nothing
  , odeEvents = mempty
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = bigV'
  , odeSolTimes = vector $ map (0.025 *) [0 .. 10]
  , odeTolerances = defaultTolerances
  }

wcL :: Double -> Double -> Double -> Double -> Double -> Double
wcL v1 v2 v3 v4 v5 = f
  where
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^(2 :: Int) + 0.25*(v1-4.0*v2+3.0*v3)^(2 :: Int)
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^(2 :: Int) + 0.25*(v2-v4)^(2 :: Int)
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^(2 :: Int) + 0.25*(3.0*v3-4.0*v4+v5)^(2 :: Int)

    c1 = 2.0e-1/((eps+s1)^(2 :: Int))
    c2 = 5.0e-1/((eps+s2)^(2 :: Int))
    c3 = 3.0e-1/((eps+s3)^(2 :: Int))

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    q1 = v1/3.0   - 7.0/6.0*v2 + 11.0/6.0*v3
    q2 =(-v2)/6.0 + 5.0/6.0*v3 + v4/3.0
    q3 = v3/3.0   + 5.0/6.0*v4 - v5/6.0

    f = (w1*q1 + w2*q2 + w3*q3)

wcR :: Double -> Double -> Double -> Double -> Double -> Double
wcR v1 v2 v3 v4 v5 = f
  where
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^(2 :: Int) + 0.25*(v1-4.0*v2+3.0*v3)^(2 :: Int)
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^(2 :: Int) + 0.25*(v2-v4)^(2 :: Int)
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^(2 :: Int) + 0.25*(3.0*v3-4.0*v4+v5)^(2 :: Int)

    c1 = 3.0e-1/(eps+s1)^(2 :: Int)
    c2 = 5.0e-1/(eps+s2)^(2 :: Int)
    c3 = 2.0e-1/(eps+s3)^(2 :: Int)

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    q1 =(-v1)/6.0    + 5.0/6.0*v2 + v3/3.0
    q2 = v2/3.0      + 5.0/6.0*v3 - v4/6.0
    q3 = 11.0/6.0*v3 - 7.0/6.0*v4 + v5/3.0

    f = (w1*q1 + w2*q2 + w3*q3)

crwenoR :: Int -> Vector Double -> Vector Double
crwenoR n u = assoc (n + 1) 0.0 [(i, f i) | i <- [0 .. n]]
  where
    f i | i == 0    = 0.0

    f i | i == 1    = wcR v1 v2 v3 v4 v5
      where
        v1 = u!(n-1)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

    f i | i == n - 1 = wcR v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!1

    f i | i == n     = wcR v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!1
        v5 = u!2

    f i | otherwise  = wcR v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

crwenoL :: Int -> Vector Double -> Vector Double
crwenoL n u = assoc (n + 1) 0.0 [(i, f i) | i <- [0 .. n]]
  where
    f i | i == 0 = wcL v1 v2 v3 v4 v5
      where
        v1 = u!(n-2)
        v2 = u!(n-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

    f i | i == 1 = wcL v1 v2 v3 v4 v5
      where
        v1 = u!(n-1)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

    f i | i == n - 1 = wcL v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!1

    f i | i == n     = 0.0

    f i | otherwise  = wcL v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

mLR :: Int -> Matrix Double
mLR n = assoc (n, n) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. n - 1]
                                              , j <- [0 .. n - 1]
                                              ]
  where
    f (i, j) | i == 0 &&
               j == 0     =  1
             | i == 0 &&
               j == n - 1 = -1
             | i == j     =  1
             | i - 1 == j = -1
             | otherwise  =  0

mL :: Int -> Matrix Double
mL n = assoc (n, n + 1) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. n - 1]
                                                 , j <- [0 .. n]
                                                 ]
  where
    f (i, j) | i == 0 &&
               j == 0     =  1
             | i == 0 &&
               j == n - 1 = -1
             | i == j     =  1
             | i - 1 == j = -1
             | otherwise  =  0

mR :: Int -> Matrix Double
mR n = assoc (n, n + 1) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. n - 1]
                                                 , j <- [0 .. n]
                                                 ]
  where
    f (i, j) | i == 0 &&
               j == 0     =  0
             | i == 0 &&
               j == 1     =  1
             | i == 0 &&
               j == n     = -1
             | i == j     = -1
             | i + 1 == j =  1
             | otherwise  =  0

preRhs :: Int -> Vector Double -> Vector Double
preRhs n v = assoc n 0.0 [(i, f i) | i <- [0 .. n - 1]]
  where
    ll = (mL n) #> (crwenoL n v)
    rr = (mR n) #> (crwenoR n v)

    f i | v!i >= 0.0 = negate (v!i * ll!i) / deltaX'
        | otherwise  = negate (v!i * rr!i) / deltaX'

rhs' :: Int -> Vector Double -> Vector Double
rhs' n v = vjoin [preRhs n v, [v!0]]

rhs'' :: Int -> Vector Double -> Vector Double
rhs'' n v = assoc n 0.0 [(i, f i) | i <- [0 .. n - 1]]
  where
    ll = (mL n) #> (crwenoL n (vjoin [v, vector [0.0]]))
    rr = (mR n) #> (crwenoR n (vjoin [vector [0.0], v]))

    f i | v!i >= 0.0 = negate (v!i * ll!i) / deltaX
        | otherwise  = negate (v!i * rr!i) / deltaX

rhs :: Int -> Vector Double -> Vector Double
rhs n u = assoc n 0.0 [(i, f i) | i <- [0 .. n - 1]]
  where
    uL = crwenoL n u
    uR = crwenoR n u
    ll = (mLR n) #> (subVector 0 n uL)
    rr = (mLR n) #> (subVector 1 n uR)
    f i | u!i >= 0.0 = negate (u!i * ll!i)
        | otherwise  = negate (u!i * rr!i)
\end{code}

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{heat1e000.png}
\caption{Caption}
\label{fig:incorrect_0}
\end{figure}

\subsection{Oscillations}

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{burgers_cds.png}
\caption{Caption}
\label{fig:oscillations_0}
\end{figure}

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{burgers_weno.png}
\caption{Caption}
\label{fig:weno_0}
\end{figure}

$$
\frac{u_{i}^{(n+1)}-u_{i}^{(n)}}{\Delta t}+u_{i}^{(n)} \frac{u_{i+1}^{(n)}-u_{i-1}^{(n)}}{2 \Delta x}=0
$$

$$
V(x) \triangleq \int_{-\infty}^{x} v(\xi) d \xi
$$

$$
\begin{aligned}
\frac{1}{\Delta x_{j}} \int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}} p(\xi) d \xi &=\frac{1}{\Delta x_{j}} \int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}} P^{\prime}(\xi) d \xi=\frac{1}{\Delta x_{j}}\left(P\left(x_{j+\frac{1}{2}}\right)-P\left(x_{j-\frac{1}{2}}\right)\right) \\
&=\frac{1}{\Delta x_{j}}\left(V\left(x_{j+\frac{1}{2}}\right)-V\left(x_{j-\frac{1}{2}}\right)\right) \\
&=\frac{1}{\Delta x_{j}}\left(\int_{-\infty}^{x_{j+\frac{1}{2}}} v(\xi) d \xi-\int_{-\infty}^{x_{j-\frac{1}{2}}} v(\xi) d \xi\right) \\
&=\frac{1}{\Delta x_{j}} \int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}} v(\xi) d \xi=\bar{v}_{j}, \quad j=i-r, \ldots, i+s
\end{aligned}
$$

$$
C_{r j}=
\sum_{m=j+1}^{k}
\frac{\sum_{l=0 \atop l \neq m}^{k} \prod_{q=0 \atop q \neq m, l}^{k}(r-q+1)}
     {\prod_{l=0 \atop l \neq m}^{k}(m-l)}
$$

\begin{code}

coeff :: Integral a => Ratio a -> Ratio a -> Ratio a -> Ratio a
coeff k r j = sum [ num m / den m | m <- [j + 1 .. k] ]
  where
    den m = product [ m - l | l <- [0 .. k], l /= m]
    num m = sum [ product [ r - q + 1 | q <- [0 .. k], q /= l, q /= m]
                | l <- [0 .. k], l /= m ]

coeffs k = [ coeff k r j | r <- [-1 .. k - 1], j <- [0 .. k - 1] ]
\end{code}

Here is the eval:
\eval{coeffs 3}

$$
u(x, 0)=(\sin (\pi x))^{100}
$$


\begin{code}
myOptions :: EncodeOptions
myOptions = defaultEncodeOptions {
      encDelimiter = fromIntegral (ord ' ')
    }

main :: IO ()
main = do
  x <- sol
  BL.writeFile "simpleAdvect.txt" $ encodeWith myOptions $ map toList $ toRows x

defaultOpts' :: method -> ODEOpts method
defaultOpts' method = ODEOpts
  { maxNumSteps = 1e5
  , minStep     = 1.0e-14
  , fixedStep   = 0.0 -- 0.001
  , maxFail     = 10
  , odeMethod   = method
  , initStep    = Nothing
  , jacobianRepr = DenseJacobian
  }

defaultOpts :: method -> ODEOpts method
defaultOpts method = ODEOpts
  { maxNumSteps = 1e5
  , minStep     = 1.0e-14
  , fixedStep   = 0
  , maxFail     = 10
  , odeMethod   = method
  , initStep    = Nothing
  , jacobianRepr = SparseJacobian
                 $ SparsePattern
                 $ cmap (fromIntegral :: I -> Int8)
                 $ cmap (\x -> case x of 0 -> 0; _ -> 1)
                 $ flatten
                 $ toInt bigA
  }

instance Element Int8

emptyOdeProblem :: OdeProblem
emptyOdeProblem = OdeProblem
      { odeRhs = error "emptyOdeProblem: no odeRhs provided"
      , odeJacobian = Nothing
      , odeInitCond = error "emptyOdeProblem: no odeInitCond provided"
      , odeEvents = mempty
      , odeTimeBasedEvents = TimeEventSpec $ return $ undefined -- 1.0 / 0.0
      , odeEventHandler = nilEventHandler
      , odeMaxEvents = 100
      , odeSolTimes = error "emptyOdeProblem: no odeSolTimes provided"
      , odeTolerances = defaultTolerances
      }

nilEventHandler :: EventHandler
nilEventHandler _ _ _ = throwIO $ ErrorCall "nilEventHandler"

defaultTolerances :: Tolerances
defaultTolerances = Tolerances
  { absTolerances = Left 1.0e-5
  , relTolerance = 1.0e-10
  }
\end{code}

\section{Age-Dependent Populations}

\subsection{McKendrick / von Foerster}

[McKendrick](https://en.wikipedia.org/wiki/Anderson\_Gray\_McKendrick)
and [von Foerster](https://en.wikipedia.org/wiki/Heinz\_von\_Foerster)
independently derived a model of age-dependent population
growth.

Let $n(a,t)$ be the density of females of age $a$ at time $t$. The
number of females between ages $a$ and $a + \delta a$ are thus $n(a,
t)\delta a$. Assuming individuals are born at age $0$, we have

$$
\frac{\partial}{\partial t}(n(a, t)\delta a) =
J(a, t) - J(a + \delta a, t) - \mu(a, t)n(a, t)\delta a
$$

where $\mu(a, t)$ is the death rate density and $J(a, t)$ denotes the
rate of entry to the cohort of age $a$. Dividing by $\delta a$ we obtain

$$
\frac{\partial}{\partial t}n(a, t) =
 - \frac{J(a + \delta a, t) - J(a, t)}{\delta a} - \mu(a, t)n(a, t)
$$

which in the limit becomes

$$
\frac{\partial}{\partial t}n(a, t) =
 - \frac{\partial J(a, t)}{\partial a} - \mu(a, t)n(a, t)
$$

We can further assume that the rate of entry to a cohort is
proportional to the density of individuals times a velocity of aging
$v(a, t)$.

$$
J(a, t) = n(a, t)v(a, t)
$$

Occasionally there is some reason to assume that aging one year is
different to experiencing one year but we further assume $v = 1$.

We thus obtain

$$
\frac{\partial n(a, t)}{\partial t} + \frac{\partial n(a, t)}{\partial a} =
- \mu(a, t)n(a, t)
$$

\subsection{Gurtin / MacCamy}

To solve any PDE we need boundary and initial conditions. The number
of births at time $t$ is

$$
n(0, t) = \int_0^\infty n(a, t) m(a, N(t))\, \mathrm{d}a
$$

where $m$ is the natality aka birth-modulus and

$$
N(t) = \int_0^\infty n(a, t)\, \mathrm{d}a
$$

and we further assume that the initial condition

$$
n(a, 0) = n_0(a)
$$

for some given $n_0$.

@gurtin1974non focus on the situation where

$$
m(a, N(t)) = \beta(N)e^{-\alpha a}
$$

and we can also assume that the birth rate of Tribbles decreases
exponentially with age and further that Tribbles can live
forever. @gurtin1974non then transform the PDE to obtain a pair of
linked ODEs which can then be solved numerically.

Of course, we know what happens in the Enterprise and rather than
continue with this example, let us turn our attention to the more
serious subject of Malaria.

\section{Malaria}

I realise now that I went a bit overboard with references. Hopefully
they don't interrupt the flow too much.

The World Health Organisation (WHO) estimated that in 2015 there were
214 million new cases of malaria resulting in 438,000 deaths (source:
[Wikipedia](https://en.wikipedia.org/wiki/Malaria\#Epidemiology)).

The lifecycle of the plasmodium parasite that causes malaria is
extremely ingenious. @Thibodeaux2011 model the human segment of the
[plasmodium
lifecycle](https://en.wikipedia.org/wiki/Malaria#Life_cycle) and
further propose a way of determing an optimal treatment for an
infected individual. @hall2013pharmacokinetic also model the effect of
an anti-malarial. Let us content ourselves with reproducing part of
the paper by @Thibodeaux2011.

At one part of its sojourn in humans, plasmodium infects erythrocytes
aka red blood cells. These latter contain haemoglobin (aka
hemoglobin).  The process by which red blood cells are produced,
[Erythropoiesis](https://en.wikipedia.org/wiki/Erythropoiesis), is
modulated in a feedback loop by
[erythropoietin](https://en.wikipedia.org/wiki/Erythropoietin). The
plasmodium parasite severely disrupts this process. Presumably the
resulting loss of haemoglobin is one reason that an infected
individual feels ill.

As can be seen in the overview by @Torbett2009,
the full feedback loop is complex. So as not to lose ourselves in the
details and following @Thibodeaux2011 and @BELAIR1995317, we consider
a model with two compartments.

* Precursors: prototype erythrocytes developing in the bone marrow
  with $p(\mu, t)$ being the density of such cells of age $\mu$ at
  time $t$.

* Erythrocytes: mature red blood cells circulating in the blood with
  $m(\mu, t)$ being the density of such cells of age $\nu$ at time
  $t$.

$$
\begin{aligned}
\frac{\partial p(\mu, t)}{\partial t} + g(E(t))\frac{\partial p(\mu, t)}{\partial \mu} &=
\sigma(\mu, t, E(t))p(\mu, t) & 0 < \mu < \mu_F & & 0 < t < T \\
\frac{\partial m(\nu, t)}{\partial t} + \phantom{g(E(t))}\frac{\partial m(\nu, t)}{\partial \nu} &=
-\gamma(\nu, t, M(t))m(\nu, t) & 0 < \nu < \nu_F & &  0 < t < T
\end{aligned}
$$

where $\sigma(\mu, t, E(t))$ is the birth rate of precursors and
$\gamma(\nu, t, M(t))$ is the death rate of erythrocytes, $g(E(t))$ is
the maturation rate of precursors and where

$$
M(t) = \int_0^{\nu_F} p(\nu, t) \,\mathrm{d}\nu
$$

As boundary conditions, we have that the number of precursors maturing
must equal the production of number of erythrocytes

$$
m(0, t) = g(E(t))p(\mu_F, t)
$$

and the production of the of the number of precursors depends on the
level of erythropoietin

$$
g(E(t))p(0, t) = \phi(t)E(t)
$$

where $\phi(t)$ is some proportionality function.

As initial conditions, we have

$$
\begin{aligned}
p(\mu, 0) &= p_0(\mu) \\
m(\nu, 0) &= m_0(\nu)
\end{aligned}
$$

We can further model the erythropoietin dynamics as

$$
\frac{\mathrm{d}E(t)}{\mathrm{d}t} = f(M(t), t) - a_E(P(t))E(t)
$$

where $f$ is the feedback function from the kidneys and the decay
rate, $a_E$ depends on the total precursor population $P(t)$
(@sawyer1987binding) although this often is taken to be a constant and
I would feel more comfortable with a more recent citation and where

$$
P(t) = \int_0^{\mu_F} p(\mu, t) \,\mathrm{d}\mu
$$

As initial condition we have

$$
E(0) = E_0
$$

A Finite Difference Attempt
---------------------------

Let us try solving the above model using a finite difference scheme
observing that we currently have no basis for whether it has a
solution and whether the finite difference scheme approximates such a
solution! We follow @Thibodeaux2011 who give a proof of convergence
presumably with some conditions; any failure of the scheme is entirely
mine.

Divide up the age and time ranges, $[0, \mu_F]$, $[0, \nu_F]$ and $[0, T]$
into equal sub-intervals,
$[\mu_i, \mu_{i+1}]$, $[\nu_j, \nu_{j+1}]$ and $[t_k, t_{k+1}]$
where

$$
\begin{aligned}
\mu_i &= i\Delta\mu & & \mathrm{for} & i = 1 \ldots n_1 \\
\nu_j &= j\Delta\nu & & \mathrm{for} & j = 1 \ldots n_2 \\
t_k   &= k\Delta t  & & \mathrm{for} & k = 1 \ldots K
\end{aligned}
$$

where $\Delta\mu = \mu_F / n_1$, $\Delta\nu = \nu_F / n_2$ and $\Delta t = T / K$.

Denoting $p(\mu_i, t_k) = p_i^k$ and similarly we obtain

$$
\begin{aligned}
\frac{p_i^{k+1} - p_i^k}{\Delta t} + g^k\frac{p_i^{k+1} - p_{i-1}^{k+1}}{\Delta\mu} &= \sigma_i^k p_i^{k+1} \\
\frac{m_j^{k+1} - m_j^k}{\Delta t} + \phantom{g^k}\frac{m_j^{k+1} - m_{j-1}^{k+1}}{\Delta\mu} &= -\gamma_j^k m_j^{k+1}
\end{aligned}
$$

and

$$
\begin{aligned}
\frac{E^{k+1} - E^k}{\Delta t} &= f^k - a_E^k E^{k+1} \\
g^k p_0^{k+1} &= \phi^k E^k \\
m_0^{k+1}     &= g^k m_{n_1}^{k+1}
\end{aligned}
$$

Re-arranging we get

$$
\begin{aligned}
-g^k\frac{\Delta t}{\Delta \mu}p_{i-1}^{k+1} +
\bigg(1 + g^k\frac{\Delta t}{\Delta \mu} - \Delta t \sigma_i^k\bigg)p_i^{k+1} &=
p_i^k \\
\frac{\Delta t}{\Delta \mu}m_{j-1}^{k+1} +
\bigg(1 + \frac{\Delta t}{\Delta \mu} + \Delta t \gamma_j^k\bigg)m_j^{k+1} &=
m_j^k
\end{aligned}
$$

Writing

$$
\begin{aligned}
d_{1,i}^k &= 1 + g^k\frac{\Delta t}{\Delta \mu} - \Delta t \sigma_i^k \\
d_{2,i}^k &= 1 + \frac{\Delta t}{\Delta \nu} - \Delta t \gamma_i^k
\end{aligned}
$$

We can express the above in matrix form

$$
\begin{bmatrix}
g^k & 0 & 0 & \ldots & 0 & 0 \\
-g^k\frac{\Delta t}{\Delta \mu} & d_{1,1}^k & 0 & \ldots & 0 & 0\\
0 & -g^k\frac{\Delta t}{\Delta \mu} & d_{1,2}^k & \ldots & 0 & 0 \\
\ldots & \ldots & \ldots & \ldots & \ldots & \ldots \\
0 & 0 & 0 & \ldots &\ -g^k\frac{\Delta t}{\Delta \mu} & d_{1,n_1}^k \\
\end{bmatrix}
\begin{bmatrix}
p_0^{k+1} \\
p_1^{k+1} \\
p_2^{k+1} \\
\ldots \\
p_{n_1}^{k+1}
\end{bmatrix}
=
\begin{bmatrix}
\phi^k E^k \\
p_1^k \\
p_2^k \\
\ldots \\
p_{n_1}^k \\
\end{bmatrix}
$$

$$
\begin{bmatrix}
1 & 0 & 0 & \ldots & 0 & 0 \\
-\frac{\Delta t}{\Delta \mu} & d_{2,1}^k & 0 & \ldots & 0 & 0\\
0 & -\frac{\Delta t}{\Delta \mu} & d_{2,2}^k & \ldots & 0 & 0 \\
\ldots & \ldots & \ldots & \ldots & \ldots & \ldots \\
0 & 0 & 0 & \ldots &\ -\frac{\Delta t}{\Delta \mu} & d_{2,n_1}^k \\
\end{bmatrix}
\begin{bmatrix}
m_0^{k+1} \\
m_1^{k+1} \\
m_2^{k+1} \\
\ldots \\
m_{n_2}^{k+1}
\end{bmatrix}
=
\begin{bmatrix}
g^k p_{n_1}^{k+1} \\
m_1^k \\
m_2^k \\
\ldots \\
m_{n_1}^k \\
\end{bmatrix}
$$

Finally we can write

$$
E^{k+1} = \frac{E^k + \Delta t f^k}{1 + a_E^k\Delta T}
$$


References
==========

\end{document}
