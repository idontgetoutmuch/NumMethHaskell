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
{-# LANGUAGE TypeFamilies      #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE TypeOperators     #-}
{-# LANGUAGE NumDecimals       #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}

{-# OPTIONS_GHC -Wall          #-}

module WENO where

import Prelude as P

import           Numeric.LinearAlgebra
import           Numeric.Sundials

import           Control.Exception
import           Data.Coerce
import           Katip
import           Katip.Monadic
import           System.IO (stderr)
import           GHC.Int
import           Foreign.C.Types

import           Data.Csv
import           Data.Char
import qualified Data.ByteString.Lazy as BL
import           Data.ByteString.Lazy (putStr, writeFile)
import           Prelude hiding (putStr, writeFile)

import           Control.Monad.Writer

import           Data.Ratio

import qualified Data.Vector as V

import           Data.List (transpose)

import qualified Language.R as R
import           Language.R (R)
import           Language.R.QQ

import           Frames
import           Frames.CSV
import           Control.Lens((^.))
import           System.IO.Unsafe (unsafePerformIO)
import qualified Data.Foldable as F
\end{code}
%endif

\begin{code}
tableTypes "NoisyObs" "lynx_hare_df.csv"

loadNoisyObs :: IO (Frame NoisyObs)
loadNoisyObs = inCoreAoS (readTable "lynx_hare_df.csv")

predPreyObs :: [((Int, Double), (Int, Double))]
predPreyObs = unsafePerformIO $
          do xs <- loadNoisyObs
             let us = F.toList $ fmap (\u -> ((u ^. year), (u ^. hare))) xs
             let vs = F.toList $ fmap (\u -> ((u ^. year), (u ^. lynx))) xs
             return $ zip us vs

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
\bar{v}_{i} \equiv \frac{1}{\Delta x_{i}} \int_{x_{i-\frac{1}{2}}}^{x_{i+\frac{1}{2}}} v(\xi) d \xi
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

We assume a uniform grid TBD. A $p$-point stencil at $i$ is defined as

$$
S(i)_{p}^{q}=\left\{x_{i-p+q}, \ldots, x_{i+q-1}\right\}
$$

where $q \in\{1, \ldots, p\}$. The explicit dependence on $i$ bloats
 the notation without adding clarity is dropped for now.

Rather than work in full generality let us work with some specific stencils

$$
\begin{aligned}
S_{5}^{3} &= \left\{x_{i-2}, x_{i-1}, x_{i}, x_{i+1}, x_{i+2}\right\} \\
S_{1}^{3} &= \left\{x_{i-2}, x_{i-1}, x_{i}\right\} \\
S_{2}^{3} &= \left\{x_{i-1}, x_{i}, x_{i+1}\right\} \\
S_{3}^{3} &= \left\{x_{i}, x_{i+1}, x_{i+2}\right\}
\end{aligned}
$$

$$
q_{q}(x)=\sum_{x_{k} \in S_{p}^{q}} y\left(x_{k}\right) \ell_{k}(x)
$$

$$
\ell_{k}(x)=\prod_{x_{j} \in S_{p}^{q}} \frac{x-x_{j}}{x_{k}-x_{j}}
$$

$$
\begin{aligned}
q_{1}(x) &= \sum_{x_{k} \in S_{1}^{3}} y\left(x_{k}\right) \ell_{k}(x) \\
q_{2}(x) &= \sum_{x_{k} \in S_{2}^{3}} y\left(x_{k}\right) \ell_{k}(x) \\
q_{3}(x) &= \sum_{x_{k} \in S_{3}^{3}} y\left(x_{k}\right) \ell_{k}(x)
\end{aligned}
$$

We can write the unique order $4$ polynomial over $S_{5}^{3}$ as

$$
p(x)=\sum_{m=-p+q+\ell}^{q} \gamma_{m}(x) q_{m}(x)
$$

in our case this is

$$
\begin{aligned}
p(x) &= \frac{(x - x_{i+1})(x - x_{i+2})}{\Gamma_1} q_1 \\
     &+ \frac{(x - x_{i-2})(x - x_{i+2})}{\Gamma_2} q_2 \\
     &+ \frac{(x - x_{i-2})(x - x_{i-1})}{\Gamma_3} q_3 \\
\end{aligned}
$$

Possibly

$$
\sum_{m=-p+q+\ell}^{q} \gamma_{m}(x)=1
$$

which should become (not sure about this)

$$
\begin{aligned}
\frac{1}{\Gamma_1} + \frac{1}{\Gamma_2} &= 1 \\
\frac{1}{\Gamma_1} + \frac{1}{\Gamma_2} + \frac{1}{\Gamma_3} &= 1 \\
\frac{1}{\Gamma_2} + \frac{1}{\Gamma_3} &= 1
\end{aligned}
$$

$$
\gamma = \frac{(x - x_{i+1})(x - x_{i+2})}{a} +
         \frac{(x - x_{i-2})(x - x_{i+2})}{b} +
         \frac{(x - x_{i-1})(x - x_{i-2})}{c}
$$

Substituting $x = x_{i-2}$ gives

$$
\frac{(x_{i-2} - x_{i+1})(x_{i-2} - x_{i+2})}{a} = 1
$$

and thus

$$
{a} = {(x_{i-2} - x_{i+1})(x_{i-2} - x_{i+2})}
$$

Substituting $x = x_{i-1}$ gives

$$
\frac{(x_{i-1} - x_{i-2})(x_{i-1} - x_{i+2})}{b} +
\frac{(x_{i-1} - x_{i+1})(x_{i-1} - x_{i+2})}{(x_{i-2} - x_{i+1})(x_{i-2} - x_{i+2})} = 1
$$

and thus

$$
b = \frac{(x_{i-1} - x_{i-2})(x_{i-1} - x_{i+2})}
         {1 - \frac{(x_{i-1} - x_{i+1})(x_{i-1} - x_{i+2})}{(x_{i-2} - x_{i+1})*(x_{i-2} - x_{i+2})}}
$$

Substituting $x = x_{i}$ gives

$$
1 = \frac{(x_{i} - x_{i+1})(x_{i} - x_{i+2})}{a} +
    \frac{(x_{i} - x_{i-2})(x_{i} - x_{i+2})}{b} +
    \frac{(x_{i} - x_{i-1})(x_{i} - x_{i-2})}{c}
$$

which on some simplification results in

$$
c = (x_{i-1} - x_{i+2})(x_{i-2} - x_{i+2})
$$

Of course we could have substituted $x = x_{i+2}$ and given ourselves
  an easier calculation but this gives a good check that the
  derivation is correct.

\begin{code}
coeff :: Integral a => Ratio a -> Ratio a -> Ratio a -> Ratio a
coeff k x j = sum [ num m / den m | m <- [j + 1 .. k] ]
  where
    den m = product [ m - l | l <- [0 .. k], l /= m]
    num m = sum [ product [ x - q + 1 | q <- [0 .. k], q /= l, q /= m]
                | l <- [0 .. k], l /= m ]

coeffs :: Integral a => Ratio a -> [Ratio a]
coeffs k = [ coeff k x j | x <- [-1 .. k - 1], j <- [0 .. k - 1] ]
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

data Rate a = Rate { theta1 :: a, theta2 :: a, theta3 :: a, theta4 :: a }
  deriving Show

meanRate :: Floating a => Rate a
meanRate = Rate { theta1 = 0.5, theta2 = 0.025, theta3 = 0.8, theta4 = 0.025 }

dzdt :: (Show a, Floating a) => Rate a -> a -> [a] -> [a]
dzdt x _t [u, v] = [ (alpha - beta * v) * u
                   , (-gamma + delta * u) * v
                   ]
  where
    alpha = theta1 x
    beta = theta2 x
    delta = theta4 x
    gamma = theta3 x
dzdt _x _t xs = error $ show xs

lotkaVolterra :: Double -> Double -> OdeProblem
lotkaVolterra h l = emptyOdeProblem
  { odeRhs = odeRhsPure $ \t x -> fromList (dzdt meanRate t (toList x))
  , odeJacobian = Nothing
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = vector [h, l]
  , odeSolTimes = vector us
  , odeTolerances = defaultTolerances
  }

us :: [Double]
us = map (* 0.05) $ map fromIntegral ([0 .. 400] :: [Int])

vs :: [Double]
vs = map (+ 1900) us

sol''' :: Double -> Double -> IO (Matrix Double)
sol''' h l = do
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  logEnv <- registerScribe "stdout" handleScribe defaultScribeSettings =<< initLogEnv "namespace" "devel"
  w <- runKatipT logEnv $ solve (defaultOpts' BOGACKI_SHAMPINE_4_2_3) (lotkaVolterra h l)
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

initPop :: (Double, Double)
initPop = (snd $ fst $ head predPreyObs, snd $ snd $ head $ predPreyObs)

main1 :: IO ()
main1 = do
  let hs, ks, ts :: [Double]
      ts = map fromIntegral $ map fst $ map fst predPreyObs
      hs = map snd $ map fst predPreyObs
      ks = map snd $ map snd predPreyObs
  y <- sol''' (hs!!0) (ks!!0)
  let rs = transpose $ toLists y
  R.runRegion $ do
    _ <- [r| library(ggplot2) |]
    c0 <- [r| ggplot() + ggtitle("Hares and Lynxes") + xlab("Year") + ylab("Animals (000s)") |]
    _  <- foldM (\c (n, rr, tt) -> do
                        [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
                c0 ([("Hares", hs, ts), ("Lynxes", ks, ts),
                     ("Predicted Hares", rs!!0, vs), ("Predicted Lynxes", rs!!1, vs)] :: [(String, [Double], [Double])])
    _ <- [r| ggsave(filename="diagrams/HudsonBay.png") |]
    return ()

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

instance Element Int8

emptyOdeProblem :: OdeProblem
emptyOdeProblem = OdeProblem
      { odeRhs = error "emptyOdeProblem: no odeRhs provided"
      , odeJacobian = Nothing
      , odeInitCond = error "emptyOdeProblem: no odeInitCond provided"
      , odeEventDirections = V.empty
      , odeEventConditions = EventConditionsHaskell V.empty
      , odeTimeBasedEvents = TimeEventSpec $ return $ 1.0 / 0.0
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
