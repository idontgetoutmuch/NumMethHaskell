% Haskell for Numerics?
% Dominic Steinitz
% 2nd June 2017

> {-# OPTIONS_GHC -Wall #-}

> import           Numeric.Sundials.ARKode.ODE
> import           Numeric.LinearAlgebra
> import           Data.Csv
> import           Data.Char
> import           Data.ByteString.Lazy (putStr, writeFile)
> import           Prelude hiding (putStr, writeFile)

With one spatial dimension we have:

$$
u_{t}=k u_{x x} + f
$$

initial condition $u(0, x)=0$

Dirichlet boundary conditions

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
\dot{u_2} \\
\dot{u_3} \\
\dot{u_4}
\end{bmatrix}
=
\begin{bmatrix}
0 &  0 &  0 &  0 & 0 \\
1 & -2 &  1 &  0 & 0 \\
0 &  1 & -2 &  1 & 0 \\
0 &  0 &  1 & -2 & 1 \\
0 &  0 &  0 &  0 & 0
\end{bmatrix}
\begin{bmatrix}
u_0 \\
u_1 \\
u_2 \\
u_3 \\
u_4
\end{bmatrix}
+
\begin{bmatrix}
f_0 \\
f_1 \\
f_2 \\
f_3 \\
f_4
\end{bmatrix}
$$

where $f_j \triangleq f(t, x_j)$.

spatial mesh size

> bigN :: Int
> bigN = 201

heat conductivity

> k :: Double
> k = 0.5

mesh spacing

> deltaX :: Double
> deltaX = 1.0 / (fromIntegral bigN - 1)
> c1, c2 :: Double
> c1 = k / deltaX / deltaX
> c2 = (-2.0) * k / deltaX / deltaX

initial time

> t0 :: Double
> t0 = 0.0

final time

> tf :: Double
> tf =1.0

total number of output times

> bigNt :: Int
> bigNt = 10

relative tolerance

> rtol :: Double
> rtol = 1.0e-6

absolute tolerance

> atol :: Double
> atol = 1.0e-10

> bigA :: Matrix Double
> bigA = assoc (bigN, bigN) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. bigN - 1]
>                                                    , j <- [0 .. bigN - 1]
>                         ]
>  where
>    f (i, j) | i     == 0        = 0.0    -- left boundary condition
>             | i     == bigN - 1 = 0.0    -- right boundary condition
>             | i     == j        = c2
>             | i - 1 == j        = c1
>             | i + 1 == j        = c1
>             | otherwise         = 0.0

> b :: Vector Double
> b = assoc bigN 0.0 [ (iSource, 0.01 / deltaX) ]
>   where
>     iSource = bigN `div` 2

> bigU0 :: Vector Double
> bigU0 = assoc bigN 0.0 []
>
> deltaT :: Double
> deltaT = (tf - t0) / (fromIntegral bigNt)

> sol :: Matrix Double
> sol = odeSolveV SDIRK_5_3_4' Nothing rtol atol (const bigU') bigU0 (vector (map (deltaT *) [0 .. 10]))
>   where
>     bigU' bigU = bigA #> bigU + b

> myOptions :: EncodeOptions
> myOptions = defaultEncodeOptions {
>       encDelimiter = fromIntegral (ord ' ')
>     }

> main :: IO ()
> main = do
>   writeFile "heat1E.txt" $ encodeWith myOptions $ map toList $ toRows sol

