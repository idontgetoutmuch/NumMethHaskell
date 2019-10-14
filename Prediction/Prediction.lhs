> {-# OPTIONS_GHC -Wall #-}

> {-# LANGUAGE DataKinds #-}
> {-# LANGUAGE QuasiQuotes #-}
> {-# LANGUAGE LambdaCase #-}
> {-# LANGUAGE GADTs #-}
> {-# LANGUAGE ScopedTypeVariables #-}
> {-# LANGUAGE OverloadedStrings #-}

> import           Numeric.Sundials.ARKode.ODE
> import           Numeric.LinearAlgebra

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

Notice that there is nothing stochastic about the biology but we express our uncertainty about the parameter
by making it a time-varying stochastic variable which says the further we go into the future the less certain we are about it.

We are going to turn this into a Fokker-Planck equation which we can then solve using e.g. the method of lines. But before turning to Fokker-Planck, let's show that we can indeed solve a diffusion equation using the method of lines.

Let us solve the heat equation over the unit square to some arbitrary point in the future.


```haskell
x1, a, x2 :: Double
x1 = 0
a = 1.0
x2 = a

y1, y2 :: Double
y1 = 0.0
y2 = 1.0

bigT :: Double
bigT = 1000.0
```


```haskell
n :: Int
n = 5

dx :: Double
dx = a / (fromIntegral n + 1)

dy :: Double
dy = a / (fromIntegral n + 1)

beta = 1.0e-5
s = beta / (dx * dx)
```

We discretize the space variables. I probably ought to show how this is done.


```haskell
bigA :: Matrix Double
bigA = assoc (n, n) 0.0 [((i, j), f (i, j)) | i <- [0 .. n - 1]
                                           , j <- [i-1..i+1]
                                           , j `elem` [0..n-1]]
 where
   f (i, j) | i     == j = (-2.0) * s -- diagonal
            | i - 1 == j = s          -- below diagonal
            | i + 1 == j = s          -- above diagonal

bigAA1 :: Matrix Double
bigAA1 = assoc (n * n, n * n) 0.0 [((i, j), f (i, j)) | i <- [0 .. n * n - 1]
                                                     , j <- [i - n, i,  i + n]
                                                     , j `elem` [0 .. n * n -1]]
 where
   f (i, j) | i     == j = (-2.0)
            | i - n == j = 1.0
            | i + n == j = 1.0
            | otherwise = error $ show (i, j)

bigAA2 :: Matrix Double
bigAA2 = diagBlock (replicate n bigA)
 where
   bigA :: Matrix Double
   bigA = assoc (n, n) 0.0 [((i, j), f (i, j)) | i <- [0 .. n - 1]
                                               , j <- [i-1..i+1]
                                               , j `elem` [0..n-1]]
     where
       f (i, j) | i     == j = (-2.0)
                | i - 1 == j = 1.0
                | i + 1 == j = 1.0

bigAA :: Matrix Double
bigAA = cmap (*s) $ bigAA1 + bigAA2
```


```haskell

```

Set up the initial conditions
-----------------------------

Since we know that a solution is the normal distribution the
variance of which gets larger as time increases, we can start with a normal distribution.


```haskell
bb :: Vector Double
bb = assoc (n * n) 0.0 []

x = linspace n (x1 + dx, x2 - dx)
y = linspace n (y1 + dy, y2 - dy)
```


```haskell
sigma :: Double
sigma = 0.05

bigUU0 = vector $
        map (\(u, v) -> exp (-((u - mux)^2 + (v - muy)^2) / (2 * sigma^2)) / (2 * pi * sigma^2))
            [(u, v) | u <- toList x, v <- toList y]
 where
   mux = (x2 - x1) / 2
   muy = (y2 - y1) / 2
```

Now we can solve our system with a rather solver and solver parameters.


```haskell
sol = odeSolveV SDIRK_5_3_4' Nothing 1.0e-6 1.0e-10 (const bigUU') bigUU0 (vector [0, bigT])
  where
    bigUU' bigUU = bigAA #> bigUU + bb
```


```haskell
myTraceT = line (aes & GP.x .~ fst
                     & GP.y .~ snd) (zip (toList x) (toList middle))
 where
   middle = toRows (reshape n (toRows soll !! 1)) !! (n `div` 2)
```


```haskell

```


```haskell
myTrace  = line (aes & GP.x .~ fst
                     & GP.y .~ snd) (zip (toList x) (toList middle))
 where
   middle = toRows (reshape n bigUU0) !! (n `div` 2)
```


```haskell
main = do
 T.writeFile "test0.html" $ renderText $ doctypehtml_ $ do
   head_ $ do meta_ [charset_ "utf-8"]
              plotlyCDN
   body_ $ toHtml $ plotly "myDiv" [myTrace]

 T.writeFile "testT.html" $ renderText $ doctypehtml_ $ do
   head_ $ do meta_ [charset_ "utf-8"]
              plotlyCDN
   body_ $ toHtml $ plotly "myDiv" [myTraceT]
```


```haskell
main
```


    


For now you will have to open the html files yourself. I think I know how to fix this but one step at a time.
