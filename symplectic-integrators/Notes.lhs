% Data Sources
% Dominic Steinitz
% 18th February 2018


Introduction
============

These are some very hasty notes on Runge-Kutta methods and IRK2 in
particular. I make no apologies for missing lots of details.

Some Uncomprehensive Theory
===========================

In general, an implicit Runge-Kutta method is given by

$$
y_{n+1} = y_n + h \sum_{i=1}^s b_i k_i
$$

where

$$
k_i = f\left( t_n + c_i h,\ y_{n} + h \sum_{j=1}^s a_{ij} k_j \right), \quad i = 1, \ldots, s
$$

and

$$
\sum_{j=1}^{i-1} a_{ij} = c_i \text{ for } i=2, \ldots, s
$$

Traditionally this is written as a Butcher tableau:

$$
\begin{array}{c|cccc}
c_1    & a_{11} & a_{12}& \dots & a_{1s}\\
c_2    & a_{21} & a_{22}& \dots & a_{2s}\\
\vdots & \vdots & \vdots& \ddots& \vdots\\
c_s    & a_{s1} & a_{s2}& \dots & a_{ss} \\
\hline
       & b_1    & b_2   & \dots & b_s\\
       & b^*_1  & b^*_2 & \dots & b^*_s\\
\end{array}
$$

and even more succintly as

$$
\begin{array}{c|c}
\mathbf{c}& A\\
\hline
          & \mathbf{b^T} \\
\end{array}
$$

For a Gauß-Legendre method we set the values of the $c_i$ to the zeros
of the shifted Legendre polynomials.

An explicit expression for the shifted Legendre polynomials is given by

$$
\tilde{P}_n(x) = (-1)^n \sum_{k=0}^n \binom{n}{k} \binom{n+k}{k} (-x)^k
$$

The first few shifted Legendre polynomials are:

$$
\begin{array}{r|r}
n & \tilde{P}_n(x) \\
\hline
0 & 1 \\
1 & 2x-1 \\
2 & 6x^2-6x+1 \\
3 & 20x^3-30x^2+12x-1 \\
4 & 70x^4-140x^3+90x^2-20x+1 \\
5 & 252x^5 -630x^4 +560x^3 - 210 x^2 + 30 x -1
\end{array}
$$

Setting

$$ q(t) \triangleq \prod_{j=0}^s (t - c_j) \quad {\text and} \quad q_l
\triangleq \frac{q(t)}{t - c_l}, \, l = 1,2, \ldots, s $$

then the co-location method gives

$$
a_{j,i} \triangleq \int_0^{c_j} \frac{q_i(\tau)}{q_i(c_i)}\,{\mathrm d}\tau
$$

$$
b_{j} \triangleq \int_0^{1} \frac{q_i(\tau)}{q_i(c_i)}\,{\mathrm d}\tau
$$

For $s = 1$ we have $2x - 1 = 0$ and thus $c_1 = 1/2$ and the Butcher tableau is

$$
\begin{array}{c|c}
1/2    & 1/2\\
\hline
       & 1 \\
\end{array}
$$

that is, the implicit RK2 method aka the mid-point method.

For $s = 2$ we have $6x^2-6x+1 = 0$ and thus $c_1 = \frac{1}{2} -
\frac{1}{2\sqrt{3}}$ and $c_2 = \frac{1}{2} + \frac{1}{2\sqrt{3}}$ and
the Butcher tableau is

$$
\begin{array}{c|cc}
\frac{1}{2} - \frac{1}{2\sqrt{3}} & \frac{1}{4} & \frac{1}{4} - \frac{1}{2\sqrt{3}}\\
\frac{1}{2} + \frac{1}{2\sqrt{3}} & \frac{1}{4} + \frac{1}{2\sqrt{3}} & \frac{1}{4}\\
\hline
       & \frac{1}{2} & \frac{1}{2}\\
\end{array}
$$

that is, the implicit RK4 method.


Implicit RK2
============

Explicitly

$$
y_{n+1} = y_n + h b_1 k_1
$$

where

$$
k_1 = f\left( t_n + c_1 h,\ y_{n} + h a_{11} k_1 \right)
$$

Substituting in the values from the tableau, we have

$$
y_{n+1} = y_n + h k_1
$$

where

$$
k_1 = f\left( t_n + \frac{1}{2} h,\ y_{n} + h \frac{1}{2} k_1 \right)
$$

and further inlining and substitution gives

$$
y_{n+1}=y_n+hf(t+\frac{h}{2},\frac{1}{2}(y_n+y_{n+1}))
$$

which can be recognised as the mid-point method.

Implementation
==============

A common package for solving ODEs is
[gsl](https://www.gnu.org/software/gsl/) which Haskell interfaces via
the [hmatrix-gsl](https://hackage.haskell.org/package/hmatrix-gsl)
package. Some of gsl is coded with the conditional compilation flag
*DEBUG* e.g. in
[msbdf.c](https://github.com/idontgetoutmuch/gsl/blob/953e81673caa21d690e5d72594bc4cd2a60ba311/ode-initval2/msbdf.c#L581)
but sadly not in the simpler methods maybe because they were made part
of the package some years earlier. We can add our own of course;
here's a
[link](https://github.com/idontgetoutmuch/gsl/commit/c2035977d65cd804169ff3370da6723cf879be75)
for reference.

Let's see how the Implicit Runge-Kutta Order 2 method does on the
following system taken from the [gsl
documenation](https://www.gnu.org/software/gsl/manual/html_node/ODE-Example-programs.html)

$$
\frac{{\mathrm d}^2u}{{\mathrm d} t^2} + \mu \frac{{\mathrm d}^u}{{\mathrm d} t} (u^2 - 1) + u = 0
$$

which can be re-written as

$$
\begin{aligned}
\frac{{\mathrm d}y_0}{{\mathrm d} t} &= y_1 \\
\frac{{\mathrm d}y_1}{{\mathrm d} t} &= -y_0 - \mu y_1 (y_0 y_0 - 1)
\end{aligned}
$$

but replacing *gsl_odeiv2_step_rk8pd* with *gsl_odeiv2_step_rk2imp*.

Here's the first few steps


```
rk2imp_apply: t=0.00000e+00, h=1.00000e-06, y:1.00000e+00 0.00000e+00 
-- evaluate jacobian
(  2,  2)[0,0]:                      0
(  2,  2)[0,1]:                      1
(  2,  2)[1,0]:                     -1
(  2,  2)[1,1]:                     -0
YZ:1.00000e+00 -5.00000e-07 
-- error estimates: 0.00000e+00 8.35739e-20 
rk2imp_apply: t=1.00000e-06, h=5.00000e-06, y:1.00000e+00 -1.00000e-06 
-- evaluate jacobian
(  2,  2)[0,0]:                      0
(  2,  2)[0,1]:                      1
(  2,  2)[1,0]:   -0.99997999999999998
(  2,  2)[1,1]: 1.00008890058234101e-11
YZ:1.00000e+00 -3.50000e-06 
-- error estimates: 1.48030e-16 1.04162e-17 
rk2imp_apply: t=6.00000e-06, h=2.50000e-05, y:1.00000e+00 -6.00000e-06 
-- evaluate jacobian
(  2,  2)[0,0]:                      0
(  2,  2)[0,1]:                      1
(  2,  2)[1,0]:  -0.999880000000002878
(  2,  2)[1,1]: 3.59998697518904009e-10
YZ:1.00000e+00 -1.85000e-05 
-- error estimates: 1.48030e-16 1.30208e-15 
rk2imp_apply: t=3.10000e-05, h=1.25000e-04, y:1.00000e+00 -3.10000e-05 
-- evaluate jacobian
(  2,  2)[0,0]:                      0
(  2,  2)[0,1]:                      1
(  2,  2)[1,0]:  -0.999380000000403723
(  2,  2)[1,1]: 9.6099972424212865e-09
YZ:1.00000e+00 -9.35000e-05 
-- error estimates: 0.00000e+00 1.62760e-13 
rk2imp_apply: t=1.56000e-04, h=6.25000e-04, y:1.00000e+00 -1.56000e-04 
-- evaluate jacobian
(  2,  2)[0,0]:                      0
(  2,  2)[0,1]:                      1
(  2,  2)[1,0]:  -0.996880000051409643
(  2,  2)[1,1]: 2.4335999548874554e-07
YZ:1.00000e+00 -4.68500e-04 
-- error estimates: 1.55431e-14 2.03450e-11 
```

Let's see if we can reproduce this in a fast and loose way in Haskell

To make our code easier to write and read let's lift some arithmetic
operators to act on lists (we should really use the
[Naperian](https://hackage.haskell.org/package/Naperian) package).

> module RK2Imp where

> instance Num a => Num [a] where
>   (+) = zipWith (+)
>   (*) = zipWith (*)

The actual algorithm is almost trivial

> rk2Step :: Fractional a => a -> a -> [a] -> [a] -> [a]
> rk2Step h t y_n0 y_n1 = y_n0 + (repeat h) * f (t + 0.5 * h) ((repeat 0.5) * (y_n0 + y_n1))

The example Van der Pol oscillator with the same parameter and initial
conditions as in the gsl example.

> f :: Fractional b => a -> [b] -> [b]
> f t [y0, y1] = [y1, -y0 - mu * y1 * (y0 * y0 - 1)]
>   where
>     mu = 10.0

> y_init :: [Double]
> y_init = [1.0, 0.0]

We can use co-recursion to find the fixed point of the Runge-Kutta
equation arbitrarily choosing the 10th iteration!

> y_1s, y_2s, y_3s :: [[Double]]
> y_1s = y_init : map (rk2Step 1.0e-6 0.0 y_init) y_1s

> nC :: Int
> nC = 10

> y_1, y_2, y_3 :: [Double]
> y_1 = y_1s !! nC

This gives

    [ghci]
    y_1

which is not too far from the value given by gsl.

```
y:1.00000e+00 -1.00000e-06 
```

Getting the next time and step size from the gsl *DEBUG* information
we can continue

> y_2s = y_1 : map (rk2Step 5.0e-6 1.0e-6 y_1) y_2s
> 
> y_2 = y_2s !! nC

    [ghci]
    y_2

gsl gives

```
y:1.00000e+00 -6.00000e-06
```

> y_3s = y_2 : map (rk2Step 2.5e-5 6.0e-6 y_2) y_3s
> y_3 = y_3s !! nC

One more time

    [ghci]
    y_3

And gsl gives

```
y:1.00000e+00 -3.10000e-05
```
