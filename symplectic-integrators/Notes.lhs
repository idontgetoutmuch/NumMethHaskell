% Data Sources
% Dominic Steinitz
% 18th February 2018

Introduction
============

For the midpoint method you approximate 

\begin{align}
\frac{y(t+h)-y(t)}h&\approx\dot y(t+\frac h2)
\\
&=f(t+\frac{h}{2},y(t+\frac{h}{2}))
\\
&\approx f(t+\frac h2,\frac12(y(t+h)+y(t))).
\end{align}

All approximations with an error $O(h^2)$

This gives the method
$$
y_{k+1}=y_k+hf(t+\frac{h}{2},\frac{1}{2}(y_k+y_{k+1}))
$$

This can be decomposed into first one step of the implicit Euler
method and then one of the explicit Euler method, both with half the
step size

\begin{align}
\text{implicit Euler: }&&y_{k+\frac12}&=y_k+\frac{h}{2} f(t+\frac{h}{2},y_{k+\frac12})\\
\text{explicit Euler: }&&y_{k+1}&=y_{k+\frac12}+\frac h2 f(t+\frac h2,y_{k+\frac12})
\end{align}

where $(t+\frac h2,y_{k+\frac12})$ is the midpoint of the segment
connecting the iteration points with
$y_{k+\frac12}=\frac12(y_k+y_{k+1})$

Other
=====

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

