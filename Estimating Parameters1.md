
# Estimating Parameters in Chaotic Systems

## Introduction

Suppose we obseve a dynamical system and we wish to estimate its
parameters. One way of doing this is to use a Monte Carlo
method. But if the system is chaotic then this approach is highly
unlikely to work as the small changes in the proposals for the
parameters will result in large changes in the path.

Here's a [well-known chaotic dynamical system](https://en.wikipedia.org/wiki/Lorenz_system):

\begin{align}
\frac{\mathrm{d}x}{\mathrm{d}t} &= \alpha (y - x), \\
\frac{\mathrm{d}y}{\mathrm{d}t} &= x (\rho - z) - y, \\
\frac{\mathrm{d}z}{\mathrm{d}t} &= x y - \beta z.
\end{align}

And here is its representation in [libbi](http://libbi.org) together with some
instructions to generate noisy values of the $x$ variable (we use
the noisy values later).


```R
library(readr)
model_file_name <- "LorenzGenerate.bi"
writeLines(read_file(model_file_name))
```

    model Lorenz {
      const rho = 45.92
      const beta = 4.0
      const alpha = 16.0
    
      state X, Y, Z
      obs X_obs
    
      sub initial {
        X ~ log_normal(log(1.0), 0.00002)
        Y ~ log_normal(log(1.0), 0.00002)
        Z ~ log_normal(log(1.0), 0.00002)
      }
    
      sub transition(delta = 0.0001) {
        ode {
          dX/dt = alpha * (Y - X)
          dY/dt = X * (rho - Z) - Y
          dZ/dt = X * Y - beta * Z
        }
      }
    
      sub observation {
        X_obs ~ normal(X, 0.2)
      }
    }
    



```R
library('rbi')
library(ggplot2)

Lorenz <- bi_model(model_file_name)

T <- 10.0
nObs <- 100
init_parameters <- list(X = 1, Y = 1, Z = 1)
```


```R
synthetic_dataset <- bi_generate_dataset(end_time=T, model=Lorenz,
                                         init=init_parameters,
                                         noutputs = nObs)
```


```R
synthetic_data <- bi_read(synthetic_dataset)
synthetic_df <- as.data.frame(synthetic_data)
tail(synthetic_df)
```




<table>
<thead><tr><th></th><th scope=col>X.time</th><th scope=col>X.value</th><th scope=col>Y.time</th><th scope=col>Y.value</th><th scope=col>Z.time</th><th scope=col>Z.value</th><th scope=col>X_obs.time</th><th scope=col>X_obs.value</th><th scope=col>clock</th></tr></thead>
<tbody>
	<tr><th scope=row>96</th><td> 9.5      </td><td>-14.764658</td><td> 9.5      </td><td>-19.536008</td><td> 9.5      </td><td>39.76409  </td><td> 9.5      </td><td>-15.044416</td><td>96309     </td></tr>
	<tr><th scope=row>97</th><td> 9.6      </td><td>-17.611354</td><td> 9.6      </td><td>-14.413804</td><td> 9.6      </td><td>53.75645  </td><td> 9.6      </td><td>-17.397395</td><td>96309     </td></tr>
	<tr><th scope=row>98</th><td> 9.7      </td><td> -9.657817</td><td> 9.7      </td><td> -5.782971</td><td> 9.7      </td><td>45.82000  </td><td> 9.7      </td><td>-10.099852</td><td>96309     </td></tr>
	<tr><th scope=row>99</th><td> 9.8      </td><td> -8.019435</td><td> 9.8      </td><td> -9.629767</td><td> 9.8      </td><td>35.46636  </td><td> 9.8      </td><td> -8.105424</td><td>96309     </td></tr>
	<tr><th scope=row>100</th><td> 9.9      </td><td>-14.220221</td><td> 9.9      </td><td>-19.872234</td><td> 9.9      </td><td>37.34805  </td><td> 9.9      </td><td>-14.511350</td><td>96309     </td></tr>
	<tr><th scope=row>101</th><td>10.0      </td><td>-18.793750</td><td>10.0      </td><td>-15.760650</td><td>10.0      </td><td>55.08816  </td><td>10.0      </td><td>-18.686090</td><td>96309     </td></tr>
</tbody>
</table>





```R
p <- ggplot(synthetic_df, aes(X.time)) +
    geom_path(aes(y = X.value, colour="alpha 16.0")) +
    theme(legend.position="bottom") +
    ggtitle("Lorenz") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Time") +
    ylab("X Value")
ggsave(filename = "diagrams/xpath.svg", plot = p)
```

    Saving 10 x 5 in image


![title](diagrams/xpath.svg)

We can check that the system becomes chaotic, in the sense that small
changes in initial conditions lead to qualitatively different
behaviour.


```R
path0 <- ggplot() +
    theme(legend.position="bottom") +
    ggtitle("Lorenz") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Time") +
    ylab("Value")


set.seed(42)

T <- 20.0

for (i in c("red", "blue", "green")) {
    init_parameters <- list(X = 1 + rnorm(1,0.0,0.01),
                            Y = 1 + rnorm(1,0.0,0.01),
                            Z = 1 + rnorm(1,0.0,0.01))

    synthetic_dataset <- bi_generate_dataset(end_time=T, model=Lorenz,
                                             init=init_parameters,
                                             noutputs = nObs)

    synthetic_data <- bi_read(synthetic_dataset)
    synthetic_df <- as.data.frame(synthetic_data)

    path0 <- path0 +
        geom_line(data = synthetic_df, aes(x = X.time, y = X.value), color = i)
}
```


```R
ggsave(filename = "diagrams/xpath4.svg", plot = path0)
```

    Saving 10 x 5 in image


![title](diagrams/xpath4.svg)

## Parameters as State

Alternatively we can model the parameter as part of the state and
assume that it undergoes Brownian Motion. This seems reasonable: the
further we go into the future, the less certain we are about its
value. An improvement might be to model it as an Ornstein-Uhlenbeck
process which is mean-reverting - after all we don't expect the
parameter to take on arbitrarily large or small values but
let's learn to walk before we learn to run.

Since we expect our parameter to positive let's model it as geometric
Brownian Motion.

\begin{equation}
\mathrm{d}\alpha = \alpha\sigma\mathrm{d}W_t
\end{equation}

By [Itô's lemma](https://en.wikipedia.org/wiki/Itô%27s_lemma#Geometric_Brownian_motion),
we have

\begin{equation}
\mathrm{d}(\log \alpha) = -\frac{\sigma^2}{2}\mathrm{d}t + \sigma\mathrm{d}W_t
\end{equation}

We can model this in libbi as:


```R
model_file_name <- "LorenzState.bi"
writeLines(read_file(model_file_name))
```

    model LorenzState {
      const rho = 45.92
      const beta = 4
    
      const h         = 0.1;    // time step
      const delta_abs = 1.0e-3; // absolute error tolerance
      const delta_rel = 1.0e-6; // relative error tolerance
    
      state X, Y, Z
      state ln_alpha
    
      param mu, sigma
    
      noise w
    
      obs X_obs
    
      sub parameter {
        mu ~ uniform(12.0, 20.0)
        sigma ~ uniform(0.0, 0.5)
      }
    
      sub proposal_parameter {
         mu ~ truncated_gaussian(mu, 0.02, 12.0, 20.0);
         sigma ~ truncated_gaussian(sigma, 0.01, 0.0, 0.5);
       }
    
      sub initial {
        X ~ log_normal(log(1.0), 0.2)
        Y ~ log_normal(log(1.0), 0.2)
        Z ~ log_normal(log(1.0), 0.2)
        ln_alpha ~ gaussian(log(mu), sigma)
      }
    
      sub transition(delta = h) {
        w ~ normal (0.0, sqrt(h))
        ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg =  'RK4(3)') {
          dX/dt = exp(ln_alpha) * (Y - X)
          dY/dt = X * (rho - Z) - Y
          dZ/dt = X * Y - beta * Z
          dln_alpha/dt = -sigma * sigma / 2 - sigma * w / h
        }
      }
    
      sub observation {
        X_obs ~ normal(X, 0.2)
      }
    }
    


We can then run the model.


```R
LorenzState <- bi_model(model_file_name)

bi_state_model <- libbi(model=LorenzState)
bi_state <- filter(bi_state_model,
                   nparticles = 8192,
                   nthreads = 1,
                   end_time = T,
                   obs = synthetic_dataset,
                   init = init_parameters,
                   ess_rel = 1,
                   sample_obs = TRUE)

bi_file_summary(bi_state$output_file_name)
bi_state
summary(bi_state)
```

    File /private/var/folders/h1/bkwn2ct12hvb6__dzrk5gwx80000gn/T/RtmpQ7zb3k/LorenzState1176a1bafc6c5/LorenzState_output1176a1d4a68c2.nc (NC_FORMAT_NETCDF4):
    
         13 variables (excluding dimension variables):
            double time[nr]   (Contiguous storage)  
            double w[np,nr]   (Contiguous storage)  
            double X[np,nr]   (Contiguous storage)  
            double Y[np,nr]   (Contiguous storage)  
            double Z[np,nr]   (Contiguous storage)  
            double ln_alpha[np,nr]   (Contiguous storage)  
            double mu[]   (Contiguous storage)  
            double sigma[]   (Contiguous storage)  
            8 byte int clock[]   (Contiguous storage)  
            int ancestor[np,nr]   (Contiguous storage)  
            double logweight[np,nr]   (Contiguous storage)  
            double loglikelihood[]   (Contiguous storage)  
            double X_obs[np,nr]   (Contiguous storage)  
    
         2 dimensions:
            np  Size:8192
            nr  Size:101
    
        3 global attributes:
            libbi_schema: ParticleFilter
            libbi_schema_version: 1
            libbi_version: 1.4.1





<table>
<thead><tr><th></th><th scope=col>Min.</th><th scope=col>1st Qu.</th><th scope=col>Median</th><th scope=col>Mean</th><th scope=col>3rd Qu.</th><th scope=col>Max.</th></tr></thead>
<tbody>
	<tr><th scope=row>mu</th><td>14.1060017</td><td>14.1060017</td><td>14.1060017</td><td>14.1060017</td><td>14.1060017</td><td>14.1060017</td></tr>
	<tr><th scope=row>sigma</th><td> 0.2652258</td><td> 0.2652258</td><td> 0.2652258</td><td> 0.2652258</td><td> 0.2652258</td><td> 0.2652258</td></tr>
</tbody>
</table>




And then after some manipulation, draw the path of the "state" parameter.


```R
output <- bi_read(bi_state)
logw <- xtabs(value ~ time + np, data = output$logweight, addNA = TRUE)
X <- output$X$value
Y <- output$Y$value
Z <- output$Z$value
A <- output$ln_alpha$value
```


```R
log2normw <- function(lw){
  w <- exp(lw - max(lw))
  return(w / sum(w))
}
```


```R
w = t(apply(X=logw, MARGIN=1, FUN=log2normw))
Xmeans = apply(X = X*w, MARGIN=1, FUN=sum)
Ymeans = apply(X = X*w, MARGIN=1, FUN=sum)
Zmeans = apply(X = Z*w, MARGIN=1, FUN=sum)
Ameans = apply(X = A*w, MARGIN=1, FUN=sum)

```


```R
synthetic_data <- bi_read(synthetic_dataset)
X_original <- synthetic_data$X$value
Y_original <- synthetic_data$Y$value
Z_original <- synthetic_data$Z$value

```


```R
synthetic_df <- as.data.frame(synthetic_data)
synthetic_df$Xmeans <- Xmeans
synthetic_df$Ymeans <- Ymeans
synthetic_df$Zmeans <- Zmeans
synthetic_df$Ameans <- Ameans
```


```R
pAmeans <- ggplot(synthetic_df, aes(X.time)) +
           geom_path(aes(y = exp(Ameans), colour="Ameans")) +
           theme(legend.position="bottom") +
           ggtitle("Lorenz") +
           theme(plot.title = element_text(hjust = 0.5)) +
           ylim(0.0, max(exp(Ameans))) +
           xlab("Time") +
           ylab("Value")

```

    simpleWarning in if (class(res$value) == "help_files_with_topic") {: the condition has length > 1 and only the first element will be used
    
    



```R
ggsave(filename = "diagrams/xpath5.svg", plot = pAmeans)
```

    Saving 10 x 5 in image


![title](diagrams/xpath5.svg)

We can try inferring the parameter for the different chaotic solutions.


```R
dataset_list <- list()
parameters_list <- list()
```


```R
for (i in c(1,2,3)) {
    init_parameters <- list(X = 1 + rnorm(1,0.0,0.01),
                            Y = 1 + rnorm(1,0.0,0.01),
                            Z = 1 + rnorm(1,0.0,0.01))

    parameters_list[[i]] <- init_parameters
    synthetic_dataset <- bi_generate_dataset(end_time=T, model=Lorenz,
                                             init=init_parameters,
                                             noutputs = nObs)

    dataset_list[[i]] <- synthetic_dataset
}
```


```R
X_list <- list()
Y_list <- list()
Z_list <- list()
A_list <- list()
```


```R
 for (i in c(1,2,3)) {
    bi_state <- filter(bi_state_model, nparticles = 8192, nthreads = 1, end_time = T, obs = dataset_list[[i]], init = parameters_list[[i]], ess_rel = 1, sample_obs = TRUE)
    output <- bi_read(bi_state)
    logw <- xtabs(value ~ time + np, data = output$logweight, addNA = TRUE)
    w = t(apply(X=logw, MARGIN=1, FUN=log2normw))
    X <- output$X$value
    Y <- output$Y$value
    Z <- output$Z$value
    A <- output$ln_alpha$value
    X_list[[i]] = apply(X = X*w, MARGIN=1, FUN=sum)
    Y_list[[i]] = apply(X = X*w, MARGIN=1, FUN=sum)
    Z_list[[i]] = apply(X = Z*w, MARGIN=1, FUN=sum)
    A_list[[i]] = apply(X = A*w, MARGIN=1, FUN=sum)
}

```


```R
path2 <- ggplot() +
    theme(legend.position="bottom") +
    ggtitle("Lorenz") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Time") +
    ylab("Value")

```

    simpleWarning in if (class(res$value) == "help_files_with_topic") {: the condition has length > 1 and only the first element will be used
    
    



```R
for (i in c(1,2,3)) {
        synthetic_data <- bi_read(dataset_list[[i]])
        synthetic_df <- as.data.frame(synthetic_data)
        synthetic_df$Ameans <- exp(A_list[[i]])
        path2 <- path2 + geom_line(data = synthetic_df,
                                   aes(x = X.time, y = Ameans), color = "blue")
}
```


```R
ggsave(filename = "diagrams/xpath7.svg", plot = path2)
```

    Saving 10 x 5 in image


![title](diagrams/xpath7.svg)

And if we take the mean of the tails then we get pretty decent
estimates of the parameter.


```R
x <- list()
```


```R
for (i in c(1:3)) {
        x[[i]] <- tail(exp(A_list[[i]]), n = 50)
}
```


```R
for (i in 1:3) print(mean(x[[i]]))
```

    [1] 15.92264
    [1] 16.00402
    [1] 16.19887



```R
for (i in 1:3) print(sd(x[[i]]))
```

    [1] 1.170018
    [1] 0.6070526
    [1] 0.763862


This notebook can be downloaded from
