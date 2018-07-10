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
