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
