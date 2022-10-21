coordinate_correction <- function(x_tilde, mu, mu_tilde, Sigma_inv, phi, S, step_size = 0.1, max_step = NULL){
  if(is.null(max_step)){
    max_step <- 1/step_size
  }
  p <- length(x_tilde)
  x_history <- x_tilde
  phi_history <- phi
  mu_tilde_history <- mu_tilde
  d <- numeric(length = p)
  if(p != length(S)){
    iter <- 0
    while((max(c(phi[S],0),na.rm = TRUE)*(1+length(S)/p) > max(c(phi[-S],0), na.rm = TRUE))&(iter <= max_step)){
      iter <- iter + 1
      correction_step <- (x_tilde[S] - mu_tilde[S])*step_size
      d[S] <- d[S] + correction_step
      x_tilde[S] <- x_tilde[S] - correction_step
      phi <- shapley_base(x_tilde, mu_tilde, Sigma_inv, inverted = TRUE)
      x_history <- rbind(x_history, x_tilde)
      phi_history <- rbind(phi_history, phi)
      mu_tilde_history <- rbind(mu_tilde_history, mu_tilde)
      # if(p == length(S)){
      #   x_tilde <- mu
      #   print("mean")
      #   break
      # }
    }
  } else{
    x_tilde <- mu
    x_history <- x_tilde
  }
  return(list("x" = x_tilde, "d" = d, "phi" = phi,
              "x_history" = x_history, "phi_history" = phi_history, "mu_tilde_history" = mu_tilde_history))
}
