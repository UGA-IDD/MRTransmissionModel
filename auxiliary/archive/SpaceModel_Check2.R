## Conclusion
## waifw.tmp%*%state.tmp is the same thing as what c code is doing

phi <- as.double(tran@waifw)*((state[tran@i.inds[this.loc],]^tran@exponent)/denom)
#phi <- (tran@waifw%*%(state[tran@i.inds[this.loc],]^tran@exponent))/denom

waifw.tmp <- tran@waifw[1:10, 1:10]
state.tmp <- ((state[tran@i.inds[this.loc],]^tran@exponent)/denom)[1:10]

waifw.tmp%*%state.tmp
as.double(waifw.tmp)*state.tmp


waifw=tran@waifw
state
s_indsR = tran@s.inds[1:300]
i_indsR = tran@i.inds[1:300]
exponentR = tran@exponent
denomR = denom
tran_matrix = tran@age.surv.matrix
n_age_class = length(i_indsR)

# Do not need phi outside of the function
phi <- numeric(n_age_class)

for (i in 1:n_age_class) {
  phi[i] <- 0.0

  for (j in 1:n_age_class) {
    phi[i] <- phi[i] + waifw[i + (j - 1) * n_age_class] *
      (state[i_indsR[j]] ^ exponentR) / denomR
  }
}

phi0 <- phi

for (i in 1:n_age_class) {
  phi[i] <- 1 - exp(-phi[i])
}


phin <- as.vector(phi[,1])
phi0/phin

