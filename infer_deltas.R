require(rstan)

load('all_rates_hierarchical.Rdata')

model_code <- "data {
  int<lower=1> N;
  int<lower=1> D;
  real p_z[D,N];
}
transformed data {
  real e[D,N];
  for (d in 1:D) {
  for (i in 1:N) {
    e[d,i] = 0;
    if (p_z[d,i]<=.5) {
      e[d,i] += p_z[d,i];
    }
    else {
      e[d,i] += 1-p_z[d,i];
    }
    if (1-p_z[d,i]<=.5) {
      e[d,i] += 1-p_z[d,i];
    }
    else {
      e[d,i] += p_z[d,i];
    }
    if (e[d,i] == 0) {
      e[d,i] += 1e-10;
    }
    if (e[d,i] == 1) {
      e[d,i] -= 1e-10;
    }
  }
  }
}
parameters {
  real<lower=0> lambda;
  real<lower=0> a[D];
  real<lower=0> b[D];
}
model {
  lambda ~ gamma(1,1);
  a ~ exponential(lambda);
  b ~ exponential(lambda);
  for (d in 1:D) {
    for (i in 1:N) {
      e[d,i] ~ beta(a[d],b[d]);
    }
  }
}
generated quantities {
  real delta_diff[D,D];
  for (i in 1:D) {
    for (j in 1:D) {
      delta_diff[i,j] = (b[i]/a[i])-(b[j]/a[j]);
    }
  }
}"

fit.list.delta <- list()
for (i in 1:length(fit.list)) {
  p_z = colMeans(extract(fit.list[[i]])$z)[,37:71]
  data.list <- list(D=nrow(p_z),N=ncol(p_z),p_z=p_z)
  fit <- stan(model_code=model_code,data=data.list,chains=3)
  fit.list.delta[[i]] <- fit
  fit.full <- sflist2stanfit(fit.list.delta)
  save.image('all_deltas.Rdata')
}