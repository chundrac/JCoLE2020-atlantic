require(loo)

load('all_rates_no_pooling.Rdata')

loo.nh <- c()

for (i in 1:10) {
  loo_i <- loo(fit.list[[i]])
  loo.nh <- c(loo.nh,loo_i$estimates[2])
}

loo.nh.mean <- loo(fit.full)$estimates[2]

load('all_rates_pooling.Rdata')

loo.cp <- c()

for (i in 1:10) {
  loo_i <- loo(fit.list[[i]])
  loo.cp <- c(loo.cp,loo_i$estimates[2])
}

loo.cp.mean <- loo(fit.full)$estimates[2]

load('all_rates_hierarchical.Rdata')

loo.h <- c()

for (i in 1:10) {
  loo_i <- loo(fit.list[[i]])
  loo.h <- c(loo.h,loo_i$estimates[2])
}

loo.h.mean <- loo(fit.full)$estimates[2]

#(loo.nh-loo.h)/(loo.nh-loo.cp)

(loo.nh.mean-loo.h.mean)/(loo.nh.mean-loo.cp.mean)