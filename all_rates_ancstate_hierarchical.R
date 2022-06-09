require(phytools)
require(rstan)

model_code = "functions {
  matrix evprob(real z, real alpha, real beta) {
    matrix[2,2] P;
	P[1,1] = (beta/(alpha+beta)) + (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
	P[1,2] = (alpha/(alpha+beta)) - (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
	P[2,1] = (beta/(alpha+beta)) - (beta/(alpha+beta)*exp(-(alpha+beta)*z));
	P[2,2] = (alpha/(alpha+beta)) + (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    return P;
  }
  real pruning(int N, int B, int[] child, int[] parent, real[] brlen, matrix tiplik, real alpha, real beta) {
    matrix[N,2] lambda;                  //likelihoods at tips+nodes
    vector[2] pi;                         //stationary probability
    lambda = log(tiplik);
    for (b in 1:B) {
      matrix[2,2] P = evprob(brlen[b], alpha, beta); //via matrix exponentiation
      for (d in 1:2) {
        lambda[parent[b],d] += log(dot_product(P[d],exp(lambda[child[b]])));
      }
    }
  	pi[1] = log(beta) - log(alpha+beta) + lambda[parent[B],1];
  	pi[2] = log(alpha) - log(alpha+beta) + lambda[parent[B],2];
    return(log_sum_exp(pi));
  }
  vector anc_state_rng(int N, int T, int B, int[] child, int[] parent, real[] brlen, matrix tiplik, real alpha, real beta) {
    matrix[N,2] lambda;                  //likelihoods at tips+nodes
    vector[2] pi;                         //stationary probability
    int z[N];
    lambda = log(tiplik);
    for (b in 1:B) {
      matrix[2,2] P = evprob(brlen[b], alpha, beta); //via matrix exponentiation
      for (d in 1:2) {
        lambda[parent[b],d] += log(dot_product(P[d],exp(lambda[child[b]])));
      }
    }
    pi[1] = log(beta) - log(alpha+beta) + lambda[parent[B],1];
  	pi[2] = log(alpha) - log(alpha+beta) + lambda[parent[B],2];
  	pi = exp(pi-log_sum_exp(pi));
  	z[parent[B]] = bernoulli_rng(pi[2]);
  	for (b in 0:(B-1)) {
  	  vector[2] d;
  	  matrix[2,2] P = evprob(brlen[B-b], alpha, beta);
  	  d = to_vector(lambda[child[B-b],] + log(P[z[parent[B-b]]+1,]));
  	  d = exp(d-log_sum_exp(d));
  	  z[child[B-b]] = bernoulli_rng(d[2]);
  	}
  	return(to_vector(z));
  }
}
data {
  int<lower=1> T;
  int<lower=1> N;
  int<lower=1> B; //number of branches
  int<lower=1> D;                       //number of features
  int<lower=1> child[B];                //child of each branch
  int<lower=1> parent[B];               //parent of each branch
  real<lower=0> brlen[B];               //length of each branch
  matrix[N,D*2] tiplik;     //likelihoods for data at tips in tree
}
parameters {
  real<lower=0> lambda;
  real<lower=0> mu;
  real<lower=0> alpha[D];
  real<lower=0> beta[D];
}
model {
  lambda ~ exponential(1);
  mu ~ exponential(1);
  alpha ~ gamma(lambda,lambda);
  beta ~ gamma(mu,mu);
  for (d in 1:D) {
    target += pruning(N,B,child,parent,brlen,tiplik[,((2*d)-1):(2*d)],alpha[d],beta[d]);
  }
}
generated quantities {
  real pi[D];
  real r[D];
  real contrasts_pi[D,D];
  real contrasts_r[D,D];
  vector[N] z[D];
  for (d in 1:D) {
    pi[d] = alpha[d]/(alpha[d]+beta[d]);
    r[d] = alpha[d]+beta[d];
  }
  for (i in 1:D) {
    for (j in 1:D) {
      contrasts_pi[i,j] = pi[i]-pi[j];
      contrasts_r[i,j] = r[i]-r[j];
    }
  }
  for (d in 1:D) {
    z[d] = anc_state_rng(N,T,B,child,parent,brlen,tiplik[,((2*d)-1):(2*d)],alpha[d],beta[d]);
  }
}"

trees <- read.nexus('trees.nex')

atlantic.data <- read.csv('../data_raw/Atlantic.csv')
rownames(atlantic.data) <- atlantic.data$Glottocode
atlantic.data <- atlantic.data[,c('ART','DEM','ADJ','CARD','PRO','VERB','PREF')]

set.seed(1234)
inds <- sample(1:length(trees),100)

fit.list <- list()

for (t in inds) {
  tree <- trees[[t]]
  tree <- reorder.phylo(tree,'pruningwise')
  atlantic.data.curr <- atlantic.data[tree$tip.label,]
  bin.states <- NULL
  for (d in 1:ncol(atlantic.data)) {
    bin.states.d <- to.matrix(as.character(atlantic.data.curr[,d]),seq=c('0','1'))
    bin.states.d[rowSums(bin.states.d)==0,] <- c(1,1)
    bin.states <- cbind(bin.states,bin.states.d)
  }
  bin.states <- rbind(bin.states,matrix(1,nrow=tree$Nnode,ncol=ncol(bin.states)))
  parent <- tree$edge[,1]
  child <- tree$edge[,2]
  b.lens <- tree$edge.length/1000
  N <- length(unique(c(parent,child)))
  T <- length(child[which(!child %in% parent)])
  tip.lik <- bin.states
  data.list <- list(N=N,
                    T=T,
                    B=length(parent),
                    brlen=b.lens,
                    child=child,
                    parent=parent,
                    tiplik=tip.lik,
                    D=ncol(tip.lik)/2)
  fit <- stan(model_code=model_code,data=data.list,chains=3,thin=10)
  fit.list <- append(fit.list,fit)
  fit.full <- sflist2stanfit(fit.list)
  save.image('all_rates_hierarchical.Rdata')
}