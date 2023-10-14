functions { // function for icar
    real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2) {
        return -0.5 * dot_self(phi[node1] - phi[node2]);
    }
}
data {
    int<lower=0> N;
    int<lower=0> k; // # of independent variable
    int<lower=0> N_edges;
    array[N_edges] int<lower=1, upper=N> node1;
    array[N_edges] int<lower=1, upper=N> node2;
    array[N] int<lower=0> Y; // dependent variable which is the # of absentees
    matrix[N, k] X; // independent variables (expenditure & income)
    vector<lower=0>[N] E; // estimated number of expected absentees
}
transformed data {
    vector[N] log_offset = log(E); //log transform expected value (offset)
}
parameters {
    real alpha; // intercept which is the overall risk (overall risk for the entire area)
    vector[k] beta; // coefficients for independent variables
    real<lower=0> sigma; // SD
    vector[N] phi; // spatial random effect
}
model { // likelihood function and priors 
    phi ~ icar_normal(N, node1, node2); // prior for the spatial random effect which is a function of # of prefectures & adjancy matrix
    Y ~ poisson_log(log_offset + alpha + X*beta + phi*sigma); // likelihood function
    alpha ~ normal(0.0, 1.0); // prior for intercept   (weak/uninformative prior) always use normal, alpha is the overall risk
    beta ~ normal(0.0, 1.0); // prior for coefficient (weak/uninformative prior) always use normal
    sigma ~ normal(0.0, 1.0); // prior for SD          (weak/uninformative prior)
    sum(phi) ~ normal(0, 0.001*N); // phi should sum to zero
}
generated quantities {
    vector[N] eta = alpha + X*beta  + phi*sigma; // relative risk for areas 
    vector[N] mu = exp(eta); // areas-specific relative risk ratios
}
