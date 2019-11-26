function LL = truVM_lik(par,x,a,b)

LL = -sum(log(truncatedVMpdf(x,0,par,a,b)));