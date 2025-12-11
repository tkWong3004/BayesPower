####
# This script is used to evaluate the absolute error of numeric integration
# across different statistical tests. Note that two-proportion tests are not included,
# as an exact method is applied for those cases.
# Data  : 2025/11/28
# Author: TK Wong
#### loading the required functions
source("functions.r")
# These functions are the modified versions of the back-end codes
# for the calculation of power. Instead of showing the power,
# absolute errors are shown.
#### t-test ####
# specification of analysis prior and sample size
df=49
location=0
scale=1
dff=1
t.cri=seq(0.1,5,.01)
model=c("t-distribution","Normal","NLP")
hypothesis ="!="
n1=50
r=1

abs_error_t1 = array(NA, dim=c(3,length(t.cri)))
abs_error_t2 = array(NA, dim=c(3,length(t.cri)))

for( ii in 1:length(model)){
for( i in 1:length(t.cri)){
  t=c(-t.cri[i],t.cri[i]) # specifying which critical value is used for power calculation
  abs_error_t1[ii,i]=t1_TPE_abs_error(t, df, model[ii], location, scale, dff)
  abs_error_t2[ii,i]=t2_TPE_abs_error(t,n1,r,model[ii] ,location ,scale,dff , hypothesis )
}}
max(abs_error_t1,abs_error_t2)
#### correlation ####
#specification of analysis prior and sample size
r.cri=seq(0.01,.5,.01)
n=50
k=1
alpha=beta=2
h0=location=0
scale=.01
dff=1
model = c("d_beta","beta","NLP")
r_TPE_abs_error(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
abs_error_r = array(NA, dim=c(3,length(r.cri)))

for( ii in 1:length(model)){
  for( i in 1:length(r.cri)){
    r=c(-r.cri[i],r.cri[i]) # specifying which critical value is used for power calculation
    abs_error_r[ii,i]=r_TPE_abs_error(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model[ii])
  }}
max(abs_error_r)
#### f-test ####
# specification of analysis prior, sample size and the number
# of predictors in the reduced model (p) and full model (k).

f.cri = seq(0.1,5,.01)
k=10
p=2
q= k-p
n=50
m= n-p
model = c("effectsize", "Moment")
dff=3
rscale=.1
f_m=sqrt(.1)
abs_error_f = array(NA, dim=c(2,length(f.cri)))
for( ii in 1:length(model)){
  for( i in 1:length(f.cri)){
    f=f.cri[i] # specifying which critical value is used for power calculation
    abs_error_f[ii,i]=F_TPE_abs_error(f,q,m,dff,rscale,f_m,model[ii])

  }}
max(abs_error_f)
#### one-proportion ####
# specification of analysis prior and sample size
n=50
h0=location=.5
alpha=beta=1
scale=.1
hypothesis = "!="
model = c("beta","Moment")
abs_error_bin = array(NA, dim=c(2,25))

lower= 24:0
upper = 26:50
for( ii in 1:length(model)){
  for( i in 1:(25)){
    x=c(lower[i],upper[i]) # specifying which critical value is used for power calculation
    abs_error_bin[ii,i]=bin_TPE_abs_error(x,n,h0,alpha,beta,location,scale,model[ii],hypothesis)


  }}
max(abs_error_bin)
#### general ####
# the maximum absolute error across all tests is obtained by:
max(abs_error_t1,abs_error_t2,abs_error_r,abs_error_f,abs_error_t1,abs_error_bin)
0.0001215392
# the results when excluding the one for f-test with effect size analysis prior
max(abs_error_t1,abs_error_t2,abs_error_r,abs_error_f[2,],abs_error_t1,abs_error_bin)
9.992642e-05
9.992642 * 10^(-5)
