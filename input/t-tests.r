#Mode of the function
mode_bf = 1 
# 1: sample size determination
# 0: design analysis for a fixed N
N = 509 # sample size in case of one sample t-test
N1=N2=50 # sample size per group for independent t-tests
r = 1 # ratio of sample size N2/N1

hypothesis = "!="  # the direction of H1"!="  , ">" or "<"
e = c(-.2,.2) # bound for equivalence tests
# Analysis prior 
model = "Cauchy"
location = 0
scale = 1
dff=1
de_an_prior = 1  
# 1: analysis and design priors are the same.
# 0: they are different.

#Design prior
model_d = "Point"
location_d =0.2
scale_d = 0
dff_d=1

# power analysis 
D = 3    # desired value of BF
target = .8 # targeted power 
alpha = .05 # targeted alpha


t1_Table(D,target,model,location,scale,dff, hypothesis,
                     model_d,location_d,scale_d,dff_d, de_an_prior,N, mode_bf ,alpha )
t1_prior_plot(D ,target,model ,location ,scale,dff , hypothesis,model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior)
bf10_t1 (D ,df, target,model ,location ,scale,dff  , hypothesis )
Power_t1(D,model,location,scale,dff, hypothesis,
                   model_d,location_d,scale_d,dff_d, de_an_prior,N)


# supporting functions
t1_BF10(t,df,model ,location,scale,dff , hypothesis )
t1_BF10_bound (D, df,model ,location ,scale,dff , hypothesis)
t1_BF01_bound (D, df,model ,location ,scale,dff , hypothesis)
t1_TPE(t , df,model,location ,scale,dff , hypothesis)
t1_FNE(t , df,model,location ,scale,dff , hypothesis)
t1_TNE(t , df,model,location ,scale,dff , hypothesis)
t1_FPE(t , df,model,location ,scale,dff , hypothesis)
t1_N_finder(D,target,model,location,scale,dff, hypothesis ,
                      model_d,location_d,scale_d,dff_d,de_an_prior,alpha)
