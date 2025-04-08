#Mode of the function
mode_bf = 1 
# 1: sample size determination
# 0: design analysis for a fixed N
n = 509 

e = c(-.2,.2) # bound for equivalence tests
# Analysis prior 
model = "beta"  # "Moment"
alpha = 1
beta = 1
location = .5
scale = 1

de_an_prior = 1  
# 1: analysis and design priors are the same.
# 0: they are different.

#Design prior
model_d = "beta"  # "Moment"
alpha_d = 1
beta_d = 1
location_d = .5
scale_d = 1


# power analysis 
D = 3    # desired value of BF
target = .8 # targeted power 
FP = .05 # targeted alpha

bin_table(D,target,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,N, mode_bf,FP)
bin_e_table(D,target,alpha,beta,location,scale,model,hypothesis,
          alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,N, mode_bf,FP,e)
