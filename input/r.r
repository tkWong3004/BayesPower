#Mode of the function
mode_bf = 1 
# 1: sample size determination
# 0: design analysis for a fixed N
N = 509 

hypothesis = "!="  # the direction of H1"!="  , ">" or "<"
e = c(-.2,.2) # bound for equivalence tests
# Analysis prior 
model = "beta"  # prior "d_beta" "NLP"
k = 1 
alpha = 1
beta = 1
h0=location = 0
scale=1
dff=1
de_an_prior = 1  
# 1: analysis and design priors are the same.
# 0: they are different.

#Design prior
model_d = "beta"    # "d_beta" "NLP" "Point"
location_d = 0
k_d = 1
alpha_d = 1
beta_d = 1
scale_d = 1
dff_d = 1

# power analysis 
D = 3    # desired value of BF
target = .8 # targeted power 
FP = .05 # targeted alpha

r_table(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
         location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,FP )

re_table(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
         location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,FP ,e)

