#Mode of the function
mode_bf = 1 
# 1: sample size determination
# 0: design analysis for a fixed N
n = 509 

e = .2 # bound for equivalence tests
# Analysis prior 
model = "Moment"  # "tdis"
p = 1
k = 3
dff = 10
rscale = 1
f_m = .1

de_an_prior = 1  
# 1: analysis and design priors are the same.
# 0: they are different.

#Design prior
model_d= "Moment"  # "tdis" "Point"
dff_d = 3
rscale_d = 1
f_m_d = .5


# power analysis 
D = 3    # desired value of BF
target = .8 # targeted power 
FP = .05 # targeted alpha

f_table(D,target,p,k,dff,rscale,f_m,model,
                  dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, mode_bf,FP )
fe_table(D,target,p,k,dff,rscale,f_m,model,
        dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, mode_bf,FP,e )
