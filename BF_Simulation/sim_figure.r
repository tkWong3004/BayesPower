####
# This scripts are used to conduct simulation to verify
# our numeric methods for power calculation, and are used to generate the plots
# in the README file of the same folder.
# Data  : 2025/11/28
# Author: TK Wong
####
# the functions to conduct the simulations are in the file "sim_functions.r"
source("sim_functions.r")
#### t-test ####

# input
iter =10000 # number of iteration

D = 3  # bound of evidence
model = "t-distribution"   # "Normal" or "t-distribution" or "NLP"
location = -.23
scale = .2
dff = 1
hypothesis = "!="  # "!="   ">"   "<"
e=   c(-0.36, 0.36)
# e for equivalence testing
# length have to 2 when hypothesis = "!="
# one should be >0 another one should be <0
# e should <0 when hypothesis = "<"
# e should >0 when hypothesis = ">"
# set e to be NULL if testing a point null of 0

de_an_prior= 0    # whether design and analysis priors are the same
model_d="Normal" # "Normal" or "t-distribution" or "NLP" or "Point"
location_d = -.23
scale_d = .2
dff_d = 1

r =1         # ratio of sample size
t1_T2 =F    # TRUE: one sample / paired t-test FALSE : two sample t-test

png(filename="example_1.png")
t_ver(iter,D,  model, location, scale, dff, hypothesis,
      model_d, location_d, scale_d, dff_d, de_an_prior, r,t1_T2,e)
dev.off()
# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


#### correlation ####
# input
iter = 10000 # number of iteration

D =3 # bound of evidence
model = "beta"   # "d_beta" or "beta" or "NLP"
k =1
alpha=1
beta=1
h0=0
scale=1
hypothesis = ">"  # "!="   ">"   "<"
e= NULL
# Let it NULL when testing a point -null
# e for equivalence testing
# length have to 2 when hypothesis = "!="
# one should be >0 another one should be <0
# e should <0 when hypothesis = "<"
# e should >0 when hypothesis = ">"
de_an_prior=1     # whether design and analysis priors are the same
model_d = "Point" ##  # "d_beta" or "beta" or "NLP" "Point"
location_d=.3
k_d=1
alpha_d=1
beta_d=1
scale_d=1


png(filename="example_2.png")

r_ver(iter, D,model,k, alpha, beta,h0,scale, hypothesis ,model_d,
      location_d,k_d, alpha_d, beta_d,scale_d,de_an_prior,e)
dev.off()

# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


####  f-test ####

### !!! computationally demanding simulation !!! ###
# input
iter = 10000 # number of iteration

D=3
p=3       # number of predictor from the reduced model
k=4       # number of predictor from the full model

# the analysis prior
model="effectsize"    # "Moment" "effectsize"
dff=3
rscale=.18
f_m=.1
e = NULL


#the design prior
de_an_prior=  0   # design and analysis prior being the same :1 ,not the same :2
model_d="Point"      # "Moment" "effectsize" "Point"
dff_d=3                   # must be > 3 if "Moment"
rscale_d=1
f_m_d= .10


png(filename="example_3.png")
f_ver(iter,D,p,k,dff,rscale,f_m,model,
      dff_d,rscale_d,f_m_d,model_d,de_an_prior,e)
dev.off()
# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


#### one proportion ####
# input
iter = 10000 # number of iteration

D =3 # bound of evidence
model = "Moment"   #  "beta" or "Moment"

alpha=1
beta=1
h0=location=0.5
scale=1
hypothesis = ">"  # "!="   ">"   "<"
e= NULL
# Let it NULL when testing a point -null
# e for equivalence testing
# length have to be 2 when hypothesis = "!="
# one should be >0 another one should be <0
# e should <0 when hypothesis = "<"
# e should >0 when hypothesis = ">"
de_an_prior=0     # whether design and analysis priors are the same
model_d = "Point" ##  # "beta" or "Moment" or "Point"
location_d=.55
alpha_d=1
beta_d=1
scale_d=1

png(filename="example_4.png")
bin_ver(iter,h0, D,model,  alpha, beta,location,scale, hypothesis ,model_d,
        location_d,k_d, alpha_d, beta_d,scale_d,de_an_prior,e)
dev.off()

# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


#### two proportion ####
# input
iter = 10000 # number of iteration
# beta prior under h0
a0=b0=1
### beta prior under h1
#p1
a1=156
b1=339
#p2
a2=151
b2=339

# if analysis prior and design prior are the same under h1
model1=model2="same"
# else specify the things below


###design prior under h1

# p1
model1="beta"    #"Point" or "beta"
da1=156
db1=339
dp1=.5
#p2
model2="beta"    #"Point" or "beta"
da2=151
db2=339
dp2=.5

png(filename="example_5.png")
p2_ver(iter,D, a0, b0, a1, b1, a2, b2, model1,da1,db1,dp1,model2,da2,db2,dp2)
dev.off()

# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


