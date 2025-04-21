`.sourceCpp_3_DLLInfo` <- dyn.load('/home/jorge/MyGitHub/BayesPower_alpha/tmp_cache/sourceCpp-x86_64-pc-linux-gnu-1.0.14/sourcecpp_d14d530e46c2/sourceCpp_4.so')

pnct <- Rcpp:::sourceCppFunction(function(x, df, ncp = as.numeric( c(0.0)), lower = TRUE) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_pnct')

rm(`.sourceCpp_3_DLLInfo`)
