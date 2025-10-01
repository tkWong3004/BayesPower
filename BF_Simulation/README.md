# Simulation Scripts

This folder contains R scripts for running simulations used to verify our method.

## Method

The simulation process consists of three steps:

1. **Simulating cumulative probabilities**  
   A vector of cumulative probabilities of length *i* is generated. These probabilities are then transformed into effect sizes based on the design prior.

2. **Generating data**  
   For each simulated effect size, data are generated. For example, given a simulated delta, a corresponding *t*-value is drawn using the `rt` function.

3. **Computing the probabilities**  
   For a specified analysis prior, the probabilities that BF > *k* under the null hypothesis and alternative hypothesis are computed based on the simulations.

Simulations are conducted for sample sizes ranging from 100 to 1000 with the increments of 100.

## Results

We ran simulations under various conditions (different design and analysis priors).  
The results show that the simulated probabilities of compelling and misleading evidence align closely with those obtained using our method.

Users can replicate or extend these simulations using the scripts in this folder.  
Each script loads the required packages and defines custom functions for simulation. Toward the end of the script, users may adjust the input parameters as needed and run the simulation by executing the final lines of code.
