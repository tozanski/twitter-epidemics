# Twitter-epidemics

+ In this repository, we share Julia code for modeling Twitter networks as Bolobas-Jensen-Riordian graphs and simulating epidemic processes in them.

## The Jupyter notebooks 

+ Plots-Higgs.ipynb is the main Julia notebook which generates all the necessary plots.
+ Running the notebook require Julia kernel for Jupyter. 
+ The most basic requirements are `IJulia` and `PyPlot` packages. 
+ Please make sure that the rest of the julia files shared in the repository are visible.

## Dataset
+ We use the Higgs Twitter dataset as a sample dataset here. The details can be found here: https://snap.stanford.edu/data/higgs-twitter.html. 
+ The formatted data (suitable to be loaded by the Julia notebook) is shared as higgs.jld2 in the "Data" folder. 

## Output 
+ As output, the code generates the estimated indegree, outdegree distributions with different BJR models and the comparison of the simulated and estimated values of different metrics such as outbreak probability, reach, prevalence and resistance.

## Contact
  + For any query, please contact Tomasz Ożański at tomaszoz[at]gmail[dot]com.
