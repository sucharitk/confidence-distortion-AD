# Formation of global confidence and its distortion in anxious-depression

  #################### #################### #################### 
####################

All code and data for reproducing analyses and figures for:

Katyal, Sucharit, Quentin J. Huys, Raymond J. Dolan, and Stephen M. Fleming. 2023. “How Underconfidence Is Maintained in Anxiety and Depression.” PsyArXiv. May 24. 
https://psyarxiv.com/qcg92/


  #################### #################### #################### 
####################

Model-free analyses along with the generation of figures for the paper was done in R v4.3.0 in RStudio
(see: code/r)

analysis and figures for Exp1:                    code/r/mbsAnalyses_Exp1.R

analysis and figures for Exp2:                    code/r/mbsAnalyses_Exp2.R

figures for the model simulations and parameters: code/r/mbsAnalyses_model.R

  ####################

Model was fitted using Matlab (v 2022b) + Jags toolbox (v 3.4.1) (see code/matlab)

Before running the model, you need to first ensure JAGS (an MCMC language similar to BUGS) is installed on your machine. See here for further details: http://mcmc-jags.sourceforge.net/

Note that there are compatibility issues between matjags and newer version of JAGS (e.g., 4.X). You will need to install JAGS 3.4.0 rather than the latest version. 

model fitting code can be called from the main script:  code/matlab/mbsMainScript.m

for model simulations (and recovery), see script: code/matlab/mbsSimulateExperiment.m

  ####################

NOTE: Before performing analyses for Exp 2, you need to unzip the two zipped .csv files in data/exp2

  #################### #################### #################### 
####################

Pre-registration doc for Exp 2 can be found at 
https://osf.io/2pq6y/


  #################### #################### #################### 
####################

stimuli directory contains javascript stimuli used in Exp 1 and Exp 2

  #################### #################### #################### 
  
  License

This code is being released with a permissive open-source license. Please feel free to use or adapt the code as long as you follow the terms of the license enumerated below. If you use the model in a publication, we ask that you cite the following paper:

Katyal, S., Huys, Q. J., Dolan, R. J., & Fleming, S. M. (2023, May 24). How underconfidence is maintained in anxiety and depression. https://doi.org/10.31234/osf.io/qcg92

Copyright (c) 2023, Sucharit Katyal

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

