# How Underconfidence Is Maintained in Anxiety and Depression

 
  #################### #################### #################### 
####################

All code and data for reproducing analyses and figures for:

Sucharit Katyal, Quentin J. Huys, Raymond J. Dolan, and Stephen M. Fleming. 2023. “How Underconfidence Is Maintained in Anxiety and Depression.” PsyArXiv. May 24. 
https://psyarxiv.com/qcg92/


 
  #################### #################### #################### 
####################

Model-free analyses along with the generation of figures for the paper was done in R v4.3.0 in RStudio
(see: _code/r_). For reproducing the analyses and figures, step through the following code:

    analysis and figures for Exp1:                    code/r/mbsAnalyses_Exp1.R
    
    analysis and figures for Exp2:                    code/r/mbsAnalyses_Exp2.R
    
    figures for the model simulations and parameters: code/r/mbsAnalyses_model.R

NOTE: Before performing analyses for Exp 2, you need to unzip the two zipped .csv files in _data/exp2_

  #################### ####################

Model was fitted using Matlab (v 2022b) + Jags toolbox (v 3.4.1) (see _code/matlab_)

To run the model / reproduce the model analyses, perform the following steps:

1) Ensure that JAGS (an MCMC language similar to BUGS) is installed on your machine. See here for further details: http://mcmc-jags.sourceforge.net/. Note that there are compatibility issues between matjags (MATLAB interface for JAGS) and newer version of JAGS (e.g., 4.X). You will need to install JAGS 3.4.0 rather than the latest version (https://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/).

2) The matlab code to fit the model is available in the script:  _code/matlab/mbsMainScript.m_ 
The code contains separate cells for fitting Exp 1 and Exp 2. For each Exp, you would need to re-evaluate the cell for each model variant by changing the _modelName_ variable to 'model_N' (where N is between 0-10 for Exp 1 and 0-11 for Exp 2). The model numbers corresponding to the model variants are the same as in the manuscript and are also specfied in the matlab script. Currently, the script specifies to run the model for 2000 iterations that takes a long time to run (>1 hour), but the model should give the same results with fewer iterations (for example, you can use 50 or 100 iterations to speeed up processing, which should run in 5-10 minutes depending on your machine). After each successful run of the model, there will be a plot of the posterior distributions of the relevant model parameters. The fitted model will be returned in a field named _fitZ_ or _fitRegZ_ inside the returned data structure.

3) After running all the model variants, you can save the model using the following matlab cell (note that this will overwrite the previous model fits). You can then reproduce the model figures from the paper by using the same R script as above: mbsAnalyses_model.R
   

For running model simulations (and recovery), see script: _code/matlab/mbsSimulateExperiment.m_

  #################### #################### #################### 
####################

To run the models on your own data with _self-performance estimates_ + _local confidence_ + _accuracy_ [+ _feedback_ (optional)]), use the following function (function help contains all necessary inputs for the model):

    code/matlab/model/fit_globalSPE.m

If you need help with running the model on your data or if you need help customising the model to your dataset, contact: sucharit.katyal [@] gmail [dot] com

  #################### #################### #################### 
####################

Pre-registration doc for Exp 2 can be found at 
https://osf.io/2pq6y/


  #################### #################### #################### 
####################

_stimuli_ directory contains javascript stimuli used in Exp 1 and Exp 2

  #################### #################### #################### 
  
  License

This code is being released with a permissive open-source license. Please feel free to use or adapt the code as long as you follow the terms of the license enumerated below. If you use the model in a publication, we ask that you cite the following paper:

Katyal, S., Huys, Q. J., Dolan, R. J., & Fleming, S. M. (2023, May 24). How underconfidence is maintained in anxiety and depression. https://doi.org/10.31234/osf.io/qcg92

Copyright (c) 2023, Sucharit Katyal

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

