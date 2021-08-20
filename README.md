# SEIR_SINDY_MPC

This tutorial article (LINK follows soon) provides the tools to apply data driven system identification methods to a model predictive control (MPC) problem.
The method was first introduced by Kaiser et al.:
  - github: https://github.com/eurika-kaiser/SINDY-MPC
  - paper: https://arxiv.org/abs/1711.05501

We apply the recent sparse identification of nonlinear dynamics (SINDy) to MPC and provide a tutorial on the control of the spread of an infectious disease using an identified SEIR model. 

In the tutorial article (LINK), we also review the SINDy and MPC methods, and demonstrate the tutorial with the main code snippets to run the SINDy-MPC algorithm.

To run the code: Code/main.m
  - select the system identification method: DMD or SINDy (line 16)
  - select if the MPC is run with or without constraint (line 20)
  
We encourage the user to play around with the code, choose different forcing functions to train the model, change the SEIR and SINDy model parameters, add compartments to the SEIR model, run different prediction and control horizons for the MPC etc.

