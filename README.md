# SEIR_SINDY_MPC

With this tutorial we want to provide the tools to apply data driven system identification methods to a model predictive control (MPC) problem.
The method was first introduced by Kaiser et al.:
  - github: https://github.com/eurika-kaiser/SINDY-MPC
  - paper: https://arxiv.org/abs/1711.05501

We apply the recent sparse identification of nonlinear dynamics (SINDy) to MPC and provide a tutorial on the control of the spread of an infectious disease using an identified SEIR model. 

In the tutorial article LINK, we also review the SINDy and MPC methods, and demonstrate the tutorial with the main code snippets to run the SINDy-MPC algorithm.

To run the code: Code/main.m
  - select the system identification method: DMD or SINDy
  - select if the MPC is run with or without constraint
  
We encourage the user to play around with the choose, choose different forcing functions to train the model, change the model parameters, add compartments to the SEIR model, run different prediction and horizon lengths for the MPC etc.

