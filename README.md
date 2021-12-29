<p align="center">
  <img src="https://github.com/IngeborgGjerde/nitsche-method-for-navier-stokes-with-slip/blob/master/potential-flow.png">
</p>

# Solving the Navier Stokes equations with slip using a Nitsche method

This repository contains source codes and an environment used to produce
results of the manuscript **Nitsche's method for Navier-Stokes equations with slip boundary conditions**
by Ingeborg Gjerde (Simula Research Laboratory) and L. Ridgway Scott (University of Chigaco)

## Dependencies
- FEniCS 2019.1.0 (with python3)


## Usage
The numerical results used in the article can be obtained by running the script run_convergence_test.py as follows

```
~/ns-with-slip$ python run_convergence_test.py --nrefs 4 --beta=-2.0
```

with 4 mesh refinements, friction parameter beta=-2.0 and the remaining variables initialized to a default value. 

The user can specify the following variables:
-- nrefs (int), number of mesh refinements, default = 4
-- cyl_refinement (int), whether the mesh should be refined around the cylinder, 1=Yes and 0=No, default=0
-- beta (float), friction parameter, default = -2.0 in which case one can use the potential flow solution
-- nu (float), viscosity, default=1.0
-- degree (int), polynomial degree, default=2
-- ngammas (int), gives a range of n different Nitsche stabilization parameters to try, default=1
-- normal_choice (str), gives the normal type, options are projected_normal and discrete_normal


Running e.g. with all default values, 
```
~/ns-with-slip$ python run_convergence_test.py
```
produces the text file `u_error_beta-2.0_uniform.txt` located in `/home/ns-with-slip/results` which gives a latex table containing the errors and convergence rates for the H1-error of u-u_h. Tables are produced for p-p_h in the L2 norm, u-u_h in the H1 norm and u-u_h in the norm induced by a. 
 

## Limitations
- The code assumes a cylinder radius R=1.0 and cylinder center in [0.0, 0.0]

## Citing
If you want to use this code for your own work, please cite our article: https://www.ams.org/journals/mcom/0000-000-00/S0025-5718-2021-03682-0/home.html
