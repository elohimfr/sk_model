# sk_model
Phase space of the SK model.

This program maps the SK model phase space using the following parameters that have to be set by the user before compiling:
 * total number of spins (N);
 * thermal average dimension (tdim);
 * number of interaction configurations (conf_num);
 * number of equilibration sweeps (thermal);
 * mean interaction array (mu);
 * array of the standard deviation of the interactions (sd). 

The mean interaction array is defined by the minimum value (mu_min), maximum value (mu_max) and the step (mu_step). The size of the mean interaction array
is given by 1+(mu_max-mu_min)/mu_step.
The array of the standard deviation of the interactions is define by the minimum value (sd_min), the maximum value (sd_max) and the step (sd_step).
The size of the array of the standard deviation of the interactions is given by 1+(max_sd-min_sd)/sd_step.

The program computes: 
 * spin glass susceptibility (Xsg);
 * uniform susceptibility (Xuni);
 * spin glass order parameter (q);
 * magnetization (m);
 * specific heat (c).

The outputs of the program are .txt files. Each file is the result of the simulation of one pair of (mu, sd), and it contains one line with seven values in this 
order: mu, sd, Xsg, Xuni, q, m, c.
The program has a built-in naive parallelization for the mu and sd arrays, i.e., after compiling the code, if the user runs the same program on 10 different instances,
for example, the total time will decrease 10 times. So the user might be interested in running an array of jobs, each job running the same program for a given set of parameters.
 
