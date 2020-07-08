# skmodel.c
Phase space of the SK model.

Please cite the following paper when you use this code.

[Ezaki T, Fonseca dos Reis E, Watanabe T, Sakaki M, Masuda N. Closer to critical resting-state neural dynamics in individuals with higher fluid intelligence. Commun Biol 3:1 (2020).](https://www.nature.com/articles/s42003-020-0774-y)

### How to use

1) Download skmodel.c

2) Download [mt19937ar.c](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html) in the same folder.

3) Set the parameters

The following parameters have to be set by the user in the code before compiling.

Parameters:
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

4) Compile and run

The program uses GSL libraries. Then, the user must link the libraries (-lgsl -lgslcblas -lm). For intance, if using gcc, type in the terminal:
```
gcc -Wall skmodel.c -o skmodel_program -lgsl -lgslcblas -lm
```

Then, run the program:
```
./skmodel_program
```

5) Output

The program computes: 
 * spin glass susceptibility (Xsg);
 * uniform susceptibility (Xuni);
 * spin glass order parameter (q);
 * magnetization (m);
 * specific heat (c).

The outputs of the program are .txt files. Each file is the result of the simulation of one pair of (mu, sd), and it contains one line with seven values in this 
order: mu, sd, Xsg, Xuni, q, m, c.

6) Parallelizing

The program has a built-in naive parallelization for the mu and sd arrays, i.e., after compiling the code, if the user runs the same program on 10 different instances,
for example, the total time will decrease 10 times. So the user might be interested in running an array of jobs, each job running the same program for a given set of parameters.
 


