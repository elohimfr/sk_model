# skmodel.c

This program simulates the Sherrington-Kirkpatrick (SK) model of spin glasses without external field ($`h=0`$) and $kT=1$. In this case, the Hamiltonian of the SK model is given by,

<a href="https://www.codecogs.com/eqnedit.php?latex=H&space;=&space;-\sum_{i<j}&space;J_{ij}&space;S_i&space;S_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H&space;=&space;-\sum_{i<j}&space;J_{ij}&space;S_i&space;S_j" title="H = -\sum_{i<j} J_{ij} S_i S_j" /></a>

The spin variables are assumed to be of the Ising type (<a href="https://www.codecogs.com/eqnedit.php?latex=S_i&space;=&space;\pm&space;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S_i&space;=&space;\pm&space;1" title="S_i = \pm 1" /></a>). The interaction $J_{ij}$ is a quenched variable independently distributed according to a Gaussian distribution given by
$$ P(J_{ij}) = \frac{1}{\sqrt{2 \pi \sigma^2}} \exp{\left\{ - \frac{(J_{ij} - \mu)^2}{2 \sigma^2} \right\}},$$
where $\mu$ and $\sigma$ are the mean and the standard deviation of $J_{ij}$, respectively.

The program uses the single-spin flip Metropolis Monte Carlo algorithm. So at each *time step* (or Monte Carlo step) one spin is chosen uniformly at random and the spin flip is accepted or rejected according to the acceptance rates of the Metropolis algorithm. For a system of $N$ spins, a *sweep* of the system consists of $N$ time steps. Before computing the variables of the model, a *transient period* of a given number of sweeps is performed and discarded.

To compute the variables of the SK model, because it is a spin glass model, two averages have to be taken. The *time average* (a.k.a. thermal averge) and the *configurational average*. The time average is the average over a given number of sweeps for a given fixed configuration of interactions $\mathbf{J} \equiv \{J_{ij}\}$. Therefore, the time average samples the variable every sweep, but not every time step. The configurational average is the averge taken over a given number of different $J_{ij}$ configurations. For example, if we allow $L$ sweeps for the time average and $M$ configurations for the configurational average, a variable $A$ is calculated as
$$ [\langle A \rangle] = \frac{1}{M} \sum_{\alpha = 1}^M \frac{1}{L} \sum_{i=1}^L A_i^{(\alpha)},$$
where $\langle \ldots \rangle$ denotes the time average and $[\ldots]$ denotes the configurational average.

Please cite the following paper when you use this code.

[Ezaki T, Fonseca dos Reis E, Watanabe T, Sakaki M, Masuda N. Closer to critical resting-state neural dynamics in individuals with higher fluid intelligence. Commun Biol 3:1 (2020).](https://www.nature.com/articles/s42003-020-0774-y)

### How to use

1) Download skmodel.c

2) Download [mt19937ar.c](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html) in the same folder

3) Set the parameters

The following parameters have to be set by the user in the code before compiling.

Parameters:
 * Total number of spins (N);
 * Number of time steps to calculate the averages given a fixed configuration of $J_{ij}$ (tdim);
 * Number of $J_{ij}$ configurations for the configurational average (conf_num);
 * Number of sweeps in the transient period (thermal);
 * Array with the mean values of $J_{ij}$ (mu). The user has to set the minimum value (mu_min), maximum value (mu_max) and the step (mu_step). The size of the array with the mean values of $J_{ij}$ is given by 1+(mu_max-mu_min)/mu_step.
 * Array with the standard deviation values of $J_{ij}$ (sd). The user has to set the minimum value (sd_min), the maximum value (sd_max) and the step (sd_step). The size of the array with the standard deviation values of $J_{ij}$ is given by 1+(max_sd-min_sd)/sd_step.

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

The program outputs the following measures: 
 * spin glass susceptibility (Xsg);
 * uniform susceptibility (Xuni);
 * spin glass order parameter (q);
 * magnetization (m);
 * specific heat (c).

The outputs of the program are multiple .txt files, one for each pair of (mu, sd), and each .txt file contains one line with seven values in this 
order: mu, sd, Xsg, Xuni, q, m, c. For example, if the user sets *n* values for the mu array and *m* values for the sd array, the program outputs *nm* .txt files. The files are named as: file_0_0.txt, file_0_1.txt, ..., file_1_0.txt, file_1_1.txt... The reason for generating multiple .txt files is for parallelizing the program as explained next.

6) Parallelizing

The program has a built-in naive parallelization for the mu and sd arrays, i.e., after compiling the code, if the user runs the same program on 10 different instances,
for example, the total time will decrease 10 times. So the user might be interested in running an array of jobs, each job running the same program for a given set of parameters.
 
7) Gather the outputs

* Download the jupyter notebook sk_heatmap.ipynb in the same folder where the outputs are
* Open the notebook
* Set the parameters accordingly, as indicated in the notebook
* Run the python code
* A figure for each measure is produced


