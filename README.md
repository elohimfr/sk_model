# skmodel.c

This program simulates the Sherrington-Kirkpatrick (SK) model of spin glasses without external field (<a href="https://www.codecogs.com/eqnedit.php?latex=h=0" target="_blank"><img src="https://latex.codecogs.com/svg.latex?h=0" title="h=0" /></a>) and, under <a href="https://www.codecogs.com/eqnedit.php?latex=kT=1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?kT=1" title="kT=1" /></a>. In this case, the Hamiltonian of the SK model is given by,

<a href="https://www.codecogs.com/eqnedit.php?latex=H&space;=&space;-\sum_{\substack{i,\,j&space;=&space;1,\,&space;\ldots,\,&space;N&space;\\&space;i<j}}&space;J_{ij}&space;S_i&space;S_j," target="_blank"><img src="https://latex.codecogs.com/svg.latex?H&space;=&space;-\sum_{\substack{i,\,j&space;=&space;1,\,&space;\ldots,\,&space;N&space;\\&space;i<j}}&space;J_{ij}&space;S_i&space;S_j," title="H = -\sum_{\substack{i,\,j = 1,\, \ldots,\, N \\ i<j}} J_{ij} S_i S_j," /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=S_i=1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?S_i=1" title="S_i=1" /></a> or <a href="https://www.codecogs.com/eqnedit.php?latex=-1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?-1" title="-1" /></a>.
The interaction <a href="https://www.codecogs.com/eqnedit.php?latex=J_{ij}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?J_{ij}" title="J_{ij}" /></a> is independently distributed according to a Gaussian distribution given by

<a href="https://www.codecogs.com/eqnedit.php?latex=P(J_{ij})&space;=&space;\frac{1}{\sqrt{2&space;\pi&space;\sigma^2}}&space;\exp{\left\{&space;-&space;\frac{(J_{ij}&space;-&space;\mu)^2}{2&space;\sigma^2}&space;\right\}}," target="_blank"><img src="https://latex.codecogs.com/svg.latex?P(J_{ij})&space;=&space;\frac{1}{\sqrt{2&space;\pi&space;\sigma^2}}&space;\exp{\left\{&space;-&space;\frac{(J_{ij}&space;-&space;\mu)^2}{2&space;\sigma^2}&space;\right\}}," title="P(J_{ij}) = \frac{1}{\sqrt{2 \pi \sigma^2}} \exp{\left\{ - \frac{(J_{ij} - \mu)^2}{2 \sigma^2} \right\}}," /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\mu" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mu" title="\mu" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\sigma" title="\sigma" /></a> are the mean and the standard deviation of <a href="https://www.codecogs.com/eqnedit.php?latex=J_{ij}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?J_{ij}" title="J_{ij}" /></a>, respectively.

The program uses the single-spin flip Metropolis Monte Carlo algorithm.
So, at each *time step* (i.e., Monte Carlo step) we choose one spin uniformly at random and accept or reject the proposed spin flip according to the acceptance rates of the Metropolis algorithm. For a system of *N* spins, a *sweep* of the system consists of *N* time steps. Before computing the observables of the model, a *transient period* of a given number of sweeps is performed and discarded.

To compute the observables of the SK model, we distinguish between the *time average* (a.k.a. thermal averge) and the *configurational average*.
The time average is the average over a given number of sweeps for a given fixed configuration of interactions <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{J}&space;\equiv&space;\{J_{ij}\}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{J}&space;\equiv&space;\{J_{ij}\}" title="\mathbf{J} \equiv \{J_{ij}\}" /></a>.
We sample the observables once every sweep, but not every time step, and average the samples over the given number of sweeps.
The configurational average is the averge taken over a given number of different configurations <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{J}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{J}" title="\mathbf{J}" /></a>.
For example, if we allow *L* sweeps for the time average and *M* configurations for the configurational average, an observable *A* is calculated as

<a href="https://www.codecogs.com/eqnedit.php?latex=[\langle&space;A&space;\rangle]&space;=&space;\frac{1}{M}&space;\sum_{\alpha&space;=&space;1}^M&space;\frac{1}{L}&space;\sum_{i=1}^L&space;A_i^{(\alpha)}," target="_blank"><img src="https://latex.codecogs.com/svg.latex?[\langle&space;A&space;\rangle]&space;=&space;\frac{1}{M}&space;\sum_{\alpha&space;=&space;1}^M&space;\frac{1}{L}&space;\sum_{i=1}^L&space;A_i^{(\alpha)}," title="[\langle A \rangle] = \frac{1}{M} \sum_{\alpha = 1}^M \frac{1}{L} \sum_{i=1}^L A_i^{(\alpha)}," /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\langle&space;\ldots&space;\rangle" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\langle&space;\ldots&space;\rangle" title="\langle \ldots \rangle" /></a> denotes the time average, <a href="https://www.codecogs.com/eqnedit.php?latex=[\ldots]" target="_blank"><img src="https://latex.codecogs.com/svg.latex?[\ldots]" title="[\ldots]" /></a> denotes the configurational average, and <a href="https://www.codecogs.com/eqnedit.php?latex=A_i^{(\alpha)}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?A_i^{(\alpha)}" title="A_i^{(\alpha)}" /></a> is the sample	of *A* in the *i*th sweep and <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\alpha" title="\alpha" /></a>th configuration.

Please cite the following paper when you use this code.

[Ezaki T, Fonseca dos Reis E, Watanabe T, Sakaki M, Masuda N. Closer to critical resting-state neural dynamics in individuals with higher fluid intelligence. Communications Biology 3:1 (2020).](https://www.nature.com/articles/s42003-020-0774-y)

### How to use

1) Download skmodel.c

2) Download [mt19937ar.c](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html) in the same folder

3) Set the parameters

The following parameters have to be set by the user in the code before compiling.

Parameters:
 * Number of spins (N);
 * *L*, i.e., number of sweeps for the time average given a fixed configuration of <a href="https://www.codecogs.com/eqnedit.php?latex=J_{ij}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?J_{ij}" title="J_{ij}" /></a> (tdim);
 * *M*, i.e., number of <a href="https://www.codecogs.com/eqnedit.php?latex=J_{ij}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?J_{ij}" title="J_{ij}" /></a> configurations for the configurational average (conf_num);
 * Number of sweeps in the transient period (thermal);
 * The minimum value of <a href="https://www.codecogs.com/eqnedit.php?latex=\mu" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mu" title="\mu" /></a> (mu_min), maximum value of <a href="https://www.codecogs.com/eqnedit.php?latex=\mu" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mu" title="\mu" /></a> (mu_max), and the step (mu_step) between the maximum and minimum values. Therefore, the number of values of <a href="https://www.codecogs.com/eqnedit.php?latex=\mu" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mu" title="\mu" /></a> the are scanned is 1+(mu_max-mu_min)/mu_step.
 * The minimum value of <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\sigma" title="\sigma" /></a> (sd_min), maximum value of <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\sigma" title="\sigma" /></a> (sd_max), and the step (sd_step) between the maximum and minimum values of <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\sigma" title="\sigma" /></a>.
Therefore, the number of values of  <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\sigma" title="\sigma" /></a> the are scanned is 1+(sd_max-sd_min)/sd_step.

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

The program outputs the following observables:
 * spin glass susceptibility (Xsg),
 * uniform susceptibility (Xuni),
 * spin glass order parameter (q),
 * magnetization (m),
 * specific heat (c).

The outputs of the program are multiple .txt files, one for each pair of (mu, sd).
Each .txt file contains only one line with seven values in this order: mu, sd, Xsg, Xuni, q, m, c.
If the user sets *n* values for the mu array and *m* values for the sd array, the program outputs *nm* .txt files. The files are named as: file_0_0.txt, file_0_1.txt, ..., file_1_0.txt, file_1_1.txt... The reason for generating multiple .txt files is for parallelizing the program as explained next.

6) Parallelizing

The program has a built-in naive parallelization for scanning the values of <a href="https://www.codecogs.com/eqnedit.php?latex=\mu" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mu" title="\mu" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\sigma" title="\sigma" /></a>, i.e., after compiling the code, if the user runs the same program on 10 different cores, for example, the total computation time will decrease 10 times.
 
7) Gather the outputs

* Download the jupyter notebook sk_heatmap.ipynb in the same folder where the outputs are.
* Open the notebook.
* Set the parameters accordingly, as indicated in the notebook.
* Run the python code.
* A figure for each measure is produced.


