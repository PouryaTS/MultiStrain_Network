# Multi-strain competition on network structure
This repository provides the source code for the following study: Pourya Toranj Simin, Juliana Taube, Elisabeta Vergu, Shweta Bansal, Lulla Opatowski, Chiara Poletto. "Contact structure and population immunity shape the selective advantage of emerging variants"
## Simulation on network 
This code simulated the spread and circulation of one, two, or three pathogens on a network structure using the SIR compartmental model.  
The main parameters of the model are:
| Parameter       |        Description       | values          | 
|-----------------|:-------------------------|:---------------:|
|$$N$$            | Number of nodes          | $10^4$      |
|$$\mu_1$$, $$\mu_2$$, $$\mu_3$$        |Recovery rate for the first, second and third strain            |  $$\mu_1 = \frac{1}{7(days)}$$,  $$\mu_2=\frac{\mu_1}{\tau_2}$$, $$\mu_3=\frac{\mu_1}{\tau_3}$$    | 
|$$\beta_1$$, $$\beta_2$$, $$\beta_3$$       |Transmission rate for the first, second and third strain       |  $$\beta_1=0.0185$$, $$\beta_2=r_2 \beta_1$$, $$\beta_3=r_3 \beta_1$$     | 
|$$\tau_2$$, $$\tau_3$$        |Scaling factor of the recovery rate of the second and third strain with respect to the first strain |  $$\tau_2=1$$, $$\tau_3=1$$     | 
|$$r_2$$, $$r_3$$     |Ratio of transmission rate of the second and third strain with respect to the first strain|    |
|$$\sigma_2$$      |Cross-immunity between the second variant and the first variant ($\sigma_2=0$ full cross-immunity, $\sigma_2=1$ no cross-immunity)|  $$[0,1]$$  |
|$$\sigma_3$$      |Cross-immunity between the third variant and the first and second variants ($\sigma_3=0$ full cross-immunity, $\sigma_3=1$ no cross-immunity)|  $$[0,1]$$  |
|$$R_1(t_2)$$      |Immunity level of the first variant (cumulative number of the infected) when the second variant is introduced into the population |  $$[0,1]$$  |
|$$\delta t$$      |Time of the emergence of the third variant since the emergence of the second variant in days | $$\ge0$$ |



### How to use
To run the executable, you need:

*	A configuration file specifying the model parameters,
*	A CSV edge list file describing the network,
*	A label that will be used to name the output files

 On Linux, you can execute the program as follows:
 
`./MultiStrainSIRonNet.exe $path_to_configuration_file $path_to_network_edgelist_file $OutputFileLabel`

### 🧩 Configuration file

Below is an example of a configuration file (`ModelConfig_2strains.txt`) used to set the model parameters.

```txt
beta1   0.0175
mu1     0.143
tau2    1
tau3    1
r2      1.5 1.6 0.25
r3      1
sigma2  0.05
sigma3  0 0.1 0.2
R1t2    0.1 0.11 0.1
deltat  0 1 3
Iinit   50 50 0
itr     100
```

The meaning of each parameter is explained in the following table:


| **Parameter** | **Example value(s)** | **Description** |
|----------------|----------------------|------------------|
| `beta1` | `0.0175` | Transmission rate of the first variant. |
| `mu1` | `0.143` | Recovery rate (inverse of the infectious period) for the first variant. |
| `tau2`, `tau3` | `1`, `1` | Scaling factors of the recovery rate for variants 2 and 3, relative to the baseline variant. |
| `r2` | `1.5 1.6 0.25` | Range of values for parameter $$r_2$$ (start = 1.5, end = 1.6, step = 0.25). |
| `r3` | `1` | Fixed value of parameter $$r_3$$. |
| `sigma2` | `0.05` | cross immunity of the second variant  $$\sigma_2$$. |
| `sigma3` | `0 0.1 0.2` | cross immunity of the third variant $$\sigma_3$$ (start = 0, end = 0.1, step = 0.2). |
| `R1t2` | `0.1 0.11 0.1` | Range of values for parameter $$R_1(t_2)$$ |
| `deltat` | `0 1 3` | Range of time delays (start = 0, end = 1, step = 3). |
| `Iinit` | `50 50 0` | Initial number of infected individuals for the first, second, and third variants, respectively. |
| `itr` | `100` | Number of stochastic iterations to perform. |

**💡 Notes:**

* Parameters with three numbers (like r2, sigma3, R1t2, deltat) represent vectors:
[start, end, step]. The program iterates through all values in that range.

* This example of a configuration file is designed for a 2-strain scenario where there is no third strain. You can modify these parameters to explore different epidemic scenarios.

