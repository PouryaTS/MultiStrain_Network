# Multi-strain competition on network structure
This repository provides the source code for the following study: Pourya Toranj Simin, Juliana Taube, Elisabeta Vergu, Shweta Bansal, Lulla Opatowski, Chiara Poletto. "Contact structure and population immunity shape the selective advantage of emerging variants"
## Simulation on network 
This code simulated the spread and circulation of one, two, or three pathogens on a network structure using the SIR compartmental model.  
The main parameters of the model are:
| Parameter       |        Description       | values          | 
|-----------------|:-------------------------|:---------------:|
|$$N$$            | Number of nodes          | $10^4$      |
|$$\mu_1$$       |recovery rate for the first strain            |  $$\mu_1 = \frac{1}{7(days)}$$   | 
|$$\mu_2$$       |recovery rate for the second strain            |  $$\mu_2=\frac{\mu_1}{\tau}$$      | 
|$$\beta_1$$     |transmission rate for the first strain        |  $$\beta_1=0.0185$$      | 
|$$\beta_2$$     |transmission rate for the second strain        |  $$\beta_2=r_\beta \beta_1$$      | 
|$$\tau$$        |scaling factor of the recovery rate of the two strains |  1  | 
|$$r_\beta$$     |ratio of transmission rate of the two strains |  $$1.5$$  |
|$$\sigma$$      |cross-immunity ($\sigma=0$ full cross-immunity, $\sigma=1$ no interaction)|  $$[0,1]$$  |



### How to use
To run the executable, you need:

*	A configuration file specifying the model parameters,
*	A CSV edge list file describing the network,
*	A label that will be used to name the output files

 On Linux, you can execute the program as follows:
 
`./MultiStrainSIRonNet.exe $path_to_configuration_file $path_to_network_edgelist_file $OutputFileLabel`

