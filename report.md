# Tow-pathogen competition on network structure

A SIR model is used to capture the behaviour of spreading two pathogens. The model includes 8 compartments regarding the state of individuals with respect to the disease. An individual can be in one of the *Susceptible*, *Infected* or *Recovered* states with respect to each pathogens. We assume a recoverd individual from one infection have a partial immunity against the other one and it's susceptibility reduced by a factor $\sigma$.

The compartmental model is illustrated bellow:

![Image of ISR](figs/SIRModel.png)

## Simulation on networks

In the first step, we assume that the network is homogeneous and we simulate the SIR model of spreading two pathogens on random network of Edos Reiny.

<span style="color: red;">- To do: Upgrade codes to run the simulation for any desire network</span>

For each set of parameter we simulate 1000 stochastic realization and at each time step of simulation we monitor the current state of the system by recording the number of individuals in each compartment and the the number of new incidence rate for each pathogen. Also at each time step, every transmition event is recorded. The simulation continues until no infected individual remains.  The parameters of the model and their values are reported in the table 1.

| Parameter       |        Description       | values          | 
|-----------------|:-------------------------|:---------------:|
|$$N$$            | Number of nodes          | $10^4$      |
|$$\bar k$$       | Mean degree of network   | 5    |
|$$R_0^f = \bar k \frac{\beta_f} {\gamma_f} $$ | reproductive number for the fast strain | 2 |
|$$R_0^s = \bar k \frac{\beta_s} {\gamma_s} $$ | reproductive number for the fast strain | $R_0^s = r R_0^f$ | 
|$$\mu_f$$       |recovery rate for the fast strain            |  0.6   | 
|$$\mu_s$$       |recovery rate for the slow strain            |  $$\mu_s^{-1}=\tau \mu_f^{-1}$$      | 
|$$\beta_f$$     |transmission rate for the fast strain        |  $$\beta_f=R_0^f \mu_f$$      | 
|$$\beta_s$$     |transmission rate for the slow strain        |  $$\beta_s=R_0^s \mu_s$$      | 
|$$\tau$$        |timescale separation between the two pathogens |  1.5  | 
|$$r$$           |ratio of basic reproductive numbers of the two pathogens |  $$[0.6,1.6]$$  |
|$$\sigma$$      |cross-immunity ($\sigma=0$ full cross-immunity, $\sigma=1$ no interaction)|  $$[0,1]$$  |



## Results

In the following the final attack rate for each pathogen and the total attack rate is ilustrated as function of cross-immunity $\sigma$ and the ratio of reprodiction numbers $r$.
![Image of ISR](figs/Attackrate_f_heatmap.png "Title") ![Image of ISR](figs/Attackrate_s_heatmap.png)  

![Image of ISR](figs/Attackrate_B_heatmap.png) ![Image of ISR](figs/Ratio_sf_heatmap.png ) 