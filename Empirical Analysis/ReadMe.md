#  Empirical Data Description

This folder contains the aggregated data used for the **empirical analysis** of the Alpha variant’s spread across U.S. states.

## Data File

The file **`USCeoffDataforRegression_alpha_period.csv`** includes the following fields:

- **`EstimatedCoeff`** — The estimated **selection coefficient** of the Alpha variant for each U.S. state.  
- **`IM`** — The estimated **immunity level** at the time of the Alpha variant’s emergence in each state.  
- **`MD`** and **`SD`** — Represent the **average over time** of (i) the **mean** number of non-household contacts and (ii) the **standard deviation** of the number of contacts, respectively, during the Alpha variant period.


##  Analysis Script

The script **`ridge_regression.R`** performs a **ridge regression analysis** on the empirical dataset to quantify associations between the Alpha variant’s selection coefficient and contact patterns or immunity levels.
