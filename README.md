# R-INLA-Code-Solea-HBSTM

Identifying persistent abundance areas for common sole in the northern Iberian waters

In this repository we can find the R-INLA code from the two final Bayesian spatio-temporal models fitted in the manuscript named above.

- Exploratory analysis code not included, Protocol from Zuur et al. (2016) was applied.
- Selection of spatio temporal structures from Iosu Paradinas can be found [here](https://github.com/FranIzquierdo/R-INLA-Code-Hake-HBSTM/blob/main/0_Spatio_temporal_structure_comparison.R).

1) On the first script, we can find the model fitting process for the catch per unit effort (CPUE) survey data of common sole (Kg), where the response variables are the probability of ocurrence (Binomial) and conditional to presence abundance (Gamma):

- Code to project WGS84 raster files and data locations in UTM included.
- Useful distance function to select the pc.prior for the spatial effect and mesh.
- Code for plotting model hyperparameters and yearly spatio-temporal effect maps included.

2) On the second script, we have the model fitting process for the length of common sole (cm), where the response variable is Gaussian:

- Useful function to select the pc.prior for the spatial effect and mesh.
- Code for plotting model hyperparameters and spatio-temporal effects included.
- Code for preparing prediction data, prediction stack and prediction spatio-temporal maps. Prediction is carried out through the common inla stack.obs/stack.pred approach.


References:

Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Elias T. Krainski, Virgilio Gómez-Rubio, Haakon Bakka, Amanda Lenzi, Daniela Castro-Camilo, Daniel Simpson, Finn Lindgren and Håvard Rue. CRC Press/Taylor and Francis Group, 2019. (https://becarioprecario.bitbucket.io/spde-gitbook/ch-stapp.html#sec:hgst)

Paradinas, I., Conesa, D., López-Quílez, A., and Bellido, J. M. (2017). Spatio-temporal model structures with shared components for semi-continuous species distribution modelling. Spatial Statistics 22, 434–450

Zuur, A. F., & Ieno, E. N. (2016). A protocol for conducting and presenting results of regression‐type analyses. Methods in Ecology and Evolution, 7(6), 636-645.

Zuur, A. F., Ieno, E. N., & Saveliev, A. A. (2017). Spatial, Temporal and Spatial-Temporal Ecological Data Analysis with R-INLA. Highland Statistics Ltd, 1.

