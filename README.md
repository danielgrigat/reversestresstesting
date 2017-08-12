# reversestresstesting
Underlying STOXX data and code to run experiments for "Reverse stress testing interbank networks"

Data source: STOXX_banks.xlsx

Code implemented in MATLAB.

0. Function to load data from .xlsx file is import_stoxx1.m

1. Function to construct STOXX network, i.e. the RAS algorithm: fitness_weights.m<br />
	requires fitness_model_param.m to fit the density parameter.

2. Main function for framework is: reverse_fun1_SREP.m.<br />For example to run the framework for the full year 2015, with \lambda max in the range from 0.1 to 1 insteps of 0.1, for targeted losses = 1, time horizon = 20 enter:
[E, E_node, IPR, loss_frac, h, adj, u, flag] = reverse_fun1_SREP(0.1:0.1:1, 1, 20, 2, 0, '2015');

3. Function to run policy experiment: reverse_policy_SREP.m<br /> Note that this function requires that the main function has been run. For this reason I provide the main files of the reverse stress testing framework results for 2014, 2015, 2016, in the folders of the same name, that can serve as inputs.<br /> For example to run the policy for the full year 2015, target losses = 0.1, for the K_i-based policy with a 5% equity increase, enter: [R1] = reverse_policy_SREP(1,0.1,2,0.05,1,2015)
