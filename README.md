# reversestresstesting
Underlying STOXX data and code to run experiments for "Reverse stress testing interbank networks"

Data source: STOXX_banks.xlsx

Code implemented in MATLAB.

0. Function to load data from .xlsx file is import_stoxx1.m

1. Function to construct STOXX network, i.e. the RAS algorithm: fitness_weights.m
	requires fitness_model_param.m to fit the density parameter.

2. Main function for framework is: reverse_fun1_SREP.m

3. Function to run policy experiment: reverse_policy_SREP.m 
	note that this function requires that the main function has been run. For this reason I provide the main files of the reverse stress testing framework results for 2014, 2015, 2016, in the folders of the same name, that can serve as inputs.
