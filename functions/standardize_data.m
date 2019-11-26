function output_standardized = standardize_data(input_data)

output_standardized = (input_data - mean(input_data))/std(input_data);