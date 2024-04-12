#!/usr/bin/env Rscript

# Read the cake types from the input file
cake_types <- scan("input.csv", what = "character", sep = ",")

# Split the cake types into individual strings
cake_types <- unlist(strsplit(cake_types, ","))

# Generate the numbers for each cake type
numbers <- rpois(length(cake_types), 10)

# Combine the cake types and numbers into a data frame
cake_data <- data.frame(cake_type = cake_types, number = numbers)

# Write the data frame to a CSV file
write.csv(cake_data, file = "cake_data.csv", row.names = FALSE)
