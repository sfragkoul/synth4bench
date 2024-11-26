source("R/libraries.R")


# folder = "D:/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/coverage_test/300_30_10"
# caller = "Mutect2"
# 
# 
# snvs_tp = fread(paste0(folder, "/Merged_", caller, "_snvs_TP.tsv"))
# snvs_fp = fread(paste0(folder, "/Merged_", caller, "_snvs_FP.tsv"))
# snvs_fn = fread(paste0(folder, "/Merged_", caller, "_snvs_FN.tsv"))

folder1 = "D:/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/coverage_test/"
caller = "Mutect2"
#TP
cov300 = fread(paste0(folder1, "300_30_10","/Merged_", caller, "_snvs_TP.tsv"))
cov300= cov300[!is.na(cov300$`Mutect2 DP`),]
cov300$delta = cov300$`Mutect2 AF` - cov300$`Ground Truth AF`
cov300$delta_dp = cov300$`Mutect2 DP` - cov300$`Ground Truth DP`
cov300$DP_percentage = cov300$`Mutect2 DP`/cov300$`Ground Truth DP`

cov700 = fread(paste0(folder1, "700_70_10","/Merged_", caller, "_snvs_TP.tsv"))
cov700= cov700[!is.na(cov700$`Mutect2 DP`),]
cov700$delta = cov700$`Mutect2 AF` - cov700$`Ground Truth AF`
cov700$delta_dp = cov700$`Mutect2 DP` - cov700$`Ground Truth DP`
cov700$DP_percentage = cov700$`Mutect2 DP`/cov700$`Ground Truth DP`

cov1000 = fread(paste0(folder1, "1000_100_10","/Merged_", caller, "_snvs_TP.tsv"))
cov1000= cov1000[!is.na(cov1000$`Mutect2 DP`),]
cov1000$delta = cov1000$`Mutect2 AF` - cov1000$`Ground Truth AF`
cov1000$delta_dp = cov1000$`Mutect2 DP` - cov1000$`Ground Truth DP`
cov1000$DP_percentage = cov1000$`Mutect2 DP`/cov1000$`Ground Truth DP`

cov3000 = fread(paste0(folder1, "3000_300_10","/Merged_", caller, "_snvs_TP.tsv"))
cov3000= cov3000[!is.na(cov3000$`Mutect2 DP`),]
cov3000$delta = cov3000$`Mutect2 AF` - cov3000$`Ground Truth AF`
cov3000$delta_dp = cov3000$`Mutect2 DP` - cov3000$`Ground Truth DP`
cov3000$DP_percentage = cov3000$`Mutect2 DP`/cov3000$`Ground Truth DP`

cov5000 = fread(paste0(folder1, "5000_500_10","/Merged_", caller, "_snvs_TP.tsv"))
cov5000= cov5000[!is.na(cov5000$`Mutect2 DP`),]
cov5000$delta = cov5000$`Mutect2 AF` - cov5000$`Ground Truth AF`
cov5000$delta_dp = cov5000$`Mutect2 DP` - cov5000$`Ground Truth DP`
cov5000$DP_percentage = cov5000$`Mutect2 DP`/cov5000$`Ground Truth DP`


summary(cov300$DP_percentage)
summary(cov700$DP_percentage)
summary(cov1000$DP_percentage)
summary(cov3000$DP_percentage)






# List of datasets with corresponding coverage values
datasets <- list(
    "300" = cov300,
    "700" = cov700,
    "1000" = cov1000,
    "3000" = cov3000,
    "5000" = cov5000
)

# Initialize a data frame to store results
results <- data.frame(
    Coverage = numeric(),
    Metric = character(),
    Mean = numeric(),
    Sigma = numeric(),
    Count = numeric()
    
)

# Calculate mean, standard deviation, and mutation count for each dataset
for (coverage in names(datasets)) {
    data <- datasets[[coverage]]
    
    # Calculate metrics
    mean_af <- mean(data$`Mutect2 AF`, na.rm = TRUE)
    sigma_af <- sd(data$`Mutect2 AF`, na.rm = TRUE)
    
    mean_dp <- mean(data$`Mutect2 DP`, na.rm = TRUE)
    sigma_dp <- sd(data$`Mutect2 DP`, na.rm = TRUE)
    
    mutation_count <- nrow(data)  # Number of mutations
    
    mean_delta <-  mean(data$delta, na.rm = TRUE)
    sigma_delta <- sd(data$delta, na.rm = TRUE)
    
    mean_delta_dp <-  mean(data$delta_dp, na.rm = TRUE)
    sigma_delta_dp <- sd(data$delta_dp, na.rm = TRUE)
    
    # Append results for AF
    results <- rbind(results, data.frame(
        Coverage = as.numeric(coverage),
        Metric = "AF",
        Mean = mean_af,
        Sigma = sigma_af,
        Count = NA  # Count is not relevant for AF
    ))
    
    # Append results for DP
    results <- rbind(results, data.frame(
        Coverage = as.numeric(coverage),
        Metric = "DP",
        Mean = mean_dp,
        Sigma = sigma_dp,
        Count = NA  # Count is not relevant for DP
    ))
    
    # Append results for Mutation count
    results <- rbind(results, data.frame(
        Coverage = as.numeric(coverage),
        Metric = "Mutations",
        Mean = mutation_count,
        Sigma = NA,  # Standard deviation is not relevant for counts
        Count = mutation_count
    ))
    
    
    # Append results for delta
    results <- rbind(results, data.frame(
        Coverage = as.numeric(coverage),
        Metric = "Delta AF",
        Mean = mean_delta,
        Sigma = sigma_delta,
        Count = NA  # Count is not relevant for AF
    ))
    
    # Append results for delta DP
    results <- rbind(results, data.frame(
        Coverage = as.numeric(coverage),
        Metric = "Delta DP",
        Mean = mean_delta_dp,
        Sigma = sigma_delta_dp,
        Count = NA  # Count is not relevant for AF
    ))
}
rm(cov1000, cov300, cov3000, cov700, cov5000)


# Load necessary libraries
library(ggplot2)

# Filter results for Delta AF and Delta DP
delta_results <- subset(results, Metric %in% c("Delta AF", "Delta DP"))

# Create the plot
ggplot(delta_results, aes(x = Coverage, y = Mean, color = Metric, group = Metric)) +
    geom_point(size = 3) +                          # Points for each metric
    geom_line() +                                   # Lines to connect the points
    scale_x_log10() +                               # Log scale for coverage
    labs(
        title = "Delta AF and Delta DP Across Different Coverage Scenarios",
        x = "Coverage (log scale)",
        y = "Mean Delta Value",
        color = "Metric"
    ) +
    theme_minimal()

# Create separate plots for Delta AF and Delta DP using facet_wrap
ggplot(delta_results, aes(x = Coverage, y = Mean, color = Metric, group = Metric)) +
    geom_point(size = 3) +                          # Points for each metric
    geom_line() +                                   # Lines to connect the points
    scale_x_log10() +                               # Log scale for coverage
    facet_wrap(~ Metric, scales = "free_y") +       # Facet by Metric (Delta AF and Delta DP)
    labs(
        title = "Delta AF and Delta DP Across Different Coverage Scenarios",
        x = "Coverage (log scale)",
        y = "Mean Delta Value"
    ) +
    theme_minimal()
