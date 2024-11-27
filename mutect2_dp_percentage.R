source("R/libraries.R")

#coverage
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



#read_length
folder2 = "D:/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/read_length/"
caller = "Mutect2"

#TP
read50 = fread(paste0(folder2, "1000_50","/Merged_", caller, "_snvs_TP.tsv"))
read50= read50[!is.na(read50$`Mutect2 DP`),]
read50$delta = read50$`Mutect2 AF` - read50$`Ground Truth AF`
read50$delta_dp = read50$`Mutect2 DP` - read50$`Ground Truth DP`
read50$DP_percentage = read50$`Mutect2 DP`/read50$`Ground Truth DP`

read75 = fread(paste0(folder2, "1000_75","/Merged_", caller, "_snvs_TP.tsv"))
read75= read75[!is.na(read75$`Mutect2 DP`),]
read75$delta = read75$`Mutect2 AF` - read75$`Ground Truth AF`
read75$delta_dp = read75$`Mutect2 DP` - read75$`Ground Truth DP`
read75$DP_percentage = read75$`Mutect2 DP`/read75$`Ground Truth DP`

read100 = fread(paste0(folder2, "1000_100","/Merged_", caller, "_snvs_TP.tsv"))
read100= read100[!is.na(read100$`Mutect2 DP`),]
read100$delta = read100$`Mutect2 AF` - read100$`Ground Truth AF`
read100$delta_dp = read100$`Mutect2 DP` - read100$`Ground Truth DP`
read100$DP_percentage = read100$`Mutect2 DP`/read100$`Ground Truth DP`

read150 = fread(paste0(folder1, "1000_100_10","/Merged_", caller, "_snvs_TP.tsv"))
read150= read150[!is.na(read150$`Mutect2 DP`),]
read150$delta = read150$`Mutect2 AF` - read150$`Ground Truth AF`
read150$delta_dp = read150$`Mutect2 DP` - read150$`Ground Truth DP`
read150$DP_percentage = read150$`Mutect2 DP`/read150$`Ground Truth DP`

read300 = fread(paste0(folder2, "1000_300","/Merged_", caller, "_snvs_TP.tsv"))
read300= read300[!is.na(read300$`Mutect2 DP`),]
read300$delta = read300$`Mutect2 AF` - read300$`Ground Truth AF`
read300$delta_dp = read300$`Mutect2 DP` - read300$`Ground Truth DP`
read300$DP_percentage = read300$`Mutect2 DP`/read300$`Ground Truth DP`


summary(read50$DP_percentage)
summary(read75$DP_percentage)
summary(read100$DP_percentage)
summary(read150$DP_percentage)


# Combine coverage data into one data.table with a 'Coverage' column
cov_data <- rbind(
    cov300[, .(DP_percentage, Coverage = "300")],
    cov700[, .(DP_percentage, Coverage = "700")],
    cov1000[, .(DP_percentage, Coverage = "1000")],
    cov3000[, .(DP_percentage, Coverage = "3000")],
    cov5000[, .(DP_percentage, Coverage = "5000")]
)
cov_data[, Coverage := factor(Coverage, levels = c("300", "700", "1000", "3000", "5000"))]  # Ensure ascending order

# Calculate mean DP_percentage for coverage
coverage_means <- cov_data[, .(Mean_DP = mean(DP_percentage, na.rm = TRUE)), by = Coverage]

# Combine read length data into one data.table with a 'ReadLength' column
read_data <- rbind(
    read50[, .(DP_percentage, ReadLength = "50")],
    read75[, .(DP_percentage, ReadLength = "75")],
    read100[, .(DP_percentage, ReadLength = "100")],
    read150[, .(DP_percentage, ReadLength = "150")],
    read300[, .(DP_percentage, ReadLength = "300")]
)
read_data[, ReadLength := factor(ReadLength, levels = c("50", "75", "100", "150", "300"))]  # Ensure ascending order

# Calculate mean DP_percentage for read length
read_means <- read_data[, .(Mean_DP = mean(DP_percentage, na.rm = TRUE)), by = ReadLength]

# Merge data to add mean values for consistent coloring
cov_data <- merge(cov_data, coverage_means, by = "Coverage")
read_data <- merge(read_data, read_means, by = "ReadLength")

# Find the global range of means across both datasets
global_min <- min(c(coverage_means$Mean_DP, read_means$Mean_DP))
global_max <- max(c(coverage_means$Mean_DP, read_means$Mean_DP))

# Define a common color scale based on global min and max
color_scale <- scale_fill_gradient(
    low = "red",
    high = "blue",
    limits = c(global_min, global_max)  # Ensures consistent color mapping
)

# Plot for Coverage using mean-based color
coverage_plot <- ggplot(cov_data, aes(x = Coverage, y = DP_percentage, fill = Mean_DP)) +
    geom_boxplot(width = 0.4) +
    color_scale +
    labs( x = "Varying Coverage", y = "DP Percentage") +
    theme_minimal() +
    ylim(0, 1.00)

# Plot for Read Length using mean-based color
read_length_plot <- ggplot(read_data, aes(x = ReadLength, y = DP_percentage, fill = Mean_DP)) +
    geom_boxplot(width = 0.4) +
    color_scale +
    labs(x = "Varying Read Length", y = "DP Percentage") +
    theme_minimal() +
    ylim(0, 1.00)

# Add centered title to the combined plot
combined_plot <- coverage_plot / read_length_plot +
    plot_annotation(
        title = "DP Percentage utilized by Mutect2",
        theme = theme(plot.title = element_text(hjust = 0.5))
    )

# Print the combined plot
print(combined_plot)
