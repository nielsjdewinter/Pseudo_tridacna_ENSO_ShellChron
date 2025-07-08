# Script to post process ShellChron results and avoid negative time jumps and missed years due to low growth rate estimates
# by Niels J. de Winter
# Project "Pseudo_tridacna_ShellChron_Han_Tao.Rproj"
# Example for specimen YX1 (similar routine applied to specimen A5)

require(tidyverse)
require(ggpubr)

# Load data
dat <- read.csv("YX1/YX1.csv") # Original d18O data with year markers
age_model_results <- read.csv("YX1/Age_model_results.csv") # Results in terms of age (DOY)
Growth_rate_results <- read.csv("YX1/Growth_rate_results.csv") # Results in terms of growth rate (um/day)

# Combine relevant data into one dataframe
merged_dat <- select(dat, D, d18Oc, YEARMARKER) |>
mutate(
    mean.day = age_model_results$mean.day,
    mean.GR = Growth_rate_results$mean.GR / 1000,
    year = cumsum(dat$YEARMARKER) # Add year labels
)

# Method 1
# Construct new age-distance model using relative growth rate per datapoint

merged_dat$day_new <- NA # Add column to keep track of new age values
Days_tot <- 365 # Assume no growth stops (365 days per year recorded, the growth slowdowns are covered by the increments with slow growth rate < 1 um/day)

# Loop through growth years
for(year in 1:(max(merged_dat$year) - 1)){
    # Find distance steps
    D_diff <- diff(
        merged_dat$D[
            which(merged_dat$YEARMARKER == 1 & merged_dat$year == year) :
            which(merged_dat$YEARMARKER == 1 & merged_dat$year == year + 1)
        ]
    )
    GR <- merged_dat$mean.GR[which(merged_dat$year == year)] # Find growth rates per step
    GR[which(GR < 1)] <- 1 # Cap minimum growth rate to prevent very small GR numbers (and very large numbers of days) for intervals with very limited growth
    Times <- D_diff / GR # Calculate time steps
    Times_days <- Times * Days_tot / sum(Times) # Convert time steps to total number of days
    Times_days_cum <- cumsum(Times_days) # Find cumulative time steps
    merged_dat$day_new[which(merged_dat$year == year)] <- Times_days_cum + 365 * (year - 1) # Add whole years if needed and add new times to data
}

# Export results
write.csv(merged_dat, "YX1/Realigned_age_model_results.csv", row.names = FALSE)

# Plot resulting distance-time relationship
Age_model_plot <- ggplot(merged_dat) +
geom_point( # Plot old results
    aes(
        x = D / 1000,
        y = mean.day,
        color = "Original ShellChron outcome"
    ),
    alpha = 0.3
) +
geom_line( # Plot eraligned results
    aes(
        x = D / 1000,
        y = day_new,
        color = "Realigned ShellChron outcome"
    ),
    linewidth = 2
) +
geom_vline( # Plot position of year markers
    xintercept = merged_dat$D[which(merged_dat$YEARMARKER == 1)] / 1000
) +
scale_x_continuous("Distance (mm)") +
scale_y_continuous(
    "Age (days)",
    sec.axis = sec_axis(
        trans = ~. / 365,
        name = "Age (years)"
    )
) +
ggtitle("Age model realignment for specimen YX1") +
theme_minimal()

# Add original isotope plot
d18O_plot <- ggplot(merged_dat) +
geom_point( # Plot old results
    aes(
        x = D / 1000,
        y = d18Oc,
        color = "Original d18O values"
    )
) +
geom_line( # Plot old results
    aes(
        x = D / 1000,
        y = d18Oc,
        color = "Original d18O values"
    )
) +
geom_vline( # Plot position of year markers
    xintercept = merged_dat$D[which(merged_dat$YEARMARKER == 1)] / 1000
) +
scale_x_continuous("Distance (mm)") +
scale_y_continuous("d18O (per mille)") +
theme_minimal()

# Combine plots
Age_model_update_combined <- ggarrange(
    Age_model_plot,
    d18O_plot,
    common.legend = TRUE,
    ncol = 1,
    heights = c(2, 1),
    align = "v"
)
