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

# Combine relevant data into one dataframe and replace NAs with 0
merged_dat <- select(dat, D, d18Oc, YEARMARKER) |>
mutate(
    mean.day = replace(age_model_results$mean.day, is.na(age_model_results$mean.day), 0),
    se.day = replace(age_model_results$se.day, is.na(age_model_results$se.day), 0),
    mean.GR = replace(Growth_rate_results$mean.GR, is.na(Growth_rate_results$mean.GR), 0) / 1000,
    se.GR = replace(Growth_rate_results$se.GR, is.na(Growth_rate_results$se.GR), 0) / 1000,
    year = cumsum(dat$YEARMARKER) # Add year labels
)

# Construct new age-distance model using relative growth rate per datapoint
Days_tot <- 365 # Assume no growth stops (365 days per year recorded, the growth slowdowns are covered by the increments with slow growth rate < 1 um/day)

# Simulate 3000 iterations of time-distance model sampling growth rates and man ages from normal uncertainty distribution.
N = 3000

# Create matrix to store all simulated age models
simulated_age_models <- matrix(NA, nrow = nrow(merged_dat), ncol = N)

for(i in 1:N){
    print(paste0("Simulating age model iteration ", i, " of ", N)) # Progress indicator
    # Loop through growth years
    for(year in 1:(max(merged_dat$year) - 1)){
        # Find distance steps
        D_diff <- diff(
            merged_dat$D[
                which(merged_dat$YEARMARKER == 1 & merged_dat$year == year) :
                which(merged_dat$YEARMARKER == 1 & merged_dat$year == year + 1)
            ]
        )
        # Sample growth rates from normal distribution, but prevent negative growth rates
        GR_sampled <- rnorm(
            n = length(D_diff),
            mean = merged_dat$mean.GR[which(merged_dat$year == year)],
            sd = merged_dat$se.GR[which(merged_dat$year == year)]
        )
        GR_sampled[which(GR_sampled < 0)] <- 0
        GR_sampled[which((GR_sampled < 0.001) & (GR_sampled > 0))] <- 0.001 # Cap minimum non-zero growth rate to prevent very small GR numbers (and very large numbers of days) for intervals with very limited growth
        first_sample_time <- rnorm( # Sample first time point from normal distribution
            n = 1,
            mean = merged_dat$mean.day[1],
            sd = merged_dat$se.day[1]
        )
        Times <- D_diff / GR_sampled # Calculate time steps
        Times[which(Times == Inf)] <- 0 # Handle growth stops
        Times_days <- first_sample_time + Times * Days_tot / (sum(Times) + first_sample_time) # Convert time steps to total number of days
        Times_days_cum <- cumsum(Times_days) # Find cumulative time steps
        if(year == 1){
            simulated_age_models[which(merged_dat$year == year), i] <- Times_days_cum
        } else {
            simulated_age_models[which(merged_dat$year == year), i] <- Times_days_cum + 365 * (year - 1)
        }
    }
    # Fill in last year
    last_year <- max(merged_dat$year) # Isolate last year
    Times_last <- diff( # Find steps in time domain for last year
        merged_dat$mean.day[
            (which(merged_dat$YEARMARKER == 1 & merged_dat$year == last_year) - 1) :
            nrow(merged_dat)
        ]
    )
    simulated_age_models[which(merged_dat$year == last_year), i] <-
        simulated_age_models[which(merged_dat$YEARMARKER == 1 & merged_dat$year == last_year) - 1, i] +
        cumsum(Times_last) # Add time values for last year based on original time steps
}

# Calculate median and 95% confidence intervals of simulated age models
merged_dat$day_new <- apply(simulated_age_models, 1, mean)
merged_dat$day_new_lower <- apply(simulated_age_models, 1, quantile, probs = 0.025)
merged_dat$day_new_upper <- apply(simulated_age_models, 1, quantile, probs = 0.975)

# Export results
write.csv(merged_dat, "YX1/Realigned_age_model_results_MC.csv", row.names = FALSE)
write.csv(simulated_age_models, "YX1/Realigned_age_model_simulations.csv", row.names = FALSE)

# Plot resulting distance-time relationship
Age_model_plot <- ggplot(merged_dat) +
geom_ribbon( # Plot uncertainty ribbon
    aes(
        x = D,
        ymin = day_new_lower,
        ymax = day_new_upper,
        fill = "Realigned ShellChron outcome"
    ),
    alpha = 0.5
) +
geom_line( # Plot realigned results
    aes(
        x = D,
        y = day_new,
        color = "Realigned ShellChron outcome"
    ),
    linewidth = 2
) +
geom_point( # Plot old results
    aes(
        x = D,
        y = mean.day,
        color = "Original ShellChron outcome"
    ),
    alpha = 0.3
) +
geom_vline( # Plot position of year markers
    xintercept = merged_dat$D[which(merged_dat$YEARMARKER == 1)]
) +
scale_x_continuous("Distance (mm)") +
scale_y_continuous(
    "Age (days)",
    sec.axis = sec_axis(
        trans = ~. / 365,
        name = "Age (years)"
    )
) +
scale_color_manual(
    name = "Age model",
    values = c(
        "Realigned ShellChron outcome" = "blue",
        "Original ShellChron outcome" = "red"
    )
) +
scale_fill_manual(
    name = "Age model",
    values = c(
        "Realigned ShellChron outcome" = "lightblue"
    )
) +
ggtitle("Age model realignment for specimen YX1") +
theme_minimal()

# Add original isotope plot
d18O_plot <- ggplot(merged_dat) +
geom_point( # Plot old results
    aes(
        x = D,
        y = d18Oc,
        color = "Original d18O values"
    )
) +
geom_line( # Plot old results
    aes(
        x = D,
        y = d18Oc,
        color = "Original d18O values"
    )
) +
geom_vline( # Plot position of year markers
    xintercept = merged_dat$D[which(merged_dat$YEARMARKER == 1)]
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
