# Test ShellChron on pseudo-tridacna records
# Data by Han Tao
# ShellChron script by Niels J. de Winter
# Project "Pseudo_tridacna_ShellChron_Han_Tao"

# devtools::install_github("nielsjdewinter/ShellChron") # Install ShellChron from GitHub

# Load packages
require(ShellChron)
require(signal)
require(tidyverse)

# ------------------------------------------------------------------------------
# Dataset YX1
dat <- read.csv("YX1/YX1.csv") |> # Load data
mutate( # Multiply D values by 1000
    D = D * 1000,
    D_err = D_err * 1000
)

# Run ShellChron model normally
wrap_function(
    input_from_file = FALSE, # Should input be read from a file?
    object_name = dat, # Name of object with input (only if input_from_file = FALSE)
    transfer_function = "GrossmanKu86", # Set transfer function of the record, default is Kim and O'Neil 1997.
    t_int = 1, # Set time interval in days
    T_per = 365, # Set annual time period in days (default = 365)
    d18Ow = -1, # Set d18Ow value or vector (default = constant year-round at 0 VSMOW). Alternative options are either one value (assumed constant year-round) or a vector with length T_per / t_int and interval t_int specifying d18Ow evolution through one year.
    t_maxtemp = 182.5, # Define the day of the year at which temperature is heighest. Default = Assume that the day of maximum temperature is helfway through the year
    SCEUApar = c(1, 25, 10000, 5, 0.01, 0.01), # Set parameters for SCEUA optimization (iniflg, ngs, maxn, kstop, pcento, peps)
    sinfit = TRUE, # Apply sinusoidal fitting to guess initial parameters for SCEUA optimization? (TRUE/FALSE)
    MC = 1000, # Number of MC simulations to include measurement error into error analysis. Default = 1000 (if MC = 0, error on D and d18O measurements not considered)
    plot = TRUE, # Should intermediate plots be given to track progress? WARNING: plotting makes the script much slower, especially for long datasets.
    plot_export = TRUE, # Should a plot of the results be saved as PDF?
    export_raw = TRUE, # Should the results of all individual model runs be exported as CSV files?
    export_path = "C:/Users/niels/Dropbox/Research/Side projects/Pseudo_tridacna_ShellChron_Han_Tao/YX1"
)
# DONE

# ------------------------------------------------------------------------------
# Dataset A5
dat <- read.csv("A5/A5.csv") |> # Load data
mutate( # Multiply D values by 1000
    D = D * 1000,
    D_err = D_err * 1000
)

# Break down to catch errors
# ------------------------------------------------------------------------------
# Set wrapper parameters
input_from_file = FALSE # Should input be read from a file?
object_name = dat # Name of object with input (only if input_from_file = FALSE)
transfer_function = "GrossmanKu86" # Set transfer function of the record, default is Kim and O'Neil 1997.
t_int = 1 # Set time interval in days
T_per = 365 # Set annual time period in days (default = 365)
d18Ow = -1 # Set d18Ow value or vector (default = constant year-round at 0 VSMOW). Alternative options are either one value (assumed constant year-round) or a vector with length T_per / t_int and interval t_int specifying d18Ow evolution through one year.
t_maxtemp = 182.5 # Define the day of the year at which temperature is heighest. Default = Assume that the day of maximum temperature is helfway through the year
SCEUApar = c(1, 25, 10000, 5, 0.01, 0.01) # Set parameters for SCEUA optimization (iniflg, ngs, maxn, kstop, pcento, peps)
sinfit = TRUE # Apply sinusoidal fitting to guess initial parameters for SCEUA optimization? (TRUE/FALSE)
MC = 1000 # Number of MC simulations to include measurement error into error analysis. Default = 1000 (if MC = 0, error on D and d18O measurements not considered)
plot = TRUE # Should intermediate plots be given to track progress? WARNING: plotting makes the script much slower, especially for long datasets.
plot_export = TRUE # Should a plot of the results be saved as PDF?
export_raw = TRUE # Should the results of all individual model runs be exported as CSV files?
export_path = "C:/Users/niels/Dropbox/Research/Side projects/Pseudo_tridacna_ShellChron_Han_Tao/A5"

# ------------------------------------------------------------------------------
# STEP 1: Import data
if(input_from_file){
    setwd(path)
    importlist <- data_import(file_name)
}else{
    importlist <- data_import_object(object_name)
}
if(length(importlist) != 2){ # Catch errors in the input data
    return("ERROR: Input data does not match the default input data format")
}
dat <- importlist[[1]]
dynwindow <- importlist[[2]]
G_per <- T_per # Period of growth rate sinusoid should equal that of the temperature sinusoid (which is given)

# ------------------------------------------------------------------------------
# STEP 2: Run the model
d18Oc <- d18Oc_err <- Omod <- NULL # Predefine variables to circumvent global variable binding error

# Prepare data arrays for storage of modeling results
resultarray <- array( # Create array to contain all modeling results of overlapping windows
    rep(as.matrix(cbind(dat, matrix(NA, ncol = length(dynwindow$x), nrow = length(dat$D)))), 9), # Replicate matrix of dat + extra columns to contain all variables
    dim = c(length(dat$D), length(dynwindow$x) + length(dat[1, ]), 9) # Six variables, being: Modeled d18O, residuals, Day of the Year, Growth between datapoints, Instantaneous growth rate at datapoint and Modeled temperature
)

parmat <- matrix(NA, nrow = 7, ncol = length(dynwindow$x)) # Matrix in which to store the modeling parameters
colnames(parmat) <- dynwindow$x
rownames(parmat) <- c("T_amp", "T_pha", "T_av", "G_amp", "G_pha", "G_av", "G_skw")

# Prepare plot to show model progress
if(plot == TRUE){
    dev.new()
    fitplot <- ggplot2::ggplot(dat, ggplot2::aes(D, d18Oc)) + # Create a plot showing the fit of each window on the original data, start with a plot of the original data
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = d18Oc - d18Oc_err, # Add error bars on model result (1 SD)
        ymax = d18Oc + d18Oc_err),
        width = dat$D_err) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = D - D_err,
        xmax = D + D_err),
        height = 0.05) +
    ggplot2::ggtitle("Plot of d18Oc values vs. depth. Black = data, Red = model, Errorbars = 1 SD")
    plot(fitplot)
}

# Estimate growth rate variability and round up to nearest higher magnitude of 10 for conservative boundary
GRavest <- max(diff(dat$D[dat$YEARMARKER == 1])) / 365 # Estimate maximum growth rate from yearmarkers
GRavmax <- 10 ^ (ceiling(log(GRavest, 10))) # Round up to nearest higher magnitude of 10

# Find tailored range of temperatures from data
d18Oc_range <- range(dat$d18Oc) # Find d18Oc range in data
if(transfer_function == "KimONeil97"){ # Find temperature range (to be superseded with inverse d18O_model function in later updates)
    T_range <- 18.03 * 1000 / (log((d18Oc_range - (0.97002 * rev(range(d18Ow)) - 29.98)) / 1000 + 1) * 1000 + 32.42) - 273.15 # Use Kim and O'Neil (1997) with conversion between VSMOW and VPDB by Brand et al. (2014)
}else if(transfer_function == "GrossmanKu86"){
    T_range <-  20.6 - 4.34 * (d18Oc_range - rev(range(d18Ow)) - 0.2) # Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
}else{
    print("ERROR: Supplied transfer function is not recognized")
}
T_max <- max(T_range)
T_amp_max <- 2 * abs(diff(T_range))

# Collate lower boundaries of parameters
parl <- c(
    T_amp = 0, # Minimum T amplitude in degrees C
    T_pha = 0, # Minimum phase in days
    T_av = -4, # Minimum average T in degrees C
    G_amp = 0, # Minimum seasonal GR range in um/d
    G_pha = 0, # Minimum GR phase in days
    G_av = -1 * GRavmax, # Minimum average GR in um/d.
    G_skw = 0 # Minimum skew factor
)

# Collate upper boundaries of parameters
paru <- c(
    T_amp = round(T_amp_max + 0.5, 0), # Maximum T amplitude in degrees C
    T_pha = 365, # Maximum phase in days
    T_av = round(T_max + 0.5, 0), # Maximum average T in degrees C
    G_amp = 2 * GRavmax, # Maximum seasonal GR range in um/d
    G_pha = 365, # Maximum GR phase in days
    G_av = GRavmax, # Maximum average GR in um/d (based on conservative boundaries of YEARMARKER indicators)
    G_skw = 100 # Maximum skew factor 
)

# Set parameters for SCEUA optimization
iniflg = SCEUApar[1] # Flag for initial parameter array (default = 1; included)
ngs = SCEUApar[2] # Number of complexes (sub-populations, default = 25)
maxn = SCEUApar[3] # Maximum number of function evaluations allowed during optimization (default = 10000)
kstop = SCEUApar[4] # Maximum number of evolution loops before convergency (default = 5)
pcento = SCEUApar[5] # Percentage change allowed in function value criterion before stop (default = 0.01)
peps = SCEUApar[6] # Convergence level for parameter set (difference between parameters required for stop; default = 0.01)

# Run the model on all windows
for(i in 1:length(dynwindow$x)){ # Loop over shell record
    print(paste("Processing Datawindow:", i, "of", length(dynwindow$x))) # Keep track of progress
    
    # Isolate year of d18O data based on window data
    Dsam <- dat$D[dynwindow$x[i]:(dynwindow$x[i] + dynwindow$y[i] - 1)]
    Osam <- dat$d18Oc[dynwindow$x[i]:(dynwindow$x[i] + dynwindow$y[i] - 1)]
    if(MC > 0){
        D_err <- dat$D_err[dynwindow$x[i]:(dynwindow$x[i] + dynwindow$y[i] - 1)] # Optional: include error on D
        O_err <- dat$d18Oc_err[dynwindow$x[i]:(dynwindow$x[i] + dynwindow$y[i] - 1)] # Optional: include error on d18Oc
    }else{
        D_err <- rep(0, dynwindow$y[1])
        O_err <- rep(0, dynwindow$y[1])
        MC <- 0
    }
    if(i == 1){ # Estimate starting parameters for first data window
        if(sinfit){ # If sinusoidal fitting is enabled
            sinlist <- sinreg(Dsam, Osam) # Run sinusoidal regression to find initial parameter values
            # Estimate starting parameters from regression results
            O_av_start <- sinlist[[1]][1] # Export starting value for annual d18O average
            O_amp_start <- sinlist[[1]][2] # Export starting value for d18O amplitude
            O_pha_start <- sinlist[[1]][4] %% sinlist[[1]][3] # Estimate position (in depth of the first peak in d18O)
            O_per_start <- sinlist[[1]][3] # Export starting value for period in distance domain
        }else{
            O_av_start <- mean(Osam) # Estimate starting value for annual d18O average by mean of d18O in record
            O_amp_start <- diff(range(Osam)) / 2 # Estimate starting value for d18O amplitude by half the difference between minimum and maximum d18Oc
            O_per_start <- diff(range(Dsam)) # Estimate starting period as thickness of isolated year
            O_pha_start <- 0.25 * O_per_start # Set starting phase to one quarter of a cycle if it cannot be estimated from sinusoidal regression 
        }

        if(transfer_function == "KimONeil97"){
            T_av_start <- 18.03 * 1000 / (1000 * log((O_av_start - (0.97002 * mean(d18Ow) - 29.98)) / 1000 + 1) + 32.42) - 273.15  # Estimate mean temperature. Use Kim and O'Neil (1997) with conversion between VSMOW and VPDB by Brand et al. (2014)
            T_amp_start <- 18.03 * 1000 / (1000 * log((O_av_start - O_amp_start - (0.97002 * mean(d18Ow) - 29.98)) / 1000 + 1) + 32.42) - 273.15 - T_av_start # Estimate temperature amplitude. Use Kim and O'Neil (1997) with conversion between VSMOW and VPDB by Brand et al. (2014)
        }else if(transfer_function == "GrossmanKu86"){
            T_av_start <- 20.6 - 4.34 * (O_av_start - mean(d18Ow) - 0.2) # Estimate mean temperature. Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
            T_amp_start <- 20.6 - 4.34 * (O_av_start - O_amp_start - mean(d18Ow) - 0.2) - T_av_start # Estimate mean temperature. Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
        }else{
            print("ERROR: Supplied transfer function is not recognized")
        }

        O_peak <- O_pha_start + Dsam[1] # Find position of d18O peak in distance domain
        T_pha_start <- ((O_pha_start - 0.5 * O_per_start) %% O_per_start) / O_per_start * T_per # Estimate position of first peak in temperature (low in d18O) relative to annual cycle (days)
        G_av_start <- O_per_start / G_per # Estimate average growth rate in distance/day

        # Collate starting parameters
        par0 <- c(
            T_amp = T_amp_start,
            T_pha = T_pha_start,
            T_av = T_av_start,

            G_amp = G_av_start / 2, # Start by estimating growth rate changes by half the average
            G_pha = T_pha_start, # Start by estimating that the peak in growth rate coincides with the peak in temperature
            G_av = G_av_start,
            G_skw = 50 # Start with a no skew
        )
    }else{
        par0 <- par1 # For remaining data windows, the starting parameters are based on the modelled parameters of the previous window
    }

    years <- 3 # Set default number of years to 3

    invisible(capture.output( # Suppress the details on converging SCEUA
        sceua_list <- rtop::sceua(growth_model,
            par0,
            T_per = T_per,
            G_per = G_per,
            years = years,
            t_int = t_int,
            transfer_function = transfer_function,
            d18Ow = d18Ow,
            Dsam = Dsam,
            Osam = Osam,
            t_maxtemp = t_maxtemp,
            parl,
            paru,
            maxn,
            kstop,
            pcento,
            ngs,
            iniflg = iniflg,
            peps = peps,
            implicit = function(pars){sum(pars[4]/2, pars[6]) < 1} # Make sure that the cumulative GR curve is not located below 0 (if G_av < - G_amp / 2)
            )
        ))

    par1 <- sceua_list[[1]] # Export parameters of final model
    names(par1) <- names(par0)

    result <- growth_model(par1, T_per, G_per, years, t_int, transfer_function, d18Ow, Dsam, Osam, t_maxtemp, plot = FALSE, MC, D_err, O_err, return = "result") # Calculate the end result of the best fit
    
    if(plot == TRUE){
        fitplot <- fitplot + # Add the new model fit to the plot to track progress of the model
            ggplot2::geom_point(data = as.data.frame(result), ggplot2::aes(Dsam, Omod), colour = "red") +
            ggplot2::geom_line(data = as.data.frame(result), ggplot2::aes(Dsam, Omod), colour = "red") 
        print(fitplot)
    }
    
    resultarray[, i + length(dat[1, ]), ] <- rbind(matrix(NA, nrow = i - 1, ncol = 9), result[, 3:11], matrix(NA, nrow = length(dat$D) - i - length(result[,1]) + 1, ncol = 9)) # Add results to result array
    parmat[, i] <- par1 # Add parameters to parameter array
}
# Continued without issues

# Provide names for all dimensions in the result array
dimnames(resultarray) <- list(
    paste("sample", 1:length(resultarray[, 1, 3])),
    c(colnames(dat), paste("window", 1:length(dynwindow$x))),
    c("Modeled_d18O", "d18O_residuals", "Time_of_year", "Instantaneous_growth_rate", "Modeled temperature", "Modeled_d18O_SD", "Time_of_Year_SD", "Instantaneous_growth_rate_SD", "Modeled_temperature_SD")
)

colnames(parmat) <- paste("window", 1:length(parmat[1,]))

# STEP 3: Align model resultis s to cumulative timescale
print("Calculating cumulative day of the year results...")
suppressWarnings(resultarray[, , 3] <- cumulative_day(resultarray, TRUE, TRUE, export_path)) # Calculate cumulative day of the year for all model runs and replace matrix in result array

# STEP 4: Order and export results and statistics
export_results(export_path, dat, resultarray, parmat, MC, dynwindow, plot, plot_export, export_raw) # Export results of model