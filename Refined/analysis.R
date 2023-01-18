# Setting working directory
setwd('D:/SAVRA/Refined/')

# Cleaning up the environment before running the script
rm(list = ls(all = TRUE))
shell("cls")

# Importing Modules
  # Taken from HCI-STATS
library(ez)
library(ggpubr)
library(dplyr)
library(car)
library(stats)
  # Taken from R companion
if(!require(psych)){install.packages("psych")}
if(!require(FSA)){install.packages("FSA")}
if(!require(lattice)){install.packages("lattice")}
if(!require(coin)){install.packages("coin")}
if(!require(PMCMRplus)){install.packages("PMCMRplus")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(DescTools)){install.packages("DescTools")}

# Defining function to run for data analysis

is_normal <- function(df, metric_name) {
  m_stp <- df[df$TPT == 'STP' & df$AA == "Audio Assistance Disabled",][, metric_name]
  m_stpa <- df[df$TPT == 'STP' & df$AA == "Audio Assistance Enabled",][, metric_name]
  m_rsvp <- df[df$TPT == 'RSVP' & df$AA == "Audio Assistance Disabled",][, metric_name]
  m_rsvpa <- df[df$TPT == 'RSVP' & df$AA == "Audio Assistance Disabled",][, metric_name]
  return(shapiro.test(m_stp)$p.value > 0.05 && shapiro.test(m_stpa)$p.value > 0.05 && shapiro.test(m_rsvp)$p.value > 0.05 && shapiro.test(m_rsvpa)$p.value > 0.05)
}

is_sphericity <- function(df, metric_name) {
  return(leveneTest(df[, metric_name] ~ df[, 'TPT'])$`Pr(>F)` > 0.05 && leveneTest(df[, metric_name] ~ df[, 'AA'])$`Pr(>F)` > 0.05)
}

run_test <- function(df, metric_name) {
  # NOTE: Considers two trials from each participants as isolated data points
  # Converting the IV, DV columns to factors
  df$Participant_ID <- as.factor(df$Participant_ID)
  df$TPT <- as.factor(df$TPT)
  df$AA <- as.factor(df$AA)
  
  # Checking if all the ANOVA requirements are met
  # NOTE: Checking for sphericity is not required as each factor does not have more than 2 levels
  
  satisfy_reqs = is_normal(df, metric_name) ## && check_levene(df, metric_name)
  if (satisfy_reqs) {
    print(paste("Running two-way repeated measures ANOVA for:", toString(metric_name)))
    # Run two-way repeated measures ANOVA (two within-subjects test)
    # This is a temporary fix for the fucking function
    df$metric_name = df[, names(df) == metric_name]
    # Running the test
    options(contrasts=c("contr.sum", "contr.poly")) # No clue what this does
    output <- ezANOVA(data=df, dv=(metric_name), wid=.(Participant_ID), within=.(TPT, AA), type=3, detailed=TRUE)
    print(output)
    
    # Checking if any of the p-values are less than 0.05
    is_significant <- FALSE
    p_values <- as.list(output[[1]][2:4,'p']) # Extracting p-values of last 3 rows.
    for (p in p_values) {
      is_significant = is_significant || p < 0.05
    }
    if (is_significant) {
      # Do post-hoc analysis
      # Running Tukey's Test for all the four comparisons
      
      # STP vs RSVP (with audio assistance disabled)
      TPTAADdf <- subset(df[df$AA == 'Audio Assistance Disabled',], select = c('TPT', 'metric_name'))
      print(TukeyHSD(aov(TPTAADdf$metric_name ~ TPTAADdf$TPT, data=TPTAADdf)))
      
      # STP vs RSVP (with audio assistance enabled)
      TPTAAEdf <- subset(df[df$AA == 'Audio Assistance Enabled',], select = c('TPT', 'metric_name'))
      print(TukeyHSD(aov(TPTAAEdf$metric_name ~ TPTAAEdf$TPT, data=TPTAAEdf)))
      
      # STP vs STPA
      AASTPdf <- subset(df[df$TPT == 'STP',], select = c('AA', 'metric_name'))
      print(TukeyHSD(aov(AASTPdf$metric_name ~ AASTPdf$AA, data=AASTPdf)))
      
      # RSVP vs RSVPA
      AARSVPdf <- subset(df[df$TPT == 'RSVP',], select = c('AA', 'metric_name'))
      print(TukeyHSD(aov(AARSVPdf$metric_name ~ AARSVPdf$AA, data=AARSVPdf)))
      
    } # else: Do nothing.
    
    # Now, check if the p-value is less than
  } else {
    print(paste("Running one-factor repeated measures friedman test for:", toString(metric_name)))
    # Running Friedman test
    # Friedman currently does not have an implementation for two-way repeated measures. Hence, running for each factor separately
    
    # Running test for text presentation technique with audio assistance disabled
    TPTAADdf <- cbind(STP=df[df$TPT=='STP' & df$AA=='Audio Assistance Disabled',][, metric_name], RSVP=df[df$TPT=='RSVP' & df$AA=='Audio Assistance Disabled',][, metric_name])
    TPTAADtest <- friedman.test(TPTAADdf)
    print(TPTAADtest)
    
    # Running test for text presentation technique with audio assistance enabled
    TPTAAEdf <- cbind(STP=df[df$TPT=='STP' & df$AA=='Audio Assistance Enabled',][, metric_name], RSVP=df[df$TPT=='RSVP' & df$AA=='Audio Assistance Enabled',][, metric_name])
    TPTAAEtest <- friedman.test(TPTAAEdf)
    print(TPTAAEtest)
    
    # Running test for audio assistance in STP
    AASTPdf <- cbind(AAD=df[df$AA=='Audio Assistance Disabled' & df$TPT=='STP',][, metric_name], AAE=df[df$AA=='Audio Assistance Enabled' & df$TPT=='STP',][, metric_name])
    AASTPtest <- friedman.test(AASTPdf)
    print(AASTPtest)
    
    # Running test for audio assistance in RSVP
    AARSVPdf <- cbind(AAD=df[df$AA=='Audio Assistance Disabled' & df$TPT=='RSVP',][, metric_name], AAE=df[df$AA=='Audio Assistance Enabled' & df$TPT=='RSVP',][, metric_name])
    AARSVPtest <- friedman.test(AARSVPdf)
    print(AARSVPtest)
    
    # Now, check if any of the comparisons are significant
    is_significant = TPTAADtest[[3]] < 0.05 || TPTAAEtest[[3]] < 0.05 || AASTPtest[[3]] < 0.05 || AARSVPtest[[3]] < 0.05
    if (is_significant) {
      print("Found significant difference in Freidman test. Running Wilcoxon Signed Rank Test")
      # Do post-hoc analysis, running the Wilcoxon signed rank test
      
      # Test: RSVP vs STP (with audio assistance disabled)
      TPTAADph <- subset(df[df$AA=='Audio Assistance Disabled',], select = c('TPT', metric_name))
      res <- pairwise.wilcox.test(TPTAADph[, metric_name], TPTAADph$TPT, p.adj = "bonferroni", exact=F, paired=T)
      if(res[[3]] < 0.05) {
        # Calculate the medians
        medianRSVP = median(TPTAADph[TPTAADph$TPT == 'RSVP',][, metric_name])
        medianSTP = median(TPTAADph[TPTAADph$TPT == 'STP',][, metric_name])
        if (medianRSVP < medianSTP) {
          print(paste(metric_name, ": RSVP < STP", res[[3]]))
        } else {
          print(paste(metric_name, ": RSVP > STP", res[[3]]))
        }
      }

      # Test: RSVP vs STP (with audio assistance enabled)
      TPTAAEph <- subset(df[df$AA=='Audio Assistance Enabled',], select = c('TPT', metric_name))
      res <- pairwise.wilcox.test(TPTAAEph[, metric_name], TPTAAEph$TPT, p.adj = "bonferroni", exact=F, paired=T)
      if(res[[3]] < 0.05) {
        # Calculate the medians
        medianRSVP = median(TPTAAEph[TPTAAEph$TPT == 'RSVP',][, metric_name])
        medianSTP = median(TPTAAEph[TPTAAEph$TPT == 'STP',][, metric_name])
        if (medianRSVP < medianSTP) {
          print(paste(metric_name, ": RSVPA < STPA", res[[3]]))
        } else {
          print(paste(metric_name, ": RSVPA > STPA", res[[3]]))
        }
      }

      # Test: RSVP vs RSVPA
      AARSVPph <- subset(df[df$TPT=='RSVP',], select = c('AA', metric_name))
      res <- pairwise.wilcox.test(AARSVPph[, metric_name], AARSVPph$AA, p.adj = "bonferroni", exact=F, paired=T)
      if(res[[3]] < 0.05) {
        # Calculate the medians
        medianRSVP = median(AARSVPph[AARSVPph$AA == 'Audio Assistance Disabled',][, metric_name])
        medianRSVPA = median(AARSVPph[AARSVPph$AA == 'Audio Assistance Enabled',][, metric_name])
        if (medianRSVP < medianRSVPA) {
          print(paste(metric_name, ": RSVP < RSVPA", res[[3]]))
        } else {
          print(paste(metric_name, ": RSVP > RSVPA", res[[3]]))
        }
      }
      
      # Test: STP vs STPA
      AASTPph <- subset(df[df$TPT=='STP',], select = c('AA', metric_name))
      res <- pairwise.wilcox.test(AASTPph[, metric_name], AASTPph$AA, p.adj = "bonferroni", exact=F, paired=T)
      if(res[[3]] < 0.05) {
        # Calculate the medians
        medianSTP = median(AASTPph[AASTPph$AA == 'Audio Assistance Disabled',][, metric_name])
        medianSTPA = median(AASTPph[AASTPph$AA == 'Audio Assistance Enabled',][, metric_name])
        if (medianSTP < medianSTPA) {
          print(paste(metric_name, ": RSVP < RSVPA", res[[3]]))
        } else {
          print(paste(metric_name, ": RSVP > RSVPA", res[[3]]))
        }
      }
    } else {
      print("Found no significant difference in Freidman test.")
    }
  }
}

# Importing all the csv data
# Reaction Times
RT <- read.csv(file = 'ReactionTimes.csv')
# Average Lane Position Offset
ALPO <- read.csv(file = 'ALPO.csv') # Not reported in the findings
# Standard Deviation in Lane Position
SDLP <- read.csv(file = 'SDLP.csv')
# Steering Wheel Reversal Rate
SWRR <- read.csv(file = 'SWRR.csv')
# Maximum Steering Wheel Rotation
MSWR <- read.csv(file = 'MSWR.csv')
# Standard Deviation in Steering Wheel Angle
SDSWA <- read.csv(file = 'SDSWA.csv')
# Maximum Braking Input
MBI <- read.csv(file = 'MBI.csv')
# Average Brake Actuation
ABA <- read.csv(file = 'ABA.csv')

# Run the tests based on their properties.
run_test(SWRR, 'SWRR')