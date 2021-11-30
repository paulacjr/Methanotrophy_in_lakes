########################################################
# Script for producing the figures in Reis et al. 2021 #
#   Role of methanotrophy in the C cycling of lakes    #
########################################################

library(tidyverse)
library(reshape2)
library(cowplot)
library(grid)
library(gridExtra)
library(rLakeAnalyzer)
library(ggrepel)


theme_set(theme_classic()+theme(panel.grid = element_blank()))


#### Load data ####
df <- read.csv("TidyData/Incubation_data.csv")

df2 <- read.csv("TidyData/MOB_DAPI_counts.csv")
  
dfp <- read.csv("TidyData/Profiles_data.csv")


##### Fig 1A: Boxplot MOB and other bacteria (DAPI) cell abundance #####

sub1a <- subset(df, select = c("AlphaMOB_DAPIsubset_cellsmL",
                                "GammaMOB_DAPIsubset_cellsmL",
                                "Total_bacteria_DAPI_cellsmL"))

# Calculate abundance of non-MOB cells (other_bac)
sub1a$other_bac <- sub1a$Total_bacteria_DAPI_cellsmL - sub1a$AlphaMOB_DAPIsubset_cellsmL - sub1a$GammaMOB_DAPIsubset_cellsmL
sub1a <- sub1a[,-3]

msub1a <- melt(sub1a)
msub1a$variable <- as.character(msub1a$variable)
msub1a$variable[msub1a$variable == "AlphaMOB_DAPIsubset_cellsmL"] <- "Alpha-MOB"
msub1a$variable[msub1a$variable == "GammaMOB_DAPIsubset_cellsmL"] <- "Gamma-MOB"
msub1a$variable[msub1a$variable == "other_bac"] <- "DAPI"

aov.out <- aov(value ~ variable, data = na.omit(msub1a))
summary(aov.out) # p < 0.0001
TukeyHSD(aov.out) # Gamma and Alpha-MOB are not different; Other bact are different from both MOB

msub1a$variable <- factor(msub1a$variable, levels = c("DAPI", "Alpha-MOB", "Gamma-MOB"))

(fig1a <- ggplot(data = na.omit(msub1a), aes(x=variable, y=value))+
    geom_boxplot()+ 
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    labs(y=expression(paste("Cell abundance (cells mL "^-1*")")))+
    theme(axis.title.x = element_blank(), 
          axis.text = element_text(size=12, colour = "black"),
          axis.title.y = element_text(size=14), 
          legend.position = "none")+
    annotate("text", x = 1, y = 17000000, label = "a", size = 5)+
    annotate("text", x = 2, y = 110000, label = "b", size = 5)+
    annotate("text", x = 3, y = 1000000, label = "b", size = 5)
)

#ggsave("Figures/fig1a.tiff", fig1a, dpi=300, width=4.1, height=3.6)

# Summary stats
msub1a %>% 
  na.omit() %>% 
  group_by(variable) %>%
  summarise(mean(value), median(value), n = n())


##### Fig 1B: Boxplot MOB and other bacteria cell size (volume) #####

sub1b <- subset(df, select = c("AlphaMOB_mean_cell_area_DAPIsubset_um2",
                                "GammaMOB_mean_cell_area_DAPIsubset_um2",
                                "DAPI_mean_cell_area_um2",
                                "CH4_O2_ratio"))

colnames(sub1b) <- c("Alpha-MOB", "Gamma-MOB","DAPI", "CH4_O2_ratio")
msub1b <- melt(sub1b, id = "CH4_O2_ratio")
msub1b$variable <- factor(msub1b$variable, levels = c("DAPI", "Alpha-MOB", "Gamma-MOB"))

msub1b$cell_radius <- sqrt(msub1b$value/pi) # circle area = pi*r2, so r = sqrt(area/pi) ("value" is the measured area of cells in um2)
msub1b$cell_vol <- 4/3 * pi * (msub1b$cell_radius)^3 # cell volume assuming spherical shape of cells

aov.out <- aov(cell_vol ~ variable, data = na.omit(msub1b))
summary(aov.out) # p < 0.0001
TukeyHSD(aov.out) # Gamma-MOB vol is different from Alpha-MOB and Other bacteria. Alpha and other bacteria are not different.

(fig1b <- ggplot(data=na.omit(msub1b), aes(x=variable, y=cell_vol))+
    geom_boxplot()+
    labs(x="", y=expression(paste("Cell size (" , mu, "m"^"3"*")")))+
    theme(axis.text = element_text(size = 12, colour = 'black'),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_blank(),
          legend.position = "none")+
    annotate("text", x = 1, y = 0.65, label = "a", size = 5)+
    annotate("text", x = 2, y = 0.75, label = "a", size = 5)+
    annotate("text", x = 3, y = 1.9, label = "b", size = 5)
)

#ggsave("Figures/fig1b_vol.tiff", fig1b, dpi = 300, width = 4.1, height = 3.6)

# Summary stats
msub1b$cell_diameter <- msub1b$cell_radius * 2
msub1b %>%
  na.omit() %>%
  group_by(variable) %>%
  summarise(mean(cell_diameter), mean(cell_vol), n=n())

# MOB biomass in incubations
mob_biomass_percent <- ((df$AlphaMOB_DAPIsubset_biomass_ugCL + df$GammaMOB_DAPIsubset_biomass_ugCL) / df$Total_bact_biomass_ugCL) * 100
mean(na.omit(mob_biomass_percent)); median(na.omit(mob_biomass_percent)); max(na.omit(mob_biomass_percent)) # max: 17.32663

alpha_biomass_percent <- (df$AlphaMOB_DAPIsubset_biomass_ugCL / df$Total_bact_biomass_ugCL) * 100
mean(na.omit(alpha_biomass_percent)); max(na.omit(alpha_biomass_percent))

gamma_biomass_percent <- (df$GammaMOB_DAPIsubset_biomass_ugCL / df$Total_bact_biomass_ugCL) * 100
mean(na.omit(gamma_biomass_percent)); max(na.omit(gamma_biomass_percent))

##### Fig 2A: Total C consumption by methanotrophy and heterotrophy by CH4:O2 molar ratio #####

## Calculate the filtration factor (how much the filtration reduced the DAPI biomass) to correct heterotrophic respiration:

# Add column 'Layer' to df2 with lake layer code
df2$Label <- as.character(df2$Label)
t <- strsplit(df2$Label, " ", fixed = T)
for(i in 1:length(t)){
  df2[i, "Layer"] <- t[[i]][1]
}

# Add column 'initial.final' to df2 with initial/final info
for(i in 1:length(t)){
  df2[i, "initial.final"] <- t[[i]][2]
}

subdf2 <- subset(df2, initial.final == 'i',
                select = c(ï..Lake,
                            SampleID, 
                            Layer,
                            Slide,
                            MOB_biomass_ugCL,
                            Total_DAPI_biomass_ugCL))

colnames(subdf2)[1] <- "Lake"
head(subdf2)

# Create column 'Treatment' (filtered/ unfiltered)
s <- strsplit(subdf2$Layer, "_", fixed = T)
for(i in 1:length(s)){
  subdf2[i, "Treatment"] <- s[[i]][2]
}
subdf2$Treatment[subdf2$Treatment == "filt"] <- "Filtered"
subdf2$Treatment[is.na(subdf2$Treatment)] <- "Unfiltered"

# Rename layers
subdf2$Layer <- as.character(subdf2$Layer)
subdf2$Layer[subdf2$Layer == "Epi" | subdf2$Layer == "Epi_filt"] <- "Epilimnion"
subdf2$Layer[subdf2$Layer == "Meta" | subdf2$Layer == "Meta_filt"] <- "Metalimnion"
subdf2$Layer[subdf2$Layer == "Hypo" | subdf2$Layer == "Hypo_filt"] <- "Oxic hypolimnion"
head(subdf2); tail(subdf2)

# Create replicates ID
r <- strsplit(subdf2$SampleID, "-", fixed = T)
for(i in 1:length(r)){
  subdf2[i, "ID"] <- r[[i]][2]
}

# Sum MOB biomass (gamma + alpha) and subtract from total DAPI
subdf2$MOB_biomass_ugCL <- ifelse(is.na(subdf2$MOB_biomass_ugCL), 0, subdf2$MOB_biomass_ugCL)
subdf2 <- subdf2 %>% 
  group_by(ID) %>% 
  mutate(total_MOB_biomass = sum(MOB_biomass_ugCL)) %>% 
  mutate(DAPI_biomass2 = Total_DAPI_biomass_ugCL - total_MOB_biomass)

subdf2 <- subset(subdf2, select = c(Lake, Layer, Treatment, ID, DAPI_biomass2, total_MOB_biomass))

mean <- subdf2 %>% 
  group_by(Lake, Layer, Treatment) %>% 
  summarise(mean_DAPI_biomass2 = mean(DAPI_biomass2))
gm <- data.frame("Geai", "Metalimnion", "Filtered", NA)
names(gm) <- c("Lake","Layer","Treatment","mean_DAPI_biomass2")
mean <- rbind(mean,gm)
goh <- data.frame("Geai", "Oxic hypolimnion", "Filtered", NA)
names(goh) <- c("Lake","Layer","Treatment","mean_DAPI_biomass2")
mean <- rbind(mean,goh)

mean <- mean %>%
  group_by(Lake, Layer) %>% 
  mutate(ratio = mean_DAPI_biomass2[Treatment == "Unfiltered"] / mean_DAPI_biomass2[Treatment == "Filtered"])

mean$ratio[is.na(mean$ratio)] <- mean$ratio[mean$Lake == "Geai"][1] # use Geai's epi ratio for meta and hypo

## Multiply heterotrophic respiration by the filtration factor (dataframe 'mean' produced above):

mean$LakeLayer <- paste0(mean$Lake, mean$Layer)
mean <- mean[,-c(1:4)] # keep only the cols we need
mean <- distinct(mean) # remove duplicate lines
df$LakeLayer <- paste0(df$Lake, df$Layer)
df <- left_join(df, mean, by = "LakeLayer")
df$Total_C_cons_HB_ugCLd <- df$BP_ugCLd + ((df$BR_rate_ugCLd*-1) * df$ratio)

sub2a <- subset(df, select = c(Lake,
                               Layer,
                               CH4_cons_rate_ugCLd,
                               Total_C_cons_HB_ugCLd,
                               CH4_O2_ratio,
                               DO_uM,
                               log10_CH4_uM))

sub2a_ll <- sub2a %>%
  na.omit() %>%
  group_by(Lake, Layer)
sub2a_ll <- sub2a_ll[,-c(3,5:7)] # for Fig 3 below

names(sub2a)[3] <-  "Methanotrophy"
names(sub2a)[4] <-  "Heterotrophy"

msub2a <- melt(sub2a, id=c("Lake", "Layer", "CH4_O2_ratio", "DO_uM", "log10_CH4_uM"))
head(msub2a); tail(msub2a)

msub2a$variable <- factor(msub2a$variable, c("Heterotrophy", "Methanotrophy"))

(fig2a <- ggplot(data=na.omit(msub2a), aes(x=CH4_O2_ratio, y=value))+
    geom_point(aes(color=variable))+
    scale_color_manual(values = c("#CC6677","#44AA99"))+
    geom_smooth(method='loess', se = T, aes(color=variable, fill=variable), alpha = 0.2)+
    scale_fill_manual(values = c("#CC6677","#44AA99"))+
    scale_y_log10(breaks = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))+
    scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = c(0.001, 0.01, 0.1, 1, 10, 100))+
    labs(y=expression(paste("Total C consumption (", mu,"gC L"^"-1"*" d"^"-1"*")")),
         x=expression(paste("CH"[4], ":O"[2], " molar ratio")))+
    theme(axis.text.y = element_text(size = 10, colour = 'black'), 
          axis.text.x = element_text(size = 9, colour = 'black'),
          axis.title = element_text(size = 11),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.position = c(0.75, 0.2))
)

# Summary stats
msub2a %>% 
  na.omit() %>% 
  group_by(variable) %>% 
  summarise(mean(value), median(value), min(value), max(value), n=n()) 
sub2a %>% 
  mutate(ratio = Methanotrophy/Heterotrophy) %>% 
  summarise(max(na.omit(ratio)))
sub2a %>% 
  summarise(max(CH4_O2_ratio)) 
sub2a %>% 
  summarise(min(na.omit(Methanotrophy)), max(na.omit(Methanotrophy)), mean(na.omit(Methanotrophy)))
sub2a %>% 
  summarise(min(na.omit(Heterotrophy)), max(na.omit(Heterotrophy)), mean(na.omit(Heterotrophy)))

# Get standard error of Methanotrophy loess curve
list <- ggplot_build(fig2a)
fit <- ggplot_build(fig2a)$data[[2]]
fit_mob <- subset(fit, fit$colour == "#44AA99")
fit_mob$x_notlog <- 10^(fit_mob$x); fit_mob$x_notlog <- round(fit_mob$x_notlog, digits = 2) # x is CH4:O2 ratio
fit_mob$se[fit_mob$x_notlog == 0.20] # 0.09759951 = 0.1


##### Fig 2B: Ratio MCC:HCC by CH4:O2 molar ratio #####

# Get predicted values from loess model of Fig 2A
# define the models
msub2a <- na.omit(msub2a)
m <- subset(msub2a, msub2a$variable == "Methanotrophy")
h <- subset(msub2a, msub2a$variable == "Heterotrophy")
model_m <- loess(log10(value) ~ log10(CH4_O2_ratio), data = m)
model_h <- loess(log10(value) ~ log10(CH4_O2_ratio), data = h)

# Predict Met and Het values for a CH4:O2 range between 0.0001 and 100 and include the standard error in the output
new_ch4o2ratio <- data.frame(CH4_O2_ratio = seq(0.0001, 100, by=0.0001))
ph <- predict(model_h, newdata = new_ch4o2ratio, allow.new.levels = T, re.form = NA)
ph <- data.frame(ph)
pm <- predict(model_m, newdata = new_ch4o2ratio, allow.new.levels = T, re.form = NA)
pm <- data.frame(pm)

# Prepare data to plot
plotdata <- data.frame(cbind(pred_h = ph$ph,
                             pred_m = pm$pm,
                             new_ch4o2ratio = new_ch4o2ratio$CH4_O2_ratio)
)

plotdata$ratio_MH <- (10^(plotdata$pred_m)) / (10^(plotdata$pred_h))

fig2b <- ggplot(data=na.omit(plotdata), aes(x=new_ch4o2ratio,y=ratio_MH)) +
  geom_point(size = 0.5) +
  scale_x_log10(limits = c(0.001, 100),
                breaks = c(0.001, 0.01, 0.1, 1, 10, 100), 
                labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
  scale_y_log10(limits = c(0.01, 10),
                breaks = c(0.01, 0.1, 1, 10, 100), 
                labels = c(0.01, 0.1, 1, 10, 100)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.6, linetype = "dashed") +
  labs(y="Ratio MCC:HCC",
       x=expression(paste("CH"[4], ":O"[2], " molar ratio"))) +
  theme(axis.text.y = element_text(size = 10, colour = 'black'),
        axis.text.x = element_text(size = 9, colour = 'black', angle = 30),
        axis.title = element_text(size = 11))


##### Fig. 2C: Boxplots C consumption by MOB and HB biomass #####

# Calculate C consumption per MOB and HB biomass (ugC/ugC/day)
df$C_cons_rate_per_MOB_biomass_ugCugCd <- df$CH4_cons_rate_ugCLd/(df$AlphaMOB_DAPIsubset_biomass_ugCL + df$GammaMOB_DAPIsubset_biomass_ugCL) # ugC per ugC biomass per day
df$C_cons_rate_per_HB_biomass_ugCugCd <- df$Total_C_cons_HB_ugCLd/(df$Total_bact_biomass_ugCL-(df$AlphaMOB_DAPIsubset_biomass_ugCL + df$GammaMOB_DAPIsubset_biomass_ugCL)) # ugC per ugC biomass per day

sub2c <- subset(df, select = c("C_cons_rate_per_MOB_biomass_ugCugCd", "C_cons_rate_per_HB_biomass_ugCugCd"))
colnames(sub2c) <- c("MOB", "HP")
msub2c <- melt(sub2c)

msub2c$variable <- factor(msub2c$variable, c("HP", "MOB"))

# Welch Two Sample t-test
t.test.out <- t.test(value ~ variable, data = na.omit(msub2c), alternative="two.sided")
t.test.out # Significative difference

(fig2c <- ggplot(na.omit(msub2c), aes(x=variable, y=value))+
    geom_boxplot(aes(fill=variable))+
    scale_fill_manual(values = c("#CC6677","#44AA99"))+
    scale_y_log10(breaks = c(0.5,5,50), labels = c(0.5,5,50))+
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10, colour = 'black'),
          axis.text.x = element_text(size = 9, colour = 'black'),
          axis.title = element_text(size = 11),
          legend.position = "none")+
    labs(y = expression(paste("C cons. per biomass ( d"^-1*")")))+
    annotate("text", x = 1, y = 5.4, label = "a", size = 4)+
    annotate("text", x = 2, y = 60, label = "b", size = 4)
)

# Summary stats
msub2c %>% 
  na.omit() %>% 
  group_by(variable) %>% 
  summarise(mean(value), median(value), n=n())

##### Combine fig2a, fig2b and fig2c #####
fig2 <- plot_grid(fig2a, fig2b, fig2c, labels = c('a', 'b', 'c'), # or fig2a
                  label_size = 12, ncol = 3, align = "h", rel_widths = c(1.1,0.7,0.4))

#ggsave("Figures/fig2_finalRevision.pdf", fig2, dpi = 300, width = 7, height = 3.3) #7.5, 3.2

##### Suppl Figure 3: Change in C consumption per MOB biomass by CH4 concentration #####
df$total_MOB_biomass <- df$AlphaMOB_DAPIsubset_biomass_ugCL + df$GammaMOB_DAPIsubset_biomass_ugCL

df_nonas <- df[!is.na(df$C_cons_rate_per_MOB_biomass_ugCugCd),]

(figS3 <- ggplot(df_nonas, aes(x=CH4_uM, y=C_cons_rate_per_MOB_biomass_ugCugCd))+
    geom_smooth(method = "loess", se = T, alpha = 0.2, color = "grey60", size = 0.5)+
    geom_point(aes(color = DO_uM))+ 
    scale_size_continuous(breaks = c(0.5, 5, 50))+
    scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = c(0.001, 0.01, 0.1, 1, 10, 100))+
    scale_y_log10(breaks = c(0.05, 0.5, 5, 50), labels = c(0.05, 0.5, 5, 50))+
    labs(x = expression(paste("CH"[4], " (", mu, "M)")), 
         y = expression(paste("C cons. per MOB biomass (", mu, "gC ", mu, "gC"^-1*" d"^-1*")")),
         color = expression(paste("O"[2], " (", mu, "M)" )),
         size = expression(paste("MOB (", mu, "gC L"^-1*")")))+
    theme(axis.text = element_text(size = 10, colour = 'black'),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 10), legend.text = element_text(size = 10))
)

# CH4 and O2 range
df %>% 
  summarise(min(CH4_uM), max(CH4_uM), min(DO_uM), max(DO_uM))


###### Fig 3: Boxplots whole-lake proportion of C consumption by Het and Met per lake #####

# Calculate methanotrophy with model from Thottathil, Reis and Prairie 2019 (DOI: 10.1007/s10533-019-00552-x)
dfp <- subset(dfp, dfp$O2_uM_cal > 0)
dfp$lnch4oxid <- 20.08 + (0.79*log(dfp$CH4_uM)) + (-5669.61/(dfp$Temp_C+273.15)) + log(exp(-0.01*dfp$O2_uM_cal) - exp(-(0.01+0.18)*dfp$O2_uM_cal)) # ln umol/L/d
dfp$ch4oxid <- exp(dfp$lnch4oxid) # umol/L/d
dfp$ch4oxid <- dfp$ch4oxid * 12 # ugC/L/d which is equal to mgC/m3/d
min(dfp$ch4oxid, na.rm = T); max(dfp$ch4oxid, na.rm = T); median(dfp$ch4oxid, na.rm = T); hist(log(dfp$ch4oxid))

# Determine Epi-, Meta- and Hypolimnion layers of profiles

# Split dataframe of profiles by LakeDate
dfp <- dfp %>% unite(LakeDate, Lake, Date, sep = "/")
dfp_split <- split(dfp, dfp$LakeDate) 

# Calculate depth of metalimnia top and bottom (meta_top and meta_bottom will be NA in mixed profiles)

meta <- data.frame(LakeDate = NA, meta_top = NA, meta_bottom = NA)

for (i in 1:length(dfp_split)) {
  meta[i, "LakeDate"] <- names(dfp_split[i])
  meta[i, "meta_top"] <- meta.depths(dfp_split[[i]]$Temp_C, dfp_split[[i]]$Depth_m, slope = 0.1, seasonal = F)[1]
  meta[i, "meta_bottom"] <- meta.depths(dfp_split[[i]]$Temp_C, dfp_split[[i]]$Depth_m, slope = 0.1, seasonal = F)[2]
}

dfp <- dfp %>% left_join(meta, by = "LakeDate")

# Create epi, meta and hypo classification (column Layer)
dfp <- dfp %>% 
  mutate(Layer = ifelse(Depth_m < meta_top, paste("Epilimnion"), NA))
dfp <- dfp %>% 
  mutate(Layer = ifelse(Depth_m < meta_bottom & Depth_m > meta_top, paste("Metalimnion"), Layer))
dfp <- dfp %>% 
  mutate(Layer = ifelse(Depth_m > meta_bottom, paste("Oxic hypolimnion"), Layer))

# Join profiles data and heterotrophic C consumption (dfp and sub2a_ll)
dfp <- dfp %>% 
  separate(LakeDate, into = c("Lake", "Date"), sep = "/")

# Define manually layers for Triton and en Coeur so that we have the respective HB metabolism aligned in the joined table
dfp$Layer[dfp$Lake == "Triton" & dfp$Depth_m > 1.0] <- "Oxic hypolimnion"
dfp$Layer[dfp$Lake == "Triton" & dfp$Depth_m <= 1.0] <- "Epilimnion"

dfp$Layer[dfp$Lake == "en Coeur" & dfp$Depth_m > 4.5] <- "Oxic hypolimnion"
dfp$Layer[dfp$Lake == "en Coeur" & dfp$Depth_m <= 4.5] <- "Epilimnion"

dfp <- dfp %>% 
  left_join(sub2a_ll, by = c("Lake", "Layer"))

# Calculate total C processing by volume of water column layer (volume weighted rates)

# Methanotrophy
dfp$ch4oxid_vol <- dfp$ch4oxid * dfp$vol_m3 # unit of ch4oxid_vol: mgC/d

dfp$Date <- as.character(dfp$Date)
dfp$Lake <- as.character(dfp$Lake)

dfp_nonas <- dfp[!is.na(dfp$ch4oxid_vol),]
unique((dfp_nonas$Date[dfp_nonas$Lake == "Geai"]))

dfp2 <- dfp_nonas %>% 
  group_by(Lake, Date) %>%
  mutate(ch4oxid_vol_sum = sum(ch4oxid_vol))
unique((dfp2$Date[dfp2$Lake == "Geai"]))

# Heterotrophy
dfp2$heterot_vol <- dfp2$Total_C_cons_HB_ugCLd * dfp2$vol_m3 # unit of heterot_vol: mgC/d (ugC/L/d = mgC/m3/d)
dfp2 <- dfp2 %>%
  group_by(Lake, Date) %>% 
  mutate(heterot_vol_sum = sum(heterot_vol))
dfp2$totalC <- dfp2$ch4oxid_vol_sum + dfp2$heterot_vol_sum
dfp2$percent_ch4oxid <- dfp2$ch4oxid_vol_sum/dfp2$totalC * 100
dfp2$percent_heterot <- dfp2$heterot_vol_sum/dfp2$totalC * 100
dfp2$percent_ch4oxid[dfp2$Lake == "Geai"]

# Plot
sub3 <- subset(dfp2, select = c(Lake, Date, percent_ch4oxid, percent_heterot))
colnames(sub3) <- c("Lake", "Date", "Met.", "Het.")

msub3 <- melt(sub3, id = c("Lake", "Date"))

msub3$Lake <- factor(msub3$Lake, c("Triton", "en Coeur", "Morency", "Croche", "Cromwell", "Geai"))
msub3$variable <- factor(msub3$variable, c("Het.", "Met."))

(fig3 <- ggplot(na.omit(msub3), aes(x = variable, y = value))+
    facet_grid(cols = vars(Lake))+
    geom_boxplot(aes(fill=variable))+ 
    geom_point(aes(fill=variable), shape=21, size=1)+
    scale_fill_manual(values = c("#CC6677","#44AA99"))+
    scale_color_manual(values = c("#CC6677","#44AA99"))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 8),
          axis.text = element_text(size = 10, color = "black"),
          legend.position = "none")+
    labs(y="Prop. of C consumed (%)")
)

#ggsave("Figures/fig3.tiff", fig3, dpi=300, width = 6.6, height = 2.6)

# Summary stats
max(na.omit(sub3$Met.)); min(na.omit(sub3$Met.))
max(na.omit(sub3$Met.[sub3$Lake == "en Coeur"]))
max(na.omit(sub3$Met.[sub3$Lake == "Triton"]))

##### Fig 4A: Linear relationship between % Methanotrophy and DOC ####

msub4 <- subset(msub3, variable == "Met.")

msub4$doc_mgL <- NA
msub4$doc_mgL[msub4$Lake == "en Coeur"] <- 4.2
msub4$doc_mgL[msub4$Lake == "Morency"] <- 3.4
msub4$doc_mgL[msub4$Lake == "Croche"] <- 4.1
msub4$doc_mgL[msub4$Lake == "Cromwell"] <- 5.2
msub4$doc_mgL[msub4$Lake == "Geai"] <- 8.2
msub4$doc_mgL[msub4$Lake == "Triton"] <- 4.6

msub4 <- msub4 %>%
  na.omit() %>%
  group_by(Lake) %>%
  mutate(mean_value = mean(value)) # Mean % methanotrophy in each lake during summer

msub4 <- distinct(msub4)

msub4$strat <- ifelse(msub4$Lake == "Triton" | msub4$Lake == "en Coeur", paste("no hypolimnion"), paste("with hypolimnion"))
msub4$strat <- factor(msub4$strat, c("with hypolimnion", "no hypolimnion"))

lm_doc <- lm(value ~ doc_mgL, data = subset(msub4, Lake != "Triton" & Lake != "en Coeur"))
summary(lm_doc) 

msub4$LakeLabel <- NA
msub4$LakeLabel[1] <- "Triton"
msub4$LakeLabel[7] <- "en Coeur"
msub4$LakeLabel[12] <- "Morency"
msub4$LakeLabel[18] <- "Geai"
msub4$LakeLabel[31] <- "Cromwell"
msub4$LakeLabel[35] <- "Croche"

(fig4a <- ggplot(msub4, aes(x=doc_mgL, y=mean_value))+
    geom_point(aes(shape=strat), size = 2)+
    scale_shape_manual(values=c(19, 1))+
    theme(axis.text = element_text(size = 10, color="black"),
          axis.title = element_text(size = 10),
          legend.title = element_blank(),
          legend.position = c(0.65, 0.25),
          legend.text = element_text(size = 12))+
    geom_smooth(data=subset(msub4, Lake != "Triton" & Lake != "en Coeur"), aes(x=doc_mgL, y=value), method = "lm", color = "grey30", size=0.5, se=F)+
    labs(x=expression(paste("DOC (mg L"^-1*")")), y=" ")+
    annotate("text", x=4.6, y=79, label = expression(paste("R"^2*"=0.60, p<0.0001")), size = 4)+ # change R2 and p accordingly
    geom_text_repel(aes(label = LakeLabel), size = 3)
)


##### Fig 4B: Linear relationship between % Methanotrophy and Kd PAR ####

msub4$kd_par <- NA
msub4$kd_par[msub4$Lake == "en Coeur"] <- 0.83
msub4$kd_par[msub4$Lake == "Morency"] <- 0.5
msub4$kd_par[msub4$Lake == "Croche"] <- 0.80
msub4$kd_par[msub4$Lake == "Cromwell"] <- 1.29
msub4$kd_par[msub4$Lake == "Geai"] <- 1.73
msub4$kd_par[msub4$Lake == "Triton"] <- 0.88

lm_kd <- lm(value ~ kd_par, data = subset(msub4, Lake != "Triton" & Lake != "en Coeur"))
summary(lm_kd)

(fig4b <- ggplot(msub4, aes(x=kd_par, y=mean_value))+
    geom_point(aes(shape=strat), size = 2)+
    scale_shape_manual(values=c(19, 1))+
    theme(axis.text = element_text(size = 10, color="black"),
          axis.title = element_text(size = 10),
          legend.position = 'none',
          legend.title = element_blank())+ 
    geom_smooth(data=subset(msub4, Lake != "Triton" & Lake != "en Coeur"), aes(x=kd_par, y=value), method = "lm", color="grey30", size = 0.5, se=F)+
    labs(y=" ", x=expression(paste("Kd PAR (m"^-1*")")))+
    annotate("text", x=0.8, y=75, label = expression(paste("R"^2*"=0.62, p<0.0001")), size = 4)+ # change R2 and p accordingly
    geom_text_repel(aes(label = LakeLabel), size = 3) 
)


##### Combine fig4a and fig4b #####
(fig4_grid <- plot_grid(fig4a + theme(legend.position="none"), 
                        fig4b, 
                        labels = c('a', 'b'), 
                        label_size = 12, 
                        align = "h"))

# Make common x axis label (x.grob)
y.grob <- textGrob("% of C cons. by methanotrophy", rot = 90, vjust = 1)

# Put graphs and x.grob together
fig4_grid_g <- grid.arrange(arrangeGrob(fig4_grid, left = y.grob))

# Make common legend
legend <- get_legend(fig4a +
                       guides(color = guide_legend(nrow=1)) +
                       theme(legend.position = "top"))

# Put all together
(fig4_grid_final <- plot_grid(legend, fig4_grid_g, ncol = 1, rel_heights = c(.1,1)))
#ggsave("Figures/fig4_grid.pdf", fig4_grid_final, dpi=300, width = 7.75, height = 3.5)
#ggsave("Figures/fig4_finalRevision.pdf", fig4_grid_final, dpi=300, width = 7, height = 3.4)


##### Suppl Figure 5: Relationships between DOC, cDOM, and kdPAR #####
msub4$cdom <- NA
msub4$cdom[msub4$Lake == "Morency"] <- 0.69
msub4$cdom[msub4$Lake == "Croche"] <- 0.92
msub4$cdom[msub4$Lake == "Cromwell"] <- 1.38
msub4$cdom[msub4$Lake == "Geai"] <- 4.26
msub4$cdom[msub4$Lake == "en Coeur"] <- 0.46
msub4$cdom[msub4$Lake == "Triton"] <- 1.04

lm_cdomdoc <- lm(doc_mgL ~ cdom, data = msub4)
summary(lm_cdomdoc)

(figS5a <- ggplot(msub4, aes(x=doc_mgL, y=cdom))+
    geom_point(aes(shape=strat), size = 2)+
    theme(axis.text = element_text(size = 10, color="black"),
          axis.title = element_text(size = 10),
          legend.position = 'none',
          legend.title = element_blank())+ 
    geom_smooth(data=msub4, aes(x=doc_mgL, y=cdom), method = "lm", color="grey30", size = 0.5)+
    labs(y=expression(paste("cDOM (m"^-1*")")), x=expression(paste("DOC (mg L"^-1*")")))+
    annotate("text", x=4.5, y=4.5, label = expression(paste("R"^2*"=0.97, p<0.0001")), size = 4)+ 
    geom_text(aes(label = Lake), size = 2.5, hjust=0.3, vjust=-0.5)
)

lm_cdomkd <- lm(kd_par ~ cdom, data = msub4)
summary(lm_cdomkd)

(figS5b <- ggplot(msub4, aes(x=cdom, y=kd_par))+
    geom_point(aes(shape=strat), size = 2)+
    theme(axis.text = element_text(size = 10, color="black"),
          axis.title = element_text(size = 10),
          legend.position = 'none',
          legend.title = element_blank())+ 
    geom_smooth(data=msub4, aes(x=cdom, y=kd_par), method = "lm", color="grey30", size = 0.5)+
    labs(y=expression(paste("Kd PAR (m"^-1*")")), x=expression(paste("cDOM (m"^-1*")")))+
    annotate("text", x=1.4, y=1.75, label = expression(paste("R"^2*"=0.88, p<0.0001")), size = 4)+
    geom_text(aes(label = Lake), size = 2.5, hjust=0.3, vjust=-0.5)
)

##### Combine figsS5A and figsS5B #####
figS5_grid <- plot_grid(figS5a, figS5b, labels = c('A', 'B'), label_size = 12, ncol = 2, align = "h")
figS5_grid <- plot_grid(legend, figS5_grid, ncol = 1, rel_heights = c(.1,1))


