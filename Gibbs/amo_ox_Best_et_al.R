# Script for MB to calculate energetics of ammonia oxidation at Frasassi, provided by HSA, modified by MB and DSJ
# Script modified from https://www.caitlincasar.com/post/thermodynamic_modeling/
# https://github.com/CaitlinCasar/Casar2020_DeMMO_MineralHostedBiofilms
# https://onlinelibrary.wiley.com/doi/full/10.1111/gbi.12391?casa_token=rUV7ebksKqMAAAAA%3A7-KRCT5Dv_T91jI3Kg0UED90NenzSLNvbWHMfSt-CG3V_qRfHBlQNrSNm33BvAutc_zRR_ViehLCGxci3A

# Install packages
#install.packages("CHNOSZ")
#install.packages("tidyverse")
library(CHNOSZ)
library(readxl)
library(tidyverse)

# set working directory
#setwd()

# Set units to joules 
E.units('J')
# Set temperature units to C
T.units('C')

# Set parameters
# Set temperature, pH, RT ranges
temperature_range = seq(0,100,by=10)
pH_range = seq(0,13,by=1)

# Import activities
# Here, activities are listed as log10(activity). I set these values based on approximate in situ concentrations in the cave air
# Log activities are used because that is the output of speciation calculations (ASU Worm portal - let's chat about this)
# This spreadsheet can be modified to add different sites and other reactants/products
# Water and minerals have activities of 1 don't affect the Q term (generally) so they are not listed in the activities sheet
log_activities = read_xlsx('logactivities.xlsx')%>%
  uncount(length(pH_range)) # replicate activities for length of pH units

# Add columns for pH range
log_activities = log_activities %>%
  add_column(`H+` = -(rep(pH_range, nrow(log_activities)/14)))%>%
  add_column(pH = rep(pH_range, nrow(log_activities)/14)) 

# Convert to activity and pivot to long format
activities = log_activities %>%
  pivot_longer(`NH3`:`H+`, names_to = 'react_prod',values_to = 'activity') %>% # pivot from wide to long for joining
  mutate(activity = 10^activity)

# Import reactions
reactions = read_xlsx('reactions.xlsx')

rxn_text = reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c('.value','set'),
               names_pattern = '(.+).(.+)') %>% # pivot from wide to long
  unite('react_prod', reactant:product,na.rm = TRUE, remove = F)%>% # unite the reactants and products into one columns
  filter(!react_prod == '')%>% # remove any rows with missing react_prod values
  group_by(rxn.name) %>%
  reframe(subcrt(react_prod, coeff, state)$reaction)%>%
  group_by(rxn.name) %>%
  nest()

rxn_text = rxn_text %>%
  add_column(text = map(rxn_text$data, describe.reaction))

reactions = reactions %>%
  left_join(rxn_text)

# Calculate ∆Gº (kJ/mol)
G0 = reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c('.value','set'),
               names_pattern = '(.+).(.+)')%>% # pivot from wide to long dataframe
  unite('react_prod', reactant:product,na.rm = TRUE, remove = F)%>% # unite the reactants and products into one columns
  filter(!react_prod == '')%>% # remove any rows with missing react_prod values
  group_by(rxn.name)%>% # group by reaction name for calculations
  reframe(G0 = subcrt(react_prod, coeff, state, T = temperature_range)$out$G/1000)%>% # Calculate ∆G0 with the subcrt function from CHNOSZ
  add_column(temperature = rep(temperature_range, nrow(reactions)))%>% # Add tempreature column
  add_column(RT = rep((temperature_range + 273.15)*.008314, nrow(reactions)))%>% # Add RT value (gas constant * temperature)
  expand_grid(pH = pH_range) # Add pH values

# Calculate logQ
logQ = reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c('.value','set'),
               names_pattern = '(.+).(.+)')%>% # pivot from wide to long
  unite('react_prod', reactant:product, na.rm = TRUE, remove = F)%>% #unite the reactant and product columns into one column called react_prod
  left_join(activities, multiple = 'all')%>% # join with activities data
  filter(!is.na(activity))%>% # remove any activities with NA value
  mutate(Q = if_else(!is.na(reactant), activity^coeff, activity^coeff))%>%
  group_by(rxn.name, Site, pH)%>% #group on the reaction number and site 
  summarise(logQ = log(prod(Q))) #calculate logQ 

# Calculate ∆G!
deltaG <- G0 %>%
  left_join(logQ, multiple = 'all') %>% #join the G0 and logQ tables 
  left_join(reactions %>% select(rxn.name, e.transfer, text)) %>% #add the reaction number, number of electrons transferred, and minerals from each reaction 
  mutate(deltaG = G0 + RT*logQ) %>% #calculate deltaG for each reaction at each site 
  mutate(deltaG.e.transfer = deltaG/e.transfer)

#Loop through sites and reactions and plot
reaction_names = reactions$rxn.name
sites = unique(log_activities$Site)

# calculate pK of NH3/NH4+
pK_NH3 = subcrt(c("NH3", "H+", "NH4+"), c(-1, -1, 1), T = temperature_range)$out$logK

# Axis labels
DGlab <- axis.label("DGr", prefix = "k")
DG.per.e = expression(paste(Delta,'G'[r],' kJ/mol e'^'-'))

# Plot! If you want separate png files, use this
for(reaction in reaction_names){
  for(site in sites){
    png(file=paste0(reaction, ".png"), height=6, width=8,units='in',res=300)
    data = deltaG %>%
      filter(rxn.name == reaction & Site == site)%>%
      select(pH, temperature, deltaG.e.transfer) %>%
      pivot_wider(names_from = pH, values_from = deltaG.e.transfer) %>%
      column_to_rownames(var = "temperature") %>%
      as.matrix()
    rxn_text = unlist(deltaG %>%
                        filter(rxn.name == reaction & Site == site)%>%
                        select(text)%>%
                        distinct())
    filled.contour(temperature_range,pH_range,data, xlab = axis.label("T"),
                   color.palette = ifelse(getRversion() >= "3.6.0", function(n) hcl.colors(n), topo.colors),
                   #plot.title = {par(cex.main=1);title(main = reaction)},
                   key.title = {par(cex.main=.8);title(main=DG.per.e)},
                   plot.axes = {
                     legend("topleft", legend = rxn_text$text, bty = "n", inset = c(0,.85))
                     lines(temperature_range, pK_NH3, lty = 2)
                     text(85, 8.6, expr.species("NH3"))
                     text(85, 6.9, expr.species("NH4+"))
                     axis(1)
                     axis(2)
                     title(ylab="pH")})
    dev.off()
  }
}

# Plot! If you want a single pdf, use this
pdf(file = 'pH_temp_gradient_plots.pdf', height = 8, width = 11, onefile = TRUE)
for(reaction in reaction_names){
  for(site in sites){
    data = deltaG %>%
      filter(rxn.name == reaction & Site == site)%>%
      select(pH, temperature, deltaG.e.transfer) %>%
      pivot_wider(names_from = pH, values_from = deltaG.e.transfer) %>%
      column_to_rownames(var = "temperature") %>%
      as.matrix()
    rxn_text = unlist(deltaG %>%
                        filter(rxn.name == reaction & Site == site)%>%
                        select(text)%>%
                        distinct())
    filled.contour(temperature_range,pH_range, data, xlab = axis.label("T"),
                   color.palette = ifelse(getRversion() >= "3.6.0", function(n) hcl.colors(n), topo.colors),
                   #plot.title = {par(cex.main=1);title(main = reaction)},
                   key.title = {par(cex.main=.8);title(main=DG.per.e)},
                   plot.axes = {
                     legend("topleft", legend = rxn_text$text, bty = "n", inset = c(0,.85))
                     lines(temperature_range, pK_NH3, lty = 2)
                     text(85, 8.6, expr.species("NH3"))
                     text(85, 6.9, expr.species("NH4+"))
                     axis(1)
                     axis(2)
                     title(ylab="pH")})
  }
}
dev.off()

#### Ammonia concentration ####
NH3 = seq(-12,-1)
O2 = seq(-12,-1)

nh3_o2 = expand.grid(NH3 = seq(-12,-1,by=1), O2 = seq(-12,-1,by=1))
gradient = seq(1:144)

log_activities = read_xlsx('logactivities.xlsx')%>%
  uncount(nrow(nh3_o2))%>%
  mutate(NH3 = nh3_o2$NH3)%>%
  mutate(`NH4+` = nh3_o2$NH3)%>%
  mutate(O2 = nh3_o2$O2)%>%
  add_column(`H+` = -7)%>%
  add_column(grad = gradient)

# Convert to activity and pivot to long format
activities = log_activities %>%
  pivot_longer(`NH3`:`H+`, names_to = 'react_prod',values_to = 'activity') %>% # pivot from wide to long for joining
  mutate(activity = 10^activity)

# Calculate ∆Gº (kJ/mol) at 25ºC
G0 = reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c('.value','set'),
               names_pattern = '(.+).(.+)')%>% # pivot from wide to long dataframe
  unite('react_prod', reactant:product,na.rm = TRUE, remove = F)%>% # unite the reactants and products into one columns
  filter(!react_prod == '')%>% # remove any rows with missing react_prod values
  group_by(rxn.name)%>% # group by reaction name for calculations
  reframe(G0 = subcrt(react_prod, coeff, state, T = 25)$out$G/1000)%>% # Calculate ∆G0 with the subcrt function from CHNOSZ at 25ºC
  add_column(RT = rep((25 + 273.15)*.008314, nrow(reactions))) # Add RT value (gas constant * temperature)

# Calculate logQ
logQ = reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c('.value','set'),
               names_pattern = '(.+).(.+)')%>% # pivot from wide to long
  unite('react_prod', reactant:product, na.rm = TRUE, remove = F)%>% #unite the reactant and product columns into one column called react_prod
  left_join(activities, multiple = 'all')%>% # join with activities data
  filter(!is.na(activity))%>% # remove any activities with NA value
  mutate(Q = if_else(!is.na(reactant), activity^coeff, activity^coeff))%>%
  group_by(rxn.name, Site, grad)%>% #group on the reaction number and site 
  summarise(logQ = log(prod(Q))) #calculate logQ 

# Calculate ∆G!
deltaG <- G0 %>%
  left_join(logQ, multiple = 'all') %>% #join the G0 and logQ tables 
  left_join(reactions %>% select(rxn.name, e.transfer, text)) %>% #add the reaction number, number of electrons transferred, and minerals from each reaction 
  left_join(log_activities)%>%
  mutate(NH3 = 10^NH3)%>%
  mutate(NH3 = 10^`NH4+`)%>%
  mutate(O2 = 10^O2) %>%
  mutate(deltaG = G0 + RT*logQ) %>% #calculate deltaG for each reaction at each site 
  mutate(deltaG.e.transfer = deltaG/e.transfer)

# Define x and y axis labels - plug in to "ylab" and "xlab" in the plot below
xlabel = expression(paste('log',italic('a'), ' NH'[3]))
ylabel = expression(paste('log',italic('a'), ' O'[2]))

# Plot! If you want separate png files, use this
for(reaction in reaction_names){
  for(site in sites){
    png(file=paste0(reaction, "_gradient.png"), height=6, width=8,units='in',res=300)
    data = deltaG %>%
      filter(rxn.name == reaction & Site == site)%>%
      select(NH3, O2, deltaG.e.transfer)%>%
      pivot_wider(names_from = NH3, values_from = deltaG.e.transfer)%>%
      column_to_rownames(var = "O2") %>%
      as.matrix()
    rxn_text = unlist(deltaG %>%
                        filter(rxn.name == reaction & Site == site)%>%
                        select(text)%>%
                        distinct())
    filled.contour(NH3,pH, data,
                   color.palette = ifelse(getRversion() >= "3.6.0", function(n) hcl.colors(n), topo.colors),
                   key.title = {par(cex.main=.8);title(main=DG.per.e)},
                   plot.axes = {
                     legend("topleft", legend = rxn_text$text, bty = "n", inset = c(0,.85))
                     axis(1)
                     axis(2)
                     title(ylab=ylabel,xlab = xlabel)})
    dev.off()
  }
}

# Plot! If you want a single pdf, use this
pdf(file = 'ammonia_concentration_gradient_plots.pdf', height = 8, width = 11, onefile = TRUE)
for(reaction in reaction_names){
  for(site in sites){
    data = deltaG %>%
      filter(rxn.name == reaction & Site == site)%>%
      select(NH3, O2, deltaG.e.transfer)%>%
      pivot_wider(names_from = NH3, values_from = deltaG.e.transfer)%>%
      column_to_rownames(var = "O2") %>%
      as.matrix()
    rxn_text = unlist(deltaG %>%
                        filter(rxn.name == reaction & Site == site)%>%
                        select(text)%>%
                        distinct())
    filled.contour(NH3,O2, data,
                   color.palette = ifelse(getRversion() >= "3.6.0", function(n) hcl.colors(n), topo.colors),
                   key.title = {par(cex.main=.8);title(main=DG.per.e)},
                   plot.axes = {
                     legend("topleft", legend = rxn_text$text, bty = "n", inset = c(0,.85))
                     axis(1)
                     axis(2)
                     title(ylab=ylabel,xlab = xlabel)})
  }
}

dev.off()

#### Switch axes for pH and temperature ####
# Plot! If you want a single pdf, use this
pdf(file = 'pH_temp_gradient_plots.pdf', height = 8, width = 11, onefile = TRUE)
for(reaction in reaction_names){
  for(site in sites){
    data = deltaG %>%
      filter(rxn.name == reaction & Site == site)%>%
      select(pH, temperature, deltaG.e.transfer) %>%
      pivot_wider(names_from = temperature, values_from = deltaG.e.transfer) %>% # switched names_from to temperature
      column_to_rownames(var = "pH") %>% # switched var = temp to var = pH
      as.matrix()
    rxn_text = unlist(deltaG %>%
                        filter(rxn.name == reaction & Site == site)%>%
                        select(text)%>%
                        distinct())
    filled.contour(pH_range, temperature_range, data, xlab = axis.label("pH"), # switched order of pH_range and temperature_range, xlab = pH
                   color.palette = ifelse(getRversion() >= "3.6.0", function(n) hcl.colors(n), topo.colors),
                   #plot.title = {par(cex.main=1);title(main = reaction)},
                   key.title = {par(cex.main=.8);title(main=DG.per.e)},
                   plot.axes = {
                     legend("topleft", legend = rxn_text$text, bty = "n", inset = c(0,.85))
                     lines(pK_NH3, temperature_range, lty = 2) # Switched order of pK_NH3 and temperature_range
                     text(8.6, 85,expr.species("NH3"))
                     text(6.9,85,  expr.species("NH4+"))
                     axis(1)
                     axis(2)
                     title(ylab="T (ºC)")})
  }
}
dev.off()


#### NH3 vs pH (pH on Y axis) ####
NH3 = seq(-12,-1)

nh3_pH = expand.grid(NH3 = seq(-12,-1,by=1), pH = pH_range)
gradient = seq(1:168)

log_activities = read_xlsx('logactivities.xlsx')%>%
  uncount(nrow(nh3_pH))%>%
  mutate(NH3 = nh3_pH$NH3)%>%
  mutate(`NH4+` = nh3_pH$NH3)%>%
  mutate(`H+` = -nh3_pH$pH)%>%
  mutate(pH = nh3_pH$pH)%>%
  add_column(grad = gradient)

# Convert to activity and pivot to long format
activities = log_activities %>%
  pivot_longer(`NH3`:`H+`, names_to = 'react_prod',values_to = 'activity') %>% # pivot from wide to long for joining
  mutate(activity = 10^activity)

# Calculate ∆Gº (kJ/mol) at 25ºC
G0 = reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c('.value','set'),
               names_pattern = '(.+).(.+)')%>% # pivot from wide to long dataframe
  unite('react_prod', reactant:product,na.rm = TRUE, remove = F)%>% # unite the reactants and products into one columns
  filter(!react_prod == '')%>% # remove any rows with missing react_prod values
  group_by(rxn.name)%>% # group by reaction name for calculations
  reframe(G0 = subcrt(react_prod, coeff, state, T = 25)$out$G/1000)%>% # Calculate ∆G0 with the subcrt function from CHNOSZ at 25ºC
  add_column(RT = rep((25 + 273.15)*.008314, nrow(reactions))) # Add RT value (gas constant * temperature)

# Calculate logQ
logQ = reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c('.value','set'),
               names_pattern = '(.+).(.+)')%>% # pivot from wide to long
  unite('react_prod', reactant:product, na.rm = TRUE, remove = F)%>% #unite the reactant and product columns into one column called react_prod
  left_join(activities, multiple = 'all')%>% # join with activities data
  filter(!is.na(activity))%>% # remove any activities with NA value
  mutate(Q = if_else(!is.na(reactant), activity^coeff, activity^coeff))%>%
  group_by(rxn.name, Site, grad)%>% #group on the reaction number and site 
  summarise(logQ = log(prod(Q))) #calculate logQ 

# Calculate ∆G!
deltaG <- G0 %>%
  left_join(logQ, multiple = 'all') %>% #join the G0 and logQ tables 
  left_join(reactions %>% select(rxn.name, e.transfer, text)) %>% #add the reaction number, number of electrons transferred, and minerals from each reaction 
  left_join(log_activities)%>%
  mutate(NH3 = 10^NH3)%>%
  mutate(NH3 = 10^`NH4+`)%>%
  mutate(deltaG = G0 + RT*logQ) %>% #calculate deltaG for each reaction at each site 
  mutate(deltaG.e.transfer = deltaG/e.transfer)

# Define x and y axis labels - plug in to "ylab" and "xlab" in the plot below
# Plot! If you want a single pdf, use this
pdf(file = 'ammonia_concentration_gradient_plots.pdf', height = 8, width = 11, onefile = TRUE)
for(reaction in reaction_names){
  for(site in sites){
    data = deltaG %>%
      filter(rxn.name == reaction & Site == site)%>%
      select( pH, NH3,deltaG.e.transfer)%>%
      pivot_wider(names_from = pH, values_from = deltaG.e.transfer)%>%
      column_to_rownames(var = "NH3") %>%
      as.matrix()
    rxn_text = unlist(deltaG %>%
                        filter(rxn.name == reaction & Site == site)%>%
                        select(text)%>%
                        distinct())
    filled.contour( NH3,pH_range,data,
                    color.palette = ifelse(getRversion() >= "3.6.0", function(n) hcl.colors(n), topo.colors),
                    key.title = {par(cex.main=.8);title(main=DG.per.e)},
                    plot.axes = {
                      legend("topleft", legend = rxn_text$text, bty = "n", inset = c(0,.85))
                      axis(1)
                      axis(2)
                      title(ylab="pH",xlab = xlabel)})
  }
}

dev.off()

