library(tidyverse)
library(magrittr)
library(lubridate)

###############################################################################
### 01. READ IN DATA

kelts <- read_csv("data/kelts_2016.csv")

merwin <- read_csv("data/merwintrap_2016.csv")

tangle <- read_csv("data/tanglenet_2016.csv") %>%
  filter(LIFE != "S") %>% #remove smolts
  mutate(DATE = mdy(gsub("/16", "/2016",DATE)),
         ORIGIN = if_else(MARKS == "N" & STUB_DORSAL == "N", 
                          #Natural origin are fish without BWT or stubby dorsal
                         "NOR",
                         "HOR")) %>% 
  arrange(PIT_ID, DATE) %>%   
  replace_na(list(PIT_ID = "No PIT"))

trap <- read_csv("data/Summary_of_MerwinTrapData_2016.csv") %>%
  mutate(DATE = mdy(DATE))

# Check to see what unique combinations of fish captures
# exist in tanlge net data
tangle %>% 
  mutate(COND = paste(RETENTION, MATURITY, sep = '_')) %>% 
  distinct(COND)


###############################################################################
### 02. Process Tangle Net Data

#3/8 - Event 1
#3/15 - Event 2
#3/23 - Event 3
#3/29 - Event 4
#4/5 - Event 5
#4/12 - Event 6
#4/19 - Event 7
#4/26 - Event 8
#5/3 - Event 9
#5/10 - Event 10
# Event 11 any trapping after 5/10

# Ecological States: 1 - Potential Population
#                    2 - Returning Fish
#                    3 - Kelt/Spawner
#                    4 - Trap
# Observed States:   0 - Not Observed
#                    2 - Tangle net trapped returner
#                    3 - Tangle Net trapped Kelt
#                    4 - Trapped

##Need to make four(?) input matrices
#Observed state matrix is data input and consists of observations only (No NA)
#Ecological state matrix is data input and consists of known ecological states and NA otherwise
#Availability matrix is just a helper and consists of either 1's (if fish *might* be available to recapture) and 0 otherwise
#Initial ecological state matrix, takes the ecological state matrix and provides initial values for it

##Need to decide what to do about PIT ID 770E4EC
##Identified as HOR on 3/23 then apparently recapture and identified as NOR on 4/19
# we can try to run the model both ways, first with it identified as HOR and 
#second with it identified as NOR
#For what it's worth fish was originally identified as BWT and stubby dorsal residual when first captured in 2014

tangle %>% 
  filter(PIT_ID == "770E4EC")

tangle_1 = tangle %>% 
  mutate(ORIGIN = if_else(PIT_ID == "770E4EC", "HOR", ORIGIN))

# tangle_2 = tangle %>% 
#   mutate(ORIGIN = if_else(PIT_ID == "770E4EC", "NOR", ORIGIN))


###############################################################################
#### 02A. Create observed state data for tangle net data

observed_state <- tangle_1 %>% 
  mutate(COND = if_else(MATURITY != "K", 2L, 3L), 
         Event_0 = 0L) %>% 
  select(DATE, ORIGIN, GENDER, LIFE, PIT_ID, COND, Event_0) %>% 
  spread(key = DATE, value = COND) %>% 
    rename(Event_1 = `2016-03-08`, 
           Event_2 = `2016-03-15`, 
           Event_3 = `2016-03-23`, 
           Event_4 = `2016-03-29`, 
           Event_5 = `2016-04-05`, 
           Event_6 = `2016-04-12`, 
           Event_7 = `2016-04-19`, 
           Event_8 = `2016-04-26`, 
           Event_9 = `2016-05-03`, 
           Event_10 = `2016-05-10`) %>% 
  mutate(Event_11 = NA)

# Check for duplicate fish IDs among rows,
# should only be for the two No_PIT fish
observed_state %>% 
  group_by(PIT_ID) %>% 
  summarise(n = n()) %>% 
  filter(n > 1)

# 3DD003BE8CB09 initially captured as kelt at event 6, then identified as ripe 
# turning original kelt designation into ripe
observed_state[observed_state$PIT_ID == "3DD003BE8CB09", "Event_6"] = 2L

# 3DD003BE8CB1E initially captured as kelt at event 5, then identified as ripe 
# turning original kelt designation into ripe
observed_state[observed_state$PIT_ID == "3DD003BE8CB1E", "Event_5"] = 2L

#3DD003BE8CB17 trapped on 4/6
observed_state[observed_state$PIT_ID == "3DD003BE8CB17", "Event_6"] = 4L

#3DD003BE8CB21 trapped on 4/13
observed_state[observed_state$PIT_ID == "3DD003BE8CB21", "Event_7"] = 4L

#3DD003BE8CAC6 trapped on 4/1
observed_state[observed_state$PIT_ID == "3DD003BE8CAC6", "Event_5"] = 4L

#3DD003BE8D67C trapped on 4/3
observed_state[observed_state$PIT_ID == "3DD003BE8D67C", "Event_5"] = 4L

#3D6001569A080 trapped on 5/17, but identified as residual
observed_state[observed_state$PIT_ID == "569A080", "Event_11"] = 4L

#3DD003BE8CB18 marked as MT recapture, but specifed as BROOD in data... Need to investigate
#observed_state[observed_state$PIT_ID == "3DD003BE8CB18", "Event_11"] = 5
#Looks like it appeared in the trap after being removed for broodstock... will ignore this datum.

###############################################################################
#### 02A. Create ecological state data for tangle net data
####      Should be NA if we don't know for certain what state fish was in

#Use observed states as starting point for ecological state matrix
ecological_state <- observed_state

#Replace NAs with 0 in observed data
observed_state %<>% 
  replace_na(list(Event_1  = 0L,
                  Event_2  = 0L,
                  Event_3  = 0L,
                  Event_4  = 0L,
                  Event_5  = 0L,
                  Event_6  = 0L,
                  Event_7  = 0L,
                  Event_8  = 0L,
                  Event_9  = 0L,
                  Event_10 = 0L,
                  Event_11 = 0L))

#Set initial state for all fish to "potential"
ecological_state$Event_0 = 1L

##Assign known states

#3DD003BE8CB17 tangle netted at Event_2, trapped at Event_6

ecological_state[ecological_state$PIT_ID == "3DD003BE8CB21",8:11] = 2
ecological_state[ecological_state$PIT_ID == "3DD003BE8CB21", 13:16] = 4

#3DD003BE8CAC6 tangle netted at Event_3, trapped at Event_5
ecological_state[ecological_state$PIT_ID == "3DD003BE8CAC6", 9] = 2
ecological_state[ecological_state$PIT_ID == "3DD003BE8CAC6", 11:16] = 4

#3DD003BE8D67C tangle netted at Event_4, trapped at Event_5
ecological_state[ecological_state$PIT_ID == "3DD003BE8D67C",11:16] = 4

#3D6001569A080 tangle netted at Event_5, trapped at Event_11
ecological_state[ecological_state$PIT_ID == "569A080",11:15] = 2


## Assign fish identified as kelts to always be kelts thereafter.
for(i in 1:nrow(ecological_state)){
  
  first_kelt = min(which(as.numeric(ecological_state[i, 6:16]) == 3))
  
  if(first_kelt != Inf){
    message(i) 
    ecological_state[i, (5 + first_kelt):16] = 3
  }
}

# Check for fish captured more than once, if in same state during each capture
# we can safely assign known state in-between
tangle_1 %>% 
  filter(PIT_ID != "No PIT") %>% 
  group_by(PIT_ID) %>% 
  summarise(n = n()) %>% 
  filter(n > 1)

#3DD003BE8CB09 captured twice - Event 6 and 8, ripe both times. stay ripe in betweeen
ecological_state[ecological_state$PIT_ID == "3DD003BE8CB09", 12] = 2

#3DD003BE8CB1E captured twice - Event 5 and 8, ripe both times. stay ripe in betweeen
ecological_state[ecological_state$PIT_ID == "3DD003BE8CB1E",11:12] = 2

#770E4EC captured twice, ripe both times. stay ripe in betweeen
ecological_state[ecological_state$PIT_ID == "770E4EC",9:11] = 2

###############################################################################
### 03. Process Trap Data

## Observed data

trap_observed <- trap %>% 
  # ignore trap capture prior to first tangle net event, those fish
  # were never subject to tangle netting and provide no information
  # on fish on spawning ground
  filter(DATE > mdy("03/08/16")) %>% 
  # Assign daily counts to tangle net events
  mutate(Event = cut(DATE, 
                     breaks = c(mdy("03/08/16"), mdy("03/15/16"), 
                                mdy("03/23/16"), mdy("03/29/16"), 
                                mdy("04/05/16"), mdy("04/12/16"),
                                mdy("04/19/16"), mdy("04/26/16"),
                                mdy("05/03/16"), mdy("05/10/16"),
                                mdy("06/27/16")),
                     labels = c("Event_2", "Event_3", "Event_4",
                                "Event_5", "Event_6", "Event_7",
                                "Event_8", "Event_9", "Event_10",
                                "Event_11"),
                     right = TRUE)) %>% 
  group_by(Event) %>% 
  # Sum total count by fish type for each interval between tangle net events
  # HOR_M_A means hatchery origin, male, adult life - will separate out in subsequent step
  # We don't have any residuals in trap data beyond the one known above from tangle net
  # So we treat all fish as anadramous. In future years, should have a separate count of residuals
  summarise(HOR_M_A = sum(HOR_M), 
            HOR_F_A = sum(HOR_F),
            NOR_M_A = sum(NOR_M),
            NOR_F_A = sum(NOR_F)) %>% 
  # Gather the counts into a single column
  gather(-Event, key = Type, value = Count) %>%
  # seperate out fish type information into three columns
  separate(Type, c("ORIGIN", "GENDER", "LIFE"), sep = "_", remove = T) %>% 
  # create row for each count
  uncount(Count) %>%
  # Assign unique ID to each fish and set state to 4
  mutate(PIT_ID = paste("Trap", 1:n(), sep = "_"),
         COND = 4L) %>% 
  # Back to wide format
  spread(key = Event, value = COND)

# Use observed data as starting point for ecological data
trap_ecological <- trap_observed %>% 
  mutate(Event_0 = 1L,
         Event_1 = NA) %>% 
  select(ORIGIN, GENDER, LIFE, PIT_ID, Event_0,
         Event_1, Event_2, Event_3, Event_4, 
         Event_5, Event_6, Event_7, Event_8, 
         Event_9, Event_10, Event_11)
  
# Reassing NAs to 0 in observed data
trap_observed %<>% 
  mutate(Event_0 = 0L,
         Event_1 = 0L) %>% 
  replace_na(list(Event_1 = 0L,
                  Event_2 = 0L,
                  Event_3 = 0L,
                  Event_4 = 0L,
                  Event_5 = 0L,
                  Event_6 = 0L,
                  Event_7 = 0L,
                  Event_8 = 0L,
                  Event_9 = 0L,
                  Event_10 = 0L,
                  Event_11 = 0L)) %>% 
  select(ORIGIN, GENDER, LIFE, PIT_ID, Event_0,
         Event_1, Event_2, Event_3, Event_4, 
         Event_5, Event_6, Event_7, Event_8, 
         Event_9, Event_10, Event_11)



## Assign fish trapped as always being trapped thereafter and
## assign fish to be returning for timestep immediately prior to trapping
for(i in 1:nrow(trap_ecological)){
  trap_idx = min(which(as.numeric(trap_ecological[i, 6:16]) == 4L)) + 5
  
  trap_ecological[i, trap_idx - 1] = 2L
  if(trap_idx < 16)
    trap_ecological[i, (trap_idx + 1):16] = 4L
}


###############################################################################
### 04. Create availability helper matrix
## Availability matrix is 0 if fish are known to be not available for recapture
## either because they taken for broodstock or dead or trapped

#First step, identify time at which fish become unavailable with 1
availability_state <- tangle_1 %>% 
  mutate(COND = ifelse(RETENTION %in% c("EUTHANIZED", "MORTALITY", "BROOD"), 1L, 2L)) %>% 
  select(DATE, ORIGIN, GENDER, LIFE,  PIT_ID, COND) %>% 
  spread(key = DATE, value = COND) %>% 
  mutate(Event_0 = 2) %>% 
  rename(Event_1 = `2016-03-08`, 
         Event_2 = `2016-03-15`, 
         Event_3 = `2016-03-23`, 
         Event_4 = `2016-03-29`, 
         Event_5 = `2016-04-05`, 
         Event_6 = `2016-04-12`, 
         Event_7 = `2016-04-19`, 
         Event_8 = `2016-04-26`, 
         Event_9 = `2016-05-03`, 
         Event_10 = `2016-05-10`) %>% 
  mutate(Event_11 = NA) %>% 
  select(ORIGIN, GENDER, LIFE, PIT_ID, Event_0,
         Event_1, Event_2, Event_3, Event_4, 
         Event_5, Event_6, Event_7, Event_8, 
         Event_9, Event_10, Event_11)

#3DD003BE8CB17 trapped on 4/6
availability_state[availability_state$PIT_ID == "3DD003BE8CB17", "Event_6"] = 4L

#3DD003BE8CB21 trapped on 4/13
availability_state[availability_state$PIT_ID == "3DD003BE8CB21", "Event_7"] = 4L

#3DD003BE8CAC6 trapped on 4/1
availability_state[availability_state$PIT_ID == "3DD003BE8CAC6", "Event_5"] = 4L

#3DD003BE8D67C trapped on 4/3
availability_state[availability_state$PIT_ID == "3DD003BE8D67C", "Event_5"] = 4L

#3D6001569A080 trapped on 5/17, but identified as residual
availability_state[availability_state$PIT_ID == "569A080", "Event_11"] = 4L

for (i in 1:nrow(availability_state)){
  unavail_idx <- min(which(as.numeric(availability_state[i,5:16]) %in% c(1L,4L))) + 4
  if (unavail_idx == Inf){
    availability_state[i,5:16] = 1L
  } else {
    availability_state[i,5:unavail_idx] = 1L
    if(unavail_idx < 16)
      availability_state[i,(unavail_idx + 1):16] = 0L
  }
}

trap_available <- trap_observed

for (i in 1:nrow(trap_available)){
  unavail_idx <- min(which(as.numeric(trap_available[i,5:16]) == 4L)) + 4
  message(paste(i, ",", unavail_idx))
  
  if (unavail_idx == Inf){
    trap_available[i,5:16] = 1L
  } else {
    trap_available[i,5:unavail_idx] = 1L
    if(unavail_idx < 16)
      trap_available[i,(unavail_idx + 1):16] = 0L
  }
}

###############################################################################
### 06. Bring it all together and prepare for JAGS

obs_use <- bind_rows(observed_state, trap_observed) %>% 
  mutate(origin_state = if_else(ORIGIN == "HOR", 0, 1),
         sex_state    = if_else(GENDER == "M", 0, 1),
         resid_state  = if_else(LIFE != "R", 0, 1)) %>% 
  arrange(PIT_ID) %>% 
  select(origin_state, sex_state, resid_state, 
         Event_0, Event_1, Event_2, Event_3, 
         Event_4, Event_5, Event_6, Event_7, 
         Event_8, Event_9, Event_10, Event_11)

eco_use <- bind_rows(ecological_state, trap_ecological) %>% 
  mutate(origin_state = if_else(ORIGIN == "HOR", 0, 1),
         sex_state    = if_else(GENDER == "M", 0, 1),
         resid_state  = if_else(LIFE != "R", 0, 1)) %>% 
  arrange(PIT_ID) %>% 
  select(origin_state, sex_state, resid_state, 
         Event_0, Event_1, Event_2, Event_3, 
         Event_4, Event_5, Event_6, Event_7, 
         Event_8, Event_9, Event_10, Event_11)

avail_use <- bind_rows(availability_state, trap_available) %>% 
  arrange(PIT_ID) %>% 
  select(Event_0, Event_1, Event_2, Event_3, 
         Event_4, Event_5, Event_6, Event_7, 
         Event_8, Event_9, Event_10, Event_11)


obs_pad <- data_frame(origin_state  = rep(NA, 1000),
                      sex_state  = NA,
                      resid_state  = NA,
                      Event_0 = 0,
                      Event_1 = 0,
                      Event_2 = 0,
                      Event_3 = 0,
                      Event_4 = 0,
                      Event_5 = 0,
                      Event_6 = 0,
                      Event_7 = 0,
                      Event_8 = 0,
                      Event_9 = 0,
                      Event_10 = 0,
                      Event_11 = 0)

eco_pad <- data_frame(origin_state  = rep(NA, 1000),
                      sex_state  = NA,
                      resid_state  = NA,
                      Event_0 = 1,
                      Event_1 = NA,
                      Event_2 = NA,
                      Event_3 = NA,
                      Event_4 = NA,
                      Event_5 = NA,
                      Event_6 = NA,
                      Event_7 = NA,
                      Event_8 = NA,
                      Event_9 = NA,
                      Event_10 = NA,
                      Event_11 = NA)

avail_pad <- data_frame(Event_0 = rep(1, 1000),
                        Event_1 = 1,
                        Event_2 = 1,
                        Event_3 = 1,
                        Event_4 = 1,
                        Event_5 = 1,
                        Event_6 = 1,
                        Event_7 = 1,
                        Event_8 = 1,
                        Event_9 = 1,
                        Event_10 = 1,
                        Event_11 = 1)

fish_type_JAGS <- bind_rows(eco_use , eco_pad) %>% 
                    select(origin_state, sex_state, resid_state)
eco_JAGS       <- bind_rows(eco_use , eco_pad) %>% 
                    select(-origin_state, -sex_state, -resid_state) %>% 
                    as.matrix()
obs_JAGS       <- bind_rows(obs_use, obs_pad) %>% 
                    select(-origin_state, -sex_state, -resid_state) %>% 
                    as.matrix()
avail_JAGS     <- bind_rows(avail_use, avail_pad) %>% as.matrix()         

##Rewrite 0 as 1 in observed data
obs_JAGS[obs_JAGS==0] = 1

pHOS.Input <- list(z      = eco_JAGS,
                   y      = obs_JAGS,
                   avail  = avail_JAGS,
                   origin = fish_type_JAGS$origin_state,
                   sex    = fish_type_JAGS$sex_state,
                   resid  = fish_type_JAGS$resid_state,
                   M      = nrow(eco_JAGS),
                   K      = ncol(eco_JAGS))


