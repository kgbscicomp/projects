library(tidyverse)
library(dplyr)
setwd("C:/Users/gkarthik/Desktop/Michael/python/")

# Feed in the behavioral data
data = read.csv("Wave1_MIDMerge_Anticipation_2019-09-04 copy.csv")


# select subject of variables
data_short = data %>% 
  select(scan_id, 
         Run, 
         Condition, 
         Cue.OnsetTime, 
         Anticipation.Duration.Model, 
         Feedback.OnsetTime, 
         FeedbackDuration) %>% 
  mutate(Feedback_end = Feedback.OnsetTime + FeedbackDuration, 
         CueOnset = Cue.OnsetTime, 
         subj = scan_id, 
         Anticipation = Anticipation.Duration.Model)


# create full run TR/TIME length, e.g. 407 (0.8s TRs) * 2 (runs), range .8 to 651 ever .8
time = data.frame(TR = seq(0.8,651,.8))

# take data from data_short, group by subject id (104 subjects x 100 trails each)
beh_grouped = data_short%>%
  group_by(subj)%>%
  mutate(N = n())

beh_grouped[is.na(beh_grouped)] = 99
subjs = unique(beh_grouped$subj)
#beh = beh_grouped
counter = 1
for(r in 1:length(subjs)){
  time = data.frame(TR = seq(0.8,651,.8))
  subj = subjs[r]
  subset_data = beh_grouped[beh_grouped$subj==subj,]
  subset_data$N = 1:100
  
  subset_data$Cue.OnsetTime[subset_data$N>50] = subset_data$Cue.OnsetTime[subset_data$N>50] +
                                                subset_data$Feedback_end[subset_data$N==50]

  #### Create a distance matrix to compute distance between every Cue.Onset time and TR values
  dist_matrix = replicate(length(subset_data$Cue.OnsetTime),seq(0.8,651,.8))
  
  dist_matrix = abs(subset_data$Cue.OnsetTime - t(dist_matrix))
  
  #### Find the index of the minimum absolute distance for each of the TR values 
  inds = apply(dist_matrix,2,which.min)
  
  time$N = inds
  
  #### Use the index calculated above to do the merge
  merged_data = merge(time, subset_data, by = "N", all =TRUE)
  
  
  ## Find the absolute differences between the cueonset time and TR times
  time$cueonset = merged_data$Cue.OnsetTime
  time$diff = abs(time$cueonset - time$TR)
  
  grouped_merged = time %>%
    group_by(N) %>%
    mutate(minimum_diff = min(diff)) %>%
    mutate(TRmerge = TR[diff==minimum_diff]) %>%
    ungroup() %>%
    select(N, TRmerge) 
  
  grouped_merged = unique(grouped_merged[c("N","TRmerge")])
  
  time = time %>% 
    select(N, TR)
  
  subset_merged = merge(subset_data, grouped_merged, by = "N")
  colnames(time) = c("N","TRmerge")
  
  subset_merged = merge(time, subset_merged, by = "TRmerge", all.x= TRUE)
  
  subset_merged = subset_merged[ , !(names(subset_merged) %in% c("N", "N.y", "subj"))]
  
  names(subset_merged)[names(subset_merged) == "N.x"] = "N"
  names(subset_merged)[names(subset_merged) == "TRmerge"] = "TR"
  
  # created directory output path where to save subjects data
  directory = "C:/Users/gkarthik/Desktop/Michael/python"
  # tailor this for each subj (i.e. current subject in loop)
  output_dir = file.path(directory, paste(subj,"/",sep=""))
  
  # If subject folder doesn't exist, create one
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir exists!")
  }
  # create output file for subjects' 100 trials (e.g., we save the 100 trials for a single subject at a time, per loop)
  outfile = file.path(output_dir, paste(subj,".csv",sep=""))
  
  # save CSV into the subjects' directory created above
  write.csv(subset_merged, file = outfile, row.names = FALSE)  
}