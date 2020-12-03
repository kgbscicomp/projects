library(dplyr)
library(ggplot2)
setwd("C:/Users/gkarthik/Desktop/Michael/python/")

df.allSubs_long = read.csv("output_timeseries_beh_match_long_test.csv")

## General bootstrap code [for when you calculate the bootstrap CIs for simple mean of a vector]

bootbold = function(percentchange){
  b = 5e3             # Specify the number of replicates for measures and ROI
  conf_level = 0.05   # Specify the confidence level for the upper and lower bounds
  
  indx = sample(1:length(percentchange), b*length(percentchange), replace = TRUE)
  Bx = percentchange[indx]
  dim(Bx) = c(b,length(percentchange))
  
  lcb_na = mean(percentchange) - qnorm(1 - conf_level/2)*sd(rowMeans(Bx))
  
  ucb_na = mean(percentchange) + qnorm(1 - conf_level/2)*sd(rowMeans(Bx))
  
  estimate = mean(percentchange)
  
  return(cbind(estimate,lcb_na,ucb_na))
}

test_plotttest = df.allSubs_long %>%
  filter(voi == "L_VS", Condition == "LgReward" | Condition == "Triangle") %>%
  select(tr, BOLD, Condition) 

allvals = data.frame()
for(i in 1:length(unique(test_plotttest$tr))){
  subdata = test_plotttest[test_plotttest$tr==i,]
  pval = t.test(subdata$BOLD[subdata$Condition=="LgReward"],
                subdata$BOLD[subdata$Condition=="Triangle"])$p.value
  allvals = rbind(allvals,cbind(i,pval))
}
colnames(allvals) = c("tr","pval")

merged = merge(allvals, test_plotttest, by = "tr")

test_plot = merged %>%
  group_by(tr, Condition) %>%
  mutate(mean_BOLD = bootbold(BOLD)[1],
         lcb = bootbold(BOLD)[2],
         ucb = bootbold(BOLD)[3]) %>%
  ungroup %>%
  select(tr, mean_BOLD, Condition, lcb, ucb, pval) %>%
  unique()

test_plot$label = ""
test_plot$label[test_plot$pval<0.05] = "*"
test_plot$label[test_plot$pval<0.01] = "**"
test_plot$label[test_plot$pval<0.005] = "***"


ggplot(data = test_plot, aes(x = tr, y = mean_BOLD, colour = Condition)) +
  geom_line(aes(linetype=Condition)) +
  geom_errorbar(aes(ymin=lcb, ymax=ucb))+
  scale_colour_brewer(palette = "Set1") +
  ylab("BOLD % Signal Change") +
  xlab("TR") +
  annotate("text", x = unique(test_plot$tr), 
             y = max(test_plot$ucb[test_plot$tr == unique(test_plot$tr)])+ 0.04, 
             label = test_plot$label[seq(1,28,2)])

