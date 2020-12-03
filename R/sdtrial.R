## Create a normal distribution
normdist = rnorm(100000, mean = 0, sd = 3)
hist(normdist)

## Find their mean and SD to double check we have the required parameters 
mean(normdist)
sd(normdist)

## Find the data points that lie within one SD on either side of the mean
onesdright = normdist[normdist<sd(normdist)&normdist>0]
onesdleft = normdist[normdist>-sd(normdist)&normdist<0]

## Find the percentage of data that lies within one SD on either side of the mean 
percent_ofdata_right = length(onesdright)/length(normdist)
percent_ofdata_left = length(onesdleft)/length(normdist)

print(paste0("Percent of data one SD to right of the mean: ", percent_ofdata_right))
print(paste0("Percent of data one SD to left of the mean: ", percent_ofdata_left))


####################### Do it for 2 SD now ###############################


## Find the data points that lie within one SD on either side of the mean
onesdright = normdist[normdist<2*sd(normdist)&normdist>sd(normdist)]
onesdleft = normdist[normdist>-2*sd(normdist)&normdist< -1*sd(normdist)]

## Find the percentage of data that lies within one SD on either side of the mean 
percent_ofdata_right = length(onesdright)/length(normdist)
percent_ofdata_left = length(onesdleft)/length(normdist)

print(paste0("Percent of data one SD to right of the mean: ", percent_ofdata_right))
print(paste0("Percent of data one SD to left of the mean: ", percent_ofdata_left))
