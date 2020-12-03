####### The following code works accurately only if the true labels are sorted and ordered such that
# 1's occur first and 0's occur next. Can make it work irrespective or ordering, but lazy.
# So while writing your function, make sure to sort the inputs accordingly. ###############

y = c(rep(1,10), rep(0,10))   ## Generating true labels 
yhat = runif(20, min=0, max=1) ## Generating sample yhats 

yhat_matrix = t(replicate(length(yhat),yhat))    ## Create a matrix of yhats 
tau = yhat                                        
yhat_matrix[yhat_matrix >= tau] = 1              ## In each row, values greater than tau is set to 1
yhat_matrix[yhat_matrix < tau] = 0               ## And values less than tau is set to 0. 

# Count true positives 
trp = t(yhat_matrix)[y==unique(y)[1]]                 ## Consider only the labels 1
dim(trp) = c(nrow(yhat_matrix)/2, ncol(yhat_matrix))  ## Extract values from yhat matrix for true labels 1
count_true_positives = colSums(trp)                   ## Count true positives 

# Count false positives 

fpos = t(yhat_matrix)[y==unique(y)[2]]
dim(fpos) = c(nrow(yhat_matrix)/2, ncol(yhat_matrix))
count_false_positives = colSums(fpos)

# Count false negatives 
fneg = t(yhat_matrix)[y==unique(y)[1]] 
dim(fneg) = c(nrow(yhat_matrix)/2, ncol(yhat_matrix))
count_false_negatives = colSums(fneg==0)

# Count true negatives 
tneg = t(yhat_matrix)[y==unique(y)[2]]
dim(tneg) = c(nrow(yhat_matrix)/2, ncol(yhat_matrix))
count_true_negatives = colSums(tneg==0)

# Sensitivity 
sensitivity = count_true_positives/(count_true_positives+count_false_negatives)

# Specificity 
specificity = count_true_negatives/(count_true_negatives+count_false_positives)

## Output Matrix- I have added an extra column from what the question asks (for ease of plotting). Remove this 

fpr = count_false_positives/(count_false_positives+count_true_negatives)
 
output = data.frame(cbind(yhat, 
                          count_true_positives, count_true_negatives,
                          count_false_positives, count_false_negatives,
                          sensitivity, specificity, fpr))

output = output[with(output, order(yhat)),]


plot(output$fpr,output$sensitivity, type = 'l')

####################### AUC with trapezoidal rule ########################### 

###### Can vectorize the following as well. But lazy #########
auc_trap_homework = function(x, y) {
  x_axis = rep(0,length(x)-1)
  y_axis = x_axis
  for (i in 2:length(x_axis)){
    x_axis[i-1] = abs(x[i] - x[i-1])
  }
  for (j in 1:length(x_axis)-1){
    y_axis[j] = (abs(y[j]) + abs(y[j+1]))/2
  }
  products = sum(x_axis*y_axis)
  return(products)

}

auc_trap_homework(output$fpr, output$sensitivity)

###### Compare with a ROC plot and AUC value using the ROCR and pROC packages in R ####
library(ROCR)
library(pROC)
pred = prediction( yhat, y)
perf = performance(pred,"tpr","fpr")
plot(perf)
package_auc = auc(y, yhat)
package_auc    ## For some cases this package gives weird values- Not sure why. Check with another package. 
