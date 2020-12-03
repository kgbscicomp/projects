## Creata funky random data
a = rnorm(100, mean=0, sd = 3)
b = rnorm(100, mean=1,sd = 3)
c = c(rep(1,50),rep(2,50))
data = cbind(a,b,c)
data = data.frame(data)
colnames(data) = c("a","b","dv")


## Build a model with interaction specified
model1 = lm(data = data, dv~a+b+a:b)
summary(model1)


##Create an "interaction" term the way R's inbuilt model calculates interactions
datacopy = data

## This is how R calculates interactions- it is literally multiplying of the two variables you are 
## interested in
datacopy$interact = data$a*data$b   

## Build a second model, this time with the interaction variable that you created as one of the IVs. 
model2 = lm(data = datacopy, dv~a+b+interact)
summary(model2)

## Now compare model1 and model2. They should be identical
anova(model1,model2)


## You can also check that they are identical by looking at their individual model summaries
summary(model1)
summary(model2)