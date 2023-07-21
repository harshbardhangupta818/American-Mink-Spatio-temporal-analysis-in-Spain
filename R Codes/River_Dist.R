# set the working directory
w_dir <- "/home/hbg/Desktop/Spain_american_mink_project"
setwd(w_dir)

# Read occurrence only data of American mink population in Spain
occ = read.csv("data/WATER_DIST_madrid.csv")
str(occ)


library(ggplot2)

library(minpack.lm)

# Histogram 
get_hist <- function(p) {
  d <- ggplot_build(p)$data[[1]]
  data.frame(x = d$x, xmin = d$xmin, xmax = d$xmax, y = d$y)
}

funx <- function(model_p,x) {
  a <- summary(model_p)$parameters[1]
  b <- summary(model_p)$parameters[2]
  y = Mm(x,a,b)
  return (y)
}

p<- ggplot(occ, aes(x=HubDist)) + geom_histogram(aes(y=..density..), bins=30, colour="black", fill="white") +
  geom_density(colour="red") + labs(title="River Distance histogram plot",x="Distance from River (Km)", y = "Density")

plot(p)


result<- hist(occ$HubDist)

res<- get_hist(p)

str(res)
str(result)

df1 <- data.frame(x =res$x, y = res$y) 
df2 <- data.frame(x =result$mids, y = result$counts) 
head(df1)

# Fitting exponential curve on the plot
library(minpack.lm)
G <- expression(a*(exp(b*x)))


Mm <- function(x, a, b) {eval(G)}

model= nlsLM(y~Mm(x,a,b), data = df1, start= list( a = 0 , b = 0), algorithm = "LM", control = list(maxiter = 1000))
model2 = nlsLM(y~Mm(x,a,b), data = df2, start= list( a = 0 , b = 0), algorithm = "LM", control = list(maxiter = 1000))

print(model)
print(model2)

s_ <-  seq(0, 100, 5)


pred1= funx(model,s_)
pred2= funx(model2,s_)


plot(res$x,res$y, xlim=c(0,100), xlab = "Distance From River(Km)", ylab ="Denisty")
lines(s_, pred1, lty = 3 , lwd = 2, col ="red")


plot(result$mids,result$counts, xlim=c(0,100), xlab = "Distance From River(Km)", ylab ="Frequency")
lines(s_, pred2, lty = 3 , lwd = 2, col ="red")
