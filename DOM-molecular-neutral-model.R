# Script to fit the neutral model of DOM molecules to high resolution mass spectrometry data - 6/1/2023
# She Zhixiang; Shezhixiang@mail.hfut.edu.cn
# She Zhixiang et al. Quantifying stochastic processes in shaping dissolved organic matter pool with high resolution mass spectrum.
# dom: A table for HRMS data with samples as rows and  molecules as columns.
# pool:  The assemblage of DOm molecules. DOM samples analyzed at distinct locations as local pools, while the collection of various samples as source pool.

#Make sure the following R packages are installed
require(vegan)
require(Hmisc)
require(minpack.lm)
require(stats4)

Sample_Name = "Dataset_Name" # Sample name for output

setwd("/path/to/ICR_data") # Load in data #

data<-read.csv('dom.csv',head=T, row.names=1, sep = "\t")

#Calculate the relative abundance of each molecular formula in each sample with a threshold value of 100
scale_to_100 <- function(x) {
  scaled <- (x / max(x)) * 100
  return(scaled)
}
dom <- as.data.frame(lapply(data, scale_to_100))
rownames(dom) <- rownames(data)
dom<-t(dom)

N <- mean(apply(dom, 1, sum)) #Calculate the number of molecules per sample
p.m <- apply(dom, 2, mean) 
p.m <- p.m[p.m != 0] 
p <- p.m/N #Calculate the mean relative abundance of each molecules across samples

#Calculate the occurrence frequency of each molecule across samples
dom.bi <- 1*(dom>0)
freq <- apply(dom.bi, 2, mean) 
freq <- freq[freq != 0]

#Combine
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]

d = 1/N #Calculate the limit of detection

#Fit model parameter m using Non-linear least squares (NLS)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.ci <- confint(m.fit, 'm', level=0.95)

m.fit  #Get the model parameter m value
summary(m.fit)#produce summaries of the results of the model fitting

freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE) #The neutral predicted frequency per molecule
pred.ci <- binconf(freq.pred*nrow(dom), nrow(dom), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  #Get the goodness of fit R2 of the model

fitstats <-data.frame(p,freq,freq.pred,pred.ci[,2:3])

# Distinguish between neutral and non-neutral molecules
fitstats$Type <- NA
fitstats$Type[fitstats$freq >= fitstats$Upper] <- "above"
fitstats$Type[fitstats$freq <= fitstats$Lower] <- "below"
fitstats$Type[fitstats$freq < fitstats$Upper & fitstats$freq > fitstats$Lower] <- "neutral"

# Export the neutral modeling results for HRMS data
write.csv(fitstats, paste(Sample_Name, "_neutral_result", sep = ""), quote = F, row.names = F)

#Draw the statistical map for neutral model fit result

library(ggplot2)

fitstats <- data.frame(p, freq, freq.pred, Lower = pred.ci[, 2], Upper = pred.ci[, 3])
fitstats$inter.col <- '#377eb8'  # Initialize all data points to blue
fitstats$inter.col[fitstats$freq <= fitstats$Lower] <- '#4daf4a'  # Change the color of non-neutral molecules below neutral prediction to green
fitstats$inter.col[fitstats$freq >= fitstats$Upper] <- '#ff7f00'  # Change the color of non-neutral molecules above neutral prediction to orange

ggplot(fitstats, aes(x = log10(p), y = freq)) +
  geom_point(shape = 20, color = fitstats$inter.col, alpha = 0.5, size = 2) +
  geom_line(aes(y = freq.pred), color = 'black', size = 1) +
  geom_line(aes(y = Lower), color = 'black', size = 1, linetype = 'dashed') +
  geom_line(aes(y = Upper), color = 'black', size = 1, linetype = 'dashed') +
  ylab('Frequency of Occurrence') +
  xlab('Mean Relative Abundance (log10)') +
  theme_bw()

draw.text <- function(just, i, j) {
  grid.text(paste(Sample_Name,"\n","Rsqr=",round(Rsqr,2),"\n","m=",round(coef(m.fit),2)), x=x[j], y=y[i], just=just)}

x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")

draw.text(c("centre", "bottom"), 4, 2)







