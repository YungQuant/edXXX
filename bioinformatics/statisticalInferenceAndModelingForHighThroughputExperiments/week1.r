library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

head(sampleInfo)
theDate="2005-06-27"
thisDate<-sampleInfo[sampleInfo$date==theDate,]
nrow(thisDate)

sum(sampleInfo$date=="2005-06-27")

sum(geneAnnotation$CHR=="chrY",na.rm=TRUE)

i = which(geneAnnotation$SYMBOL=="ARPC1A")
j = which(sampleInfo$date=="2005-06-10")
geneExpression[i,j]

colMedians=apply(geneExpression,2,median)
median(colMedians)

g <- factor(sampleInfo$group)
myttest <- function(e){
    x <- e[g==1]
    y <- e[g==0]
    return( t.test(x,y)$p.value )
}
pVals=apply(geneExpression,1,myttest) # "1" because we act on each row
min(pVals)
hist(pVals)

set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)
pvals <- replicate(1000,{
    control = sample(population[,1],12)
    treatment = sample(population[,1],12)
    t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)

mean(pvals<0.05)

mean(pals<0.05)

cases = rnorm(10,30,2)
controls = rnorm(10,30,2)
t.test(cases,controls)

set.seed(100)
pvals <- replicate(2000,{
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.val
})
head(pvals)
hist(pvals)
x<-seq(0,1,0.1)
y<-sapply(x,function(x){mean(pvals<x)})
plot(x,y)
abline(0,1,col="red")


set.seed(100)
B=1000
sigp<-replicate(B,{
    pvals <- replicate(20,{
        cases = rnorm(10,30,2)
        controls = rnorm(10,30,2)
        t.test(cases,controls)$p.val
    })
    sum(pvals<=0.05)
})
head(sigp)
table(sigp) ##just for illustration
mean(plessthan)
hist(sigp)
mean(sigp)

#reject null hypothesis at least once
mean(sigp>0)

