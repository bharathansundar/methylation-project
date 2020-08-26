setwd("C:/Users/{kr.pA}/Downloads")

#all work with test data of first two rows
test = read.csv("methylated data/t1.csv")
test = as.data.frame(test)

for (val in test["ID_REF"])
{
  rownames(test) <- val
}

#this method runs into problems later
test$ID_REF <- NULL

logtest <- test
logtest[, 4:7] <- log(test[4:7], exp(1))


#this is the real stuff with full augmented dataset
dataframe = read.csv("methylated data/values.csv")
typeof(dataframe)
dataframe = as.data.frame(dataframe)

for (val in dataframe["X"])
{
  rownames(dataframe) <- val
}

#dropping extra column
dataframe = dataframe[-c(1)]

logdf <- dataframe
logdf[, 1:94] <- log(dataframe[1:94], exp(1))
logdf$gene <- rownames(logdf)

#ggplot(data=logdf, aes(x=t.statistics, y=-log10(p.values), label=gene)) + geom_point() + theme_minimal() + geom_text() + geom_vline(xintercept=c(-75, 75), col="red")

library(EnhancedVolcano)
EnhancedVolcano(logdf, lab = rownames(logdf), x='t.statistics', y=('p.values'), xlim = c(-100, 100))



