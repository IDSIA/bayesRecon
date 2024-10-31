# Generate the monthly count time series for the mixCond and TDcond tests
# CHANGE THE WORKING DIRECTORY BEFORE RUNNING
rm(list=ls())
library(bayesRecon)

set.seed(42)
vals <- stats::rpois(12*10,lambda = 2)


write.table(vals,file="./Monthly-Count_ts.csv",row.names = FALSE,sep=',',
            col.names = FALSE,quote = FALSE)

