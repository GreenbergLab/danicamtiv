library(data.table)

unreg <- fread("motility/unregulated-average-values.csv")
reg <- fread("motility//regulated-average-values.csv")

library(BSDA)

tsum.test(mean.x = 369.96, s.x = 57.51, n.x = 2,
          mean.y = 309.42, s.y = 33.46, n.y = 3)
