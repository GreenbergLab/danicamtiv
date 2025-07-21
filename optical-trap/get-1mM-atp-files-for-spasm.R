library(data.table)

files <- list.files("~/lasertrapr/project_10uM-OM_10uM-atp",
                    pattern = "options.csv",
                    full.names = TRUE,
                    recursive = TRUE)

dat <- rbindlist(lapply(files, fread), fill = TRUE)

dat <- dat[include == TRUE & review == TRUE]

for(i in 1:nrow(dat)){
  print(i)
 f1 <-  list.files("/media/brent/50FABD0BFABCEE7A/trap-data",
             pattern = dat$original_filename[[i]],
             recursive = TRUE,
             full.names = TRUE)
 file.copy(from = f1, to = "~/washu/danicamtiv/optical-trap/spasm-analysis/10uM-om")
}
