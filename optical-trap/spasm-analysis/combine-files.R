library(data.table)
library(readxl)

spasm_files <- list.files("optical-trap/spasm-analysis/control",
                         pattern = "spasm",
                         recursive = TRUE,
                         full.names = TRUE)

dat <- lapply(spasm_files, read_xlsx)

writexl::write_xlsx(dat, "optical-trap/spasm-analysis/control/spasm-included.xlsx")



#########################################3

spasm_files <- list.files("optical-trap/spasm-analysis/10uM-dani",
                         pattern = "spasm-manual.xlsx",
                         recursive = TRUE,
                         full.names = TRUE)

dat <- lapply(spasm_files, read_xlsx, col_names = F)

writexl::write_xlsx(dat, "optical-trap/spasm-analysis/10uM-dani/spasm-included-manual.xlsx", col_names = F)
