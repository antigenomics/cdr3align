avian <- read.csv('avianHabitat.csv')

str(avian)
head(avian)
summary(avian)

any(!complete.cases(avian))
any(avian$PDB < 0)
any(avian$PDB > 100)

check_percent_range <- function(x) {
  any(x < 0 | x > 100)
}

names(avian)

coverage_variables <- names(avian)[-(1:4)][c(T, F)]

avian$total_coverage <- rowSums(avian[, coverage_variables])
summary(avian$total_coverage)
