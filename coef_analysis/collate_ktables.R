library(tidyverse)

df.cts = read.csv("./slicing_ktable_CTS.csv")
df.exp = read.csv("./slicing_ktable_EXP.csv")
df.cts.mto = read.csv("./slicing_ktable_CTS-MTO.csv")

df.sto = bind_rows(df.cts, df.exp)
df = df.sto %>%
  full_join(df.cts.mto %>%
              select(-expe_set), by = c("miR", "rxnID"))

write.csv(df,
          "./slicing_ktable.csv",
          quote = T,
          row.names = F)
