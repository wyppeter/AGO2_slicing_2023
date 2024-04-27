# Helper objects and functions for name resolution, information look-up, basic analyses, and more
# Run before most R scripts

library(dplyr)

difflim = 60

# Get quant ratio ( exact(FB)/approx(guide*) ) -----
# Some experiments were planned using approximate concentrations
Q_0 = 1  # default
qratio = read.csv(
  paste0(
    "./",
    "quantratiobook.csv"
    ),
  stringsAsFactors = F
  ) %>%
  pull(quantratio, name = miR)
qLookUp = function(miRname){
  # Provides values
  raw.out = qratio[as.character(miRname)]
  raw.out[is.na(raw.out)] = Q_0
  return(unname(raw.out))
}
qLookUp.bool = function(miRname){
  # Provides boolean for ones that do not have a known ratio
  raw.out = qratio[as.character(miRname)]
  return(unname(is.na(raw.out)))  # TRUE == no known correction ratio
}

# Get, convert, and parse names/IDs -----
# Includes items to:
# - Convert legacy naming scheme to official naming scheme
# - Convert ID names to descriptive names for figures
# - Parse, and include or omit detailed information in labels (for figures and specific analyses)
base.name = function(name){
  # Retrieve parent miRNA
  gsub("#.+|_.+", "", name)
}
namedict = read.csv(
  paste0(
    "./",
    "guide_name_dict.csv"
  ),
  stringsAsFactors = F
) %>%
  mutate(guidename = base.name(rxnID),
         miR.base = base.name(miR))
namedict.print = namedict %>%
  select(miR, rxnID, rxnID.explicit, rxnID.explicit.wspe) %>%
  mutate(rxnID.print = sub("\\.2", ", 2",
                       sub("\\(",  " (",
                       sub("#",    ", ",
                       sub("_",    ", ",
                       sub("mer",  "-mer",
                       sub("16bp",  "16-bp",
                      gsub("\\|",  ", ",
                           rxnID.explicit))))))))
namedict.print.pub = namedict.print %>%
  mutate(rxnID.print = gsub("flank", " flanking", gsub("prep", " prep.",
                            gsub("-", "\uad",
                                 gsub("altTarg", "alt. target",
                                      gsub("altLabel", "alt. label",
                                           rxnID.print))))))
formatexplicitID = function(idname){
  # Specifically spells out the explicit (exact substitutions) part of the names
  gsub("\\|", ", ",
       gsub("\\(", " (",
            idname
       )
  )
}
idLookUp = function(miRname, explicitID = FALSE, wspe = FALSE){
  # Convert from old naming scheme
  IDdict = namedict %>%
    pull(if_else(explicitID,
                 if_else(wspe,
                         "rxnID.explicit.wspe",
                         "rxnID.explicit"),
                 "rxnID"),
         name = miR)
  return(unname(IDdict[as.character(miRname)]))
}

# Get seqs -----
seqdict = read.csv(
  paste0(
    "./",
    "miR-seqs.csv"
  ),
  stringsAsFactors = F
) %>%
  pull(seq, name = miR)
seqLookUp = function(miRname){
  # Return sequence of guide
  return(unname(seqdict[as.character(miRname)]))
}

# Get dets -----
detLookUp = function(seqs, det = c("W7", "R10", "W17", "len")){
  # Returns boolean for W7, R10, W17, or integer len
  dets = case_when(
    det == "W7"  ~ as.integer(substr(seqs, 7,  7)  %in% c("A","U")),
    det == "R10" ~ as.integer(substr(seqs, 10, 10) %in% c("A","G")),
    det == "W17" ~ as.integer(substr(seqs, 17, 17) %in% c("A","U")),
    det == "len" ~ nchar(seqs),
    TRUE ~ NA_integer_
  )
  return(dets)
}
