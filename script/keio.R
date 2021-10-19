library(data.table)
library(stringr)
library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)

#x = "Plate_mapping_october2020/PlateMappingReads_DSS3/DSS3_fasta/test2/D3114_9.csv"



x = args[1]
print(x)
plate_name = gsub(".csv","", basename(x) )


  tag_plate <- read.csv(x, header = T, sep = ",", stringsAsFactors = FALSE)

  # Correct one of the tag , as it was written in small letter 
  tag_plate$up_constant <- gsub("BSup2","BSUp2", tag_plate$up_constant)
  
  cutseq_uniq <- unique(tag_plate$cutseq)
  # empty list to save data from loop
  best_ups <- list()
  best_downs <- list()
  best_count <- list()

  # loop via the length of cutseq_uniq, and for each tag-pair find the tag-pairs with highest count.
  for (i in 1:length(cutseq_uniq)){
    if ((i %% 500) == 0){
      message("Running ", i, " out of ", length(cutseq_uniq))
    }
    cs <- cutseq_uniq[i]
    tp <- tag_plate %>%
      filter(cutseq==cs) %>%
      group_by(cutseq, up_constant, down_constant) %>%
      tally() %>%
      arrange(desc(n))
    best_row <-  tp[1, ]
    best_ups [i] <-  best_row$up_constant
    best_downs [i] <- best_row$down_constant
    best_count [i] <- best_row$n
  }

  # create a dataframe with cutseq, best_ups, best_downs, and associated read counts.
  best_df <- data.frame(cutseq=cutseq_uniq, ups = unlist(best_ups), down = unlist(best_downs), count = unlist(best_count))
  head(best_df)
  # return (best_df)

  cutseq_count <- tag_plate %>%
    group_by(cutseq) %>%
    tally()

  # barcodeMapping <- read.delim("DSS3BlPl.txt")
  # head(barcodeMapping)

  barcodeMapping <- read.csv("DSS3BlPlpJan2021.csv", header=TRUE)
  head(barcodeMapping)



  barcode_cutseq_count <- barcodeMapping %>%
    group_by(barcode) %>%
    tally() %>%
    arrange(desc(n))

  ups <-paste0("BSUp",1:16)
  downs <- paste0("BSDn",1:24,"_rc")

  # create an empty matrix and convert to dataframe with row and column name. This one for saving count.
  Mcount <- data.frame(matrix(ncol = length(downs), nrow = length(ups)))
  colnames(Mcount) <- downs
  rownames(Mcount) <- ups

  # create an empty matrix and convert to dataframe with row and column name. This one for saving count normalized by total reads in the library i.e. 1860004
  Mcount_outof_TotalLibrary <- data.frame(matrix(ncol = length(downs), nrow = length(ups)))
  colnames(Mcount_outof_TotalLibrary) <- downs
  rownames(Mcount_outof_TotalLibrary) <- ups


  # create an empty matrix and convert to dataframe with row and column name. This one for saving count, normalized by the reads observed for each unique barcode
  Mcount_outof_cutseqCount <- data.frame(matrix(ncol = length(downs), nrow = length(ups)))
  colnames(Mcount_outof_cutseqCount) <- downs
  rownames(Mcount_outof_cutseqCount) <- ups

  # create an empty matrix and convert to dataframe with row and column name. This one for saving rb20.
  Mrb <- data.frame(matrix(ncol = length(downs), nrow = length(ups)))
  colnames(Mrb) <- downs
  rownames(Mrb) <- ups

  # create an empty matrix and convert to dataframe with row and column name. This one for saving count, normalized by the reads observed for each unique barcode
  BarCount <- data.frame(matrix(ncol = length(downs), nrow = length(ups)))
  colnames(BarCount) <- downs
  rownames(BarCount) <- ups


  # create an empty matrix and convert to dataframe with row and column name. This one for saving count, normalized by the reads observed for each unique barcode
  MBarcode <- data.frame(matrix(ncol = length(downs), nrow = length(ups)))
  colnames(MBarcode ) <- downs
  rownames(MBarcode ) <- ups


  # for each tag-paris find the best cutseq
  for (r in 1:length(ups)){
    rc <- ups[r]
    for (col in 1:length(downs)){
      cc <- downs[col]
      getind <- best_df %>% filter(ups==rc & down == cc)
      if (dim(getind)[1]< 1){
        Mcount[r,col] <- NA
        Mrb[r,col] <- NA
      }else{
	max_one <- getind %>% arrange(desc(count))
        max_one_row <-  max_one[1, ]
        Mcount[r,col] <-  max_one_row$count
        getcutseq <- as.character(max_one_row$cutseq)
        Mrb[r,col] <- getcutseq
        Mcount_outof_TotalLibrary[r,col] <- max_one_row$count / dim(tag_plate)[1]  * 100
        get_count_for_cutseq <- cutseq_count %>% filter(cutseq==getcutseq) %>% pull(n)
        Mcount_outof_cutseqCount[r,col] <- max_one_row$count / get_count_for_cutseq  * 100

        # find matching getcutseq in barcode mapping in barcode_cutseq_count
        aa = barcode_cutseq_count %>%
          filter(barcode == getcutseq)
        if (dim(aa)[1]< 1){
          BarCount[r,col] <- NA
          MBarcode[r,col] <- NA
        }else{
          BarCount[r,col] <- aa$n
          MBarcode[r,col] <- as.character(aa$barcode)
        }
      }
    }
  }

  ll_list = list(Mcount,
                 Mcount_outof_TotalLibrary,
                 Mcount_outof_cutseqCount,
                 Mrb,
                 BarCount,
                 MBarcode)

  names(ll_list) <- c("Mcount",
                      "Mcount_outof_TotalLibrary",
                      "Mcount_outof_cutseqCount",
                      "Mrb",
                      "BarCount",
                      "MBarcode")

OUTPUT = ll_list

  d1 <- as.data.frame.table(as.matrix(OUTPUT$Mcount))
  d2 <- as.data.frame.table(as.matrix(OUTPUT$Mcount_outof_TotalLibrary))
  d3 <- as.data.frame.table(as.matrix(OUTPUT$Mcount_outof_cutseqCount))
  d4 <- as.data.frame.table(as.matrix(OUTPUT$Mrb))
  d5 <- as.data.frame.table(as.matrix(OUTPUT$BarCount))
  d6 <- as.data.frame.table(as.matrix(OUTPUT$MBarcode))
  dfall = Reduce(function(x,y) merge(x,y,by=c("Var1","Var2"),all=TRUE) ,list(d1,d2,d3,d4,d5,d6))
  colnames(dfall) <- c("Ups_Tag","Downs_Tag", "Mcount","Mcount_outof_TotalLibrary",
                       "Mcount_outof_cutseqCount","Mrb","BarCount","MBarcode")
  dfall <- data.frame(Plate = plate_name, dfall)

saveasName <- paste0("results/",plate_name,"_summarized.csv")
write.csv(dfall,saveasName,row.names=FALSE)



