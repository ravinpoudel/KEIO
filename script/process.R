library(data.table)
library(stringr)
library(tidyverse)

rm(list=ls())

# this willread in *_summarized.csv files. Output from slurm script.
ALLFiles <- list.files(path = "results/", pattern = ".csv", full.names = TRUE)
dat_csv = plyr::ldply(ALLFiles, read_csv)

## remove any entry with lacking gene level barcode MBarcode
dat_csv_filt <- dat_csv %>%
  drop_na(MBarcode)
  
  
# Barcode mapping file
# output from Random Barcode TnSeq method, which includes the 20-bp barcode and insertion location in the chromosome. 
# This output was generated following https://bitbucket.org/berkeleylab/feba/src/master/
barcodeMapping <- read.delim("DSS3PlateMappingMarch2021/DSS3fpJan21.txt")
head(barcodeMapping)


### GFF
dat_csv_filt_withGeneBarcode = left_join(dat_csv_filt, barcodeMapping, by = c("MBarcode" = "barcode"))
dat_csv_filt_withGeneBarcode$BarCount <- dat_csv_filt_withGeneBarcode$nTot
gff <- read.csv("DSS3_gff.csv", header = TRUE, check.names = FALSE)




######

pos <- dat_csv_filt_withGeneBarcode$pos
scaffold <- dat_csv_filt_withGeneBarcode$scaffold
strand <- dat_csv_filt_withGeneBarcode$strand




locus_tag <- rep(NA, length(pos))
protein_product <-rep(NA, length(pos))
locus_start <-rep(NA, length(pos))
locus_end <-rep(NA, length(pos))
protein_length <-rep(NA, length(pos))
protein_name <-rep(NA, length(pos))
gff_strand <- rep(NA, length(pos))

for ( i in 1:length(pos)){
  bb = gff[(gff$Start<=pos[i] & gff$Stop >=pos[i]), ] %>%
    filter(Accession == scaffold[i])
  # we cannot use strand informtion because strand in gff means different from strand in mapping (orientation of transposon of integration)
  # %>%filter(Strand == strand[i])
  if (dim(bb)[1] < 1){
    locus_tag[i] <- "Intergenic"
    protein_product[i] <- "Intergenic"
  } else if (dim(bb)[1]==1) {
    locus_tag[i] <- paste(bb$`Locus tag`, collapse = ", ")
    protein_product[i] <-paste(bb$`Protein product`, collapse = ", ")
    locus_start[i] <- paste(bb$Start, collapse = ", ")
    locus_end[i] <- paste(bb$Stop, collapse = ", ")
    protein_length[i] <- paste(bb$Length, collapse = ", ")
    protein_name[i] <- paste(bb$`Protein Name`, collapse = ", ")
    gff_strand[i] <- paste(bb$`Strand`, collapse = ", ")
  } else {
    locus_tag[i] <- "Multiple rb for a gene"
    protein_product[i] <- "Multiple rb for a gene"
  }
}


dat_csv_filt_withGeneBarcode["locus_tag"]<- locus_tag
dat_csv_filt_withGeneBarcode["protein_product"] <- protein_product
dat_csv_filt_withGeneBarcode["locus_start"] <- locus_start
dat_csv_filt_withGeneBarcode["locus_end"] <- locus_end
dat_csv_filt_withGeneBarcode["protein_length"] <- protein_length 
dat_csv_filt_withGeneBarcode["protein_name"] <- protein_name
dat_csv_filt_withGeneBarcode["gff_strand"] <- gff_strand

dat_csv_filt_withGeneBarcode <- dat_csv_filt_withGeneBarcode %>% 
  rename(orientation_transposon_integration = strand,
         Strand = gff_strand)


dat_csv_filt_withGeneBarcode <- dat_csv_filt_withGeneBarcode %>%
  select(-n, -nTot)

write.csv(dat_csv_filt_withGeneBarcode,"dat_csv_filt_withGeneBarcode.csv")

###########################
# https://stackoverflow.com/questions/39001057/find-replace-or-map-using-a-lookup-table-in-r
library(RSQLite)
df <- dat_csv_filt_withGeneBarcode
df["id"]<- 1:dim(df)[1]

dbname <- "big_data_mapping.db" # db to create
csvname <- df 

upcodeMap <- data.frame(
  upcode = paste0("BSUp",1:16),
  upvalue = LETTERS[1:16]
)
# build db
con <- dbConnect(SQLite(), dbname)
dbWriteTable(con, name="main", value=csvname, overwrite=TRUE)
dbWriteTable(con, name="upcodeMap", upcodeMap, overwrite=TRUE)

# join the tables
dat <- dbGetQuery(con, "SELECT * FROM main JOIN upcodeMap ON main.Ups_Tag=upcodeMap.upcode")

# finish
dbDisconnect(con)
file.remove(dbname)

## similarly for downtag

dbname <- "big_data_mapping.db" # db to create
csvname <- dat

downcodeMap <- data.frame(
  downcode = paste0("BSDn",1:24,"_rc"),
  downvalue = 1:24)

# build db
con <- dbConnect(SQLite(), dbname)
dbWriteTable(con, name="main", value=csvname, overwrite=TRUE)
dbWriteTable(con, name="downcodeMap", downcodeMap, overwrite=TRUE)

# join the tables
dat2 <- dbGetQuery(con, "SELECT * FROM main JOIN downcodeMap ON main.Downs_Tag=downcodeMap.downcode")

# finish
dbDisconnect(con)
file.remove(dbname)

dat3 <-dat2 %>%
  select(-Ups_Tag, -Downs_Tag, -upcode, -downcode, -id)

dat4 <- dat3 %>% relocate(Plate, Row = upvalue, Column = downvalue)


# write.csv(dat4, "DSS3BlPlpJan2021_keio_mappnig_results.csv", row.names = FALSE)



## READ IN new gff file 
new_gff <- read.delim('NC_003911.gff', comment.char = '#')

colnames(new_gff) <- c("Accession",
                       "GenBank",
                       "Feature",
                       "Start",
                       "Stop",
                       "Strand",
                       "What1",
                       "What2",
                       "Description")

new_gff_genes_cds = new_gff %>%
  filter(Feature == "gene" | Feature == "CDS") %>%
  filter(Feature == "gene")


get_old_locus <- function(x){
  nn = tail(unlist(strsplit(x,";")), n=1) %>%
    gsub('old_locus_tag=','',.) %>%
    gsub('Name=','',.) %>%
    gsub('locus_tag=','',.)
  
  return(nn)
}

get_locus <- function(x){
  mm = head(unlist(strsplit(x,";")), n=1) %>%
    gsub('ID=','',.) %>%
    gsub('Name=','',.) %>%
    gsub('locus_tag=','',.)
  
  return(mm)
}

new_gff_genes_cds["old_locus_tag"] <- unlist(lapply(new_gff_genes_cds$Description, get_old_locus))
new_gff_genes_cds["locus_tag"] <- unlist(lapply(new_gff_genes_cds$Description, get_locus))

new_gff_genes_cds <- new_gff_genes_cds %>%
  select(Description,old_locus_tag,locus_tag)


write.csv(new_gff_genes_cds,"new_gff_genes_cds.csv")
dat4_updated <- left_join(dat4, new_gff_genes_cds, by = c("locus_tag" = "locus_tag"))


dat4_updated_col_arranged <- dat4_updated %>%
  select("Plate","Row","Column","Mcount","Mcount_outof_TotalLibrary","Mcount_outof_cutseqCount",
         "Mrb","BarCount","MBarcode","rcbarcode","scaffold","orientation_transposon_integration",
         "pos","locus_tag","old_locus_tag","protein_product","locus_start","locus_end","protein_length","protein_name",
         "Strand","Description")

write.csv(dat4_updated_col_arranged,"DSS3fpJan21_keio_mappnig_results.csv", row.names = FALSE)


###
############ Combining results from three runs#########
# Previously I was running three runs... if you have single run then you can just use output from single run-- which in this case is DSS3fpJan21_keio_mappnig_results.csv
rm(list=ls())
#run1 <- read.csv("DSS3_all_results/DSS3BlPlpJan2021_keio_mappnig_results.csv", header = TRUE)
#run2 <- read.csv("DSS3_all_results/DSS3PlateMappingMarch2021_DSS3fpJan21_keio_mappnig_results.csv", header = TRUE)

# Pick the highest Mcount outof cutseqCount, and then set a threshold for the lowest allowed Mcount outof TotalLibrary,
# unless there's only one or two mutants in the gene, then it should be included regardless of counts. 
# all <- rbind(run1, run2)

all <- read.csv("DSS3fpJan21_keio_mappnig_results.csv", header = TRUE)

all["Mcount_outof_cutseqCount_x_Mcount_outof_TotalLibrary"] <- all$Mcount_outof_TotalLibrary * all$Mcount_outof_cutseqCount
locus <- unique(all$locus_tag)
all["TotalCutSeq"] <- (all$Mcount * 100)/ all$Mcount_outof_cutseqCount
all["TotalReadsbyPlate"] <- (all$Mcount * 100)/ all$Mcount_outof_TotalLibrary
locus_unique <- locus[locus != "Intergenic"] %>% sort()


### Distance from Start Codon
all$locus_end <- as.numeric(all$locus_end)
all$locus_start <- as.numeric(all$locus_start)
all$protein_length <- as.numeric(all$protein_length)
all["DistancefromStartCodon"] <- as.numeric(ifelse(all$Strand=="+", abs(all$locus_start-all$pos), abs(all$pos-all$locus_end)))
all["idf_orf"] <- as.numeric(round(all$DistancefromStartCodon/(all$protein_length*3),2))


all <- all %>% dplyr::rename(
  `Insertion Distance as frac of ORF` = idf_orf,
  `Tn orientation` = orientation_transposon_integration,
  `Insertion Location`=pos,
  Protein = protein_product,
  `Protein Length`= protein_length,
  `Protein Name` = protein_name,
  McountTcountScore = Mcount_outof_cutseqCount_x_Mcount_outof_TotalLibrary)


all["RawScore"] <- ifelse(all$`Insertion Distance as frac of ORF` >= 0.05 & all$`Insertion Distance as frac of ORF` <=0.70, 10, 5)
all["Score"] <- all$McountTcountScore * all$RawScore

all <- all %>%
  select("Score","RawScore","McountTcountScore","Insertion Distance as frac of ORF","Plate","Row","Column","Mrb","BarCount","Mcount","MBarcode","rcbarcode",
         "scaffold","Tn orientation","Insertion Location","locus_tag","old_locus_tag","Protein","Protein Length","Protein Name",
         "Strand","DistancefromStartCodon"
  )



# write.csv(all,"all.csv")
# 
# all <- all %>%
#   relocate("Mcount","Mcount_outof_TotalLibrary","Mcount_outof_cutseqCount", .after = "Description")


best_clone <- list()

for (i in 1:length(locus_unique)){
  bylocus = all %>% filter(locus_tag == locus_unique[i]) %>% arrange(desc(Score))
  best_clone[[i]] <- bylocus[1, ]
}



best_clone_4locus <- do.call(rbind, best_clone)

write.csv(best_clone_4locus,"best_clone_4locus.csv")


