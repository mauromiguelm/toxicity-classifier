### IMPORT PACKAGES ###
setwd('\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\drug-metadata')
library(ggplot2)
library(tidyr)
library(readxl)
library(KEGGREST)
library(dplyr)

### IMPORT map file CAS to NSC ###

NSC_to_CAS  <-  read.csv('NSC_CAS_Sept2013.csv', header = F, stringsAsFactors = F)

# IMPORT drug to target file from drugbank

drug_to_target <- read.csv('drugbank_all_target_polypeptide_ids.csv', stringsAsFactors = F)

#import kegg Cancer drugs

kEGG_Cancer_drug_set = read.csv('Drugs_MoA_Kegg_Antineoplastics.csv', stringsAsFactors = F)

# ### ORGANIZE DATA ### ---------------------------------------------------

colnames(NSC_to_CAS)[1] = 'NSC' 

colnames(NSC_to_CAS)[2] = 'CAS'

# add compounds that were not originally there

NSC_to_CAS[NSC_to_CAS$CAS == '152459-95-5', 'CAS'] = '220127-57-1' #imatinib

NSC_to_CAS[NSC_to_CAS$CAS == '50-18-0', 'CAS'] = '6055-19-2' #cyclophosphamide

NSC_to_CAS[NSC_to_CAS$CAS == '37076-68-9', 'CAS'] = '17902-23-7' #tengafur

NSC_to_CAS[NSC_to_CAS$CAS == '444731-52-6', 'CAS'] = '635702-64-6' #Pazopanib

NSC_to_CAS[NSC_to_CAS$CAS == '849217-68-1', 'CAS'] = '1140909-48-3' #Cabozantinib

Cabozantinib <- data.frame(NSC = 761068, CAS = '1140909-48-3')

BPTES <- data.frame(NSC = 798303, CAS = '314045-39-1') 

CB839 <- data.frame(NSC = 783415, CAS = '1439399-58-2')

Metformin <- data.frame(NSC = 91485, CAS = '1115-70-4') 

Panzem <- data.frame(NSC = 659853, CAS = '362-07-2')

Geldanamycin <- data.frame(NSC = 330507, CAS = '75747-14-7')

YC_1 <- data.frame(NSC = 728165, CAS = '170632-47-0')

Ponatinib <- data.frame(NSC = 758487, CAS = '1114544-31-8')

Ixazomib <- data.frame(NSC = 758254, CAS = '1239908-20-3')

Osimertinib <- data.frame(NSC = 779217, CAS = '1421373-66-1')

Alectinib <- data.frame(NSC = 764040, CAS = '1256580-46-7')

Lenvatinib <- data.frame(NSC = 755980, CAS = '417716-92-8')

Ceritinib <- data.frame(NSC = 776422, CAS = '1032900-25-6')

Belinostat <- data.frame(NSC = 758774, CAS = '414864-00-9')

Ibrutinib <- data.frame(NSC = 761910, CAS = '936563-96-1')

Capecitabine <- data.frame(NSC = 712807 , CAS = '154361-50-9')

Gemcitabine <- data.frame(NSC = 613327, CAS = '122111-03-9')

pemetrexed_disodium <- data.frame(NSC = 698037, CAS = '150399-23-8')

lapatinib <- data.frame(NSC = 745750, CAS = '388082-78-8')

erlotinib <- data.frame(NSC = 718781, CAS = '183319-69-9')

everolymus <- data.frame(NSC = 733504, CAS = '159351-69-6')

temsirolimus <- data.frame(NSC = 683864, CAS = '162635-04-3')

Bosutinib <- data.frame(NSC = 765694, CAS = '380843-75-4')

Regorafenib <- data.frame(NSC = 763932, CAS = '755037-03-7')

Panobinostat <- data.frame(NSC = 761190, CAS = '404950-80-7')

Panobinostat <- data.frame(NSC = 761190, CAS = '404950-80-7')

NSC_to_CAS <- rbind(NSC_to_CAS, Panobinostat, Regorafenib, Bosutinib, everolymus, temsirolimus, erlotinib, lapatinib, pemetrexed_disodium, Gemcitabine, Capecitabine,
                    Ibrutinib, Belinostat, Ceritinib, Lenvatinib, Alectinib, Osimertinib, Ixazomib, Ponatinib, Cabozantinib, BPTES, CB839, Metformin,
                    Panzem, Geldanamycin, YC_1)

rm(Bosutinib, everolymus, temsirolimus, erlotinib, lapatinib, pemetrexed_disodium, Gemcitabine, Capecitabine,
   Ibrutinib, Belinostat, Ceritinib, Lenvatinib, Alectinib, Osimertinib, Ixazomib, Ponatinib, Cabozantinib, Regorafenib, Panobinostat, BPTES, CB839, Metformin,
   Panzem, Geldanamycin, YC_1)

# remove drug entities not relevant for humans in drug_to_target

drug_to_target = drug_to_target[drug_to_target$Species == 'Human',]

# merging drug IDs and targets prior to merging with datasets 

colnames(drug_to_target)[13] = 'DrugBank.ID'

colnames(drug_to_target)[2] = 'Target_name'

drug_to_target = drug_to_target %>% #split the target into columns
  mutate(DrugBank.ID = strsplit(DrugBank.ID,';')) %>%
  unnest(DrugBank.ID)

drug_to_target$DrugBank.ID = gsub(" ", "", drug_to_target$DrugBank.ID) # remove whitespaces after the drug_to_target strsplit


# import kegg drug metadata from API


tmp_store_DrugBankID_CAS <- list()

for (i in 1:length(kEGG_Cancer_drug_set$KEGG)){
  
  k = kEGG_Cancer_drug_set$KEGG[i]

  tmp_getDrug <- keggGet(k)

  tmp_getDrug_2 <- tmp_getDrug[[1]]$DBLINKS
  
  tmp_store_DrugBankID_CAS[i] <- list(c(ifelse(length(tmp_getDrug_2[grepl(pattern = 'CAS', x = tmp_getDrug_2)]) == 0 , 'NA', tmp_getDrug_2[grepl(pattern = 'CAS', x = tmp_getDrug_2)]), 
                                   ifelse(length(tmp_getDrug_2[grepl(pattern = 'PubChem', x = tmp_getDrug_2)]) == 0 , 'NA', tmp_getDrug_2[grepl(pattern = 'PubChem', x = tmp_getDrug_2)])))
  }

tmp_store_DrugBankID_CAS.df <- t(data.frame(tmp_store_DrugBankID_CAS))

kEGG_Cancer_drug_set <- cbind(kEGG_Cancer_drug_set, tmp_store_DrugBankID_CAS.df)

kEGG_Cancer_drug_set$`1` <- as.character(kEGG_Cancer_drug_set$`1`)

kEGG_Cancer_drug_set$`2` <- as.character(kEGG_Cancer_drug_set$`2`)

colnames(kEGG_Cancer_drug_set)[9:10] <- c('CAS','Pubchem')

kEGG_Cancer_drug_set$CAS <- gsub(pattern = "CAS: ", replacement = "", x = kEGG_Cancer_drug_set$CAS) 

kEGG_Cancer_drug_set$Pubchem <- gsub(pattern = "PubChem: ", replacement = "", x = kEGG_Cancer_drug_set$Pubchem) 

rm(k, i, tmp_getDrug, tmp_getDrug_2, tmp_store_DrugBankID_CAS, tmp_store_DrugBankID_CAS.df)

rownames(kEGG_Cancer_drug_set) <- NULL

#fix drugs IDs in kegg metadata  

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Nilotinib', 'CAS'] = '641571-10-0'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Thioguanine', 'CAS'] = '154-42-7'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Trifluridine, mixt', 'CAS'] = '70-00-8'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Panobinostat', 'CAS'] = '404950-80-7'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Afatinib', 'CAS'] = '439081-18-2'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Bosutinib', 'CAS'] = '380843-75-4'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Regorafenib', 'CAS'] = '755037-03-7'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Alectinib', 'CAS'] = '1256580-46-7'


# include experimental drugs into kegg_Cancer drug set. These drugs were not in the Kegg pannel of anticancer drugs, and needed to me included manually

experimental_drugs <- data.frame(Class =             c("OXPHOS",                            "Fatty acid biosynthesis"     ,"Molecularly targeted agent" ,"Fatty acid oxidation"                       ,"Fatty acid oxidation"                       ,"Glutamine inhibitor" , "Glutamine inhibitor" , "OXPHOS"        , "OXPHOS"    , "OXPHOS"         , "OXPHOS"      , "OXPHOS",       "Glycolysis"),
                                 Target =            c("Mitochondrial Pyruvate Transporter","ATP-citrate lyase inhibitor" ,"mTOR kinase inhibitor"      ,"Carnitine palmitoyltransferase-1 inhibitor" ,"Carnitine palmitoyltransferase-1 inhibitor" ,"Glutaminase GLS1"    , "Glutaminase GLS1"    , "ATP Synthase"  , "Complex I" , "HIF1-2"         , "HIF1A"       , "HIF1A",        "PKM2"),
                                 Generic.name =      c("UK5099"                            ,"MEDICA_16"                   ,"Rapamycin"                  ,"Etomoxir"                                   ,"Oxfenicine"                                 , "BPTES"              , "CB-839"              , "Oligomycin A"  , "Metformin" , "Panzem-2-ME2"   , "17-AAG"      , "YC-1",         "Shikonin"),  
                                 ATC =               c(NA),
                                 KEGG =              c(NA),
                                 Products..USA. =     c(NA),
                                 Products..Japan. =  c(NA),
                                 Indications..USA. =  c(NA),
                                 CAS =               c("56396-35-1"                        ,"87272-20-6"                  ,"53123-88-9"                 ,"828934-41-4"                                ,"32462-30-9"                                 ,"314045-39-1"         , "1439399-58-2"        , "579-13-5"      , "1115-70-4" , "362-07-2" , "75747-14-7"  , "170632-47-0",  "517-89-5"),
                                 Pubchem =           c(NA))



kEGG_Cancer_drug_set <- rbind(kEGG_Cancer_drug_set, experimental_drugs)

# Include whether drug is in screen or not

kEGG_Cancer_drug_set$is_screen <- 0
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Erlotinib" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Lenvatinib" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Irinotecan" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Topotecan" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Clofarabine" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Cladribine" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Mercaptopurine" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Fluorouracil" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Gemcitabine" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Decitabine" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Pemetrexed" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Methotrexate" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Docetaxel" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Paclitaxel" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Everolimus" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Temsirolimus" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Chlormethine" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Oxaliplatin" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Asparaginase" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Omacetaxine" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "BPTES" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "CB-839" ,"is_screen"] <- 0
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Oligomycin A" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Shikonin" ,"is_screen"] <- 0
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Panzem-2-ME2" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "17-AAG" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "YC-1" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Metformin" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Trametinib" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "UK5099" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "MEDICA_16" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Rapamycin" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Etomoxir" ,"is_screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Oxfenicine" ,"is_screen"] <- 1

# add CAS to the drug metadata

drug_metadata <- right_join(NSC_to_CAS, kEGG_Cancer_drug_set, by = 'CAS')

rm(list = ls()[!ls()%in% "drug_metadata"])

#  SAVE drug metadata into a readable format

setwd("\\\\d.ethz.ch/groups/biol/sysbc/sauer_1/users/Mauro/Cell_culture_data/190310_LargeScreen/clean_data")

write.csv(x = drug_metadata, file = "drug_metadata.csv")
