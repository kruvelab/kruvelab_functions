
# All functions in kruvelab

library(tidyverse)
library(rcdklibs)
library(rcdk)
library(stringr) #Anneli thinks that this also comes with tidyverse
#library(berryFunctions) #what is this?
library(rjson)
library(glue)
library(readxl)
#library(MS2Tox)
library(dplyr) #ANneli thinks this is in tidyverse
library(rlist)
library(Rdisop)
library(tibble) #Anneli thinks that this is in tidyverse
library(readr)
library(magrittr)
library(xgboost)
library(tidyselect) #is this needed?
library(stringr)
library(caret)

# By start of April:

# 1) re-check all functions - updates/new fns?
# 2) re-structure? + read_me in github

# GENERAL about functions:
# all dataframe functions --> vectorize!
# naming - lowercase and underscores (except names)!!! 
# name same for same variable everywhere!!!
# specifying input types in all function
# which libraries are actually needed?



# Similarity - Anneli
## all similarity functions unification (applicable for both FP and mass-spec? input of what it is defined in function?)
## similarites function from Harry?
# binning function - Yvonne

# descriptors - Helen
## one general function for calculating all descriptors and fingerprints
# PaDEL - calculating for individual smiles vs dataframe? - Anneli
# Eluent composition - Anneli

# MassBank - Yvonne

# modelling functions - Helen
## IE functions from Anneli & Helen (applying the models)
## linearity - Helen

# Querying information (based on Drew's getMetadata) - Drew
## one general function for querying metadata (from: SMILES; to: Standardized_smiles etc)

# MSDIAL workflow functions - Anneli and Helen
## MSDIAL -> .txt alignment file -> create .ms files -> function to put files to SIRIUS -> summary table of SIRIUS results
## + merging with suspect list - Anneli and Helen
## centering out! 
## update SIRIUS functions to be compatible with SIRIUS version 5
## read in ID from the back and remove non-numbers for when name has MS-DIAL ID nr

# progenesis - Harry

# toxicity functions from Pille's work (that is not in MS2Tox) - Yvonne & Anneli (scraping)







# **************************************
# ---- Eluent composition functions ----
# **************************************

organic_percentage = function(gradient_dataframe = tibble(),
                             ret_time = numeric(),
                             time_column_name,
                             organic_percent_column_name){
  print(time_column_name)
  print(deparse(substitute(time_column_name)))
  time = gradient_dataframe[[deparse(substitute(time_column_name))]]
  #time = gradient_dataframe[[time_column_name))]]
  print(time)
  B = gradient_dataframe[[deparse(substitute(organic_percent_column_name))]]
  ApproxFun = approxfun(x = time, 
                        y = B)
  organic = ApproxFun(ret_time)
  return(organic)
}
# organic_percentage = Vectorize(organic_percentage)

polarity_index = function(organic_percent = numeric(), 
                          organic_modifier = character()){
  polarity_index = case_when(
    organic_modifier == "MeCN" ~ (organic_percent/100)*5.1+((100-organic_percent)/100)*10.2,
    organic_modifier == "MeOH" ~ (organic_percent/100)*5.1+((100-organic_percent)/100)*10.2)
  return(polarity_index)
}
# polarity_index = Vectorize(polarity_index)

surface_tension = function(organic_percentage = numeric(),
                          organic_modifier = "MeCN"){
  surface_tension = case_when(
    organic_modifier == "MeCN" ~ 71.76-2.906*71.76*(organic_percentage/100)+(7.138*27.86+2.906*71.76-71.76)*(organic_percentage/100)^2+(27.86-7.138*27.86)*(organic_percentage/100)^3,
    organic_modifier == "MeOH" ~ 71.76-2.245*71.76*(organic_percentage/100)+(5.625*22.12+2.245*71.76-71.76)*(organic_percentage/100)^2+(22.12-5.625*22.12)*(organic_percentage/100)^3)
  return(surface_tension)
}
# surface_tension = Vectorize(surface_tension)

viscosity = function(organic_percentage = numeric(),
                      organic_modifier = "MeCN"){
  viscosity = case_when(
    organic_modifier == "MeCN" ~ (-0.000103849885417527)*organic_percentage^2+0.00435719229180079*organic_percentage+0.884232851261593,
    organic_modifier == "MeOH" ~ (-0.00035908)*organic_percentage^2+0.031972067*organic_percentage+0.90273943)
  return(viscosity)
}
# viscosity = Vectorize(viscosity)

add_mobile_phase_composition = function(data = tibble(),
                                        ret_time,
                                        organic_modifier = "MeCN",
                                        pH_aq = 7.0,
                                        gradient_file_name = character(),
                                        time_column_name,
                                        organic_percent_column_name) {
  
  gradient = read_delim(gradient_file_name,
                        col_types = cols())
  print(gradient)
  print(time_column_name)
  data = data %>%
    mutate(organic_percent = organic_percentage(gradient, 
                                                ret_time,
                                                time_column_name,
                                                organic_percent_column_name),
           viscosity = viscosity(organic_percent, 
                                 organic_modifier),
           surface_tension = surfacetension(organic_percent, 
                                            organic_modifier),
           polarity_index = polarityindex(organic_percent, 
                                          organic_modifier), 
           organic_modifier = organic_modifier,
           pH_aq = pH_aq)
  return(data)
}


# *****************************************
# ---- Functions to query information ----
# *****************************************

standardize_smiles <- function(smiles_string){
  molecule <- NULL
  tryCatch({
    molecule <- parse.smiles(smiles_string)[[1]]
  }, error = function(e){
    print(paste("Parsing failed for:",smiles_string))
    print(e)
  })
  if (!is.null(molecule)){
    std_smiles <- get.smiles(molecule,
                             flavor = smiles.flavors("Canonical"))
  } else {
    std_smiles <- NULL
  }
  return(std_smiles)
}
standardize_smiles <- Vectorize(standardize_smiles) 

smiles_to_standardinchikey <- function(smiles_string){
  standardinchikey <- NULL
  tryCatch({
    standardinchikey = webchem::cir_query(smiles_string, "stdinchikey")[[1]]
  }, error = function(e){
    standardinchikey = "NA"
  })
  return(standardinchikey)
}
smiles_to_standardinchikey <- Vectorize(smiles_to_standardinchikey)

smiles_to_standardinchi <- function(smiles_string){
  standardinchi <- NULL
  tryCatch({
    standardinchi = webchem::cir_query(smiles_string, "stdinchi")[[1]]
  }, error = function(e){
    standardinchi = "NA"
  })
  return(standardinchi)
}
smiles_to_standardinchi <- Vectorize(smiles_to_standardinchi)



# *********************************************
# ---- Functions for pretreating .ms files: CENTERING OUT! ----
# *********************************************

centroid <- function(data, 
                     mz_filter,
                     background_thershold = 0.5) {
  data <- data %>%
    filter(mass < mz_filter)
  
  noise <- data %>%
    group_by(mass) %>%
    mutate(mass_difect = mass - trunc(mass)) %>%
    ungroup() %>%
    filter(mass_difect > 0.5 & mass_difect < 0.8)
  
  background <- quantile(data$Intensity, background_thershold)
  
  data_processed <- tibble()
  
  for (i in 3:(nrow(data)-2)) {
    if (data$Intensity[i] > background){ #int > limiit
      if (data$Intensity[i-2] < data$Intensity[i-1]) { #5 punkti
        if (data$Intensity[i+2] < data$Intensity[i+1]){
          if (data$Intensity[i-1] < data$Intensity[i]) {
            if (data$Intensity[i+1] < data$Intensity[i]){
              data_processed <- data_processed %>%
                bind_rows(data[i,])
            }
          }
        }
      }
    }
  }
  return(data_processed)
}

# centering MS file 

centeringMSFile <- function(folder,
                            ms_maxinum = 20,
                            background_thershold = 0.5) {
  setwd(folder)
  filelist <- as.data.frame(dir(folder, all.files = TRUE, recursive = TRUE, pattern = "*.ms"))
  colnames(filelist) <- "filename"
  file <- filelist
  
  dir.create("centering")
  
  for(i in 1:length(file$filename)) {
    
    setwd(folder)
    print(file$filename[i])
    fileConnection <- file$filename[i]
    record <- readLines(fileConnection)
    #close(fileConnection)
    
    id <- substring(grep('>compound', record, value = TRUE, fixed = TRUE),11)
    parent_ms <- as.numeric(substring(grep('>parentmass', record, value = TRUE, fixed = TRUE),12))
    filter_ms <- parent_ms + ms_maxinum   #Mass Spectra Threshold
    ionization <- substring(grep('>ionization', record, value = TRUE, fixed = TRUE),13)
    #ms1peaks <- substring(grep('>ms1peaks', record, value = TRUE, fixed = TRUE),0)
    
    #Get ms1 and ms2 data indexes
    for (j in 1:length(record)) {
      if (record[j] == ">ms1") {
        break
      }
    }
    for (k in 1:length(record)) {
      if (record[k] == ">ms2") {
        break
      }
    }
    
    ms1_data = tibble()
    for(index in (j+1):(k-2)) {
      split_str <- strsplit(record[index], split = " ") 
      ms1_data <- ms1_data %>%
        bind_rows(tibble("mass" = split_str[[1]][1],
                         "Intensity" = split_str[[1]][2]))
    }
    
    ms1_data <- ms1_data %>%
      mutate(mass = as.numeric(mass),
             Intensity = as.numeric(Intensity))
    
    ms2_data = tibble()
    for(index in (k+1):length(record)) {
      split_str <- strsplit(record[index], split = " ") 
      ms2_data <- ms2_data %>%
        bind_rows(tibble("mass" = split_str[[1]][1],
                         "Intensity" = split_str[[1]][2]))
    }
    ms2_data <- ms2_data %>%
      mutate(mass = as.numeric(mass),
             Intensity = as.numeric(Intensity))
    
    
    ms1_data <- centroid(ms1_data, filter_ms, background_thershold)
    ms2_data <- centroid(ms2_data, filter_ms, background_thershold)
    
    ms1_file <- rbind(c(">ms1", ""), ms1_data)
    ms1_file <- rbind(c("", ""), ms1_file)
    ms1_file <- rbind(c('>ionization',ionization), ms1_file)
    ms1_file <- rbind(c('>parentmass',parent_ms), ms1_file)
    ms1_file <- rbind(c('>compound',id), ms1_file)
    
    ms2_file <- rbind(c(">ms2", ""), ms2_data)
    ms2_file <- rbind(c("", ""), ms2_file)
    
    file_toExport <- tibble()
    file_toExport <- rbind(ms1_file, ms2_file)
    
    
    names <- paste("centering_", id,".ms", sep = "")
    
    wd_name <- paste(folder, "/centering", sep = "")
    setwd(wd_name)
    
    write_delim(file_toExport, path = names, delim = " ", col_names = FALSE)
    
  }
}


# deleting .ms files without MS2 
delete_msfiles_withoutMS2 <- function(directories_list) {
  for (file in 1:length(directories_list)) {
  lines <- readLines(directories_list[file])
  for (i in 1:length(lines)) {
    if (lines[i] == ">ms2") {
      if (is.na(lines[i+1])){
        file.remove(directories_list[file])
      }
    }
  }
}
  
}


# ************************************************************
# ---- All molecular descriptors/fingerprints calculation functions ----
# ************************************************************

molecularmass <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  #calcuate molecular weight
  MW <- MolecularWeight(formula = ListFormula(formula))
  return(MW)
}

smiles_to_formula <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  return(formula)
}

isotopedistribution <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  
  # Chemical formula to isotope distribution
  data(isotopes)
  pattern<-isopattern(isotopes,
                      formula,
                      threshold=0.1,
                      plotit=FALSE,
                      charge=FALSE,
                      algo=1)
  isotopes <- as.data.frame(pattern[[1]])
  isotope_dist <- as.numeric(sum(isotopes$abundance))
  return(isotope_dist)
}

# is this computational? Is this one good to use?
fn_logP <- function(SMILES){       
  mol <- parse.smiles(SMILES)[[1]]
  convert.implicit.to.explicit(mol)
  get.tpsa(mol)
  logP <- get.xlogp(mol)
  return(logP)
}


source_python("C:/Users/HelenSepman/OneDrive - Kruvelab/IS_surrogate_IE/descriptors/Mordred_python.py")
np = import("numpy", convert = FALSE)
pandas = import("pandas", convert = FALSE)
rdkit = import("rdkit", conver = FALSE)
mordred = import("mordred", convert = FALSE)


PaDEL_original = function(standards) {
  SMILES_list = standards %>% 
    select(SMILES) %>% 
    unique() %>%
    na.omit() %>%
    mutate(Name = row_number())
  standards = standards %>%
    left_join(SMILES_list)
  
  write_delim(SMILES_list %>% select(SMILES) %>% unique(),
              "SMILES.smi",
              col_names = FALSE)
  
  command = "java -jar PaDEL-Descriptor/PaDEL-Descriptor.jar -dir" #file name where they will be calculated
  command_final = paste(command, "SMILES.smi", "-file", "descs_calc.csv", "-2d", sep =" ") #makes text for command prompt
  javaOutput = system(command_final, intern = TRUE) #goes into commant prompt
  #PaDEL saves the descriptors to the local folder
  descs = read_delim("descs_calc.csv",
                     delim = ",",
                     col_names = TRUE)
  
  descs = descs %>%
    mutate_all(~replace(., str_detect(., "Infinity"), as.numeric(0))) %>%
    group_by(Name) %>%
    mutate(Name = str_split(Name, pattern = "_")[[1]][length(str_split(Name, pattern = "_")[[1]])]) %>%
    ungroup() %>%
    mutate(Name = as.integer(Name)) 
  
  cols <- names(descs)
  descs[cols] <- lapply(descs[cols], as.numeric)
  
  
  descs = descs %>%
    left_join(standards) %>%
    select(colnames(standards), everything()) %>%
    select(-Name)
  
  return(descs)
}

# Mordred
Mordred_calculation_R = function(SMILES) {
  descriptors = py_to_r(Mordred_calculation(SMILES))
  return(descriptors)
}
unlisting <- function(listed) {
  if (is.numeric(listed[[1]])){
    new <- listed[[1]]
  } else{
    new <- NA
  }
  return(new)
} ## Ilmselt saaks selle teisiti kirjutada, nii nagu Veronikale tegime
Mordred_descriptors = function(SMILES_list) {
  descriptors = tibble()
  for (SMILES in SMILES_list$SMILES) {
    descriptor = Mordred_calculation_R(SMILES)
    descriptor_this = tibble(value = descriptor$`_values`)
    if (dim(descriptors)[1] == 0) {
      descriptors = descriptor_this
    } else {
      descriptors = descriptors %>%
        bind_cols(descriptor_this)
    }
  }
  names_Mordred = read_delim("Mordred_names.txt",
                             delim = ",",
                             col_names = FALSE) %>%
    rename(desc_name = X1)
  descriptors = descriptors %>%
    bind_cols(names_Mordred) 
  descriptors = descriptors%>%
    gather(key = "compound", value = "value", -dim(descriptors)[2]) %>%
    spread(desc_name, value) %>%
    bind_cols(tibble("SMILES" = SMILES_list$SMILES)) %>%
    select(-compound) %>%
    select(SMILES, everything())
  
  descriptors = descriptors%>%
    group_by(SMILES) %>%
    mutate_if(is.list, funs(unlisting(.)))%>%
    ungroup()
  return(descriptors)
}


# Morgan-2 
source_python("C:/Users/HelenSepman/OneDrive - Kruvelab/IS_surrogate_IE/descriptors/Morgan2_python.py")
Morgan2_calculation_R = function(SMILES) {#, radius = 2, nBits = 1024) {
  descriptors = Morgan2_calculation(SMILES)#, radius, nBits))
  return(descriptors)
}
Morgan2_descriptors = function(SMILES_list) {
  descriptors = tibble()
  for (SMILES in SMILES_list$SMILES) {
    descriptor = Morgan2_calculation_R(SMILES)
    descriptor_this = tibble(descriptor)
    if (dim(descriptors)[1] == 0) {
      descriptors = descriptor_this
    } else {
      descriptors = descriptors %>%
        bind_cols(descriptor_this)
    }
  }
  names_Morgan2 = tibble(desc_name = 1:1024) %>%
    mutate(desc_name = str_c("Morgan2_", desc_name))
  descriptors = descriptors %>%
    bind_cols(names_Morgan2) 
  descriptors = descriptors%>%
    gather(key = "compound", value = "value", -dim(descriptors)[2]) %>%
    spread(desc_name, value) %>%
    bind_cols(tibble("SMILES" = SMILES_list$SMILES)) %>%
    select(-compound) %>%
    select(SMILES, everything())
  return(descriptors)
}



# Cannot calculate first 18 descriptors as SMILES do not have 3D info
rcdk_fingerprints <- function(SMILES_list) {
  dn <- get.desc.names(type = "all")
  mols <- parse.smiles(SMILES_list$SMILES)
  allDescs <- eval.desc(mols, dn)
  allDescs <- allDescs %>%
    rownames_to_column(var = "SMILES")
  return(allDescs)
}

substructure_FP <- function(compoundslist, folderwithFPinfo){
  
  setwd(folderwithFPinfo)
  
  list <- read_delim("substructureCDK_fingerprints_names.txt", delim = ",") %>% ## Only used to get the nr of FP but should be able to get from rcdk, fix it later
    mutate(row = as.numeric(rowCDK)) %>%
    select(row)
  
  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column() %>%
    filter(rowname != "row")
  
  error_comp <- tibble()
  
  for (n in 1:length(compoundslist$SMILES)){
    #col <- compoundslist$id[[n]]
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      substr_fingerprints <- get.fingerprint(mol2,
                                             type = "substructure")
      
      table <- as.data.frame(substr_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      substr_table <- substr_table %>%
        add_row(newrow)
    }
    
    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }
  
  names(substr_table) <- c("SMILES", as.vector(list$row))
  substr_table[is.na(substr_table)] <- 0
  
  return(substr_table)
}

MACCS_FP <- function(compoundslist, folderwithFPinfo){
  
  setwd(folderwithFPinfo)
  
  list <- read_delim("MACCS_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowMAC)) %>%
    select(row)
  
  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column() %>%
    filter(rowname != "row")
  
  error_comp <- tibble()
  
  for (n in 1:length(compoundslist$SMILES)){
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      maccs_fingerprints <- get.fingerprint(mol2,
                                            type = "maccs")
      
      table <- as.data.frame(maccs_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      substr_table <- substr_table %>%
        add_row(newrow)
    }
    
    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }
  
  names(substr_table) <- c("SMILES", as.vector(list$row))
  substr_table[is.na(substr_table)] <- 0
  
  return(substr_table)
}

KlekRoth_FP <- function(compoundslist, folderwithFPinfo){
  
  setwd(folderwithFPinfo)
  
  list <- read_delim("KlekotaRoth_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowKR)) %>%
    select(row)
  
  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column() %>%
    filter(rowname != "row")
  
  error_comp <- tibble()
  
  for (n in 1:length(compoundslist$SMILES)){
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      kr_fingerprints <- get.fingerprint(mol2,
                                         type = "kr")
      
      table <- as.data.frame(kr_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      substr_table <- substr_table %>%
        add_row(newrow)
    }
    
    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }
  
  names(substr_table) <- c("SMILES", as.vector(list$row))
  substr_table[is.na(substr_table)] <- 0
  
  return(substr_table)
}

rings_FP <- function(compoundslist, folderwithFPinfo){
  
  setwd(folderwithFPinfo)
  
  list <- read_delim("ringsystem_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowRing)) %>%
    select(row)
  
  smartslistring <- as.vector(read_delim("ringsystem_fingerprints_names.txt", delim = ",")$description_ring)
  
  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column() %>%
    filter(rowname != "row")
  
  error_comp <- tibble()
  
  for (n in 1:length(compoundslist$SMILES)){
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      ring_fingerprints <- get.fingerprint(mol2,
                                           type = "substructure",
                                           substructure.pattern = smartslistring)
      
      table <- as.data.frame(ring_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      substr_table <- substr_table %>%
        add_row(newrow)
    }
    
    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }
  
  names(substr_table) <- c("SMILES", as.vector(list$row))
  substr_table[is.na(substr_table)] <- 0
  
  return(substr_table)
}

PubChem_FP <- function(compoundslist, folderwithFPinfo){
  
  setwd(folderwithFPinfo)
  
  list <- read_delim("pubchem_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowPC)) %>%
    select(row)
  
  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column() %>%
    filter(rowname != "row")
  
  error_comp <- tibble()
  
  for (n in 1:length(compoundslist$SMILES)){
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      pc_fingerprints <- get.fingerprint(mol2,
                                         type = "pubchem")
      
      table <- as.data.frame(pc_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      substr_table <- substr_table %>%
        add_row(newrow)
    }
    
    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }
  
  names(substr_table) <- c("SMILES", as.vector(list$row))
  substr_table[is.na(substr_table)] <- 0
  
  return(substr_table)
}

CustomSMARTS_FP <- function(compoundslist, folderwithFPinfo){
  
  setwd(folderwithFPinfo)
  
  list <- read_delim("custommadeSMARTS_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowCustom)) %>%
    select(row)
  
  smartslist <- as.vector(read_delim("custommadeSMARTS_fingerprints_names.txt", delim = ",")$description_cust)
  
  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column() %>%
    filter(rowname != "row")
  
  error_comp <- tibble()
  
  for (n in 1:length(compoundslist$SMILES)){
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      custom_fingerprints <- get.fingerprint(mol2,
                                             type = "substructure",
                                             substructure.pattern = smartslist)
      
      table <- as.data.frame(custom_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      substr_table <- substr_table %>%
        add_row(newrow)
    }
    
    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }
  
  names(substr_table) <- c("SMILES", as.vector(list$row))
  substr_table[is.na(substr_table)] <- 0
  
  return(substr_table)
}

SIRIUS_FP <- function(compoundslist, folderwithFPinfo){
  setwd(folderwithFPinfo)
  
  #substructure fingerprints
  list <- read_delim("substructureCDK_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowCDK)) %>%
    select(row)
  
  col_names_substr <- read_delim("substructureCDK_fingerprints_names.txt", delim = ",") %>%
    select(sir_absol_ind)
  col_names_substr <- as.vector(paste("Un", col_names_substr$sir_absol_ind, sep = ""))
  
  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  error_comp <- tibble()
  
  for (n in 1:length(compoundslist$SMILES)){
    col <- compoundslist$id[[n]]
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      substr_fingerprints <- get.fingerprint(mol2,
                                             type = "substructure")
      
      table <- as.data.frame(substr_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      substr_table <- substr_table %>%
        add_row(newrow)
    }
    
    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }
  
  substr_table[is.na(substr_table)] <- 0
  
  substr_table <- substr_table %>%
    filter(rowname != "row")
  
  names(substr_table) <- c("SMILES", col_names_substr)
  
  #MACCS fingerprints
  list <- read_delim("MACCS_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowMAC)) %>%
    select(row)
  
  col_names_maccs <- read_delim("MACCS_fingerprints_names.txt", delim = ",") %>%
    select(sir_absol_ind)
  col_names_maccs <- as.vector(paste("Un", col_names_maccs$sir_absol_ind, sep = ""))
  
  maccs_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist$SMILES)){
    col <- compoundslist$id[[n]]
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      maccs_fingerprints <- get.fingerprint(mol2,
                                            type = "maccs")
      
      table <- as.data.frame(maccs_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      maccs_table <- maccs_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  maccs_table[is.na(maccs_table)] <- 0
  
  maccs_table <- maccs_table %>%
    filter(rowname != "row")
  
  names(maccs_table) <- c("SMILES", col_names_maccs)
  
  
  #pubchem fingerprints
  list <- read_delim("pubchem_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowPC)) %>%
    select(row)
  
  col_names_pubchem <- read_delim("pubchem_fingerprints_names.txt", delim = ",") %>%
    select(sir_absol_ind)
  col_names_pubchem <- as.vector(paste("Un", col_names_pubchem$sir_absol_ind, sep = ""))
  
  pubchem_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist$SMILES)){
    col <- compoundslist$id[[n]]
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      pubchem_fingerprints <- get.fingerprint(mol2,
                                              type = "pubchem")
      
      table <- as.data.frame(pubchem_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      pubchem_table <- pubchem_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  pubchem_table[is.na(pubchem_table)] <- 0
  
  pubchem_table <- pubchem_table %>%
    filter(rowname != "row")
  
  names(pubchem_table) <- c("SMILES", col_names_pubchem)
  
  #KlekRoth fingerprints
  list <- read_delim("KlekotaRoth_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowKR)) %>%
    select(row)
  
  col_names_KlekotaRoth <- read_delim("KlekotaRoth_fingerprints_names.txt", delim = ",") %>%
    select(sir_absol_ind)
  col_names_KlekotaRoth <- as.vector(paste("Un", col_names_KlekotaRoth$sir_absol_ind, sep = ""))
  
  KlekotaRoth_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist$SMILES)){
    col <- compoundslist$id[[n]]
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      KlekotaRoth_fingerprints <- get.fingerprint(mol2,
                                                  type = "kr")
      
      table <- as.data.frame(KlekotaRoth_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      KlekotaRoth_table <- KlekotaRoth_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  KlekotaRoth_table[is.na(KlekotaRoth_table)] <- 0
  
  KlekotaRoth_table <- KlekotaRoth_table %>%
    filter(rowname != "row")
  
  names(KlekotaRoth_table) <- c("SMILES", col_names_KlekotaRoth)
  
  
  #custommadeSMARTS
  list <- read_delim("custommadeSMARTS_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowCustom)) %>%
    select(row)
  
  smartslist <- as.vector(read_delim("custommadeSMARTS_fingerprints_names.txt", delim = ",")$description_cust)
  
  col_names_custom <- read_delim("custommadeSMARTS_fingerprints_names.txt", delim = ",") %>%
    select(sir_absol_ind)
  col_names_custom <- as.vector(paste("Un", col_names_custom$sir_absol_ind, sep = ""))
  
  custom_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist$SMILES)){
    col <- compoundslist$id[[n]]
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      custom_fingerprints <- get.fingerprint(mol2,
                                             type = "substructure",
                                             substructure.pattern = smartslist)
      
      table <- as.data.frame(custom_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      custom_table <- custom_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  custom_table[is.na(custom_table)] <- 0
  
  custom_table <- custom_table %>%
    filter(rowname != "row")
  
  names(custom_table) <- c("SMILES", col_names_custom)
  
  
  #ringsystems
  list <- read_delim("ringsystem_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowRing)) %>%
    select(row)
  
  smartslistring <- as.vector(read_delim("ringsystem_fingerprints_names.txt", delim = ",")$description_ring)
  
  col_names_ring <- read_delim("ringsystem_fingerprints_names.txt", delim = ",") %>%
    select(sir_absol_ind)
  col_names_ring <- as.vector(paste("Un", col_names_ring$sir_absol_ind, sep = ""))
  
  ring_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist$SMILES)){
    col <- compoundslist$id[[n]]
    SMILES <- compoundslist$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      ring_fingerprints <- get.fingerprint(mol2,
                                           type = "substructure",
                                           substructure.pattern = smartslistring)
      
      table <- as.data.frame(ring_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      datarow <- list %>%
        left_join(table) %>%
        select(-row)
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      ring_table <- ring_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  ring_table[is.na(ring_table)] <- 0
  
  ring_table <- ring_table %>%
    filter(rowname != "row")
  
  names(ring_table) <- c("SMILES", col_names_ring)
  
  
  
  #Combining all FP together
  final_table <- substr_table %>%
    left_join(maccs_table) %>%
    left_join(pubchem_table) %>%
    left_join(KlekotaRoth_table) %>%
    left_join(custom_table) %>%
    left_join(ring_table) %>%
    
    return(final_table)
  #return(error_comp)
}

# ---- NB! filter out the FPs that SIRIUS actually does not calculate!! ----
sirius_true_FP_pos <- read_delim("C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/toxicity/csi_fingerid.tsv",
                                 delim = "\t")
sirius_true_FP_neg <- read_delim("C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/toxicity/csi_fingerid_neg.tsv",
                                 delim = "\t")


col_names_calcSirius_pos <- as.vector(paste("Un", sirius_true_FP_pos$absoluteIndex, sep = ""))
col_names_calcSirius_neg <- as.vector(paste("Un", sirius_true_FP_neg$absoluteIndex, sep = ""))
col_names_calcSmiles <- colnames(SIRIUS_descriptors)
intersection <- intersect(col_names_calcSirius_pos, col_names_calcSirius_neg)
intersection <- intersect(intersection, col_names_calcSmiles)
intersection <- append("SMILES", intersection)

SIRIUS_descriptors <- SIRIUS_descriptors[names(SIRIUS_descriptors)[names(SIRIUS_descriptors) %in% intersection] ]



# ***************************************************
# ---- Similarity calculations (dice and cosine) ----
# ***************************************************

# ---- All functions ----

# Existing binary fingerprints
fingerprints_present = function(row_of_a_dataset) {
  present = row_of_a_dataset %>% 
    select_if(function(col) mean(col) == 1) %>% 
    colnames()
  return(present)
}
Vectorize(fingerprints_present)

Dice_score = function(fingerprint1, fingerprint2) {
  TP = intersect(unlist(fingerprint1), 
                 unlist(fingerprint2))
  Dice_score = 2*length(TP)/(length(unlist(fingerprint1)) + 
                               length(unlist(fingerprint2)))
  return(Dice_score)
}
Vectorize(Dice_score)

# same length fingerprints?? what about mass spectra?
similarity = function(fingerprint1, fingerprint2, is_binary) {
  if(is_binary){
    return(Dice_score(fingerprint1, fingerprint2))
  } else {
    return(as.numeric(cosine(unlist(fingerprint1), unlist(fingerprint2))))
  }
}
Vectorize(similarity)

similarities_table <- function(list_of_molecules_descriptors){
  
  descriptors <- list_of_molecules_descriptors %>%
    select_if(is.numeric) 
  
  # binary or continuous?
  is_binary = all(apply(descriptors, 2, function(x) { all(na.omit(x) %in% 0:1)}))
  
  if(!is_binary) {
    descriptors <- descriptors %>%
      select(-c(nearZeroVar(descriptors, 
                            freqCut = 100/0))) %>%
      scale() %>%
      as_tibble()
  } 
  
  data_fingerprints = tibble()
  for (i in 1:dim(descriptors)[1]) {
    this_data = descriptors[i,] %>%
      mutate(fingerprints = case_when(
        is_binary ~ I(list(fingerprints_present(.))),
        TRUE ~ I(list(as.numeric(.)))))
    data_fingerprints = data_fingerprints %>%
      bind_rows(this_data)
  }
  
  data_fingerprints = data_fingerprints %>%
    select(fingerprints) %>%
    bind_cols(list_of_molecules_descriptors %>%
                select(SMILES))
  
  all_combinations = data_fingerprints %>%
    left_join(data_fingerprints,
              by = character()) %>%
    filter(SMILES.x != SMILES.y)
  
  all_combinations = all_combinations %>%
    group_by(SMILES.x, SMILES.y) %>%
    mutate(similarity = similarity(fingerprints.x, fingerprints.y, is_binary)) %>%
    ungroup()
  
  return(all_combinations)
}


# ********************************
# ---- Modelling functions ----
# ********************************

cleaning_descriptors = function(descs,
                                nearZeroVar_freqCut = 80/20,
                                highlyCorrelated_cutoff = 0.75) {
  
  # Removing columns with missing values
  descs = descs %>%
    select_if(~ sum(is.na(.))< 10,) %>%
    drop_na()
  
  # new_descs <- descs
  # new_descs[, sapply(new_descs, class) != "environment"]
  
  smiles = descs %>%
    select(SMILES)
  
  descs = descs %>%
    dplyr::select(-SMILES)
  
  # Checking that any of the categorical values would not have more than 80% of existing/non-existing values
  descs = descs %>%  
    select(-c(nearZeroVar(descs, 
                          freqCut = nearZeroVar_freqCut)))
  
  # Removing columns that are correlated to each other
  correlationMatrix <- cor(descs, use = "complete.obs")
  highlyCorrelated <- findCorrelation(correlationMatrix, 
                                      cutoff = highlyCorrelated_cutoff)
  
  descs <- descs %>%
    dplyr::select(-highlyCorrelated) %>%
    bind_cols(smiles)
  
  return(descs)
}


training_logIE_pred_model = function(data,
                                     folds = 5,
                                     fitControlMethod = "boot",
                                     method = "xgbTree",
                                     split = NULL,
                                     save_model_name = NULL) {
  #data = data_clean
  #split = 0.8
  set.seed(123)
  if (!is.null(split)) {
    training_set = tibble()
    test_set = tibble()
    for (data_type_this in levels(factor(data$data_type))) { # data type as variable
      data_this = data %>%
        filter(data_type == data_type_this)
      name = data_this %>% select(name) %>% unique()
      split_train_test = sample.split(name$name, SplitRatio = split)
      name = name %>% mutate(split_train_test = split_train_test)
      data_this = data_this %>%
        left_join(name)
      training_set = training_set %>%
        bind_rows(data_this %>% 
                    filter(split_train_test)) %>%
        select(-split_train_test)
      test_set = test_set %>%
        bind_rows(data_this %>% 
                    filter(!split_train_test)) %>%
        select(-split_train_test)
    }
  } else {
    training_set = data
  }
  
  folds = groupKFold(training_set$name, k = folds) 
  fitControl <- trainControl(method = fitControlMethod, index = folds)
  
  model <- train(logIE ~ ., 
                 data = training_set %>% select(-name, -data_type, -SMILES),
                 method = method,
                 trControl = fitControl)
  
  
  
  training_set <- training_set %>%
    mutate(logIE_predicted = predict(model, newdata = training_set))
  
  RMSE_training_set = rmse(training_set$logIE, training_set$logIE_predicted)
  
  ## Determining most influential descriptors ----
  variable_importance <- varImp(model)
  variable_importance <- as.data.frame(variable_importance$importance)
  
  data = list("training_set" = training_set)
  metrics = list("RMSE_training_set" = RMSE_training_set)
  
  if (!is.null(save_model_name)) {
    saveRDS(model,save_model_name)
  }
  
  if (!is.null(split)) {
    test_set <- test_set %>%
      mutate(logIE_predicted = predict(model, newdata = test_set))
    RMSE_test_set = rmse(test_set$logIE, test_set$logIE_predicted)
    
    data = list("training_set" = training_set,
                "test_set" = test_set)
    metrics = list("RMSE_training_set" = RMSE_training_set,
                   "RMSE_test_set" = RMSE_test_set)
  }
  
  model = list("model" = model,
               "data" = data,
               "metrics" = metrics,
               "variable_importance" = variable_importance)
  return(model)
  
}


# *************************************************************
# ---- Homologue series functions  (Melanie, Lara, Thomas) ----
# *************************************************************

is_homologue = function(analyte_SMILES,
                        standard_SMILES,
                        homologue_pattern_SMART) {
  analyte = get.mol2formula(parse.smiles(analyte_SMILES)[[1]])@isotopes
  standard = get.mol2formula(parse.smiles(standard_SMILES)[[1]])@isotopes
  pattern = get.mol2formula(parse.smiles(homologue_pattern_SMART)[[1]])@isotopes
  homologue_of_analyte = as_tibble(analyte) %>%
    left_join(as_tibble(pattern), by = c("isoto" = "isoto")) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate(smaller = as.numeric(number.x) - as.numeric(number.y),
           bigger = as.numeric(number.x) + as.numeric(number.y)) %>%
    select(isoto, smaller, bigger)
  comparison = homologue_of_analyte %>%
    left_join(as_tibble(standard) %>%
                mutate(number = as.numeric(number)),
              by = c("isoto" = "isoto")) %>%
    mutate_all(~replace(., is.na(.), 0))
  if(all(comparison$smaller == comparison$number)) {
    return("smaller")
  } else if  (all(comparison$bigger == comparison$number)) {
    return("bigger")
  } else {
    return(NA)
  }
}

concentration_forAnalytes_model <- function(filename_data, 
                                            filename_smiles, 
                                            filename_eluent, 
                                            pred_model,
                                            compounds_to_be_removed_as_list = c(),
                                            organic_modifier = "MeCN",
                                            pH.aq. = 7.0,
                                            NH4 = 1,
                                            additive = "ammoniumacetate",
                                            additive_concentration_mM = 2,
                                            instrument = "Orbitrap",
                                            source = "ESI") {
  
  analysis_data <- read_excel_allsheets(filename_data)
  SMILES_data <- read_SMILES(filename_smiles, compounds_to_be_removed_as_list)
  
  analysis_data = analysis_data %>%
    mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
    left_join(SMILES_data) %>%
    drop_na(SMILES) %>%
    mutate(RT = as.numeric(RT),
           area_IC = Area*IC,
           Theoretical_conc_uM = Theoretical_amt/Molecular_weight) 
  
  analysis_data_descr <- analysis_data %>%    
    group_by(SMILES, Compound) %>%
    mutate(slope = linear_regression(area_IC, Theoretical_conc_uM)$slope,
           RT = mean(RT)) %>%
    ungroup()
  
  analysis_data_descr = add_mobile_phase_composition(data = analysis_data_descr,
                                                     eluent_file_name = filename_eluent,
                                                     organic_modifier = organic_modifier,
                                                     pH.aq. = pH.aq.,
                                                     NH4 = NH4,
                                                     additive = additive,
                                                     additive_concentration_mM = additive_concentration_mM,
                                                     instrument = instrument,
                                                     source = source)
  
  
  
  analysis_data_descr <- PaDEL_original(analysis_data_descr)
  analysis_data_descr <- analysis_data_descr %>%
    mutate(logIE_predicted = predict(pred_model$model, newdata = analysis_data_descr))
  
  # lm function to find RF of suspect compound
  linMod <- lm(log10(slope) ~ logIE_predicted, data = analysis_data_descr %>%
                 drop_na(slope) %>%             # need some linearity check here?!
                 filter(Compound != "PFHpS-br"))     # Take this out in final!!
  
  
  # Plot measured vs predicted
  plot_predictedIE_slope = ggplot() +
    geom_point(data = analysis_data_descr,
               mapping = aes(logIE_predicted, log10(slope),
                             text = Compound),
               color = "black",
               alpha = 0.5,
               size = 3) +
    geom_abline(slope = summary(linMod)$coefficients[2], intercept = summary(linMod)$coefficients[1]) 
  ggplotly(plot_predictedIE_slope)
  
  
  # Find RF values from predicted IEs to all non-calibrant analytes
  analytes_concentrations <- analysis_data_descr %>%
    # filter(is.na(slope)) %>%
    mutate(slope_pred = 10^(summary(linMod)$coefficients[2]*logIE_predicted + summary(linMod)$coefficients[1])) %>%
    mutate(conc_pred = area_IC/slope_pred) %>%
    select(colnames(analysis_data), slope_pred, conc_pred)
  
  plot_predicted_theoretical_conc <- ggplot() +
    geom_point(data = analytes_concentrations %>%
                 drop_na(Theoretical_conc_uM) 
               #   group_by(Compound, Theoretical_conc_uM) %>%
               #   mutate(conc_pred = mean(conc_pred)) %>%
               #   ungroup()
               ,
               mapping = aes(Theoretical_conc_uM, conc_pred,
                             color = Compound)) +
    geom_abline(slope = 1, intercept = 0) +
    scale_y_log10(limits = c(10^-5, 10^0)) +
    scale_x_log10(limits = c(10^-5, 10^0)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_abline(slope = 1, intercept = 1) +
    geom_abline(slope = 1, intercept = -1) +
    theme(aspect.ratio = 1)
  ggplotly(plot_predicted_theoretical_conc)
  
  
  #Return list with data, plot
  data_predicted = list("plot_predictedIE_slope" = plot_predictedIE_slope,
                        "plot_predicted_theoretical_conc" = plot_predicted_theoretical_conc,
                        "data" = analytes_concentrations)
  
  return(data_predicted)
}


concentration_forAnalytes_homolog <- function(filename_data, 
                                              filename_smiles,
                                              homolog_pattern_SMILES,
                                              findHomolog_onlyForAnalytes = TRUE) {
  
  analysis_data <- read_excel_allsheets(filename_data)
  SMILES_data <- read_SMILES(filename_smiles)
  
  analysis_data = analysis_data %>%
    left_join(SMILES_data) %>%
    mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
    drop_na(SMILES) %>%
    mutate(RT = as.numeric(RT),
           area_IC = Area*IC,
           Theoretical_conc_uM = Theoretical_amt/Molecular_weight) 
  
  analysis_data_descr <- analysis_data %>%    
    group_by(SMILES, Compound) %>%
    mutate(slope = linear_regression(area_IC, Theoretical_conc_uM)$slope,
           RT = mean(RT)) %>%
    ungroup()
  
  if(findHomolog_onlyForAnalytes == TRUE) {
    # Create two datasets of analytes and calibrants - quantification based on homologous calibrants used in this analysis?
    data_training_original <- analysis_data_descr %>%
      filter(!is.na(slope)) %>%
      select(Compound, SMILES, slope) %>%
      unique()
    
    data_analytes <- analysis_data_descr %>%
      filter(is.na(slope)) %>%
      select(Compound, Filename, SMILES, area_IC) %>%
      unique()
    
    data_joined <- data_training_original %>%
      left_join(data_analytes, by = character())
    
  } else if (findHomolog_onlyForAnalytes == FALSE) {
    data_joined <- analysis_data_descr %>%
      select(Compound, SMILES, slope) %>%
      unique() %>%
      left_join(analysis_data_descr %>%
                  select(Compound, Filename, SMILES, area_IC) %>%
                  unique(), 
                by = character())
  }
  
  
  # data_joined <- data_joined %>%
  #   group_by(SMILES.x, SMILES.y) %>%
  #   mutate(pattern = is_homologue(SMILES.x, SMILES.y, homolog_pattern_SMILES)) %>% # Problems with isomeric compounds???  "F[C+2]F"
  #   ungroup() %>%
  #   mutate(conc = case_when((pattern == "smaller" | pattern == "bigger") ~ area_IC/slope))
  
  data_joined_smiles <- data_joined %>%
    select(SMILES.x, SMILES.y) %>%
    unique() %>%
    group_by(SMILES.x, SMILES.y) %>%
    mutate(pattern = is_homologue(SMILES.y, SMILES.x, homolog_pattern_SMILES)) %>%
    ungroup()
  
  data_joined <- data_joined %>%
    left_join(data_joined_smiles) %>%
    mutate(conc = case_when((pattern == "smaller" | pattern == "bigger") ~ (area_IC)/slope)) %>%
    filter(!is.na(pattern))
  
  return(data_joined)
}




# ********************************
# ---- Linearity functions ----
# ********************************

linear_regression <- function(y, x) {
  df <- tibble(x = x, y = y) %>%
    na.omit()
  y = df$y
  x = df$x
  if(length(x) == 0){
    slope = NA
    intercept = NA
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
  } else if (length(y) > 5) {
    for (i in length(y):5){
      y = y[2:length(y)]
      x = x[2:length(x)]
      slope = summary(lm(y ~ x))$coefficients[2]
      intercept = summary(lm(y ~ x))$coefficients[1]
      residuals = (y - (slope*x +intercept))/y*100
      regression_parameters <- list("slope" = slope, "intercept" = intercept)
      if (max(abs(residuals)) < 10) {
        return(regression_parameters)
        break
      }
    }
    return(regression_parameters)
  } else {
    slope = summary(lm(y ~ x))$coefficients[2]
    intercept = summary(lm(y ~ x))$coefficients[1]
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
  }
}




# *******************
# ---- SIRIUS ----
# *******************

# need to rewrite the function with more variables (ala elements as text etc)
ms_to_sirius_output = function(folder_with_ms_files,
                               folder_with_SIRIUS,
                               output_dir_name = character()) {
  setwd(folder_with_ms_files)
  filenames = dir(pattern="*.ms")
  for (file in filenames) {
    setwd(folder_with_ms_files)
    data = readLines(file)
    for (j in 1:length(data)) {
      if (!is.na(data[j]) & data[j] == ">formula ") {
        data = data[-j]
        write_lines(data, file, sep="\n")
        break
      }
    }
    for (i in 5:length(data)) {
      if (data[i] == ">ms2 " & !is.na(data[i+1])) {
        setwd(folder_with_SIRIUS)
        system(paste("cd ", folder_with_SIRIUS))
        command = str_glue(" sirius ",
                           "-i ", str_glue('"', folder_with_ms_files, "/", file, '"', sep = ""),
                           " -o ", str_glue('"', folder_with_ms_files, "/", output_dir_name, '"', sep = ""),
                           " formula -c 10",
                           " -p orbitrap",
                           " -E CHO -e NP[8]B[11]Si[9]S[12]Cl[18]Se[2]Br[10]FI[6]K[1]Na[1]As[2]",
                           # " -d PUBCHEM",
                           " --ppm-max=5",
                           " structure",
                           sep =" ") #makes text for command promp
        javaOutput = system(command, intern = TRUE) #goes into commant prompt
        print(file)
      }
    }
  }
}



# *********************************************
# ---- Pille and Harry - SIRIUS tables and toxicity ----
# *********************************************

#this code is only for users who have comTox data
#the address of the website: https://comptox.epa.gov/dashboard/
GetStandardSMILESandMSreadyForm <- function(compound_file_unstandard, comtox_file){
  unstandard_data <- read_delim(compound_file_unstandard, delim = ";") %>%
    mutate(stand_SMILES = standardize_smiles(SMILES)) %>%
    mutate(stand_inchi = smiles_to_standardinchi(stand_SMILES)) %>%
    mutate(stand_inchikey = smiles_to_standardinchikey(stand_SMILES)) %>%
    group_by(stand_inchikey) %>%
    mutate(stand_inchikey = strsplit(stand_inchikey, split = "=")[[1]][2]) %>%
    ungroup()
  
  MSready_file <- read_delim(comtox_file, delim = ";") %>%
    select(DSSTox_Substance_ID, `DSSTox_Compound_ID (MS-Ready)`, `Formula (MS-Ready)`, `SMILES (MS-Ready)`, `InChIString (MS-Ready)`, `InChIKey (MS-Ready)`)
  colnames(MSready_file)[1] <- "dtxsid"
  
  unstandard_data <- unstandard_data %>%
    left_join(MSready_file) %>%
    mutate(stand_SMILES = as.character(stand_SMILES)) %>%
    mutate(`InChIString (MS-Ready)` = as.character(`InChIString (MS-Ready)`)) %>%
    group_by(`InChIKey (MS-Ready)`) %>%
    mutate(InChIKey_SIRIUS = substring(`InChIKey (MS-Ready)`, first = 0, last = 14)) %>%
    ungroup()
  
  return(unstandard_data)
}

GetSmilesFromSIRIUS <- function(SIRIUS_file) {
  setwd(SIRIUS_file)
  subfolder <- dir(SIRIUS_file, all.files = TRUE, recursive = TRUE, pattern = "*.tsv")
  
  compiled_data <- tibble()
  
  for (filename in subfolder) {
    if (grepl("/fingerid/", filename, fixed=TRUE)){
      id <- str_split(filename, pattern = "_")[[1]][3]
      
      data_smile <- read_delim(filename, delim = "\t") %>%
        select(rank, molecularFormula, score, smiles, xlogp, tanimotoSimilarity, inchi, inchikey2D, dbflags) %>% 
        mutate(filename = filename) %>%
        mutate(id = id) %>%
        select(filename, id, everything())
      
      names(data_smile)[names(data_smile) == "smiles"] <- "SMILES"
      names(data_smile)[names(data_smile) == "inchikey2D"] <- "InChIKey"
      names(data_smile)[names(data_smile) == "score"] <- "FingerID_Score"
      names(data_smile)[names(data_smile) == "rank"] <- "Structure_rank"
      
      compiled_data <- compiled_data %>%
        bind_rows(data_smile) %>%
        unique()
    }
  }
  return(compiled_data)
}

SiriusScores <- function(SIRIUS_file, Compounds_info_file) {
  setwd(SIRIUS_file)
  subfolder <- dir(SIRIUS_file, all.files = TRUE, recursive = TRUE, pattern = "*.info")
  
  compiled_data <- tibble()
  compiled_data <- compiled_data %>%
    mutate(isotscore = "isotscore") %>%
    mutate(treescore = "treescore") %>%
    mutate(siriusscore = "siriusscore") %>%
    mutate(topCSIscore = "topCSIscore") %>%
    mutate(confidenceScore = "confidenceScore")
  
  for (filename in subfolder) {
    if (grepl("/scores/", filename, fixed=TRUE)){
      
      id <- str_split(filename, pattern = "_")[[1]][3]
      file_name <- str_split(filename, "/")
      molecularFormula <- str_split(file_name[[1]][lengths(file_name)], "_")[[1]][1]
      
      fileConnection <- file(filename)
      record <- readLines(fileConnection)
      close(fileConnection)
      
      IsotopeScore <- substring(grep('sirius.scores.IsotopeScore', record, value = TRUE, fixed = TRUE),28)
      TreeScore <- substring(grep('sirius.scores.TreeScore', record, value = TRUE, fixed = TRUE),25)
      SiriusScore <- substring(grep('sirius.scores.SiriusScore', record, value = TRUE, fixed = TRUE),27)
      ConfidenceScore <- substring(grep('fingerid.ConfidenceScore', record, value = TRUE, fixed = TRUE),27)
      TopCSIscore <- substring(grep('fingerid.blast.TopCSIScore', record, value = TRUE, fixed = TRUE),28)
      
      if (length(TopCSIscore) != 0 | length(ConfidenceScore) != 0) {
        if (length(TopCSIscore) == 0) {
          TopCSIscore <- "-99"
        }
        if (length(ConfidenceScore) == 0) {
          ConfidenceScore <- "-99"
        }
        
        filedata <- data.frame(id)
        
        filedata <- filedata %>%
          mutate(molecularFormula = molecularFormula) %>%
          mutate(isotscore = IsotopeScore) %>%
          mutate(treescore = TreeScore) %>%
          mutate(siriusscore = SiriusScore) %>%
          mutate(confidenceScore = ConfidenceScore) %>%
          mutate(topCSIscore = TopCSIscore)
        
        compiled_data <- compiled_data %>%
          bind_rows(filedata)
      }
    }
  }
  
  data_scores <- compiled_data %>%
    mutate(isotscore = as.numeric(isotscore)) %>%
    mutate(treescore = as.numeric(treescore)) %>%
    mutate(siriusscore = as.numeric(siriusscore)) %>%
    mutate(confidenceScore = as.numeric(confidenceScore)) %>%
    mutate(topCSIscore = as.numeric(topCSIscore)) %>%
    group_by(id) %>%
    arrange(desc(siriusscore)) %>%
    mutate(Formula_rank = row_number()) %>%
    ungroup()
  
  data_comp_info <- read_delim(Compounds_info_file, delim = ";") %>%
    select(id, dtxsid, parentMS, RT, ionmode)
  
  
  data_scores <- left_join(data_scores,data_comp_info) %>%
    select(id, molecularFormula, Formula_rank, dtxsid, parentMS, RT, ionmode, everything())
  return(data_scores)
}

# *********************************
# ---- Pille - MassBank codes ----
# *********************************


# MassBank table
data_path_txt = "C:/Users/pipe0902/Kruvelab/Anneli Kruve - R codes for variouse things/MassBank-data-2020.03"
setwd(data_path_txt)
subfolder <- dir(data_path_txt, all.files = TRUE, recursive = TRUE, pattern = "*.txt")

all_MassBank_parameter <- as_tibble()
for (filename in subfolder) {
  fileConnection <- file(filename)
  record <- readLines(fileConnection)
  close(fileConnection)
  all_MassBank_parameter <- all_MassBank_parameter %>%
    bind_rows(readMassBankfilesParam(record))
}


# Function to get parameters row/tabel for one txt data
readMassBankfilesParam <- function(record){ 
  
  # get basic information
  id <- substring(grep("ACCESSION:", record, value = TRUE, fixed = TRUE), 12)
  
  # get information of chemical compound
  chnames <- list()
  chnames <- as.list(substring(grep('CH$NAME:',record, value = TRUE, fixed = TRUE),10))
  chnames <- substring(grep("CH$NAME:", record, value = TRUE, fixed = TRUE), 10)
  formula <- substring(grep("CH$FORMULA:", record, value = TRUE, fixed = TRUE), 13)
  exactMass <- as.numeric(substring(grep("CH$EXACT_MASS:", record, value = TRUE, fixed = TRUE), 16))
  cas <- substring(grep('CH$LINK: CAS', record, value = TRUE, fixed = TRUE), 14)   #added 
  dtxsid <- substring(grep('CH$LINK: COMPTOX', record, value = TRUE, fixed = TRUE), 18)
  
  # writing a table with information
  
  paramlist = list(id,dtxsid,chnames,formula,exactMass,cas)
  
  summary_table <- as.tibble(t(paramlist))
  
  
  return(summary_table)
}  



