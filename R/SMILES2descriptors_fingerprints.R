#' Calculate the true (binary) SIRIUS fingerprint for given SMILES notations
#' 
#' This function reads structural patterns from files downloaded from https://github.com/boecker-lab/sirius, 
#' which correspond to SIRIUS:CSI:FingerID outputted fingerprint features. 
#' After that it generates binary SIRIUS fingerprints for all input SMILES notations.
#'
#' @author Ida Rahu
#' 
#' @import rcdk
#' @import rJava
#' @import dplyr
#' @import fingerprint
#' @import data.table
#' 
#' @param SMILES_list data.frame: Data frame containing SMILES notations
#' @param write_raw_fp logical: Choose whether to write out the file (.tsv) that contains raw 
#'                              (all possible features, 8924) SIRIUS fingerprints (default = TRUE)
#' @param mode character: Specifies which SIRIUS fingerprint features will be returned:
#'                       - 'pos': Positive ionization mode features (3878)
#'                       - 'neg': Negative ionization mode features (4072)
#'                       - 'overlapping': Features common to both ionization modes (3494)
#'                       - 'raw': Raw SIRIUS fingerprint (all 8925 features)
#'                       
#' @return data.frame: Contains SMILES and calculated fingerprint features based on the selected mode
#' @export
#' 
SMILES2SIRIUS_fp <- function(SMILES_list, write_raw_fp=T, mode='raw') {
  options(java.parameters = "- Xmx1024m")
  # Function for reading in the structural patterns (from relevant files) that
  # correspond to the SIRIUS+CSI:FingerID fingerprint features 
  pattern_file_reader <- function(file_name, split_pattern) {
    con <- file(file_name)
    lines <- readLines(con)
    close(con)
    slines <- strsplit(lines, split_pattern)
    colCount <- max(unlist(lapply(slines, length)))
    
    patterns <- data.frame(matrix(nrow=0, ncol=colCount))
    
    for (i in 1:length(slines)) {
      line <- slines[[i]]
      for (j in 1:length(line)) {
        patterns[i, j] <- line[j]
      }
    }
    return(patterns)
  }
  
  # Fingerprints that cover all SIRIUS fingerprints 
  # (files downloaded from: https://github.com/boecker-lab/sirius)
  OpenBabelFP3_names <- paste0('AbsIdx_', c(0:54))
  suppressWarnings(OpenBabelFP3_SMARTS <- unlist(pattern_file_reader(system.file('fingerprint_files', 'OpenBabel_FP3_patterns.txt', package = 'kruvelabFns'), '\t')[1], use.names=F))

  CDKsubstructure_names <- paste0('AbsIdx_', c(55:361))
  MACCS_names <- paste0('AbsIdx_', c(362:527))
  PubChem_names <- paste0('AbsIdx_', c(528:1408))
  KlekotaRoth_names <- paste0('AbsIdx_', c(1409:6268))
  
  ECFP6_names <- paste0('AbsIdx_', c(6269:8178))
  ECFP6_hashes <- unlist(read.table(system.file('fingerprint_files', 'ecfp_fp_hashes.txt', package = 'kruvelabFns'), header=F), use.names=F)
  
  custommadeSMARTS_names <- paste0('AbsIdx_', c(8179:8461))
  suppressWarnings(custommade_SMARTS <- unlist(pattern_file_reader(system.file('fingerprint_files', 'biosmarts_aka_custom_made_fps.txt', package = 'kruvelabFns'), '\n'), use.names=F))
  
  ringsystems_names <- paste0('AbsIdx_', c(8462:8924))
  ringsystems_SMARTS <- unlist(pattern_file_reader(system.file('fingerprint_files', 'ringsystem_fps.txt', package = 'kruvelabFns'), '\n'), use.names=F)
  
  
  no_columns = 1 + length(OpenBabelFP3_names) + length(CDKsubstructure_names) +
    length(MACCS_names) + length(PubChem_names) + 
    length(KlekotaRoth_names) + length(ECFP6_names) + 
    length(custommadeSMARTS_names) + length(ringsystems_names)
  
  
  final_fp_data <- data.frame(matrix(nrow=nrow(SMILES_list), ncol=no_columns))
  
  colnames(final_fp_data) <- c('SMILES', OpenBabelFP3_names, 
                               CDKsubstructure_names, MACCS_names, 
                               PubChem_names, KlekotaRoth_names, ECFP6_names, 
                               custommadeSMARTS_names, ringsystems_names) 
  final_fp_data[, 1] <- SMILES_list
  
  fp_index <- data.frame(matrix(nrow=8, ncol=3))
  colnames(fp_index) <- c('fingerprint', 'start_index', 'end_index')
  fp_index$fingerprint <- c('FP3', 'substructure', 'maccs', 'pubchem', 'kr', 'ecfp6', 'custom', 'ring')
  fp_index$start_index <- c(2, length(OpenBabelFP3_names)+2, 
                            length(OpenBabelFP3_names)+length(CDKsubstructure_names)+2,
                            length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+2,
                            length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+2,
                            length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+2,
                            length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names)+2,
                            length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names)+length(custommadeSMARTS_names)+2)
  
  fp_index$end_index <- c(1+length(OpenBabelFP3_names), 
                          1+length(OpenBabelFP3_names)+length(CDKsubstructure_names),
                          1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names),
                          1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names),
                          1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names),
                          1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names),
                          1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names)+length(custommadeSMARTS_names),
                          1+length(OpenBabelFP3_names)+length(CDKsubstructure_names)+length(MACCS_names)+length(PubChem_names)+length(KlekotaRoth_names)+length(ECFP6_names)+length(custommadeSMARTS_names)+length(ringsystems_names))
  
  i = 1
  for (SMILES in SMILES_list[c(1:nrow(SMILES_list)),1]) {
    mol <- parse.smiles(SMILES)[[1]]
    if(!is.null(mol)){
      openbabel_fingerprints <- get.fingerprint(mol, type='substructure',
                                                substructure.pattern=OpenBabelFP3_SMARTS)
      final_fp_data[i, c(fp_index$start_index[1]:fp_index$end_index[1])] <- strsplit(as.character(openbabel_fingerprints), "")[[1]]
      
      substr_fingerprints <- get.fingerprint(mol, type='substructure')
      final_fp_data[i, c(fp_index$start_index[2]:fp_index$end_index[2])] <- strsplit(as.character(substr_fingerprints), "")[[1]]
      
      maccs_fingerprints <- get.fingerprint(mol, type='maccs')
      final_fp_data[i, c(fp_index$start_index[3]:fp_index$end_index[3])] <- strsplit(as.character(maccs_fingerprints), "")[[1]]
      
      pubchem_fingerprints <- get.fingerprint(mol, type='pubchem')
      final_fp_data[i, c(fp_index$start_index[4]:fp_index$end_index[4])] <- strsplit(as.character(pubchem_fingerprints), "")[[1]]
      
      kr_fingerprints <- get.fingerprint(mol, type='kr')
      final_fp_data[i, c(fp_index$start_index[5]:fp_index$end_index[5])] <- strsplit(as.character(kr_fingerprints), "")[[1]]
      
      ecfp_fingerprints <- get.fingerprint(mol, type='circular', circular.type='ECFP6', fp.mode='count')
      for (idx in 1:length(ecfp_fingerprints@features)) {
        hash <- strsplit(as.character(ecfp_fingerprints@features[[idx]]), ':')[[1]][1]
        if (hash %in% ECFP6_hashes) {
          right_hash <- which(hash == ECFP6_hashes)
          column_name <- ECFP6_names[right_hash]
          final_fp_data[i, column_name] <- 1
        }
      }
      
      custommade_fingerprints <- get.fingerprint(mol, type="substructure",
                                                 substructure.pattern=custommade_SMARTS)
      final_fp_data[i, c(fp_index$start_index[7]:fp_index$end_index[7])] <- strsplit(as.character(custommade_fingerprints), "")[[1]]
      
      ring_fingerprints <- get.fingerprint(mol, type="substructure",
                                           substructure.pattern=ringsystems_SMARTS)
      final_fp_data[i, c(fp_index$start_index[8]:fp_index$end_index[8])] <- strsplit(as.character(ring_fingerprints), "")[[1]]
      
    }
    else {
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>% rbind(wrong_smile)
      print(SMILES)
    }
    gc()
    i = i + 1
  }
  
  final_fp_data <- final_fp_data %>% mutate_at(ECFP6_names, ~replace_na(as.numeric(.), 0))
  final_fp_data <- final_fp_data %>% mutate_at(c(3:ncol(final_fp_data)), as.numeric)
  
  if (write_raw_fp) {
    write.table(final_fp_data, 'raw_SIRIUS_fp.tsv', row.names=F, col.names=T, sep='\t', quote=F)
  }
  
  # Generating the names for features so calculated fingerprint features would 
  # match with SIRIUS+CSI:FingerID absolute index.
  positive_idxs <- pattern_file_reader(system.file('fingerprint_files', 'csi_fingerid.tsv', package = 'kruvelabFns'), '\t')
  colnames(positive_idxs) <- positive_idxs[1,]
  positive_idxs <- positive_idxs[-1, ] 
  
  negative_idxs <- pattern_file_reader(system.file('fingerprint_files', 'csi_fingerid_neg.tsv', package = 'kruvelabFns'), '\t')
  colnames(negative_idxs) <- negative_idxs[1,]
  negative_idxs <- negative_idxs[-1, ] 
  
  together_idx <- merge(positive_idxs, negative_idxs, by='absoluteIndex')
  together_idx <- together_idx[order(as.numeric(as.character(together_idx$absoluteIndex))), ]
  positive_idxs$absoluteIndex <- sub('^', 'AbsIdx_', positive_idxs$absoluteIndex)
  negative_idxs$absoluteIndex <- sub('^', 'AbsIdx_', negative_idxs$absoluteIndex)
  together_idx$absoluteIndex <- sub('^', 'AbsIdx_', together_idx$absoluteIndex)
  
  # print output reference
  cat(crayon::green("Fingerprint/descriptor files downloaded from: https://github.com/boecker-lab/sirius\n"))
  
  # Generating data frames containing fingerprint features overlapping with those 
  # outputted by SIRIUS+CSI:FingerID, based on the selected ionization mode ('pos', 'neg', 'overlapping' or 'raw').
  if (mode == 'pos') {
    positive_mode_data <- final_fp_data %>% dplyr::select(c('SMILES', positive_idxs$absoluteIndex))
    return(positive_mode_data)
  } else if (mode == 'neg') {
    negative_mode_data <- final_fp_data %>% dplyr::select(c('SMILES', negative_idxs$absoluteIndex))
    return(negative_mode_data)
  } else if (mode == 'overlapping'){
    general_mode_data <- final_fp_data %>% dplyr::select(c('SMILES', together_idx$absoluteIndex))
    return(general_mode_data)
  } else {
    return(final_fp_data)
  }
  
}
