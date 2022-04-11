# Script name 01_Functions.R as part of
# 00_SetGlobalVars_XXX.R
# 01_Functions.R
# 02_MakeSingleCellDataList.R
# 03_DensityPlots.R
# 04_MakeSummaryData.R
# 05_DynamicsPlots.R
# 06_MakeOtherPlots.R
# 07_CombineReplicates.R
# 08_FinalMeanFigures.R
# 
# Short description:
# Functions used to 
# 1) read in the single cell data,
# 2) to make a summary data frame, and
# 3) to make different kinds of plots.

# Starting data: 2018-07-19
# Last modified: 2019-09-30
# Written by Muriel Heldring

# R version 3.5.2 (2018-12-20)

# -----------------------------------------------------------------------------
# ----------------------------DEFINE FUNCTIONS---------------------------------
# -----------------------------------------------------------------------------

# 1. Function to make a condition vector based on the compounds of interest
makeCondition <- function(compounds = NULL, others = NULL){
  # If there are no compounds or specific other states that need to be loaded only,
  # an empty character is returned
  if (is.null(compounds) & is.null(others)){
    return(character())} 
  else {
    # There is some state that needs to be loaded, so the 'where' statement in taken as start
    begin <- " where "
    # If there are specific compounds, that correct string is made
    if (!is.null(compounds)){
      # TO DO: GENERALISE SUCH THAT WHITE SPACES AROUND THE COMPOUND NAME ARE IGNORED
      # cmp <- paste0("treatment in ","('\"",paste(compounds,collapse = " \"', '\"")," \"')")
      cmp <- paste0("treatment in ","('\"",paste(COMPS_OF_INTEREST,collapse = " \"', '\"")," \"') or treatment in ",
                    "('\"", paste(COMPS_OF_INTEREST,collapse = "\"', '\""),"\"')")}
    
    # If there are other specific interests, that correct string is made
    if (!is.null(others)){
      other <- others}
    
    # Make the complete character
    if (exists("cmp") & exists("other")){
      return(paste0(begin,cmp," and ",other))}
    else if (!exists("cmp")) {
      return(paste0(begin,other))}
    else {
      return(paste0(begin,cmp))}}
}

# 2. Function to load in the data
sql_load <- function(x, selection = "*", condition = NULL, file.format = list(header = TRUE, sep = "\t")){
  f <- file(x)
  df <- sqldf(paste0("select ", selection, " from f", condition), 
              file.format = file.format)
  print(paste0("Loaded file ",f))
  close(f)
  return(df)
}

# 3a. Function to read in the data and store it in a list
# Read in the data contained in the FILES. The data will be stored in a list
# with named entries in the following order: protein, experiment replicate, channel
makeDataFrameList <- function(path = character(), files = character(), proteins = character(), 
                              replicates = character(), channels = character(), condition = character()){
  # Make a list, called dfListSCD, 
  # of which each element has the name of a protein and is a list ...
  # of which each element has the name of the replicate and is a list ...
  # of which each element has the name of the channel and is a data frame.
  dfListSCD <- sapply(proteins,function(protName){
    list(sapply(replicates,function(repl){
      list(sapply(channels,function(channel){
        if(file.exists(paste0(path,list(files[grepl(protName,files) & grepl(repl,files) & grepl(channel,files)])))){
          sql_load(x = paste0(path,list(files[grepl(protName,files) & grepl(repl,files) & grepl(channel,files)])), 
                   condition = condition)}
      }))
    }))
  })
  
  # Remove entries of dfListSCD[[proteinName]][[replicate]][[channel]] for which there is no data frame
  dfListSCD <- sapply(names(dfListSCD),function(protName){
    sapply(names(dfListSCD[[protName]]),function(replName){
      dfListSCD[[protName]][[replName]][unlist(lapply(dfListSCD[[protName]][[replName]], length) != 0)]
    }, simplify = F, USE.NAMES = T)
  }, simplify = F, USE.NAMES = T)
  
  # Remove entries of dfListSCD[[proteinName]][[replicate]] for which there is no entry
  dfListSCD <- sapply(names(dfListSCD),function(protName){
    dfListSCD[[protName]][unlist(lapply(dfListSCD[[protName]],length) != 0)]
  }, simplify = F, USE.NAMES = T)
  
  # Remove entries of dfListSCD[[proteinName]] for which there is no entry
  dfListSCD <- dfListSCD[unlist(lapply(dfListSCD,length) != 0)]
  
  # Return the list of data frames
  return(dfListSCD)
}

# 3b. Function to read in the data and store it in a list for all compounds
# Read in the data contained in the FILES. The data will be stored in a list
# with named entries in the following order: protein, experiment replicate, channel
makeDataFrameListAllCompounds <- function(path = character(), files = character(), proteins = character(), 
                                          replicates = character(), channels = character()){
  dfListSCD <- list()
  # Make a list, called dfListSCD, 
  # of which each element has the name of a protein and is a list ...
  # of which each element has the name of the replicate and is a list ...
  # of which each element has the name of the channel and is a data frame.
  for (protName in proteins){
    for (repl in replicates) {
      for (channel in channels) {
        if( file.exists(paste0(path,list(files[grepl(protName,files) & grepl(repl,files) & grepl(channel,files)]))) ){
          dfListSCD[[protName]][[repl]][[channel]] <- read_tsv(file = paste0(path,list(files[grepl(protName,files) & grepl(repl,files) & grepl(channel,files)])), 
                                                            na = "NA")}
      }
    }
  }
  
  # dfListSCD <- sapply(proteins,function(protName){
  #   list(sapply(replicates,function(repl){
  #     list(sapply(channels,function(channel){
  #       if(file.exists(paste0(path,list(files[grepl(protName,files) & grepl(repl,files) & grepl(channel,files)])))){
  #         fread(paste0(path,list(files[grepl(protName,files) & grepl(repl,files) & grepl(channel,files)])))}
  #     }))
  #   }))
  # })
  
  # Remove entries of dfListSCD[[proteinName]][[replicate]][[channel]] for which there is no data frame
  dfListSCD <- sapply(names(dfListSCD),function(protName){
    sapply(names(dfListSCD[[protName]]),function(replName){
      dfListSCD[[protName]][[replName]][unlist(lapply(dfListSCD[[protName]][[replName]], length) != 0)]
    }, simplify = F, USE.NAMES = T)
  }, simplify = F, USE.NAMES = T)
  
  # Remove entries of dfListSCD[[proteinName]][[replicate]] for which there is no entry
  dfListSCD <- sapply(names(dfListSCD),function(protName){
    dfListSCD[[protName]][unlist(lapply(dfListSCD[[protName]],length) != 0)]
  }, simplify = F, USE.NAMES = T)
  
  # Remove entries of dfListSCD[[proteinName]] for which there is no entry
  dfListSCD <- dfListSCD[unlist(lapply(dfListSCD,length) != 0)]
  
  # Return the list of data frames
  return(dfListSCD)
}

# 3c. Function to read in the data and store it in a list for all compounds
# Read in the data contained in the FILES. The data will be stored in a list
# with named entries in the following order: protein, experiment replicate, channel
makeDataFrameListRelevantCompounds <- function(path = character(), files = character(), proteins = character(), 
                                          replicates = character(), channels = character(),
                                          compsOfInt = character()){
  dfListSCD <- list()
  # Make a list, called dfListSCD, 
  # of which each element has the name of a protein and is a list ...
  # of which each element has the name of the replicate and is a list ...
  # of which each element has the name of the channel and is a data frame.
  for (protName in proteins){
    for (repl in replicates) {
      for (channel in channels) {
        if( file.exists(paste0(path,list(files[grepl(protName,files) & grepl(repl,files) & grepl(channel,files)]))) ){
          dftmp <- read_tsv(file = paste0(path,list(files[grepl(protName,files) & grepl(repl,files) & grepl(channel,files)])), 
                                                               na = "NA")
          dftmp <- dftmp %>% filter(treatment %in% compsOfInt)
          dfListSCD[[protName]][[repl]][[channel]] <- dftmp
          gc()}
      }
    }
  }
  
  # Remove entries of dfListSCD[[proteinName]][[replicate]][[channel]] for which there is no data frame
  dfListSCD <- sapply(names(dfListSCD),function(protName){
    sapply(names(dfListSCD[[protName]]),function(replName){
      dfListSCD[[protName]][[replName]][unlist(lapply(dfListSCD[[protName]][[replName]], length) != 0)]
    }, simplify = F, USE.NAMES = T)
  }, simplify = F, USE.NAMES = T)
  
  # Remove entries of dfListSCD[[proteinName]][[replicate]] for which there is no entry
  dfListSCD <- sapply(names(dfListSCD),function(protName){
    dfListSCD[[protName]][unlist(lapply(dfListSCD[[protName]],length) != 0)]
  }, simplify = F, USE.NAMES = T)
  
  # Remove entries of dfListSCD[[proteinName]] for which there is no entry
  dfListSCD <- dfListSCD[unlist(lapply(dfListSCD,length) != 0)]
  
  # Return the list of data frames
  return(dfListSCD)
}

# 4. Function to remove quotes from every quoted value
removeQuotes <- function(dfListSCD){
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      lapply(reps, function(df){
        df <- data.frame(apply(df,2,function(x){
          if (grepl("\"", x[1])){
            x <- gsub("\"", "",x)}
          else {x <- x}
        }))
      })
    })
  })
  return(dfListSCD)
}

# 5. Function to remove whitespace from a column
removeWS <- function(dfListSCD){
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      lapply(reps, function(df){
        df <- data.frame(apply(df,2,function(x){
          if (grepl(" ", x[1])){
            x <- trimws(x, which = c("both"))} 
          else {x <- x}
        }))
      })
    })
  })
  return(dfListSCD)
}

# 6. Function to remove the irrelevant compounds
removeCompounds <- function(dfListSCD, proteins = character(), 
                            replicates = character(), channels = character(), compounds = character()){
  for (protName in proteins){
    for (repl in replicates) {
      for (channel in channels) {
        dfListSCD[[protName]][[repl]][[channel]] <- dfListSCD[[protName]][[repl]][[channel]][dfListSCD[[protName]][[repl]][[channel]]$treatment %in% compounds, ]
      }
    }
  }
  return(dfListSCD)
}

# 6. Function to convert columns to the desired format
convertColumns <- function(dfListSCD, channel = CHANNELS, columnNames, to = "factor"){
  channel <- match.arg(channel)
  lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      sapply(names(reps), function(chanName){
        if (chanName == channel){
          if (to == "numeric"){
            reps[[chanName]][,columnNames] <- sapply(columnNames, function(x) as.numeric(as.character(reps[[chanName]][[x]])))
            reps[[chanName]]}
          else if (to == "character"){
            reps[[chanName]][,columnNames] <- sapply(columnNames, function(x) as.character(reps[[chanName]][[x]]))
            reps[[chanName]]}
          else if (to == "factor"){
            if (length(columnNames) > 1){
              reps[[chanName]][,columnNames] <- sapply(columnNames, function(x) as.factor(as.character(reps[[chanName]][[x]])))
              reps[[chanName]]}
            else {
              reps[[chanName]][,columnNames] <- as.factor(as.character(reps[[chanName]][[columnNames]]))
              reps[[chanName]]}}
          else if (to == "numericFactor"){
            if (length(columnNames) > 1){
              reps[[chanName]][,columnNames] <- sapply(columnNames, function(x) as.factor(as.numeric(as.character(reps[[chanName]][[x]]))))
              reps[[chanName]]}
            else {
              reps[[chanName]][,columnNames] <- as.factor(as.numeric(as.character(reps[[chanName]][[columnNames]])))
              reps[[chanName]]}}} 
        else {reps[[chanName]]}
      },simplify = F, USE.NAMES = T)
    })
  })
}


# 7. Order the data, according to a columnName
orderDF <- function(dfListSCD, by){
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      lapply(reps, function(df){
        df <- df[do.call(order,as.list(df[,by])),]
      })
    })
  })
}

# 8. Add a well column
extractWell <- function(dfListSCD){
  lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      lapply(reps, function(df){
        df$well <- sapply(strsplit(as.character(df$locationID),"_"),function(x){x[1]})
        df
      })
    })
  })
}

# Filter out the right compounds
filterForCompounds <- function(dfListSCD, compounds){
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      lapply(reps, function(df){
        # Remove the rows with other compounds
        df <- df[df$treatment %in% compounds,]
      })
    })
  })
}

# Rename factor level
renameFactorLevel <- function(dfListSCD, colName, oldName, newName){
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      lapply(reps, function(df){
        if (oldName %in% df[,colName]){
          #df[,colName] <- recode(df[,colName], oldName = newName)
          expression <- paste0("df[,colName] <- recode(df[,colName], \"",oldName,"\" = \"",newName,"\")")
          eval_bare(parse_expr(expression))
          df}
        else {
          df
        }
      })
    })
  })
  return(dfListSCD)
}


# 9. Convert too low levels of GFP intensity < 0.00025 to 0.00001 instead of 0, such that with 
#    a log-transform, not all the 0 values are considered infinite
addValueToZero <- function(dfListSCD, channel = CHANNELS, columnNames){  
  # Check whether input channel is correct
  channel <- match.arg(channel)
  
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      sapply(names(reps), function(chanName){
        if (chanName == channel){
          reps[[chanName]][,columnNames] <- sapply(columnNames,function(x){
            reps[[chanName]][[x]][reps[[chanName]][[x]]<0.0001] <- 0.00001
            reps[[chanName]][[x]]})
          reps[[chanName]]}
        else {
          reps[[chanName]]}
      },simplify = F, USE.NAMES = T)
    })
  })
}

# 10. Function to make new columns with logTransformed data
logTransform <- function(dfListSCD, channel = CHANNELS, columnNames){
  # Check whether input channel is correct
  channel <- match.arg(channel)
  
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      sapply(names(reps), function(chanName){
        if (chanName == channel){
          reps[[chanName]][,paste0("log",columnNames)] <- sapply(columnNames,function(x){log(reps[[chanName]][[x]])})
          reps[[chanName]]}
        else {
          reps[[chanName]]}
      },simplify = F, USE.NAMES = T)
    })
  })
}

# 10. Function to make new columns with logTransformed data
log10Transform <- function(dfListSCD, channel = CHANNELS, columnNames){
  # Check whether input channel is correct
  channel <- match.arg(channel)
  
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      sapply(names(reps), function(chanName){
        if (chanName == channel){
          reps[[chanName]][,paste0("log10",columnNames)] <- sapply(columnNames,function(x){log10(reps[[chanName]][[x]])})
          reps[[chanName]]}
        else {
          reps[[chanName]]}
      },simplify = F, USE.NAMES = T)
    })
  })
}

# 11. Calculate geometric mean
geomMean <- function(x){exp(mean(log(x), na.rm = T))}

# 12. Calculate geometric standarddeviation
geomSD <- function(x){exp(sd(log(x), na.rm = T))}

# 11. Calculate geometric mean
geomMean10 <- function(x){exp(mean(log10(x), na.rm = T))}

# 12. Calculate geometric standarddeviation
geomSD10 <- function(x){exp(sd(log10(x), na.rm = T))}

# columnsForPropCalc = c("Cytoplasm_Intensity_MeanIntensity_image_GFP",
#                        "Cytoplasm_Intensity_IntegratedIntensity_image_GFP",
#                        "Nuclei_Intensity_MeanIntensity_image_GFP",
#                        "Nuclei_Intensity_IntegratedIntensity_image_GFP")
# refTreatment = "DMEM"
# refConcentration = numeric()
# forTreatment = "CDDP"
# sharedGroupVar = "timeID"
# 
# df <- dfListSCD[[1]][[1]][["GFP"]]
# 
# # Make subset data sets of the specified treatments
# dfRef <- df[df$treatment == refTreatment,]
# dfFor <- df[df$treatment %in% c(refTreatment,forTreatment),]
# dfRest <- df[!df$treatment %in% c(refTreatment,forTreatment),]
# 
# # And of the specified concentration
# if(!length(refConcentration) == 0){
#   dfRefTmp <- dfRef[dfRef$dose_uM == refConcentration,]} else {dfRefTmp <- dfRef}
# 
# 
# expression <- paste0("dfRefTmp %>% group_by(", paste0(sharedGroupVar,collapse = ","), ") %>% ",
#                      "mutate(", paste0("mean",refTreatment,columnsForPropCalc,
#                                        " = mean(", 
#                                        columnsForPropCalc,
#                                        ")",
#                                        collapse = ", "),", ",
#                      paste0("sd",refTreatment,columnsForPropCalc,
#                             " = sd(", 
#                             columnsForPropCalc,
#                             ")",
#                             collapse = ", "),")")
# dfRefTmp <- eval_bare(parse_expr(expression))
# names(dfRefTmp)
# # Merge (cbind) the dfFor and dfRefTmp data frames...
# dfFor <- merge(dfFor,dfRefTmp, by = sharedGroupVar, all.x = T)
# dim(dfFor)
# dim(dfRefTmp)
# 
#   parse(
#     text=paste0("ddply(dfRefTmp,.(",
#                 paste0(sharedGroupVar,collapse = ","),
#                 "), summarize,",
#                 paste0("mean",refTreatment,columnsForPropCalc,
#                        " = mean(",
#                        columnsForPropCalc,
#                        ")",
#                        collapse = ", "),", ",
#                 paste0("sd",refTreatment,columnsForPropCalc,
#                        " = sd(",
#                        columnsForPropCalc,
#                        ")",
#                        collapse = ", "),")"))


# 13. Add a boolean column to each data frame where the entry is TRUE if a
# value is larger than the response threshold and FALSE when it is smaller
addBoolResponseColumn <- function(dfListSCD, columnsForPropCalc, 
                                  refTreatment = COMPS_OF_INTEREST, 
                                  refConcentration = numeric(),
                                  forTreatment = COMPS_OF_INTEREST, 
                                  sharedGroupVar = character(),
                                  divGroupVar = character(),
                                  channel = "GFP"){
  # Are the arguments valid?
  refTreatment <- match.arg(refTreatment)
  ifelse(all(forTreatment %in% COMPS_OF_INTEREST),T,stop("Compound does not exist."))
  
  # Add the meanDMEM GFP intensities per timepoint to the data frame
  sapply(dfListSCD,function(prots){
    sapply(prots,function(reps){
      sapply(names(reps),function(chanName){
        # Make the column only for the specified channel
        if (chanName == channel){
          df <- reps[[chanName]]
          
          # Make subset data sets of the specified treatments
          dfRef <- df[df$treatment == refTreatment,]
          dfFor <- df[df$treatment %in% c(refTreatment,forTreatment),]
          dfRest <- df[!df$treatment %in% c(refTreatment,forTreatment),]
          
          # And of the specified concentration
          if(!length(refConcentration) == 0){
            dfRefTmp <- dfRef[dfRef$dose_uM == refConcentration,]}
          else {dfRefTmp <- dfRef}
          
          expression <- paste0("dfRefTmp %>% 
                               group_by(", paste0(sharedGroupVar,collapse = ","), ") %>% ",
                               "summarise(", paste0("mean",refTreatment,columnsForPropCalc,
                                                 " = mean(", 
                                                 columnsForPropCalc,
                                                 ")",
                                                 collapse = ", "),", ",
                               paste0("sd",refTreatment,columnsForPropCalc,
                                      " = sd(", 
                                      columnsForPropCalc,
                                      ")",
                                      collapse = ", "),")")
          dfRefTmp <- eval_bare(parse_expr(expression))
          
          # Merge (cbind) the dfFor and dfRefTmp data frames...
          dfFor <- merge(dfFor,dfRefTmp, by = sharedGroupVar, all.x = T)
          
          # Make new column names in the form of response_Nuclei_Intensity
          newColNames <- paste0("response_",columnsForPropCalc)
          
          # ... and add a boolean column
          dfFor[,newColNames] <- dfFor[,columnsForPropCalc] > (
            dfFor[,paste0("mean",refTreatment,columnsForPropCalc)] + (RESPONSE_THRESHOLD * 
                                                                        dfFor[,paste0("sd",refTreatment,columnsForPropCalc)]))
          
          dfTotal <- as.data.frame(rbindlist(list(dfFor,dfRef,dfRest),fill = T))
          dfTotal} 
        else {
          reps[[chanName]]}
      },USE.NAMES = T,simplify = F)
    },USE.NAMES = T,simplify = F)
  },USE.NAMES = T,simplify = F)
}

# 14. Function to check the input data of plotDensOverTime
checkPlotDensInput <- function(dfListSCD,protein,replicate,channel,treatmentName,timePoints,doses, facet = F){
  # Function that does some error handling. It checks the input for the plots and 
  # provides suitable errors.
  
  # If one of the variables is not provided, a suitable error is given
  if (is.null(dfListSCD) | is.null(protein) | is.null(replicate) | is.null(channel) | 
      is.null(treatmentName) | is.null(timePoints)){
    variables <- c("dfListSCD","protein","replicate","channel","treatmentName","timePoints","doses")
    idxNull <- which(sapply(list(dfListSCD,protein,replicate,channel,treatmentName,timePoints,doses),is.null))
    stop(paste0("Please provide a "),variables[idxNull]," variable.")}
  if (!facet){
    if (is.null(doses)) {
      stop("Please provide a doses variable.")}}
  
  # If one of the variables is not present in the list, a suitable error is given
  if (is.null(dfListSCD[[protein]]) | is.null(dfListSCD[[protein]][[replicate]]) | 
      is.null(dfListSCD[[protein]][[replicate]][[channel]])){
    stop("The protein, replicate or channel variable is not an entry of the list.")} 
  
  # If all variables are entered, the data frame can be obtained
  df <- dfListSCD[[protein]][[replicate]][[channel]]
  
  # If the treatment under treatmentName, timePoints or doses does not exist, this is reported
  # and the user gets the values it can choose from.
  if (!treatmentName %in% df$treatment){
    warning("Treatment name does not exist. You can choose one of the following treatments:\n",
            paste(levels(as.factor(df$treatment)),collapse = " or "))}
  # If treatment exists, make a local df, called dfSubset,
  # that contains the data for that specific treatment
  else {dfSubset <- df[df$treatment==treatmentName,]}
  
  if (any(!timePoints %in% as.factor(as.character(dfSubset$timeID)))){
    warning("Some or all of the timepoints don't exist. You can choose one of the following timepoints:\n",
            paste(levels(as.factor(sort(as.numeric(as.character(df[df$treatment==treatmentName,"timeID"]))))),collapse = " or "))}
  if (!facet){
    if (any(!doses %in% as.factor(as.character(dfSubset$dose_uM)))){
      warning("Some or all of the doses don't exist. You can choose one of the following doses:\n",
              paste(levels(as.factor(sort(as.numeric(as.character(df[df$treatment==treatmentName,"dose_uM"]))))),collapse = " or "))}}
}

# 15. Function to choose the variable for the density plotting
chooseSCDVariable <- function(which){
  if (which=="nuclei_integrated"){
    xvar <- "Nuclei_Intensity_IntegratedIntensity_image_GFP"
    xlab <- "Integrated intensity in nuclei"} 
  else if (which=="cytoplasm_integrated"){
    xvar <- "Cytoplasm_Intensity_IntegratedIntensity_image_GFP"
    xlab <- "Integrated intensity in cytoplasm"} 
  else if (which=="nuclei_mean") {
    xvar <- "Nuclei_Intensity_MeanIntensity_image_GFP"
    xlab <- "Mean intensity in nuclei"}
  else if (which=="cytoplasm_mean") {
    xvar <- "Cytoplasm_Intensity_MeanIntensity_image_GFP"
    xlab <- "Mean intensity in cytoplasm"}
  else if (which=="nuclei_logintegrated") {
    xvar <- "logNuclei_Intensity_IntegratedIntensity_image_GFP"
    xlab <- "Log(nuclei integrated intensity)"}
  else if (which=="cytoplasm_logintegrated") {
    xvar <- "logCytoplasm_Intensity_IntegratedIntensity_image_GFP"
    xlab <- "Log(cytoplasm integrated intensity)"}
  else if (which=="nuclei_logmean") {
    xvar <- "logNuclei_Intensity_MeanIntensity_image_GFP"
    xlab <- "Log(nuclei mean intensity)"}
  else if (which=="cytoplasm_logmean") {
    xvar <- "logCytoplasm_Intensity_MeanIntensity_image_GFP"
    xlab <- "Log(cytoplams integrated intensity)"}
  else if (which=="area") {
    xvar <- "Nuclei_AreaShape_Area"
    xlab <- "Nuclei area"}
  else if (which=="logarea") {
    xvar <- "logNuclei_AreaShape_Area"
    xlab <- "Log(nuclei area)"}
  xvals <- c(xvar,xlab)
  return(xvals)
}

# 16. 
findCutoff <- function(x, model, proba=0.5, i) {
  ## Cutoff such that Pr[drawn from bad component] == proba
  f <- function(x) {
    proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /
               (model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + 
                  model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
  }
  output <- try(uniroot(f=f, lower=-10, upper=10), silent = T)
  if (inherits(output, "try-error")) {
    return(NA)}
  else {
    return(output$root)}  # Careful with division by zero if changing lower and upper
}

# 17. 
getCutoffs <- function(dfListSCD = NULL, protein = PROTS_OF_INTEREST, replicate = REPLICATES,
                       channel = CHANNELS, treatmentName = NULL,
                       timePoints = NULL, doses = NULL,
                       which = c("nuclei_integrated","nuclei_mean",
                                 "nuclei_logintegrated","nuclei_logmean",
                                 "cytoplasm_integrated","cytoplasm_mean",
                                 "cytoplasm_logintegrated","cytoplasm_logmean",
                                 "area","logarea"),
                       proba = 0.75){
  # Check whether the input is valid
  protein <- match.arg(protein)
  replicate <- match.arg(replicate)
  channel <- match.arg(channel)
  which <- match.arg(which)
  
  # If input is checked, the dataframes can be made
  df <- dfListSCD[[protein]][[replicate]][[channel]]
  
  # Get the data of one specific treatment
  dfSubset1 <- df[df$treatment == treatmentName,]
  
  # Convert timePoints vector to character
  timePoints <- as.character(timePoints)
  
  # Get the data at one specific time point (tp)
  cutoffs <- sapply(timePoints,function(tp){
    df1 <- df[df$timeID==tp,]
    
    # Get the data at one specific dose
    sapply(doses, function(d){
      df2 <- df1[df1$dose_uM==d,]
      xvals <- chooseSCDVariable(which)
      
      # Get cutoffs
      x <- df2[,xvals[1]]
      output <- try(model <- normalmixEM(x=x, k=2),silent = T)
      if (inherits(output, "try-error")) {
        cutoff <- NA}
      else {
        indexLower <- which.min(model$mu)  # Index of component with lower mean
        cutoff <- findCutoff(x = x, model = model, proba=proba, i = indexLower)}
    }, simplify = T)
  }, simplify = T)
  return(cutoffs)}

# 18. 
makeFacetPlot <- function(df, timePoints, which, 
                          doses, protein, altTreatmentName){
  # Convert timePoints vector to numeric
  timePoints <- as.numeric(timePoints)
  
  # Get the x-variable
  xvals <- chooseSCDVariable(which)
  
  # Get the data for the specified timepoints and doses
  dfSubset1 <- df[df$timeID %in% timePoints &
                    df$dose_uM %in% doses,]
  
  gplot <- ggplot(dfSubset1, aes_string(x=xvals[1])) + 
    geom_density() + 
    labs(title=paste0("HepG2 ", protein,"-GFP treated with ", altTreatmentName),
         x=xvals[2], 
         y = "Density") + 
    theme(panel.background = element_blank(),
          plot.title=element_text(size=10))
  
  gplot + facet_grid(dose_uM ~ timeID, scales = "free_y")
}

# 19. 
grid_arrange_shared_legend <- function(ggplotList, 
                                       nrow = 1, ncol = length(ggplotList), 
                                       position = c("bottom", "right"),
                                       xlab = "x axis title",
                                       ylab = "y axis label",
                                       title = "main title",
                                       save, show,
                                       outputPath,
                                       plotWidth,
                                       plotHeight) {
  
  plots <- ggplotList
  position <- match.arg(position)
  for (i in seq_along(plots)){
    g <- ggplotGrob(plots[[i]] + theme(legend.position = position))$grobs
    try(legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]], silent=T)}
  
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- arrangeGrob(grobs = gl,
                    layout_matrix = matrix(1:length(gl),nrow = nrow,ncol = ncol))
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(gl,
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight),
                                            top = textGrob(title, gp=gpar(fontsize=14)),
                                            bottom=textGrob(xlab, gp=gpar(fontsize=12)),
                                            left=textGrob(ylab, gp=gpar(fontsize=12),rot = 90, vjust = 1)),
                     "right" = arrangeGrob(gl,
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth),
                                           top = textGrob(title, gp=gpar(fontsize=14)),
                                           bottom=textGrob(xlab, gp=gpar(fontsize=12)),
                                           left=textGrob(ylab, gp=gpar(fontsize=12),rot = 90, vjust = 1)))
  if (show) {grid.newpage()
    grid.draw(combined)}
  if (save){
    pdf(outputPath, width = plotWidth, height = plotHeight)    
    grid.draw(combined)
    dev.off()}
}

# 20. 
makeGridPlot <- function(df, dfRef, timePoints, which, 
                         doses, protein, refConcentration,
                         altTreatmentName, xlim, proba,
                         save,show,outputPath,
                         plotWidth, plotHeight){
  # Convert timePoints vector to character
  timePoints <- as.character(timePoints)
  
  # Get the x-axis lables
  xvals <- chooseSCDVariable(which)
  
  # If at least one of the doses is present in one of the timepoints, the plots can be made.
  # Otherwise the function is skipped completely
  if (any(doses %in% levels(droplevels(df[df$timeID %in% timePoints,"dose_uM"])))){
    
    # Get the data at one specific time point (tp)
    dfSubset2 <- sapply(timePoints,function(tp){
      df1 <- df[df$timeID==tp,]
      dfRef <- dfRef[dfRef$timeID==tp,]
      
      # Get the data at one specific dose
      sapply(doses, function(d){
        # But only when the dose exists
        if (d %in% levels(droplevels(df1$dose_uM))){
          df2 <- df1[df1$dose_uM==d,]
          dfRef <- dfRef[dfRef$dose_uM==refConcentration,]
          
          # Get cutoffs
          x <- df2[,xvals[1]]
          output <- try(model <- normalmixEM(x=x, k=2),silent = T)
          if (inherits(output, "try-error")) {
            cutoff <- data.frame(cutoff = NA)}
          else {
            indexLower <- which.min(model$mu)  # Index of component with lower mean
            cutoff <- data.frame(cutoff = findCutoff(x = x, model = model, proba=proba, i = indexLower))}
          
          # Make plot
          gplot <- ggplot(NULL) + 
            geom_density(data = df2, aes_string(x=xvals[1])) + 
            labs(title=paste0(d, " uM, t = ", round(df2$timeAfterExposure[1]), " hr"),
                 x="", 
                 y = "") + 
            # Mean intensity
            geom_vline(data = df2, aes_string(xintercept=paste0("mean(",xvals[1],", na.rm = T)")), 
                       color = "blue", linetype="dotted",
                       size=0.75) +
            geom_hline(aes_string(yintercept=0, color = shQuote("mean"), 
                                  linetype=shQuote("mean")), alpha = 0) +
            # Median intensity
            geom_vline(data = df2, aes_string(xintercept=paste0("median(",xvals[1],", na.rm = T)")), 
                       color = "red", linetype="dotted",
                       size=0.75) +
            geom_hline(aes_string(yintercept=0, 
                                  color = shQuote("median"), linetype=shQuote("median")),
                       alpha = 0) +
            # Mean + 2*SD
            geom_vline(data = dfRef, aes_string(xintercept=paste0("mean(",xvals[1],", na.rm = T) + 2*(sd(",xvals[1],",na.rm = T))")), 
                       color = "green", linetype="dashed",
                       size=0.75) +
            geom_hline(aes_string(yintercept=0, 
                                  color = shQuote("threshold "), linetype=shQuote("threshold ")),
                       alpha = 0)
          # Cutoff
          if (!is.na(cutoff$cutoff)){
            gplot <- gplot + geom_vline(data = cutoff, aes_string(xintercept="cutoff"), 
                                        color = "orange", linetype="dashed",
                                        size=0.75) +
              geom_hline(aes_string(yintercept=0, 
                                    color = shQuote("cutoff"), linetype=shQuote("cutoff")),
                         alpha = 0)}
          else {
            gplot <- gplot +  geom_vline(data = data.frame(cutoff = 0), aes_string(xintercept="cutoff"), 
                                         color = "orange", linetype="dashed",
                                         alpha = 0) +
              geom_hline(aes_string(yintercept=0, 
                                    color = shQuote("cutoff"), linetype=shQuote("cutoff")),
                         alpha = 0)}
          
          gplot <- gplot + theme(panel.background = element_blank(),
                                 plot.title=element_text(size=10,hjust = 0.5),
                                 legend.key = element_rect(fill = "white", colour = "white")) +
            scale_color_manual(name = "",values=c("mean"="blue", "median" = "red","threshold " = "green","cutoff" = "orange")) +
            scale_linetype_manual(name = "",values=c("mean"="dotted", "median" = "dotted","threshold " = "dashed","cutoff" = "dashed")) +
            guides(shape = guide_legend(override.aes = list(alpha = 1, size = 1)),
                   colour = guide_legend(override.aes = list(alpha = 1, size = 1)))
          
          if (length(xlim) > 0){
            gplot <- gplot + xlim(xlim)
            gplot}
          gplot
        } else {
          print(paste0("Dose ", d, " does not exist"))
          gplot <- ggplot(NULL) + geom_blank() + 
            theme(panel.background = element_blank(),
                  plot.title=element_text(size=10,hjust = 0.5),
                  legend.key = element_rect(fill = "white", colour = "white")) 
          gplot}
      }, simplify = F)
    }, simplify = T)
    
    grid_arrange_shared_legend(ggplotList = dfSubset2, 
                               nrow = length(doses), 
                               ncol = length(timePoints), 
                               position = "right",
                               xlab = xvals[2],
                               ylab = "Density",
                               title = paste0(protein," expression in HepG2 treated with ",altTreatmentName),
                               save = save,show = show,
                               outputPath = outputPath,
                               plotWidth = plotWidth, plotHeight = plotHeight)}
  else {print(paste0("None of the doses ", paste0(doses,collapse = ", "), " exist at timepoints ", 
                     paste0(timePoints,collapse = ", ")))}
}

# 21. 
plotDensOverTime <- function(dfListSCD = NULL, protein = PROTS_OF_INTEREST, replicate = REPLICATES,
                             channel = CHANNELS, treatmentName = NULL,
                             timePoints = NULL, doses = NULL, facet = F,
                             which = c("nuclei_integrated","nuclei_mean",
                                       "nuclei_logintegrated","nuclei_logmean",
                                       "cytoplasm_integrated","cytoplasm_mean",
                                       "cytoplasm_logintegrated","cytoplasm_logmean",
                                       "area","logarea"), 
                             xlim = numeric(), refTreatment = CONTROL_TREATMENTS,
                             refConcentration = numeric(),
                             proba,save = F,show = F,
                             outputPath = "~/",
                             plotWidth = 5,
                             plotHeight = 5){
  
  # Check whether the input is valid
  if(any(!refTreatment %in% CONTROL_TREATMENTS)){
    warning(paste0("Reference treatment is not one of the control treatments. Choose one of ", 
                   paste0(CONTROL_TREATMENTS, collapse = ", ")))}
  protein <- match.arg(protein)
  replicate <- match.arg(replicate)
  channel <- match.arg(channel)
  which <- match.arg(which)
  checkPlotDensInput(dfListSCD = dfListSCD, protein = protein, replicate = replicate, channel = channel,
                     treatmentName = treatmentName, timePoints = timePoints, doses = doses, facet = facet)
  
  # If all input is checked, the dataframes can be made
  df <- dfListSCD[[protein]][[replicate]][[channel]]
  
  # Get the data of one specific treatment
  dfSubset1 <- df[df$treatment == treatmentName,]
  dfRef <- df[df$treatment == refTreatment,]
  
  # Get alternative treatment name
  altTreatmentName <- ALT_COMP_NAMES[which(COMPS_OF_INTEREST == treatmentName)]
  
  # Make a grid or facet plot 
  if (facet){
    makeFacetPlot(df = dfSubset1, timePoints = timePoints, which = which, 
                  doses = doses, protein = protein, altTreatmentName = altTreatmentName)}
  else if (!facet){
    makeGridPlot(df = dfSubset1, dfRef = dfRef, timePoints = timePoints, which = which, 
                 doses = doses, protein = protein, refConcentration = refConcentration,
                 altTreatmentName = altTreatmentName, xlim = xlim, proba = proba,
                 save = save,show = show,
                 outputPath = outputPath,plotWidth = plotWidth, plotHeight = plotHeight)}
}

# plotDensOverTime <- function(dfListSCD = NULL, protein = PROTS_OF_INTEREST, replicate = REPLICATES,
#                              channel = CHANNELS, treatmentName = NULL,
#                              timePoints = NULL, doses = NULL, facet = F,
#                              which = c("nuclei_integrated","nuclei_mean",
#                                        "nuclei_logintegrated","nuclei_logmean",
#                                        "cytoplasm_integrated","cytoplasm_mean",
#                                        "cytoplasm_logintegrated","cytoplasm_logmean",
#                                        "area","logarea"), 
#                              xlim = numeric(), refTreatment = CONTROL_TREATMENTS){
#   # Check whether the input is valid
#   refTreatment <- if(any(!refTreatment %in% CONTROL_TREATMENTS)){
#     warning(paste0("Reference treatment is not one of the control treatments. Choose one of ", paste0(CONTROL_TREATMENTS, collapse = ", ")))}
#   protein <- match.arg(protein)
#   replicate <- match.arg(replicate)
#   channel <- match.arg(channel)
#   which <- match.arg(which)
#   checkPlotDensInput(dfListSCD = dfListSCD, protein = protein, replicate = replicate, channel = channel,
#                      treatmentName = treatmentName, timePoints = timePoints, doses = doses, facet = facet)
#   
#   # If all input is checked, the dataframes can be made
#   df <- dfListSCD[[protein]][[replicate]][[channel]]
#   
#   # Get the data of one specific treatment
#   dfSubset1 <- df[df$treatment == treatmentName,]
#   dfRef <- df[df$treatment == refTreatment,]
#   
#   # Get alternative treatment name
#   altTreatmentName <- ALT_COMP_NAMES[which(COMPS_OF_INTEREST == treatmentName)]
#   
#   # Make a grid or facet plot
#   if (facet){
#     # Convert timePoints vector to numeric
#     timePoints <- as.numeric(timePoints)
#     
#     # Get the x-variable
#     xvals <- chooseSCDVariable(which)
#     
#     # Get the data for the specified timepoints
#     dfSubset1 <- dfSubset1[which(dfSubset1$timeID %in% timePoints),]
#     
#     gplot <- ggplot(dfSubset1, aes_string(x=xvals[1])) + 
#       geom_density() + 
#       labs(title=paste0("HepG2 ", protein,"-GFP treated with ", altTreatmentName),
#            x=xvals[2], 
#            y = "Density") + 
#       theme(panel.background = element_blank(),
#             plot.title=element_text(size=10))
#     
#     gplot + facet_grid(dose_uM ~ timeID, scales = "free_y")}
#   else if (!facet){
#     # Convert timePoints vector to character
#     timePoints <- as.character(timePoints)
#     
#     # Get the data at one specific time point (tp)
#     dfSubset2 <- sapply(timePoints,function(tp){
#       df1 <- dfSubset1[dfSubset1$timeID==tp,]
#       # Get the data at one specific dose
#       sapply(doses, function(d){
#         df2 <- df1[df1$dose_uM==d,]
#         xvals <- chooseSCDVariable(which)
#         gplot <- ggplot(df2, aes_string(x=xvals[1])) + 
#           geom_density() + 
#           labs(title=paste0(protein,", ",d, " uM ", altTreatmentName, ", t = ", tp),
#                x=xvals[2], 
#                y = "Density") + 
#           geom_vline(aes_string(xintercept=paste0("mean(",xvals[1],", na.rm = T)")), color="blue", linetype="dashed", size=1) +
#           geom_vline(aes_string(xintercept=paste0("median(",xvals[1],", na.rm = T)")), color="red", linetype="dashed", size=1) +
#           geom_vline() +
#           theme(panel.background = element_blank(),
#                 plot.title=element_text(size=10))
#         
#         if (length(xlim) > 0){
#           gplot <- gplot + xlim(xlim)
#           gplot}
#         gplot
#           
#       }, simplify = F)
#     }, simplify = T)
#     
#     # Now make one plot including all the densityplots
#     grid.arrange(grobs = dfSubset2, ncol = length(timePoints), as.table = FALSE)}
# }

# 22. 
plotHistOverTime <- function(dfListSCD = NULL, protein = PROTS_OF_INTEREST, replicate = REPLICATES,
                             channel = CHANNELS, treatmentName = NULL,
                             timePoints = NULL, doses = NULL,
                             which = c("nuclei_integrated","nuclei_mean",
                                       "nuclei_logintegrated","nuclei_logmean",
                                       "cytoplasm_integrated","cytoplasm_mean",
                                       "cytoplasm_logintegrated","cytoplasm_logmean",
                                       "area","logarea"), 
                             lower.bound = 0,
                             upper.bound = 50,
                             step = 1){
  # Check whether the input is valid
  protein <- match.arg(protein)
  replicate <- match.arg(replicate)
  channel <- match.arg(channel)
  which <- match.arg(which)
  checkPlotDensInput(dfListSCD = dfListSCD, protein = protein, replicate = replicate, channel = channel,
                     treatmentName = treatmentName, timePoints = timePoints, doses = doses)
  
  # If all input is checked, the dataframes can be made
  df <- dfListSCD[[protein]][[replicate]][[channel]]
  
  # Get the data of one specific treatment
  dfSubset1 <- df[df$treatment == treatmentName,]
  
  # Get alternative treatment name
  altTreatmentName <- ALT_COMP_NAMES[which(COMPS_OF_INTEREST == treatmentName)]
  
  # Convert timePoints vector to character
  timePoints <- as.character(timePoints)
  
  # Get the data at one specific time point (tp)
  dfSubset2 <- sapply(timePoints,function(tp){
    df1 <- dfSubset1[dfSubset1$timeID==tp,]
    # Get the data at one specific dose
    sapply(doses, function(d){
      df2 <- df1[df1$dose_uM==d,]
      xvals <- chooseSCDVariable(which)
      gplot <- ggplot(data=df2,aes_string(x=xvals[1])) + 
        geom_histogram(breaks=seq(lower.bound, upper.bound, by = step), 
                       col="black", 
                       fill="grey", 
                       alpha = .4) + 
        labs(title=paste0(protein,", ",d, " uM ", altTreatmentName, ", t = ", tp),
             x=xvals[2], 
             y = "Density") +
        geom_vline(aes_string(xintercept=paste0("mean(",xvals[1],", na.rm = T)")), color="blue", linetype="dashed", size=1) +
        geom_vline(aes_string(xintercept=paste0("median(",xvals[1],", na.rm = T)")), color="red", linetype="dashed", size=1) +
        theme(panel.background = element_blank(),
              plot.title=element_text(size=10)) +
        xlim(lower.bound,upper.bound)
      
    }, simplify = F)
  }, simplify = T)
  
  # Now make one plot including all the densityplots
  grid.arrange(grobs = dfSubset2, ncol = length(timePoints), as.table = FALSE)
}

# # 23. Function to make summary dataframes
# # Input is a list of data frames that you want to be merged with names corresponding to channel names
# makeSummary <- function(listOfDFs, channel = CHANNELS){
#   # Check whether channel arguments exist
#   channel <- match.arg(channel)
#   
#   # Make a summary data frame
#   if (channel == CHANNELS[1]){
#     # Make a summary for the GFP data frame
#     ACD <- ddply(listOfDFs[[channel]], .(treatment,dose_uM,cell_line,plateID,
#                                          replID,timeID,timeAfterExposure),summarize,
#                  well = strsplit(as.character(locationID[1]),"_")[[1]][1],
#                  imageCountParentObj = mean(as.numeric(levels(as.factor(imageCountParentObj)))),
#                  GMCytoplasm_Mean_Intensity = geomMean(Cytoplasm_Intensity_MeanIntensity_image_GFP),
#                  GSDCytoplasm_Mean_Intensity = geomSD(Cytoplasm_Intensity_MeanIntensity_image_GFP),
#                  GMCytoplasm_Integrated_Intensity = geomMean(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),
#                  GSDCytoplasm_Integrated_Intensity = geomSD(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),
#                  GMNuclei_Mean_Intensity = geomMean(Nuclei_Intensity_MeanIntensity_image_GFP),
#                  GSDNuclei_Mean_Intensity = geomSD(Nuclei_Intensity_MeanIntensity_image_GFP),
#                  GMNuclei_Integrated_Intensity = geomMean(Nuclei_Intensity_IntegratedIntensity_image_GFP),
#                  GSDNuclei_Integrated_Intensity = geomSD(Nuclei_Intensity_IntegratedIntensity_image_GFP),
#                  Nuclei_AreaShape_Area = mean(Nuclei_AreaShape_Area),
#                  propNonRespCytoplasm_Mean_Intensity = 1 - (ifelse(!all(is.na(response_Cytoplasm_Intensity_MeanIntensity_image_GFP)),
#                                                                    sum(response_Cytoplasm_Intensity_MeanIntensity_image_GFP,na.rm = T),NA)/
#                                                               length(response_Cytoplasm_Intensity_MeanIntensity_image_GFP)),
#                  propNonRespCytoplasm_Integrated_Intensity = 1 - (ifelse(!all(is.na(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP)),
#                                                                          sum(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP,na.rm = T),NA)/
#                                                                     length(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP)),
#                  propNonRespNuclei_Mean_Intensity = 1 - (ifelse(!all(is.na(response_Nuclei_Intensity_MeanIntensity_image_GFP)),
#                                                                 sum(response_Nuclei_Intensity_MeanIntensity_image_GFP,na.rm = T),NA)/
#                                                            length(response_Nuclei_Intensity_MeanIntensity_image_GFP)),
#                  propNonRespNuclei_Integrated_Intensity = 1 - (ifelse(!all(is.na(response_Nuclei_Intensity_IntegratedIntensity_image_GFP)),
#                                                                       sum(response_Nuclei_Intensity_IntegratedIntensity_image_GFP,na.rm = T),NA)/
#                                                                  length(response_Nuclei_Intensity_IntegratedIntensity_image_GFP)))}
#   else if (channel == CHANNELS[2]){
#     # Make a summary for the AnVPI data frame
#     ACD <- ddply(listOfDFs[[channel]], .(treatment,dose_uM,cell_line,plateID,
#                                          replID,timeID,timeAfterExposure),summarize,
#                  well = strsplit(as.character(locationID[1]),"_")[[1]][1],
#                  cellCount_AnVPI = as.numeric(levels(as.factor(imageCountParentObj)))[1],
#                  propCount_PI = sum(as.numeric(PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area>0.1)/
#                                       sum(as.numeric(levels(as.factor(imageCountParentObj)))), na.rm = T),
#                  propCount_AnV = sum(as.numeric(AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area>0.1)/
#                                        sum(as.numeric(levels(as.factor(imageCountParentObj)))), na.rm = T))}
#   return(ACD)
# }

# convertColumnsDf <- function(df, columnNames, to = "factor"){
#   if (to == "numeric"){
#     reps[[chanName]][,columnNames] <- sapply(columnNames, function(x) as.numeric(as.character(reps[[chanName]][[x]])))
#     reps[[chanName]]}
#   else if (to == "character"){
#     reps[[chanName]][,columnNames] <- sapply(columnNames, function(x) as.character(reps[[chanName]][[x]]))
#     reps[[chanName]]}
#   else if (to == "factor"){
#     if (length(columnNames) > 1){
#       reps[[chanName]][,columnNames] <- sapply(columnNames, function(x) as.factor(as.character(reps[[chanName]][[x]])))
#       reps[[chanName]]}
#     else {
#       reps[[chanName]][,columnNames] <- as.factor(as.character(reps[[chanName]][[columnNames]]))
#       reps[[chanName]]}}
#   else if (to == "numericFactor"){
#     if (length(columnNames) > 1){
#       reps[[chanName]][,columnNames] <- sapply(columnNames, function(x) as.factor(as.numeric(as.character(reps[[chanName]][[x]]))))
#       reps[[chanName]]}
#     else {
#       reps[[chanName]][,columnNames] <- as.factor(as.numeric(as.character(reps[[chanName]][[columnNames]])))
#       reps[[chanName]]}}}

# 23a. Function to make summary dataframes
# Input is a list of data frames that you want to be merged with names corresponding to channel names
makeSummary <- function(listOfDFs, channel = CHANNELS){
  # Check whether channel arguments exist
  channel <- match.arg(channel)
  
  # Make a summary data frame
  if (channel == CHANNELS[1]){
    # Make a summary for the GFP data frame
    ACD <- listOfDFs[[channel]] %>% 
      group_by(treatment,dose_uM,cell_line,plateID,replID,timeID,timeAfterExposure) %>%
      summarise(well = as.factor(strsplit(as.character(locationID[1]),"_")[[1]][1]),
                imageCountParentObj = mean(as.numeric(levels(as.factor(imageCountParentObj)))),
                GMCytoplasm_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                    geomMean(Cytoplasm_Intensity_MeanIntensity_image_GFP),NA),
                GSDCytoplasm_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                     geomSD(Cytoplasm_Intensity_MeanIntensity_image_GFP),NA),
                GMCytoplasm_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                          geomMean(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),NA),
                GSDCytoplasm_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                           geomSD(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),NA),
                GMNuclei_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                 geomMean(Nuclei_Intensity_MeanIntensity_image_GFP),NA),
                GSDNuclei_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                  geomSD(Nuclei_Intensity_MeanIntensity_image_GFP),NA),
                GMNuclei_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                       geomMean(Nuclei_Intensity_IntegratedIntensity_image_GFP),NA),
                GSDNuclei_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                        geomSD(Nuclei_Intensity_IntegratedIntensity_image_GFP),NA),
                GM10Cytoplasm_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                      geomMean10(Cytoplasm_Intensity_MeanIntensity_image_GFP),NA),
                GSD10Cytoplasm_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                       geomSD10(Cytoplasm_Intensity_MeanIntensity_image_GFP),NA),
                GM10Cytoplasm_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                            geomMean10(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),NA),
                GSD10Cytoplasm_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                             geomSD10(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),NA),
                GM10Nuclei_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                   geomMean10(Nuclei_Intensity_MeanIntensity_image_GFP),NA),
                GSD10Nuclei_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                    geomSD10(Nuclei_Intensity_MeanIntensity_image_GFP),NA),
                GM10Nuclei_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                         geomMean10(Nuclei_Intensity_IntegratedIntensity_image_GFP),NA),
                GSD10Nuclei_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                          geomSD10(Nuclei_Intensity_IntegratedIntensity_image_GFP),NA),
                Nuclei_AreaShape_Area = mean(Nuclei_AreaShape_Area),
                propNonRespCytoplasm_Mean_Intensity = 1 - (ifelse(!all(is.na(response_Cytoplasm_Intensity_MeanIntensity_image_GFP)),
                                                                  sum(response_Cytoplasm_Intensity_MeanIntensity_image_GFP,na.rm = T),NA)/
                                                             length(response_Cytoplasm_Intensity_MeanIntensity_image_GFP)),
                propNonRespCytoplasm_Integrated_Intensity = 1 - (ifelse(!all(is.na(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP)),
                                                                        sum(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP,na.rm = T),NA)/
                                                                   length(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP)),
                propNonRespNuclei_Mean_Intensity = 1 - (ifelse(!all(is.na(response_Nuclei_Intensity_MeanIntensity_image_GFP)),
                                                               sum(response_Nuclei_Intensity_MeanIntensity_image_GFP,na.rm = T),NA)/
                                                          length(response_Nuclei_Intensity_MeanIntensity_image_GFP)),
                propNonRespNuclei_Integrated_Intensity = 1 - (ifelse(!all(is.na(response_Nuclei_Intensity_IntegratedIntensity_image_GFP)),
                                                                     sum(response_Nuclei_Intensity_IntegratedIntensity_image_GFP,na.rm = T),NA)/
                                                                length(response_Nuclei_Intensity_IntegratedIntensity_image_GFP)))}
  else if (channel == CHANNELS[2]){
    # Make a summary for the AnVPI data frame
    ACD <- listOfDFs[[channel]] %>% 
      group_by(treatment,dose_uM,cell_line,plateID,replID,timeID,timeAfterExposure) %>% 
      summarise(well = as.factor(strsplit(as.character(locationID[1]),"_")[[1]][1]),
                 cellCount_AnVPI = as.numeric(levels(as.factor(imageCountParentObj)))[1],
                 propCount_PI = sum(as.numeric(PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area>0.1)/
                                      sum(as.numeric(levels(as.factor(imageCountParentObj)))), na.rm = T),
                 propCount_AnV = sum(as.numeric(AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area>0.1)/
                                       sum(as.numeric(levels(as.factor(imageCountParentObj)))), na.rm = T))}
  return(ACD)
}

# 23b. Function to make summary dataframes
# Input is a list of data frames that you want to be merged with names corresponding to channel names
makeSummaryPerImage <- function(listOfDFs, channel = CHANNELS){
  # Check whether channel arguments exist
  channel <- match.arg(channel)
  
  # Make a summary data frame
  if (channel == CHANNELS[1]){
    # Make a summary for the GFP data frame
    ACD <- listOfDFs[[channel]] %>% 
      group_by(treatment,dose_uM,cell_line,plateID,replID,timeID,timeAfterExposure,well,locationID) %>%
      summarise(imageCountParentObj = mean(as.numeric(levels(as.factor(imageCountParentObj)))),
                GMCytoplasm_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                    geomMean(Cytoplasm_Intensity_MeanIntensity_image_GFP),NA),
                GSDCytoplasm_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                     geomSD(Cytoplasm_Intensity_MeanIntensity_image_GFP),NA),
                GMCytoplasm_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                          geomMean(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),NA),
                GSDCytoplasm_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                           geomSD(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),NA),
                GMNuclei_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                 geomMean(Nuclei_Intensity_MeanIntensity_image_GFP),NA),
                GSDNuclei_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                  geomSD(Nuclei_Intensity_MeanIntensity_image_GFP),NA),
                GMNuclei_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                       geomMean(Nuclei_Intensity_IntegratedIntensity_image_GFP),NA),
                GSDNuclei_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                        geomSD(Nuclei_Intensity_IntegratedIntensity_image_GFP),NA),
                GM10Cytoplasm_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                      geomMean10(Cytoplasm_Intensity_MeanIntensity_image_GFP),NA),
                GSD10Cytoplasm_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                       geomSD10(Cytoplasm_Intensity_MeanIntensity_image_GFP),NA),
                GM10Cytoplasm_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                            geomMean10(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),NA),
                GSD10Cytoplasm_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                             geomSD10(Cytoplasm_Intensity_IntegratedIntensity_image_GFP),NA),
                GM10Nuclei_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                   geomMean10(Nuclei_Intensity_MeanIntensity_image_GFP),NA),
                GSD10Nuclei_Mean_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                    geomSD10(Nuclei_Intensity_MeanIntensity_image_GFP),NA),
                GM10Nuclei_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                         geomMean10(Nuclei_Intensity_IntegratedIntensity_image_GFP),NA),
                GSD10Nuclei_Integrated_Intensity = ifelse(sum(is.na(Cytoplasm_Intensity_MeanIntensity_image_GFP))/length(Cytoplasm_Intensity_MeanIntensity_image_GFP) < 0.8, 
                                                          geomSD10(Nuclei_Intensity_IntegratedIntensity_image_GFP),NA),
                Nuclei_AreaShape_Area = mean(Nuclei_AreaShape_Area),
                propNonRespCytoplasm_Mean_Intensity = 1 - (ifelse(!all(is.na(response_Cytoplasm_Intensity_MeanIntensity_image_GFP)),
                                                                  sum(response_Cytoplasm_Intensity_MeanIntensity_image_GFP,na.rm = T),NA)/
                                                             length(response_Cytoplasm_Intensity_MeanIntensity_image_GFP)),
                propNonRespCytoplasm_Integrated_Intensity = 1 - (ifelse(!all(is.na(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP)),
                                                                        sum(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP,na.rm = T),NA)/
                                                                   length(response_Cytoplasm_Intensity_IntegratedIntensity_image_GFP)),
                propNonRespNuclei_Mean_Intensity = 1 - (ifelse(!all(is.na(response_Nuclei_Intensity_MeanIntensity_image_GFP)),
                                                               sum(response_Nuclei_Intensity_MeanIntensity_image_GFP,na.rm = T),NA)/
                                                          length(response_Nuclei_Intensity_MeanIntensity_image_GFP)),
                propNonRespNuclei_Integrated_Intensity = 1 - (ifelse(!all(is.na(response_Nuclei_Intensity_IntegratedIntensity_image_GFP)),
                                                                     sum(response_Nuclei_Intensity_IntegratedIntensity_image_GFP,na.rm = T),NA)/
                                                                length(response_Nuclei_Intensity_IntegratedIntensity_image_GFP)))}
  else if (channel == CHANNELS[2]){
    # Make a summary for the AnVPI data frame
    ACD <- listOfDFs[[channel]] %>% 
      group_by(treatment,dose_uM,cell_line,plateID,replID,timeID,timeAfterExposure,well,locationID) %>% 
      summarise(cellCount_AnVPI = as.numeric(levels(as.factor(imageCountParentObj)))[1],
                propCount_PI = sum(as.numeric(PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area>0.1)/
                                     sum(as.numeric(levels(as.factor(imageCountParentObj)))), na.rm = T),
                propCount_AnV = sum(as.numeric(AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area>0.1)/
                                      sum(as.numeric(levels(as.factor(imageCountParentObj)))), na.rm = T))}
  return(ACD)
}

# 24a. Function to combine the summary data frames of the GFP and AnVPI channels into one data frame
combineChannelSummaries <- function(summaryDFList, proteinName){
  # Function to combine two or more summary data frames stored as list
  # Check whether the list is not empty
  if (length(summaryDFList) != 0){
    summaryData <- summaryDFList[[1]]
    
    # Check whether there is more than 1 data frame
    if (length(summaryDFList) > 1){
      # Join the data frames into one
      for (dfIdx in 2:length(summaryDFList)){
        summaryData <- left_join(summaryData, summaryDFList[[dfIdx]], 
                                 by=c("treatment","dose_uM","cell_line","plateID",
                                      "replID","timeID","timeAfterExposure","well"))}}
    
    # Add a column that contains the protein name
    summaryData[["protein"]] <- proteinName

    return(summaryData)}
  else {
    print("Empty list of summary data frames")}
}

# 24b. Function to combine the summary data frames of the GFP and AnVPI channels into one data frame
combineChannelSummariesPerImage <- function(summaryDFList, proteinName){
  # Function to combine two or more summary data frames stored as list
  # Check whether the list is not empty
  if (length(summaryDFList) != 0){
    summaryData <- summaryDFList[[1]]
    
    # Check whether there is more than 1 data frame
    if (length(summaryDFList) > 1){
      # Join the data frames into one
      for (dfIdx in 2:length(summaryDFList)){
        summaryData <- left_join(summaryData, summaryDFList[[dfIdx]], 
                                 by=c("treatment","dose_uM","cell_line","plateID",
                                      "replID","timeID","timeAfterExposure","well","locationID"))}}
    
    # Add a column that contains the protein name
    summaryData[["protein"]] <- proteinName
    
    return(summaryData)}
  else {
    print("Empty list of summary data frames")}
}

# 25a. Function to join the list of GFP and AnVPI data frames
joinChannelDFs <- function(listOfDFs, channels, proteinName){
  # Function to join the two channels (GFP and AnVPI) in to one summary data frame
  # Input is a list of data frames that have to be merged with names corresonding to channel names
  # First make a summary per channel and store them in a list 
  summaryDFList <- sapply(channels,function(chanName){
    # Check whether the channel exists
    if (chanName %in% names(listOfDFs)){
      chanName = makeSummary(listOfDFs, channel = chanName)}}, simplify = F, USE.NAMES = T)
  
  # Remove entries of summaryDFList[[channel]] for which there is no entry
  summaryDFList <- summaryDFList[unlist(lapply(summaryDFList,length) != 0)]
  
  # Then paste the data frames together into one data frame according to corresponding columns
  summaryDF <- combineChannelSummaries(summaryDFList, proteinName)
  return(summaryDF)
}

# 25b. Function to join the list of GFP and AnVPI data frames
joinChannelDFsPerImage <- function(listOfDFs, channels, proteinName){
  # Function to join the two channels (GFP and AnVPI) in to one summary data frame
  # Input is a list of data frames that have to be merged with names corresonding to channel names
  # First make a summary per channel and store them in a list 
  summaryDFList <- sapply(channels,function(chanName){
    # Check whether the channel exists
    if (chanName %in% names(listOfDFs)){
      chanName = makeSummaryPerImage(listOfDFs, channel = chanName)}}, simplify = F, USE.NAMES = T)
  
  # Remove entries of summaryDFList[[channel]] for which there is no entry
  summaryDFList <- summaryDFList[unlist(lapply(summaryDFList,length) != 0)]
  
  # Then paste the data frames together into one data frame according to corresponding columns
  summaryDF <- combineChannelSummariesPerImage(summaryDFList, proteinName)
  return(summaryDF)
}

# 26a. Makes a list of the summary data files
makeDFListSumCD <- function(dfListSCD, channels){
  # This function makes a list, called dfListSumCD, 
  # of which each element has the name of a protein and is a list ...
  # of which each element has the name of the replicate and is a data frame of the summary data for that replicate.
  dfListSumCD <- sapply(names(dfListSCD),function(protName){
    sapply(dfListSCD[[protName]],function(reps){
      joinChannelDFs(reps,channels,proteinName = protName)
    },simplify = F, USE.NAMES = T)
  },simplify = F, USE.NAMES = T)
  return(dfListSumCD)
}

# 26b. Makes a list of the summary data files
makeDFListSumCDPerImage <- function(dfListSCD, channels){
  # This function makes a list, called dfListSumCD, 
  # of which each element has the name of a protein and is a list ...
  # of which each element has the name of the replicate and is a data frame of the summary data for that replicate.
  dfListSumCD <- sapply(names(dfListSCD),function(protName){
    sapply(dfListSCD[[protName]],function(reps){
      joinChannelDFsPerImage(reps,channels,proteinName = protName)
    },simplify = F, USE.NAMES = T)
  },simplify = F, USE.NAMES = T)
  return(dfListSumCD)
}

# 27a. Combines the separate replicate data frames into one data frame
combineReplicateSummaries <- function(dfListSCD, channels){
  # Function to combine two or more summary data frames stored as list in a list
  # First make the dfListSumCD list, i.e. a list of summaryDFs per replicate
  dfListSumCD <- makeDFListSumCD(dfListSCD, CHANNELS)
  
  # Then row-bind the data frames
  sapply(dfListSumCD,function(summaryDFList){
    # Check whether the list is not empty
    if (length(summaryDFList) != 0){
      summaryData <- summaryDFList[[1]]
      
      # Check whether there is more than 1 data frame
      if (length(summaryDFList) > 1){
        # Join the data frames into one
        for (dfIdx in 2:length(summaryDFList)){
          # Make sure that the columns are of the right class
          summaryData <- bind_rows(summaryData, summaryDFList[[dfIdx]])}}
      
      return(summaryData)}
    else {
      print("Empty list of summary data frames")}
  },simplify = F, USE.NAMES = T)
}

# 27b. Combines the separate replicate data frames into one data frame
combineReplicateSummariesPerImage <- function(dfListSCD, channels){
  # Function to combine two or more summary data frames stored as list in a list
  # First make the dfListSumCD list, i.e. a list of summaryDFs per replicate
  dfListSumCD <- makeDFListSumCDPerImage(dfListSCD, CHANNELS)
  
  # Then row-bind the data frames
  sapply(dfListSumCD,function(summaryDFList){
    # Check whether the list is not empty
    if (length(summaryDFList) != 0){
      summaryData <- summaryDFList[[1]]
      
      # Check whether there is more than 1 data frame
      if (length(summaryDFList) > 1){
        # Join the data frames into one
        for (dfIdx in 2:length(summaryDFList)){
          # Make sure that the columns are of the right class
          summaryData <- bind_rows(summaryData, summaryDFList[[dfIdx]])}}
      
      return(summaryData)}
    else {
      print("Empty list of summary data frames")}
  },simplify = F, USE.NAMES = T)
}

# # 28. Combines the separate protein data frames into one data frame
# combineProteinSummaries <- function(dfListSCD, channels = CHANNELS){
#   # Function to combine two or more summary data frames stored as list in a list
#   # First make the dfListSumCD list, i.e. a list of summaryDFs per replicate
#   dfListSumCD <- combineReplicateSummaries(dfListSCD, CHANNELS)
#   
#   # Then row-bind the data frames
#   # Check whether the list is not empty
#   if (length(dfListSumCD) != 0){
#     summaryData <- dfListSumCD[[1]]
#     
#     # Check whether there is more than 1 data frame
#     if (length(dfListSumCD) > 1){
#       # Join the data frames into one
#       for (dfIdx in 2:length(dfListSumCD)){
#         summaryData <- rbind.fill(summaryData, dfListSumCD[[dfIdx]])}}
#     
#     return(summaryData)}
#   else {
#     print("Empty list of summary data frames")}
# }
# df1 <- data.frame(x = seq(1,5), y = seq(2,6))
# df2 <- data.frame(x = seq(1,5), y = seq(2,6), z = seq(11,15))
# bind_rows(df1,df2)

# 28a. Combines the separate protein data frames into one data frame
combineProteinSummaries <- function(dfListSCD, channels = CHANNELS){
  # Function to combine two or more summary data frames stored as list in a list
  # First make the dfListSumCD list, i.e. a list of summaryDFs per replicate
  dfListSumCD <- combineReplicateSummaries(dfListSCD, CHANNELS)
  
  # Then row-bind the data frames
  # Check whether the list is not empty
  if (length(dfListSumCD) != 0){
    summaryData <- dfListSumCD[[1]]
    
    # Check whether there is more than 1 data frame
    if (length(dfListSumCD) > 1){
      # Join the data frames into one
      for (dfIdx in 2:length(dfListSumCD)){
        summaryData <- bind_rows(summaryData, dfListSumCD[[dfIdx]])}}
    
    return(summaryData)}
  else {
    print("Empty list of summary data frames")}
}

# 28b. Combines the separate protein data frames into one data frame
combineProteinSummariesPerImage <- function(dfListSCD, channels = CHANNELS){
  # Function to combine two or more summary data frames stored as list in a list
  # First make the dfListSumCD list, i.e. a list of summaryDFs per replicate
  dfListSumCD <- combineReplicateSummariesPerImage(dfListSCD, CHANNELS)
  
  # Then row-bind the data frames
  # Check whether the list is not empty
  if (length(dfListSumCD) != 0){
    summaryData <- dfListSumCD[[1]]
    
    # Check whether there is more than 1 data frame
    if (length(dfListSumCD) > 1){
      # Join the data frames into one
      for (dfIdx in 2:length(dfListSumCD)){
        summaryData <- bind_rows(summaryData, dfListSumCD[[dfIdx]])}}
    
    return(summaryData)}
  else {
    print("Empty list of summary data frames")}
}

# 29. Function to check the input data of plotDynOverTime
checkPlotDynInput <- function(df,protein,replID,treatmentName,doses, combine.by){
  # Function that does some error handling. It checks the input for the plots and 
  # provides suitable errors.
  
  # If one of the variables is not provided, a suitable error is given
  if (is.null(df) | is.null(treatmentName)){
    variables <- c("df","treatmentName")
    idxNull <- which(sapply(list(df,treatmentName),is.null))
    stop(paste0("Please provide a "),variables[idxNull]," variable.")}
  if (combine.by == "none"){
    if (is.null(replID) | is.null(doses)){
      variables <- c("replID","doses")
      idxNull <- which(sapply(list(replID,doses),is.null))
      stop(paste0("Please provide a "),variables[idxNull]," variable.")}
    if (length(protein) > 1){
      stop("Please provide the name of only one protein")}} 
  
  # If the treatment under treatmentName, timePoints or doses does not exist, this is reported
  # and the user gets the values it can choose from
  if (!treatmentName %in% df$treatment){
    stop(paste0(treatmentName," does not exist. You can choose one of the following treatments:\n",
         paste(levels(as.factor(df$treatment)),collapse = " or ")))}
  else if (!protein %in% df$protein){
    stop("Protein name does not exist. You can choose one of the following proteins:\n",
         paste(levels(as.factor(df$protein)),collapse = " or "))}
  
  # If the plots will not be combined, the replID and doses variables are checked
  if (combine.by == "none"){
    if (!replID %in% levels(as.factor(as.character(df[df$protein == protein,"replID"])))){
      stop("replID does not exist. You can choose one of the following replIDs:\n",
           paste(levels(as.factor(as.character(df[df$protein == protein,"replID"]))),collapse = " or "))}
    # If treatment, replicate and protein exists, make a new df, called dfSubset, 
    # that contains the data for that specific treatment
    else {dfSubset <- df[df$protein == protein & df$treatment == treatmentName & df$replID == replID,]}
    
    if (any(!doses %in% as.factor(as.character(dfSubset$dose_uM)))){
      stop("Some or all of the doses don't exist. You can choose one of the following doses:\n",
           paste(levels(as.factor(sort(as.numeric(as.character(df[df$treatment==treatmentName,"dose_uM"]))))),collapse = " or "))}}
}

# 30. Function to choose the variable for the dynamics plotting
chooseACDVariable <- function(which){
  if (which=="nuclei_integrated"){
    var <- "GM10Nuclei_Integrated_Intensity"
    varsd <- "GSD10Nuclei_Integrated_Intensity"
    lab <- "Integrated intensity in nuclei (a.u.)"}
  else if (which=="nuclei_mean") {
    var <- "GM10Nuclei_Mean_Intensity"
    varsd <- "GSD10Nuclei_Mean_Intensity"
    lab <- "Mean intensity in nuclei (a.u.)"}
  else if (which=="cytoplasm_integrated") {
    var <- "GM10Cytoplasm_Integrated_Intensity"
    varsd <- "GSD10Cytoplasm_Integrated_Intensity"
    lab <- "Integrated intensity in cytoplasm (a.u.)"}
  else if (which=="cytoplasm_mean") {
    var <- "GM10Cytoplasm_Mean_Intensity"
    varsd <- "GSD10Cytoplasm_Mean_Intensity"
    lab <- "Mean intensity in cytoplasm (a.u.)"}
  else if (which=="nuclei_integrated_norm"){
    var <- "GM10Nuclei_Integrated_Intensity_Norm"
    varsd <- "GSD10Nuclei_Integrated_Intensity_Norm"
    lab <- "Integrated intensity in nuclei (a.u.)"}
  else if (which=="nuclei_mean_norm") {
    var <- "GM10Nuclei_Mean_Intensity_Norm"
    varsd <- "GSD10Nuclei_Mean_Intensity_Norm"
    lab <- "Mean intensity in nuclei (a.u.)"}
  else if (which=="cytoplasm_integrated_norm") {
    var <- "GM10Cytoplasm_Integrated_Intensity_Norm"
    varsd <- "GSD10Cytoplasm_Integrated_Intensity_Norm"
    lab <- "Integrated intensity in cytoplasm (a.u.)"}
  else if (which=="cytoplasm_mean_norm") {
    var <- "GM10Cytoplasm_Mean_Intensity_Norm"
    varsd <- "GSD10Cytoplasm_Mean_Intensity_Norm"
    lab <- "Mean intensity in cytoplasm (a.u.)"}
  else if (which=="propPI") {
    var <- "propCount_PI"
    varsd <- "NULL"
    lab <- "Proportion of PI positive cells"}
  else if (which=="propAnV") {
    var <- "propCount_AnV"
    varsd <- "NULL"
    lab <- "Proportion of AnV positive cells"}
  else if (which=="cellCount_GFP") {
    var <- "imageCountParentObj"
    varsd <- "NULL"
    lab <- "Number of cells"}
  else if (which=="cellCount_AnVPI") {
    var <- "cellCount_AnVPI"
    varsd <- "NULL"
    lab <- "Number of cells"}
  else if (which=="relativeCellCount") {
    var <- "relativeCellCount"
    varsd <- "NULL"
    lab <- "Number of cells"}
  else if (which=="nuclei_area") {
    var <- "Nuclei_AreaShape_Area"
    varsd <- "NULL"
    lab <- "Average nuclear area"}
  vals <- c(var,varsd,lab)
  return(vals)
}

# 31. Function to make dynamic plots per protein, treatment, replicate and dose
makeSeparateDynPlots <- function(df, protein, replID,treatmentName, doses, which, 
                                 withSD, altTreatmentName, yvals, withControl,logScale){
  # If all input is checked, the dataframe can be made
  # Get the data of one specific treatment
  dfSubset1 <- df[df$protein == protein & df$treatment == treatmentName & df$replID == replID,]
  
  # Check whether yvals exist in the data frame
  if (yvals[1] %in% colnames(dfSubset1)){
    
    # Get the data at one specific dose
    dfSubset2 <- sapply(doses, function(d){
      df2 <- dfSubset1[dfSubset1$dose_uM==d,]
      gplot <- ggplot(df2, aes_string(x="timeAfterExposure", y = yvals[1])) + 
        geom_line() + 
        labs(title=paste0("HepG2 ", protein,"-GFP, ",d, " uM ",altTreatmentName),
             x = "Time (h)", 
             y = yvals[3]) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) 
      if (logScale){
        gplot <- gplot + scale_y_continuous(trans='log')} #, breaks=c(1,2.5,5,10, 20, 40),limits=c(0.5, 50)
      if (withSD == T) {
        gplot <- gplot + geom_errorbar(aes_string(ymin=paste0(yvals[1],"/",yvals[2]), ymax=paste0(yvals[1],"*",yvals[2])), width=1)}
      if (withControl == "DMSO"){
        dfSubset3 <- df[df$protein == protein & df$treatment == "DMSO" & df$replID == replID & df$dose_uM == 0.2,]
        gplot <- gplot + geom_line(data = dfSubset3, linetype = 2, color = "black")}
      else if (withControl == "DMEM"){
        dfSubset3 <- df[df$protein == protein & df$treatment == "DMEM" & df$replID == replID,]
        gplot <- gplot + geom_line(data = dfSubset3, linetype = 2, color = "black")}
      gplot
    },simplify = F)
    
    # Now make one plot including all the dynamics plots
    grid.arrange(grobs = dfSubset2, nrow = length(doses), as.table = FALSE)}
  else {
    print(paste0("Error: column ", yvals[1], " doesn't exist."))}
}

# 32. Function to make dynamic plots per protein and treatment, but combines the replicates
makeCombinedDynPlots <- function(df, treatmentName, which, withSD, combine.by, 
                                 protein, altTreatmentName, yvals, withControl,logScale){
  # If all input is checked, the dataframe can be made
  # Get the data of one specific treatment and protein
  dfSubset1 <- df[df$treatment == treatmentName, ]
  
  # Check whether yvals exist in the data frame
  if (yvals[1] %in% colnames(dfSubset1)){
    
    # If doses are combined, this is done for all three replicates separately
    if (combine.by == "doses"){
      #newcolours <- colorRamps::blue2red(length(levels(droplevels(dfSubset1$dose_uM))))
      gplot <- ggplot(dfSubset1, aes_string(x="timeAfterExposure", y=yvals[1], color = "dose_uM")) + 
        geom_line() + #scale_color_manual(name = "Dose (uM)",values = newcolours) +
        labs(title=paste0("HepG2 cells treated with ", altTreatmentName),
             x = "Time (h)", 
             y = yvals[3]) + 
        scale_color_viridis(name = "Dose (uM)",discrete=T) +
        scale_fill_viridis(discrete=T) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5),
              title =element_text(size=14)) 
      
      # Plot the data on a log-scale
      if (logScale){
        gplot <- gplot + scale_y_continuous(trans='log')} #, breaks=c(1,2.5,5,10, 20, 40),limits=c(0.5, 50)
      # Plot the SD
      if (withSD == T) {
        gplot <- gplot + geom_errorbar(aes_string(ymin=paste0(yvals[1],"/",yvals[2]), ymax=paste0(yvals[1],"*",yvals[2])), width=1)}
      
      # If DMSO or DMEM should be plotted 
      if (withControl == "DMSO"){
        dfSubset2 <- df[df$treatment == "DMSO" & df$dose_uM == 0.2, ]
        gplot <- gplot + geom_line(data = dfSubset2, color = "black", 
                                   linetype = "dashed")}
      else if (withControl == "DMEM"){
        dfSubset2 <- df[df$treatment == "DMEM" & df$dose_uM == 100, ]
        gplot <- gplot + geom_line(data = dfSubset2, color = "black", #size = 0.5, aes(size = "DMEM"), 
                                   linetype = "dashed")}
      
      # Make it a facet plot
      gplot + facet_rep_grid(replID ~ protein)}
    
    # If replicates are combined, this is done for all doses separately
    else if (combine.by == "replicates"){
      #newcolours <- colorRamps::blue2red(length(levels(droplevels(dfSubset1$replID))))
      gplot <- ggplot(dfSubset1, aes_string(x="timeAfterExposure", y=yvals[1], color = "replID")) +
        geom_line() + #scale_colour_manual(name = "Replicate",values = newcolours) +
        scale_color_viridis(name = "Replicate",discrete=T) +
        scale_fill_viridis(discrete=T) +
        labs(title=paste0("HepG2 cells treated with ", altTreatmentName),
             x = "Time (h)",
             y = yvals[3]) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5),
              title =element_text(size=14))
      
      # Plot the data on a log-scale
      if (logScale){
        gplot <- gplot + scale_y_continuous(trans='log')} #, breaks=c(1,2.5,5,10, 20, 40),limits=c(0.5, 50)
      # Plot the SD
      if (withSD) {
        gplot <- gplot + geom_errorbar(aes_string(ymin=paste0(yvals[1],"/",yvals[2]), ymax=paste0(yvals[1],"*",yvals[2])), width=1)}
      
      # Make it a facet plot
      gplot + facet_rep_grid(dose_uM ~ protein)}
    
    else if (combine.by == "proteins"){
      dfSubset2 <- dfSubset1[dfSubset1$protein == protein,]
      #newcolours <- colorRamps::blue2red(length(levels(droplevels(dfSubset1$replID))))
      gplot <- ggplot(dfSubset2, aes_string(x="timeAfterExposure", y=yvals[1], color = "dose_uM")) +
        geom_line() + #scale_colour_manual(name = "Replicate",values = newcolours) +
        scale_color_viridis(name = "Dose (uM)",discrete=T) +
        scale_fill_viridis(discrete=T) +
        labs(title=paste0("HepG2 cells treated with ", altTreatmentName),
             x = "Time (h)",
             y = yvals[3]) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5),
              title =element_text(size=14))
      
      # Plot the data on a log-scale
      if (logScale){
        gplot <- gplot + scale_y_continuous(trans='log')} #, breaks=c(1,2.5,5,10, 20, 40),limits=c(0.5, 50)
      # Plot the SD
      if (withSD) {
        gplot <- gplot + geom_errorbar(aes_string(ymin=paste0(yvals[1],"/",yvals[2]), ymax=paste0(yvals[1],"*",yvals[2])), width=1)}
      
      # # If DMSO or DMEM should be plotted 
      # if (withControl == "DMSO"){
      #   dfSubset3 <- df[df$treatment == "DMSO" & df$dose_uM == 0.2 & df$protein == protein, ]
      #   dfSubset3 <- dfSubset3[dfSubset3$well == dfSubset3$well[1],]
      #   gplot <- gplot + geom_line(data = dfSubset3, color = "black", 
      #                              linetype = "dashed")}
      # else if (withControl == "DMEM"){
      #   dfSubset3 <- df[df$treatment == "DMEM" & df$dose_uM == 100 & df$protein == protein, ]
      #   dfSubset3 <- dfSubset3[dfSubset3$well == dfSubset3$well[1],]
      #   gplot <- gplot + geom_line(data = dfSubset3, color = "black", #size = 0.5, aes(size = "DMEM"), 
      #                              linetype = "dashed")}
      
      # Make it a facet plot
      gplot + facet_rep_grid(replID ~ dose_uM)}}
  else {
    print(paste0("Error: column ", yvals[1], " doesn't exist."))}
}

# 33. Function to plot the dynamics of the integrated intensity at different timepoints and 
#     doses
plotDynOverTime <- function(df = NULL, protein = PROTS_OF_INTEREST, replID = NULL,
                            treatmentName = NULL, doses = NULL, 
                            which = c("nuclei_integrated","nuclei_mean","cytoplasm_integrated",
                                      "cytoplasm_mean","nuclei_integrated_norm","nuclei_mean_norm",
                                      "cytoplasm_integrated_norm",
                                      "cytoplasm_mean_norm", "propPI","propAnV","cellCount_GFP","cellCount_AnVPI",
                                      "relativeCellCount", "nuclei_area"), withSD = T, 
                            combine.by = c("none","replicates","doses","proteins"),
                            withControl = c("none","DMSO","DMEM"),
                            logScale = F){
  # Check whether the input is valid
  withControl <- match.arg(withControl)
  combine.by <- match.arg(combine.by)
  protein <- match.arg(protein)
  which <- match.arg(which)
  checkPlotDynInput(df = df, protein = protein, replID = replID,
                    treatmentName = treatmentName, doses = doses, combine.by = combine.by)
  
  # Get alternative treatment name
  altTreatmentName <- ALT_COMP_NAMES[which(COMPS_OF_INTEREST == treatmentName)]
  
  # Get the y-values
  yvals <- chooseACDVariable(which)
  
  if (combine.by == "none"){
    makeSeparateDynPlots(df = df, protein = protein, replID = replID,treatmentName = treatmentName, yvals = yvals,
                         altTreatmentName = altTreatmentName, doses = doses, which = which, withSD = withSD, withControl,
                         logScale)} 
  else if (combine.by == "replicates" | combine.by == "doses" | combine.by == "proteins"){
    makeCombinedDynPlots(df = df, treatmentName = treatmentName, yvals = yvals, 
                         altTreatmentName = altTreatmentName, which = which, withSD = withSD, 
                         combine.by = combine.by, protein = protein, withControl,logScale)}
}

# 34. 
doBackgroundSubtraction <- function(df, columnsForSubtraction, 
                                    refTreatment = COMPS_OF_INTEREST, 
                                    refConcentration = numeric(),
                                    forTreatment = COMPS_OF_INTEREST, 
                                    sharedGroupVar = character(),
                                    divGroupVar = character()){
  # Are the arguments valid?
  refTreatment <- match.arg(refTreatment)
  ifelse(all(forTreatment %in% COMPS_OF_INTEREST),T,stop("Compound does not exist."))
  
  # Make subset data sets of the specified treatments
  dfRef <- df[df$treatment == refTreatment,]
  dfFor <- df[df$treatment %in% c(refTreatment,forTreatment),]
  dfRest <- df[!df$treatment %in% c(refTreatment,forTreatment),]
  
  # And of the specified concentration
  if(!length(refConcentration) == 0){
    dfRef <- dfRef[dfRef$dose_uM == refConcentration,]}
  
  # Make a new dfRef data frame with the average of the concentrations
  dfRefTmp <- eval(
    parse(
      text=paste0("ddply(dfRef,.(",
                  paste0(sharedGroupVar,collapse = ","),
                  "), summarize,",
                  paste0(columnsForSubtraction,
                         " = mean(", 
                         columnsForSubtraction,
                         ")",
                         collapse = ","),")")))
  
  # Merge (cbind) the dfFor and dfRefTmp data frames...
  dfNew <- merge(dfFor,dfRefTmp, by = sharedGroupVar, all.x = T)
  # ... and order dfNew and dfFor by the same order
  dfNew <- eval(parse(text=paste0("dfNew[with(dfNew,order(treatment,",paste0(c(divGroupVar,sharedGroupVar),collapse = ","),")),]")))
  dfFor <- eval(parse(text=paste0("dfFor[with(dfFor,order(treatment,",paste0(c(divGroupVar,sharedGroupVar),collapse = ","),")),]")))
  
  # ... and do background subtraction of the specified columns
  dfFor[,columnsForSubtraction] <- dfNew[,eval(parse(text=paste0("c(",
                                                                 paste0("\"",
                                                                        columnsForSubtraction,
                                                                        ".x", 
                                                                        collapse = "\","),
                                                                 "\")")))] - 
    dfNew[,eval(parse(text=paste0("c(",
                                  paste0("\"",
                                         columnsForSubtraction,
                                         ".y", 
                                         collapse = "\","),
                                  "\")")))]
  
  dfTotal <- rbind(dfFor,dfRest)
  return(dfTotal)
}

minMaxNormalisation <- function(df, columnsForNormalisation = list(), combineCytAndNuc = T) {
  # Initiate the output data frame
  output <- data.frame()
  if (combineCytAndNuc){
    # Make the new normalized column mutations for every combination
    mutations <- c()
    nwvars <- c()
    for (combi in columnsForNormalisation){
      CytAndNuc <- paste0("c(",paste0(combi,collapse = ","),")")
      
      # Make new variable names
      nwvars <- c(nwvars,paste0(combi,"_MinMaxNorm"))
      # Write the calculation for the mutation
      mutations <- c(mutations,sapply(combi,function(col){
        paste0("(", col, " - min(", CytAndNuc, ",na.rm = T)) / (max(", CytAndNuc, ",na.rm = T) - min(", CytAndNuc, ",na.rm = T))")}))}
    
    # Make a new data frame with the new columns per set
    # Write the mutate expression
    expression <- paste0("mutate(",paste0(paste0(nwvars, " = ", mutations),collapse = ", "),")")
    
    # Mutate some columns
    df <- eval_bare(parse_expr(paste0("df %>% ", expression)))
    
    # And create the output data frame by binding the separate data frames
    output <- bind_rows(output,df)}
  
  else {
    columnsForNormalisation <- unlist(columnsForNormalisation)
    # Make new variable names
    nwvars <- paste0(columnsForNormalisation,"_MinMaxNorm")
    # Write the calculation for the mutation
    mutations <- c(sapply(columnsForNormalisation,function(col){
      paste0("(", col, " - min(", col, ",na.rm = T)) / (max(", col, ",na.rm = T) - min(", col, ",na.rm = T))")}))
    # Write the mutate expression
    expression <- paste0("mutate(",paste0(paste0(nwvars, " = ", mutations),collapse = ", "),")")
    
    # Mutate some columns
    df <- eval_bare(parse_expr(paste0(" df %>% ", expression)))
    # Create the output data frame
    output <- bind_rows(output,df)}
  output
}

minMaxNormPerPlate <- function(dfListSCD, proteins = character(), 
                               replicates = character(), channels = character(),
                               columnsForNormalisation = list(), combineCytAndNuc = T) {
  for (protName in proteins){
    for (repl in replicates) {
      for (channel in channels) {
        if( channel == "GFP" ){
          dfListSCD[[protName]][[repl]][[channel]] <- minMaxNormalisation(df = dfListSCD[[protName]][[repl]][[channel]],
                                                                          columnsForNormalisation = columnsForNormalisation,
                                                                          combineCytAndNuc = combineCytAndNuc)}
      }
    }
  }
  return(dfListSCD)
}

# 35. Min-max normalisation for the compounds and the control treatments
doMinMaxNormAllConditions <- function(df, columnsForNormalisation = list(), by = "", 
                                      normPerCompCtrlSet = T,
                                      treatmentSets = list(), combineCytAndNuc = F, bindTo = data.frame()){
  
  # Initiate the output data frame
  output <- bindTo
  if (combineCytAndNuc){
    
    # Make the new normalized column mutations for every combination
    mutations <- c()
    nwvars <- c()
    for (combi in columnsForNormalisation){
      CytAndNuc <- paste0("c(",paste0(combi,collapse = ","),")")
      
      # Make new variable names
      nwvars <- c(nwvars,paste0(combi,"_Norm"))
      # Write the calculation for the mutation
      mutations <- c(mutations,sapply(combi,function(col){
        paste0("(", col, " - min(", CytAndNuc, ",na.rm = T)) / (max(", CytAndNuc, ",na.rm = T) - min(", CytAndNuc, ",na.rm = T))")}))}
    
    if (normPerCompCtrlSet) {
      # Normalisation per set of compound with its control 
      # Make a new data frame with the new columns per set
      for (set in treatmentSets) {
        nwvars2 <- c(nwvars,"Treatment_Control")
        mutations2 <- c(mutations,paste0("paste0(set[1],","\"_\",","set[2])"))
        
        # Write the mutate expression
        expression <- paste0("mutate(",paste0(paste0(nwvars2, " = ", mutations2),collapse = ", "),")")
        
        # Filter the data frame. The normalisation is then performed for every specified set of a compound and its control.
        df <- df %>% filter(treatment %in% set) 
        # Mutate some columns
        df <- eval_bare(parse_expr(paste0(" df %>% group_by(",paste0(by, collapse = ","),") %>% ", expression)))
        
        # And create the output data frame by binding the separate data frames
        output <- bind_rows(output,df)}}
    else if (!normPerCompCtrlSet) {
      # The normalisation is performed for an entire plate, to be able to compare different compounds with each other. 
      # Write the mutate expression
      expression <- paste0("mutate(",paste0(paste0(nwvars, " = ", mutations),collapse = ", "),")")
      print(expression)
      
      # Mutate some columns
      df <- eval_bare(parse_expr(paste0(" df %>% group_by(",paste0(by, collapse = ","),") %>% ", expression)))
      
      # Create the output data frame
      output <- bind_rows(output,df)}}
  
  else if (!combineCytAndNuc) {
    if (normPerCompCtrlSet) {
      for (set in treatmentSets){
        columnsForNormalisation <- unlist(columnsForNormalisation)
        # Make new variable names
        nwvars <- c(paste0(columnsForNormalisation,"_Norm"),"Treatment_Control")
        # Write the calculation for the mutation
        mutations <- c(sapply(columnsForNormalisation,function(col){
          paste0("(", col, " - min(", col, ",na.rm = T)) / (max(", col, ",na.rm = T) - min(", col, ",na.rm = T))")}),
          paste0("paste0(set[1],","\"_\",","set[2])"))
        # Write the mutate expression
        expression <- paste0("mutate(",paste0(paste0(nwvars, " = ", mutations),collapse = ", "),")")
        
        # Filter the data frame
        df <- df %>% filter(treatment %in% set) 
        # Mutate some columns
        df <- eval_bare(parse_expr(paste0(" df %>% group_by(",paste0(by, collapse = ","),") %>% ", expression)))
        
        # Create the output data frame
        output <- bind_rows(output,df)}} 
    else if (!normPerCompCtrlSet) {
      columnsForNormalisation <- unlist(columnsForNormalisation)
      # Make new variable names
      nwvars <- c(nwvars,paste0(combi,"_Norm"))
      # Write the calculation for the mutation
      mutations <- sapply(columnsForNormalisation,function(col){
        paste0("(", col, " - min(", col, ",na.rm = T)) / (max(", col, ",na.rm = T) - min(", col, ",na.rm = T))")})
      
      # Write the mutate expression
      expression <- paste0("mutate(",paste0(paste0(nwvars, " = ", mutations),collapse = ", "),")")
      
      # Mutate some columns
      df <- eval_bare(parse_expr(paste0(" df %>% group_by(",paste0(by, collapse = ","),") %>% ", expression)))
      # Create the output data frame
      output <- bind_rows(output,df)}}
  output
} 

# 35. Min-max normalisation for the compounds alone, so not for the control treatments
doMinMaxNormCompoundConditions <- function(df,columnsForNormalisation,by, combine = F, controlTreatments = CONTROL_TREATMENTS){
  # Exclude the controls
  dfFor <- df[!df$treatment %in% controlTreatments,]
  dfRest <- df[df$treatment %in% controlTreatments,]
  
  # Do min-max normalisation for the specified columns
  if (!combine) {
    dfFor <- eval(parse(text=paste0("ddply(dfFor,.(",
                                    paste0(by,collapse = ","),
                                    "), mutate,",
                                      paste0(columnsForNormalisation,"_Norm",
                                             " = (", 
                                             columnsForNormalisation,
                                             " - min(",
                                             columnsForNormalisation,
                                             ",na.rm = T)) / (max(",
                                             columnsForNormalisation,
                                             ",na.rm = T) - min(",
                                             columnsForNormalisation,
                                             ",na.rm = T))",
                                             collapse = ","),
                                    ")")))
    return(rbind.fill(dfFor,dfRest))}
  else if (combine) {
    # Do min-max normalisation for the specified columns, but take the overall minimum and maximum
    # of the specified columns
    combi <- paste0("c(",paste0(columnsForNormalisation,collapse = ","),")")
    dfFor <- eval(parse(text=paste0("ddply(dfFor,.(",
                                    paste0(by,collapse = ","),
                                    "), mutate,",
                                      paste0(columnsForNormalisation,"_Norm",
                                             " = (", 
                                             columnsForNormalisation,
                                             " - min(",
                                             combi,
                                             ",na.rm = T)) / (max(",
                                             combi,
                                             ",na.rm = T) - min(",
                                             combi,
                                             ",na.rm = T))",
                                             collapse = ","),
                                    ")")))}
  return(rbind.fill(dfFor,dfRest))}

# 36.
plotCorrelation <- function(df,
                            treatment = COMPS_OF_INTEREST, 
                            protein = PROTS_OF_INTEREST, 
                            replicateNr, fromTime = 0,
                            dose, groupingVar,
                            x, y, 
                            xname, yname,
                            xlab = "x", ylab = "y",
                            title = NULL, save = F,
                            output_path = "~/"){
  # Check input
  treatment <- match.arg(treatment)
  protein <- match.arg(protein)
  
  # Get alternative treatment name
  altTreatmentName <- ALT_COMP_NAMES[which(COMPS_OF_INTEREST == treatment)]
  
  # Do the correlation
  correl <- cor.test(df[df$treatment == treatment &
                          df$protein == protein &
                          df$replID == replicateNr &
                          df$timeAfterExposure > fromTime &
                          df$dose_uM == dose, 
                        y],
                     df[df$treatment == treatment &
                          df$protein == protein &
                          df$replID == replicateNr &
                          df$timeAfterExposure > fromTime &
                          df$dose_uM == dose, 
                        x])
  df <- df[df$treatment == treatment &
             df$protein == protein &
             df$replID == replicateNr &
             df$timeAfterExposure > fromTime &
             df$dose_uM == dose,]
  
  # Make a title if there is none
  if (is.null(title)){
    title <- paste0(protein," replicate ", replicateNr,", ", 
                    dose," uM ", altTreatmentName)
    if (fromTime > 0){
      title <- paste(title," for t > ",fromTime," hr")}}
  
  # Make the plot
  gplot <- ggplot(df) + 
    geom_point(aes_string(x = x, 
                          y = y, 
                          group = groupingVar, color = groupingVar)) + 
    geom_smooth(aes_string(x = x, 
                           y = y), method = "lm", se = F, color = "black") +
    scale_color_gradient(name = "Time (hr)",low="blue", high="red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=16),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 18)) +
    labs(x = xlab,
         y = ylab,
         title = title) +
    annotate("text", x = min(df[,x], na.rm = T) + (max(df[,x], na.rm = T) - min(df[,x], na.rm = T)) * 0.4, 
             y = min(df[,y], na.rm = T) + (max(df[,y], na.rm = T) - min(df[,y], na.rm = T)) * 0.95, 
             label= paste0("rho = ",round(correl$estimate,2)),
             size = 6) + 
    annotate("text", x = min(df[,x], na.rm = T) + (max(df[,x], na.rm = T) - min(df[,x], na.rm = T)) * 0.4, 
             y = min(df[,y], na.rm = T) + (max(df[,y], na.rm = T) - min(df[,y], na.rm = T)) * 0.9, 
             label = paste0("p-value = ",signif(correl$p.value,digits = 2)),
             size = 6)
  
  
  if (save){
    gplot
    ggsave(filename = paste0(output_path,protein,"_rep",replicateNr,"_",dose,"uM_",treatment,
                             "_",yname,"_vs_",xname,"_from_",fromTime,"h",".pdf"),
           width = 7, height = 5)}
  else {
    gplot}
}

# 37.
rankFunction <- function(v){
  df <- data.frame(value = sort(unique(v)), rank = seq(1,length(unique(v)),1))
  
  rankVector <- c()
  for (value in v){
    rankVector <- c(rankVector,df[df$value==value,"rank"])}
  return(rankVector)
}

# 38
prepareD4M <- function(df1, absolute = T){
  # Function to prepare data for modeling
  # Prepare the data to be used in the model fitting script
  # Make a data frame that will contain the final data
  data4Model <- df1 %>% select(protein,replID,treatment,timeID,timeAfterExposure, dose_uM, dose_uMadj, timepoints, data4modelInterpol,data4modelReal)
  # Make a data frame with the cell counts that will eventually be appended to the CDDP4Model data
  if (absolute){
    df2 <- df1 %>% select(protein,replID,treatment,timeID,timeAfterExposure, dose_uM, dose_uMadj, timepoints, interpol_absoluteCellCount, imageCountParentObj)
    # Make a data frame that contains new replicate IDs
    replNamesDF <- df1 %>% group_by(protein, replID) %>% summarise(replName = replID[1])
    replNamesDF$newReplID <- seq(1,nrow(replNamesDF),1)
    replNamesDF$replName <- NULL
    # Merge the two data frame and rename the columns
    df2 <- left_join(df2,replNamesDF, by = c("protein","replID"))
    # Rename the replID column
    df2$replID <- df2$newReplID
    df2$newReplID <- NULL
    df2$replID <- factor(as.numeric(df2$replID))
    # Rename the cell count data columns
    df2 <- df2 %>% dplyr::rename(data4modelInterpol = interpol_absoluteCellCount) %>% dplyr::rename(data4modelReal = imageCountParentObj)
  }
  else {
    df2 <- df1 %>% select(protein,replID,treatment,timeID,timeAfterExposure, dose_uM, dose_uMadj, timepoints, interpol_relativeCellCount, relativeCellCount)
    # Make a data frame that contains new replicate IDs
    replNamesDF <- df1 %>% group_by(protein, replID) %>% summarise(replName = replID[1])
    replNamesDF$newReplID <- seq(1,nrow(replNamesDF),1)
    replNamesDF$replName <- NULL
    # Merge the two data frame and rename the columns
    df2 <- left_join(df2,replNamesDF, by = c("protein","replID"))
    # Rename the replID column
    df2$replID <- df2$newReplID
    df2$newReplID <- NULL
    df2$replID <- factor(as.numeric(df2$replID))
    # Rename the cell count data columns
    df2 <- df2 %>% dplyr::rename(data4modelInterpol = interpol_relativeCellCount) %>% dplyr::rename(data4modelReal = relativeCellCount)
  }
  # Get rid of the protein names
  df2$protein <- factor("N")
  # Merge the two data frames
  data4Model <- bind_rows(data4Model,df2)
  # Convert to factors
  data4Model$protein <- factor(data4Model$protein)
  data4Model$replID <- factor(data4Model$replID)
  return(data4Model)
}

# Function to make replace column names
replaceColumns <- function(dfListSCD, channel = CHANNELS, colList){
  # Check whether input channel is correct
  channel <- match.arg(channel)
  
  dfListSCD <- lapply(dfListSCD,function(prots){
    lapply(prots,function(reps){
      sapply(names(reps), function(chanName){
        if (chanName == channel){
          for (colPair in colList){
            names(reps[[chanName]])[names(reps[[chanName]]) == colPair[1]] <- colPair[2]}
          reps[[chanName]]}
        else {
          reps[[chanName]]}
      },simplify = F, USE.NAMES = T)
    })
  })
}

# Function to make replace column names
checkReplID <- function(dfListSCD, channel = CHANNELS){
  # Check whether input channel is correct
  channel <- match.arg(channel)
  
  dfListSCD <- sapply(names(dfListSCD),function(protName){
    sapply(names(dfListSCD[[protName]]),function(repName){
      sapply(names(dfListSCD[[protName]][[repName]]), function(chanName){
        if (chanName == channel){
          df <- dfListSCD[[protName]][[repName]][[chanName]]
          if (!unique(df$replID) == as.numeric(substr(repName,4,4))){
            df$replID <- as.character(substr(repName,4,4))}
          df}
        else {
          dfListSCD[[protName]][[repName]][[chanName]]}
      },simplify = F, USE.NAMES = T)
    },simplify = F, USE.NAMES = T)
  },simplify = F, USE.NAMES = T)
}

