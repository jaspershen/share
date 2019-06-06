setGeneric(name = "SXTMTmatch",
           def = function(data1,
                          data2,
                          mz.tol,
                          #rt.tol is relative
                          rt.tol = 30,
                          rt.error.type = c("relative", "abs")){
             rt.error.type <- match.arg(rt.error.type)
             #
             if (nrow(data1) == 0 | nrow(data2) == 0) {
               result <- NULL
               return(result)
             }
             # mz1 <- as.numeric(data1[, 1])
             # rt1 <- as.numeric(data1[, 2])
             info1 <- data1[,c(1,2),drop = FALSE]
             info1 <- apply(info1, 1, list)
             
             mz2 <- as.numeric(data2[, 1])
             rt2 <- as.numeric(data2[, 2])
             
             result <- pbapply::pblapply(info1, function(x) {
               temp.mz1 <- x[[1]][[1]]
               temp.rt1 <- x[[1]][[2]]
               mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
               if(rt.error.type == "relative"){
                 rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
               }else{
                 rt.error <- abs(temp.rt1 - rt2)
               }
               
               j <- which(mz.error <= mz.tol & rt.error <= rt.tol)
               if(length(j) == 0){
                 matrix(NA, ncol = 7)
               }else{
                 cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j])
               }
             })
             
             if(length(result) == 1){
               result <- cbind(1,result[[1]])
             }else{
               result <- mapply(function(x,y){list(cbind(x,y))},
                                x <- 1:length(info1),
                                y = result)
               result <- do.call(rbind, result)
             }
             
             result <- matrix(result[which(!apply(result,1,function(x) any(is.na(x)))),], ncol = 8)
             if(nrow(result) == 0) return(NULL)
             colnames(result) <-
               c("Index1",
                 "Index2",
                 "mz1",
                 "mz2",
                 "mz error",
                 "rt1",
                 "rt2",
                 "rt error")
             result <- result
           })








#formula functions
#############------------------------------------------------------------------
#' @title sumFormula
#' @description Combine metabolite and adduct as a new sum formula.
#' If there are no enough element to remove, return NA.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param adduct The adduct of metabolite.
#' @return  A sum formula.
#' @export

setGeneric(name = "sumFormula",
           def = function(formula = "C9H11NO2",
                          adduct = "M-H2O+H"){
             
             if(is.na(formula)) return(NA)
             if(is.na(adduct)) return(formula)
             if(adduct == "M+" | adduct == "M-"){
               return(formula)
             }
             
             formula1 <- splitFormula(formula)
             adduct1 <- strsplit(x = adduct, split = "\\-|\\+")[[1]][-1]
             polymer <- as.numeric(gsub(pattern = "M", replacement = "",
                                        strsplit(x = adduct, split = "\\-|\\+")[[1]][1]))
             if (is.na(polymer)) polymer <- 1
             
             plusorminus <- strsplit(x = adduct, split = "")[[1]]
             plusorminus <- grep("\\+|\\-", plusorminus, value = TRUE)
             
             formula1$number <- formula1$number * polymer
             
             adduct1 <- mapply(function(x, y){
               temp <- splitFormula(x)
               temp$number <- temp$number * ifelse(y == "+", 1, -1)
               list(temp)
             },
             x = adduct1,
             y = plusorminus)
             
             adduct1 <- do.call(rbind, adduct1)
             
             formula <- rbind(formula1, adduct1)
             rownames(formula) <- NULL
             
             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               if(any(formula$number < 0)) {
                 return(NA)
               }else{
                 formula$number[formula$number==1] <- "W"
                 formula <- paste(paste(formula$element.name, formula$number, sep = ""), collapse = "")
                 formula <- strsplit(formula, split = "")[[1]]
                 formula[formula == "W"] <- ""
                 formula <- paste(formula, collapse = "")
                 return(formula)
               }
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })
               
               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })
               
               formula <- do.call(rbind, formula)
               formula <- formula[formula[,2] != 0,]
               colnames(formula) <- c("element.name", "number")
               if(any(formula$number < 0)) {return(NA)}else{
                 formula$number[formula$number==1] <- "W"
                 formula <- paste(paste(formula$element.name, formula$number, sep = ""), collapse = "")
                 formula <- strsplit(formula, split = "")[[1]]
                 formula[formula == "W"] <- ""
                 formula <- paste(formula, collapse = "")
                 return(formula)
               }
             }
           })






#-----------------------------------------------------------------------------
#' @title splitFormula
#' @description Split a formula into element and number.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @return  A splited formula.
#' @export

setGeneric(name = "splitFormula",
           def = function(formula = "C9H11NO2"){
             temp.formula <- strsplit(formula, split = "")[[1]]
             
             number <- NULL
             for(i in 1:length(temp.formula)){
               if(length(grep("[0-9]{1}", temp.formula[i])) == 0){break}
               number[i] <- temp.formula[i]
             }
             
             if(!is.null(number)) {
               number <- as.numeric(paste(number, collapse = ""))
             }else{
               number <- 1
             }
             ##first select the Na, Cl and so on element
             idx1 <- gregexpr("[A-Z][a-z][0-9]*", formula)[[1]]
             len1 <- attributes(idx1)$match.length
             ##no double element
             if(idx1[1] == -1) {
               double.formula <- matrix(NA, ncol = 2)
               formula1 <- formula
             }else{
               double.letter.element <- NULL
               double.number <- NULL
               remove.idx <- NULL
               for (i in 1:length(idx1)) {
                 double.letter.element[i] <- substr(formula, idx1[i], idx1[i] + len1[i] - 1)
                 if(nchar(double.letter.element[i]) == 2){
                   double.number[i] <- 1
                 }else{
                   double.number[i] <- as.numeric(substr(double.letter.element[i], 3, nchar(double.letter.element[i])))
                 }
                 double.letter.element[i] <- substr(double.letter.element[i], 1, 2)
                 remove.idx <- c(remove.idx, idx1[i] : (idx1[i] + len1[i] - 1))
               }
               
               double.formula <- data.frame(double.letter.element,
                                            double.number, stringsAsFactors = FALSE)
               formula1 <- strsplit(formula, split = "")[[1]]
               formula1 <- formula1[-remove.idx]
               formula1 <- paste(formula1, collapse = "")
             }
             
             ## no one element
             if(formula1 == ""){
               one.formula <- matrix(NA, ncol = 2)
             }else{
               idx2 <- gregexpr("[A-Z][0-9]*", formula1)[[1]]
               len2 <- attributes(idx2)$match.length
               one.letter.element <- NULL
               one.number <- NULL
               for (i in 1:length(idx2)) {
                 one.letter.element[i] <- substr(formula1, idx2[i], idx2[i] + len2[i] - 1)
                 if(nchar(one.letter.element[i]) == 1){
                   one.number[i] <- 1
                 }else{
                   one.number[i] <- as.numeric(substr(one.letter.element[i], 2, nchar(one.letter.element[i])))
                 }
                 one.letter.element[i] <- substr(one.letter.element[i], 1, 1)
               }
               one.formula <- data.frame(one.letter.element, one.number,
                                         stringsAsFactors = FALSE)
             }
             
             colnames(double.formula) <- colnames(one.formula) <- c("element.name","number")
             formula <- rbind(double.formula, one.formula)
             formula <- formula[!apply(formula, 1, function(x) any(is.na(x))),]
             
             formula <- formula[order(formula$element.name),]
             formula$number <- formula$number * number
             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               return(formula)
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })
               
               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })
               
               formula <- do.call(rbind, formula)
               colnames(formula) <- c("element.name", "number")
               return(formula)
             }
           })

#-----------------------------------------------------------------------------
#' @title pasteElement
#' @description Paste formula and element.
#' Combine metabolite and adduct as a new sum formula.
#' If there are no enough element to remove, return NA.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param element The element.
#' @param mode Add or remove a module
#' @export
#' @return  A formula.


pasteElement <- function(formula = "C9H11NO2",
                         element = "H",
                         mode = c("plus", "minus")){
  
  mode <- match.arg(mode)
  formula <- splitFormula(formula = formula)
  element <- splitFormula(formula = element)
  
  
  ## mode = plus
  if(mode == "plus"){
    for (i in 1:nrow(element)){
      temp.name <- as.character(element[i,1])
      temp.number <- as.numeric(element[i,2])
      temp.idx <- match(temp.name, formula[,1])
      if(is.na(temp.idx)) {
        formula <- rbind(formula, element[i,])
      }else{
        formula[temp.idx, 2] <- formula[temp.idx, 2] + temp.number
      }
    }
  }else{
    for (i in 1:nrow(element)){
      temp.name <- as.character(element[i,1])
      temp.number <- as.numeric(element[i,2])
      temp.idx <- match(temp.name, formula[,1])
      if(is.na(temp.idx)) {
        # warning("Formula has no element in adduct!\n")
        return(NA)
      }else{
        formula[temp.idx,2] <- formula[temp.idx,2] - temp.number
        if(formula[temp.idx,2] < 0) {
          # warning("Formula has no enough element in adduct!\n")
          return(NA)
        }
      }
    }
  }
  
  ###return formula
  formula <- as.data.frame(formula)
  formula <- formula[formula[,2] != 0, , drop = FALSE]
  formula <- c(t(formula))
  formula <- gsub(pattern = " ", replacement = "", x = formula)
  formula <- formula[formula != "1"]
  formula <- paste(formula, collapse = "")
  return(formula)
}




#-----------------------------------------------------------------------------
#' @title checkElement
#' @description Check a formula can add one adduct or not.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param adduct The adduct.
#' @return  valid, return TRUE; invalid.
#' @export

setGeneric(name = "checkElement",
           def = function(formula = "C9H11NO2",
                          adduct = "M-H2O+H"){
             formula1 <- splitFormula(formula)
             adduct1 <- strsplit(x = adduct, split = "\\-|\\+")[[1]][-1]
             plusorminus <- strsplit(x = adduct, split = "")[[1]]
             plusorminus <- grep("\\+|\\-", plusorminus, value = TRUE)
             if(all(plusorminus == "+")) return(TRUE)
             
             adduct1 <- mapply(function(x, y){
               temp <- splitFormula(x)
               temp$number <- temp$number * ifelse(y == "+", 1, -1)
               list(temp)
             },
             x = adduct1,
             y = plusorminus)
             
             adduct1 <- do.call(rbind, adduct1)
             
             formula <- rbind(formula1, adduct1)
             rownames(formula) <- NULL
             
             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               if(any(formula$number < 0)) {return(FALSE)}else{return(TRUE)}
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })
               
               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })
               
               formula <- do.call(rbind, formula)
               colnames(formula) <- c("element.name", "number")
               if(any(formula$number < 0)) {return(FALSE)}else{return(TRUE)}
             }
           })


setGeneric(name = "getExactMass", 
           def = function(mz = 315.2314,
                          polarity = c("positive", "negative"),
                          colum = c("hilic", "rp"),
                          adduct = c("M+", "(2M+K)+")){
             ##load adduct table
             if(polarity == "positive" & column == "hilic"){
               data("hilic.pos", envir = environment())
               adduct.table <- hilic.pos
             }
             
             if(polarity == "positive" & column == "rp"){
               data("rp.pos", envir = environment())
               adduct.table <- rp.pos
             }
             
             if(polarity == "negative" & column == "hilic"){
               data("hilic.neg", envir = environment())
               adduct.table <- hilic.neg
             }
             
             if(polarity == "negative" & column == "rp"){
               data("rp.neg", envir = environment())
               adduct.table <- rp.neg
             }
             
             mz <- as.numeric(mz)
             remove.mz <- adduct.table$mz[match(adduct, adduct.table$adduct)]
             m.number <- stringr::str_extract(string = adduct, pattern = "[0-9]{0,1}M")
             m.number <- stringr::str_replace_all(string = m.number, pattern = "M", replacement = "")
             m.number[m.number == ""] <- 1
             m.number <- as.numeric(m.number)
             
             exact.mass <- (mz - remove.mz)/m.number
             exact.mass
             
           }
)




#'@title readMZXML
#'@description Read mzXML data.
#'@author Xiaotao Shen
#'\email{shenxt1990@@163.com}
#'@param file The vector of names of ms2 files. MS2 file must be mzXML or mzML.
#'@return Return ms2 data. This is a list, 
#'@export

setGeneric(name = "readMZXML",
           def = function(file,
                          threads = 3){
             # pbapply::pboptions(style = 1)
             cat("Reading MS2 data...\n")
             # mzxml.data.list <- pbapply::pblapply(file, ListMGF)
             ms2 <- MSnbase::readMSData(files = file, msLevel. = 2, mode = "onDisk")
             cat("Processing...\n")
             
             new.ms2 <- ProtGenerics::spectra(object = ms2)
             rm(list = c("ms2"))
             temp.fun <- function(idx, ms2){
               temp.ms2 <- ms2[[idx]]
               rm(list = c("ms2"))
               info <- data.frame(name = paste("mz", temp.ms2@precursorMz,
                                               "rt", temp.ms2@rt, sep = ""),
                                  "mz" = temp.ms2@precursorMz, 
                                  "rt" = temp.ms2@rt, 
                                  "file" = file[temp.ms2@fromFile],
                                  stringsAsFactors = FALSE)
               duplicated.name <- unique(info$name[duplicated(info$name)])
               if(length(duplicated.name) > 0){
                 lapply(duplicated.name, function(x){
                   info$name[which(info$name == x)] <- paste(x, c(1:sum(info$name == x)), sep = "_")
                 })
               }
               
               rownames(info) <- NULL
               spec <- data.frame("mz" = temp.ms2@mz,
                                  "intensity" = temp.ms2@intensity, 
                                  stringsAsFactors = FALSE)
               list(info = info, spec = spec)
             }
             
             new.ms2 <- BiocParallel::bplapply(X = c(1:length(new.ms2)), 
                                               FUN = temp.fun,
                                               ms2 = new.ms2,
                                               BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                                 progressbar = TRUE))
             new.ms2 <- new.ms2
           })






