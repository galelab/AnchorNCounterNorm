#' Load nCounter files
#'
#' This code contain qc functions from Michael Love's NCounter analysis
#' 
imagingQC <- function(rcc) {


    #### INPUT: rcc - input from rcc
    #### OUTPUT: flag for imaging quality

    fovRatio <- as.numeric(rcc$Lane_Attributes[3]) / as.numeric(rcc$Lane_Attributes[2])
    if (!(fovRatio > .75)) {
        return("Flag")
    }
    if (fovRatio > .75) {
        return("No flag")
    }
}

bindingDensityQC <- function(rcc, low, high) {


    #### INPUT: rcc - input from rcc
    ####         low, high - the lower and upper limits for binding density
    #### OUTPUT: flag for binding density

    bd <- as.numeric(rcc$Lane_Attributes[6])
    if (!(bd < high & bd > low)) {
        return("Flag")
    }
    if (bd < high & bd > low) {
        return("No flag")
    }
}

limitOfDetectionQC <- function(rcc, numSD = 0) {

    #### INPUT: rcc - input from rcc
    ####         numSD - number of standard deviations to calibrate the LOD
    #### OUTPUT: flag for limit of detection

    counts <- rcc$Code_Summary
    posE <- as.numeric(counts$Count[counts$Name == "POS_E"])
    negControls <- as.numeric(counts$Count[grepl("NEG", counts$Name)])
    if (!(posE > mean(negControls) + numSD * sd(negControls))) {
        return("Flag")
    }
    if (posE > mean(negControls) + numSD * sd(negControls)) {
        return("No flag")
    }
}

positiveLinQC <- function(rcc) {

    #### INPUT: rcc - input from rcc
    #### OUTPUT: flag for linearity for positive controls


    counts <- rcc$Code_Summary
    posControls <- as.numeric(counts$Count[grepl("POS_", counts$Name)])
    known <- c(128, 128 / 4, 128 / 16, 128 / 64, 128 / 256, 128 / (256 * 4))
    r2 <- summary(lm(sort(posControls) ~ sort(known)))$r.squared
    if (!(r2 > .95) | is.na(r2)) {
        return("Flag")
    }
    if (r2 > .95) {
        return("No flag")
    }
}