##########################################################################################################################
###################### MixInfect2 - Identify mixed samples from VCF file #################################################
##########################################################################################################################

# Setup CRAN mirror once
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Helper to load or install packages quietly
safe_library <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Load required packages
safe_library("mclust")
safe_library("stringr")
safe_library("optparse")
safe_library("foreach")
safe_library("doParallel")

options(stringsAsFactors = FALSE)

MixInfect2 <- function(VCFfile, prefix = "output", maskFile = NULL, useFilter = TRUE, 
                       minQual = 20, LowCov = 10, minDepth = 5, 
                       popFreq_threshold = 1, SNPwindow = 100, n_threads = 4) {

  cl <- makeCluster(n_threads)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  vcf <- read.table(VCFfile, stringsAsFactors = FALSE)
  cat("Total variants read from VCF:", nrow(vcf), "\n")
  header_input <- as.matrix(read.table(VCFfile, comment.char = " ", sep = "\n", stringsAsFactors = FALSE))
  end_head <- which(grepl("#CHROM", header_input) == TRUE)
  header <- as.data.frame(header_input[1:(end_head - 1)])
  names <- unlist(strsplit(header_input[end_head], "\t"))
  colnames(vcf) <- names
  format <- which(names == "FORMAT")
  head_start <- names[1:format]
  names <- names[(format + 1):length(names)]
  rm(header_input)

  ind <- lapply(1:nrow(vcf), function(i) {
    length(unlist(strsplit(vcf[i, 4], ""))) > 1 || length(unlist(strsplit(unlist(strsplit(vcf[i, 5], ","))[1], ""))) > 1
  })
  vcf <- vcf[which(ind == FALSE), ]
  cat("After removing indels:", nrow(vcf), "\n")

  if (useFilter) {
    vcf <- vcf[vcf[, which(head_start == "FILTER")] == "PASS", ]
    cat("After FILTER == 'PASS':", nrow(vcf), "\n")
  }

  vcf <- vcf[vcf[, which(head_start == "QUAL")] > minQual, ]
  cat("After QUAL >", minQual, ":", nrow(vcf), "\n")

  vcf <- vcf[grep("*", vcf[, which(head_start == "ALT")], fixed = TRUE, invert = TRUE), ]
  cat("After removing spanning variants (*):", nrow(vcf), "\n")

  if (!is.null(maskFile)) {
    maskFile <- read.csv(maskFile)
    res1 <- unlist(apply(maskFile, 1, function(row) seq(row[1], row[2])))
    vcf <- vcf[!vcf[, 2] %in% res1, ]
    cat("After masking regions:", nrow(vcf), "\n")
  }

  if (nrow(vcf) == 0) {
    stop("VCF is empty after filtering. Adjust filters or check input.")
  }

  AD <- which(unlist(stringr::str_split(vcf[1, format], ":")) == 'AD')
  GT <- which(unlist(stringr::str_split(vcf[1, format], ":")) == 'GT')
  DP <- which(unlist(stringr::str_split(vcf[1, format], ":")) == 'DP')

  GT_mat <- matrix(ncol = length((format + 1):ncol(vcf)), nrow = nrow(vcf))
  AD_mat <- matrix(ncol = length((format + 1):ncol(vcf)), nrow = nrow(vcf))
  DP_mat <- matrix(ncol = length((format + 1):ncol(vcf)), nrow = nrow(vcf))

  for (i in 1:ncol(GT_mat)) {
    parsed <- stringr::str_split(vcf[, i + format], ":")

    GT_mat[, i] <- sapply(parsed, function(x) {
      if (is.null(x) || length(x) < GT) NA else x[[GT]]
    })

    DP_mat[, i] <- sapply(parsed, function(x) {
      if (is.null(x) || length(x) < DP) "0" else x[[DP]]
    })

    AD_mat[, i] <- sapply(parsed, function(x) {
      if (is.null(x) || length(x) < AD) NA else x[[AD]]
    })
  }

  DP_mat <- apply(DP_mat, 2, function(x) as.numeric(ifelse(is.na(x) | x == ".", 0, x)))
  GT_mat[DP_mat < LowCov] <- "?"

  outfile <- data.frame(
    SampleName = names,
    Mix.Non.mix = rep("Non-mix", length(names)),
    hSNPs = NA_integer_,
    Total.SNPs = NA_integer_,
    Proportion.hSNPs_totalSNPs = NA_real_,
    No.strains = 1,
    Major.strain.proportion = NA_real_,
    stringsAsFactors = FALSE
  )

  mixed_calls <- c("0/1", "0/2", "0/3", "1/2", "1/3", "2/3", "0|1", "0|2", "0|3", "1|2", "1|3", "2|3")
  alt_calls <- c("1/1", "2/2", "3/3", "1|1", "2|2", "3|3")

  for (col in 1:ncol(AD_mat)) {
    ADmix <- stringr::str_split(AD_mat[, col], ",")
    for (m in seq_along(ADmix)) {
      AD_site <- as.numeric(unlist(ADmix[m]))
      AD_site <- AD_site[AD_site != 0]
      if (length(AD_site) > 1) {
        AD_site <- sort(AD_site, decreasing = TRUE)
        if (AD_site[2] >= minDepth) {
          GT_mat[m, col] <- "0/1"
        }
      }
    }
  }

  keep <- apply(GT_mat, 1, function(row) any(row %in% c(mixed_calls, alt_calls)))
  GT_mat <- GT_mat[keep, , drop = FALSE]
  vcf <- vcf[keep, , drop = FALSE]

  propMix <- apply(GT_mat, 1, function(row) sum(row %in% mixed_calls) / length(row))
  propMix <- which(propMix > popFreq_threshold)
  if (length(propMix) > 0) {
    GT_mat <- GT_mat[-propMix, ]
    vcf <- vcf[-propMix, ]
  }

  mixes <- matrix(0, ncol = ncol(GT_mat), nrow = 4)
  for (i in 1:ncol(GT_mat)) {
    mixes[1, i] <- sum(GT_mat[, i] %in% mixed_calls)
    mixes[2, i] <- sum(GT_mat[, i] %in% alt_calls)
    mixes[3, i] <- mixes[1, i] + mixes[2, i]
    mixes[4, i] <- (mixes[1, i] / mixes[3, i]) * 100
  }

  outfile$hSNPs <- mixes[1, ]
  outfile$Total.SNPs <- mixes[3, ]
  outfile$Proportion.hSNPs_totalSNPs <- round(mixes[4, ], 2)

  mixnames <- outfile$SampleName[which(outfile[, 5] > 1.5 & outfile[, 3] > 10)]
  if (length(mixnames) > 0) {
    mix_GT <- GT_mat[, outfile[, 5] > 1.5 & outfile[, 3] > 10, drop = FALSE]
    mix_VCF <- vcf[, which(outfile[, 5] > 1.5 & outfile[, 3] > 10) + format, drop = FALSE]
    positions <- vcf[, 2]
    BICvalues <- data.frame(Sample = mixnames, G2 = NA_real_, G4 = NA_real_, G6 = NA_real_)

    results <- foreach(i = 1:nrow(BICvalues),
                       .packages = c("mclust", "stringr")) %dopar% {
      samplemix_sites <- mix_VCF[which(mix_GT[, i] %in% mixed_calls), i]
      samplemix_AD <- sapply(stringr::str_split(samplemix_sites, ":"), function(x) {
        if (is.null(x) || length(x) == 0 || all(is.na(x))) {
          NA
        } else if (length(x) >= AD) {
          x[[AD]]
        } else {
          NA
        }
      })
      samplePos <- positions[which(mix_GT[, i] %in% mixed_calls)]
      samplemaj_prop <- samplemin_prop <- finalPos <- numeric()

      for (k in seq_along(samplemix_AD)) {
        sampleAD <- as.numeric(unlist(stringr::str_split(samplemix_AD[k], ",")))
        sampleAD <- sampleAD[sampleAD != 0]
        sampleAD <- sort(sampleAD, decreasing = TRUE)
        if (length(sampleAD) == 2) {
          samplemaj_prop <- c(samplemaj_prop, sampleAD[1] / sum(sampleAD))
          samplemin_prop <- c(samplemin_prop, sampleAD[2] / sum(sampleAD))
          finalPos <- c(finalPos, samplePos[k])
        }
      }

      distances <- diff(finalPos)
      group_indices <- c(1, cumsum(distances >= SNPwindow) + 1)
      samplemaj_prop <- tapply(samplemaj_prop, group_indices, median)
      samplemin_prop <- tapply(samplemin_prop, group_indices, median)
      b <- c(samplemaj_prop, samplemin_prop)

      bic_values <- c(NA_real_, NA_real_, NA_real_)
      mix_status <- 'Non-mix'
      no_strains <- 1
      major_strain_proportion <- NA_real_

      if (length(b) > LowCov) {
        a <- mclustBIC(b, G = c(2, 4, 6), verbose = FALSE)[, 2]
        if (length(a) == 3) bic_values <- a

        d <- Mclust(b, G = 2, verbose = FALSE)
        if (!is.na(bic_values[1]) && bic_values[1] >= 20) {
          mix_status <- 'Mix'
          no_strains <- 2
          major_strain_proportion <- max(d$parameters$mean)
        } else if (!is.na(bic_values[3]) && bic_values[3] >= 20) {
          d <- Mclust(b, G = 6, verbose = FALSE)
          mix_status <- 'Mix'
          no_strains <- 3
          means <- sort(d$parameters$mean, decreasing = TRUE)
          major_strain_proportion <- if (sum(means[5:6]) < 0.5) means[3] else means[4]
        }
      }

      return(c(i, bic_values, mix_status, no_strains, major_strain_proportion))
    }

    for (res in results) {
      i <- as.numeric(res[1])
      BICvalues[i, 2:4] <- as.numeric(res[2:4])
      ind <- which(BICvalues$Sample[i] == outfile$SampleName)
      outfile[ind, 2] <- res[5]
      outfile[ind, 6] <- as.numeric(res[6])
      outfile[ind, 7] <- round(as.numeric(res[7]), 2)
    }

    write.csv(BICvalues, paste0(prefix, "_BICvalues.csv"), row.names = FALSE)
    write.csv(outfile, paste0(prefix, "_MixSampleSummary.csv"), row.names = FALSE)
  } else {
    message("No mixed infection detected.")
  }

  return(outfile)
}

option_list <- list(
  make_option(c("--VCFfile"), type = "character", default = NULL),
  make_option(c("--prefix"), type = "character", default = "output"),
  make_option(c("--maskFile"), type = "character", default = NULL),
  make_option(c("--useFilter"), type = "logical", default = TRUE),
  make_option(c("--minQual"), type = "numeric", default = 20),
  make_option(c("--LowCov"), type = "numeric", default = 10),
  make_option(c("--popFreq_threshold"), type = "numeric", default = 1),
  make_option(c("--minDepth"), type = "integer", default = 5),
  make_option(c("--SNPwindow"), type = "numeric", default = 100),
  make_option(c("--n_threads"), type = "numeric", default = 4)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$VCFfile)) {
  print_help(opt_parser)
  stop("The VCF file must be specified.", call. = FALSE)
}

MixInfect2(opt$VCFfile, opt$prefix, opt$maskFile, opt$useFilter, 
           opt$minQual, opt$LowCov, opt$minDepth, 
           opt$popFreq_threshold, opt$SNPwindow, opt$n_threads)
