

## Rename thresholds
renameEpitopes <- function(thresholds){

  thresholds$model <- plyr::revalue(thresholds$model, c(PKYVKQNTLKLAT_cdr3b = "InfA_HA_PKY",
                                                        GILGFVFTL_cdr3b     = "InfA_M1_GIL",
                                                        NLVPMVATV_cdr3b     = "CMV_p65_NLV",
                                                        TPRVTGGGAM_cdr3b    = "CMV_p65_TPR",
                                                        RAKFKQLL_cdr3b      = "EBV_BZLF1_RAF",
                                                        IPSINVHHY_cdr3b     = "CMV_p65_IPS",
                                                        YVLDHLIVV_cdr3b     = "EBV_BRLF1_YVL",
                                                        RPRGEVRFL_cdr3b     = "HSV2_VP22_RPR",
                                                        GLCTLVAML_cdr3b     = "EBV_BMLF1_GLC"
                                                      ))
  return(thresholds)

}


## Load thresholds used for TCRGP predictions
thresholds_ver1 <- read.delim("results/tcrgp/threshold_table.txt", stringsAsFactors = F)
thresholds_ver2 <- read.delim("results/tcrgp/threshold_table.csv", stringsAsFactors = F, sep = ',')
thresholds_ver2 <- thresholds_ver2[,-6]

colnames(thresholds_ver1) <- colnames(thresholds_ver2)
thresholds <- rbind(thresholds_ver1, thresholds_ver2)
thresholds <- thresholds %>% renameEpitopes()
