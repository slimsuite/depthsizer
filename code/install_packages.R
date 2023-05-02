##################################################################################
### INSTALL_PACKAGES: Installs R packages needed for SLiMSuite R scripts ~~~~~ ###
### VERSION: 0.1.0                                                       ~~~~~ ###
### LAST EDIT: 02/05/23                                                  ~~~~~ ###
### AUTHORS: Richard Edwards 2023                                        ~~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au / GitHub @slimsuite                 ~~~~~ ###
##################################################################################

# This script is for installing the R packages required by SLiMSuite R scripts.

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial version for ChromSyn and DepthSizer.
version = "v0.1.0"

####################################### ::: FUNCTIONS ::: ##########################################
### ~ logWrite function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite <- function(logstr){
  writeLines(paste0("[",date(),"] ",logstr),con=stdout())
}
logWrite(paste("#RCODE install_packages.R:",version))
logWrite(paste("#PATH Running from:",getwd()))

####################################### ::: INSTALL ::: ############################################
#i# Check and install packages
for(rpack in c("tidyverse","RColorBrewer","gtools","ggstatsplot","writexl","parallel")){
  if(! rpack %in% installed.packages()[,"Package"]){
    install.packages(rpack, repos = "https://cloud.r-project.org/")
  }
}
#i# Report installation status
for(rpack in c("tidyverse","RColorBrewer","gtools","ggstatsplot","writexl","parallel")){
  logWrite(paste('#RPACK Package',rpack,'installed:',rpack %in% installed.packages()[,"Package"]))
}

logWrite("#RCODE install_packages.R finished.")


