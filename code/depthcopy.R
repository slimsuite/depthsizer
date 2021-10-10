########################################################
### DepthCopy SC depth functions               ~~~~~ ###
### VERSION: 0.5.0                             ~~~~~ ###
### LAST EDIT: 07/10/21                        ~~~~~ ###
### AUTHORS: Richard Edwards 2021              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for DepthCopy and DepthSizer single copy depth analyses

################# ::: HISTORY ::: ######################
# v0.0.0 : Initial version based on DepthCopy.Rmd data exploration
# v0.1.0 : Added PNG plot outputs.
# v0.2.0 : Added window-based region calculations.
# v0.3.0 : Added katfile=FILE, katself=FILE, homfile=FILE, and chromcheck=LIST
# v0.3.1 : Fixed up kmer/homology plots and added pngdir=PATH [pngplots]
# v0.3.2 : Updated log output for DepthSizer and DepthKopy python programs.
# v0.4.0 : Added ggstatsplot and BUSCO v5 table parsing.
# v0.5.0 : Added seqstats output for full-length assembly sequences.

################# ::: USAGE ::: ######################
# Example use:
# Rscript depmode.R depfile=FILE [busco=FILE] [scdepth=INT] [regfile=FILE] [reghead=LIST] [gfftype=LIST] [chromcheck=LIST] [katfile=FILE] [katself=FILE] [homfile=FILE] 
# : depfile=FILE = Full depth data table (one value per base, pseudo-fasta format)
# : busco=FILE = Full BUSCO table. If a number, will interpret as scdepth value.
# : regfile=FILE = Optional delimited file or GFF to calculate depths for
# : reghead=LIST = Optional SeqName,Start,End or GFF feature type. Defaults to SeqName,Start,End or "gene"
# : gfftype=LIST = List of GFF types for GFF input ["gene"]
# : chromcheck=LIST = List of chromosomes (or min. chrom size) for individual window analysis
# : seqstats=T/F = Whether to output statistics for each sequence
# : katfile=FILE = high accuracy read KAT kmer frequencies for plotting 
# : katself=FILE = assembly vs self KAT kmer frequencies for plotting 
# : homfile=FILE = assembly self-mapping homology depth file for plotting 
# : pngdir=PATH = output path for PNG graphics

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}depthcopy.R {1}'.format(rdir, optionstr)).readlines()

################# ::: OUTPUTS ::: ######################
#!# List the outputs here

################# ::: TODO ::: ######################
#!# Add parallel processing of regions if speedup needed.
#?# Add additional statistics, including kmer frequencies. (Another tool?)
#!# Add processing of all individual sequences as extra output.

################# ::: SETUP ::: ######################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
settings = list(adjust=12,scdepth=0,busco="",depfile="",regfile="",threads=1,
                winsize=0,winstep=1,debug=FALSE,chromcheck="None",
                katfile='None', katself='None', homfile='None',
                seqstats=FALSE,
                pngwidth=1200,pngheight=900,pointsize=16,pngdir="pngplots",
                basefile="depthcopy",
                gfftype="gene",reghead="SeqName,Start,End",outlog=stdout())
argvec = commandArgs(TRUE)
if(FALSE){
  argvec = strsplit("winsize=100000 depfile=Mqui_v1.fastmp regfile=Mqui_v1.easy.lca_genes.gff busco=full_table_Mqui_v1.tsv basefile=test chromcheck=1e6 chromcheck=0 katfile=Mqui_v1.kat-hifi-counts.cvg katself=Mqui_v1.kat-counts.cvg homfile=Mqui_v1.self.fasthom",split=" ")[[1]]
}
for(cmd in argvec){
  cmdv = strsplit(cmd,'=',TRUE)[[1]]
  if(length(cmdv) > 1){
    settings[[cmdv[1]]] = cmdv[2]
  }else{
    settings[[cmdv[1]]] = TRUE    
  }
}
for(cmd in c("adjust","pngwidth","pngheight","pointsize","threads","winsize")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
for(cmd in c("scdepth","winstep")){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
for(cmd in c("reghead","gfftype","chromcheck")){
  settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
}
for(cmd in c("debug","seqstats")){
  settings[[cmd]] = as.logical(settings[[cmd]])
}

oldwarn <- getOption("warn")

if(settings$debug){
  writeLines(argvec)
}else{
  options(warn = -1)
}

logWrite <- function(logstr){
  writeLines(paste0("[",date(),"] ",logstr),con=settings$outlog)
}
for(cmd in names(settings)[order(names(settings))]){
  logWrite(paste("CMD:",cmd,"=",paste0(settings[[cmd]],collapse=",")))
}
settings$buscofile = settings$busco
dir.create(settings$pngdir, showWarnings = FALSE)

### ~ Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
settings$ggstatsplot = "ggstatsplot" %in% installed.packages()[,"Package"]
if(! "tidyverse" %in% installed.packages()[,"Package"]){
  install.packages("tidyverse")
}
library(tidyverse)
if(settings$ggstatsplot){
  library(ggstatsplot)
}

################# ::: FUNCTIONS ::: ######################

### ~ Load Depth File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load lengths of depths or kmers in SeqName=Vector list
#i# deplist = depthList(filename)
depthList <- function(filename){
  deplist = list()
  indata = readLines(filename)
  i = 1
  logWrite(paste0("Loading depth data from ",filename,"..."))
  while(i < length(indata)){
    cat(".", file = stderr())
    seqname = indata[i]
    seqname = substr(seqname,2,nchar(seqname))
    deplist[[seqname]] = as.integer(strsplit(indata[i+1]," ")[[1]])
    i = i + 2
  }
  cat("\n", file = stderr())
  logWrite(paste("#DEPTH Depth data for",length(deplist),"sequences loaded from",filename))
  return(deplist)  
}



### ~ Load BUSCO File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load BUSCO table into tibble
#i# buscodb = buscoTable(filename)
#i# BUSCOv3 headings:
v3head = c("BuscoID","Status","Contig","Start","End","Score","Length")
#i# BUSCOv5 headings:
v5head = c("BuscoID","Status","Contig","Start","End","Strand","Score","Length","OrthoDBURL","Description")
buscoTable <- function(filename){
  buscodb = read.table(filename,fill=TRUE,row.names = NULL,sep="\t")
  if(ncol(buscodb) > length(v3head)){
    logWrite(paste("#BUSCOV BUSCO v5 format"))
    colnames(buscodb) = v5head
  }else{
    logWrite(paste("#BUSCOV BUSCO v3 format"))
    colnames(buscodb) = v3head
  }
  buscodb$Contig = as.character(buscodb$Contig)
  buscodb = buscodb[buscodb$Status == "Complete",]
  #logWrite(paste(nrow(buscodb),"Complete BUSCO genes loaded from",filename))
  logWrite(paste('#BUSCO',nrow(buscodb),"Complete BUSCO genes loaded from",filename))
  return(buscodb)
}
#i# dupdb = buscoDupTable(filename)
buscoDupTable <- function(filename){
  buscodb = read.table(filename,fill=TRUE,row.names = NULL,sep="\t")
  if(ncol(buscodb) > length(v3head)){
    #logWrite(paste("#BUSCOV BUSCO v5 format"))
    colnames(buscodb) = v5head
  }else{
    #logWrite(paste("#BUSCOV BUSCO v3 format"))
    colnames(buscodb) = v3head
  }
  buscodb$Contig = as.character(buscodb$Contig)
  buscodb = buscodb[buscodb$Status == "Duplicated",]
  #logWrite(paste(nrow(buscodb),"Duplicated BUSCO genes loaded from",filename))
  logWrite(paste('#BUSCO',nrow(buscodb),"Duplicated BUSCO genes loaded from",filename))
  return(buscodb)
}

### ~ Load Region File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into Region file
#i# regdb = regionTable(filename,delimit="\t",uniqreg=FALSE)
regionTable <- function(filename,delimit="\t",uniqreg=FALSE){
  regdb = read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL)
  if(uniqreg){
    logWrite(paste(nrow(regdb),"lines loaded from",filename))
    #?# Reduce to unique SeqName, Start, End
    regdb = regdb %>% select(headers) %>% distinct()
    #logWrite(paste(nrow(regdb),"unique regions loaded from",filename))
    logWrite(paste('#REGION',nrow(regdb),"unique regions loaded from",filename))
  }else{
    #logWrite(paste(nrow(regdb),"regions loaded from",filename))
    logWrite(paste('#REGION',nrow(regdb),"regions loaded from",filename))
  }
  return(regdb)
}

#i# Load GFF
gffIDFromAtt <- function(attstr){
  for(attel in strsplit(attstr,";",TRUE)[[1]]){
    if(startsWith(attel,"ID=") | startsWith(attel,"id=")){
      return(strsplit(attel,"=",TRUE)[[1]][2])
    }
  }
  return("")
}
#i# gffdb = gffTable(filename,gfftype="gene")
# > filename = GFF file
# > gfftype = vector of allowed feature types. * will keep everything.
gffTable <- function(filename,gfftype="gene"){
  gffdb = read.table(filename,sep="\t",fill=TRUE,header=FALSE,row.names = NULL)
  colnames(gffdb) = c('SeqName', 'Source', 'FType', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes')
  if(! "*" %in% gfftype){
    gffdb = gffdb %>% filter(FType==gfftype)
  }
  gffdb$ID = ""
  for(i in 1:nrow(gffdb)){
    gffdb$ID[i] = gffIDFromAtt(gffdb$Attributes[i]) 
  }
  gffdb = gffdb %>% select(SeqName, Source, FType, Start, End, Strand, ID, Attributes)
  #logWrite(paste(nrow(gffdb),"regions loaded from",filename))
  logWrite(paste('#GFF',nrow(gffdb),"GFF features loaded from",filename))
  return(gffdb)
}

### ~ Tidy tables for output ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Sets dp etc. for output fields in a table then returns for output
tidyTable <- function(regdb){
  ifields = c("Start","End",reghead[2],reghead[3],"ModeX")
  dp2fields = c("MeanX","MedX","DensX","DensK","HomPC")
  dp3fields = c("SelfK","CN","CIsyst","CIrand")
  #i# Integers
  for(colname in ifields){
    if(colname %in% colnames(regdb)){
      regdb[[colname]] = as.integer(regdb[[colname]])
    }
  }
  #i# 2 d.p.
  for(colname in dp2fields){
    if(colname %in% colnames(regdb)){
      regdb[[colname]] = round(regdb[[colname]],2)
    }
  }
  #i# 3 d.p.
  for(colname in dp3fields){
    if(colname %in% colnames(regdb)){
      regdb[[colname]] = round(regdb[[colname]],2)
    }
  }
  return(regdb)
}
 

### ~ Calculate SC Depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Straight mode of depth vector
getmode <- function(v) {
  uniqv <- unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}
#i# Density plotting function
densModePlot <- function(bdens,centre,scdepth,plotmain){
  
  plot(bdens,xlim=c(0,scdepth*2), lwd = 2,col="blue",
       main=plotmain,ylab="Depth Density",xlab="XDepth",
       ylim=c(0,max(bdens$y)))
  
  abline(v=centre,col="grey")
  text(centre-2,0,round(centre,2),adj=c(1,0),col="grey")
  
  dmode=bdens$x[bdens$y == max(bdens$y)]
  abline(v=dmode,col="steelblue")
  text(dmode+2,0,round(dmode,2),adj=0,col="steelblue")
  abline(h=max(bdens$y),col="steelblue")
  
  abline(v=scdepth,col="red")
  text(scdepth+2,max(bdens$y)/2,round(scdepth,2),adj=0,col="red")
  
}
#i# Density-based mode calculation
densModeZoom <- function(depvec,adjust=16,plotbase=NA,plotmain=NA){
  depvec = depvec[!is.na(depvec)]
  if(length(depvec) < 1){
    if(settings$debug){
      logWrite("DepVec length <1")
    }
    return(0)
  }
  centre = getmode(depvec)
  if(centre <= 0){
    if(settings$debug){
      logWrite("Centre <= 0")
    }
    return(0)
  }
  #newvec = depvec[depvec <= (centre*4) & depvec >= centre/4]
  newvec = depvec[depvec <= max(1000,centre*4)]
  if(length(newvec) < 1){
    if(settings$debug){
      logWrite("NewVec length <1")
    }
    return(0)
  }
  n = 2048
  while(max(newvec)*5 > n){
    n = n * 2
  }
  depdens = density(newvec,n=n,adjust=adjust)
  scdepth = depdens$x[depdens$y == max(depdens$y)]
  if(! is.na(plotbase)){
    bdens = density(newvec,n=n)
    ## --- ##
    pngfile = paste0(plotbase,".raw.png")
    pngfile = file.path(settings$pngdir,pngfile)
    png(pngfile,width=settings$pngwidth,height=settings$pngheight,pointsize=settings$pointsize)
    densModePlot(bdens,centre,scdepth,paste(plotmain,"raw density profile"))
    dev.off()
    logWrite(paste(plotmain,"raw density profile output to:",pngfile))
    ## --- ##
    pngfile = paste0(plotbase,".scdepth.png")
    pngfile = file.path(settings$pngdir,pngfile)
    png(pngfile,width=settings$pngwidth,height=settings$pngheight,pointsize=settings$pointsize)
    densModePlot(depdens,centre,scdepth,paste(plotmain,"smoothed density profile"))
    dev.off()
    logWrite(paste(plotmain,"adjusted density profile output to:",pngfile))
    ## --- ##
  }

  return(scdepth)
}

### ~ Predict CN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Extract depth profiles from BUSCO Complete data
buscoDepData <- function(buscofile,deplist){
  buscodb = buscoTable(buscofile)
  buscovec = list()  # List of BuscoID: depvec
  depvec = c()       # Total depth vector
  # Extract depth vectors
  for(i in 1:nrow(buscodb)){
    seqname = buscodb$Contig[i]
    posi = buscodb$Start[i]
    posj = buscodb$End[i]
    seqdep = deplist[[seqname]][posi:posj]
    buscovec[[buscodb$BuscoID[i]]] = seqdep
    if(settings$debug){
      writeLines(paste(i,buscodb$BuscoID[i],seqname,posi,posj,length(seqdep)))
    }
    depvec = c(depvec,seqdep)
    cat(".", file = stderr())
  }
  cat("\n", file = stderr())
  writeLines("Calculating BUSCO SC Depth...")
  print(summary(depvec))
  scdepth = densModeZoom(depvec,adjust=settings$adjust,plotbase=settings$basefile,plotmain="BUSCO")
  #!# Generate PNG BUSCO depth profile
  logWrite(paste("#SCDPETH BUSCO SC Depth:",round(scdepth,2)))
  ### ~ Save SC Depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  if(scdepth > 0){
    outfile = paste(settings$depfile,"scdepth",sep=".",collapse=".")
    logWrite(paste(nrow(buscodb),"BUSCO Complete SC Depth output to",outfile))
    write(scdepth,outfile)
  }else{
    logWrite(paste(nrow(buscodb),"BUSCO Complete SC Depth calculation failed!"))
    quit("no",1)
  }
  return( list(depvec=depvec, buscodb=buscodb, buscovec=buscovec, scdepth=scdepth) )
}

#i# Combine BUSCO data with depth profiles
#># buscodat = buscoDepCalc(buscodat)
buscoDepCalc <- function(buscodat){
  depvec = buscodat$depvec
  scdepth = buscodat$scdepth
  buscodb = buscodat$buscodb
  buscovec = buscodat$buscovec
  # Setup new fields. Some of these are just for testing
  buscodb$MeanX = 0
  buscodb$MedX = 0  
  buscodb$ModeX = 0  
  buscodb$DensX = 0  
  buscodb$CN = 0
  buscodb$CIsyst = 0.0
  buscodb$CIrand = 0.0
  # Calculate gene stats
  for(i in 1:nrow(buscodb)){
    #logWrite(buscodb$BuscoID[i])
    seqdep = buscovec[[buscodb$BuscoID[i]]]
    if(length(seqdep) >= 1){
      buscodb$MeanX[i] = mean(seqdep)
      buscodb$ModeX[i] = getmode(seqdep)
      buscodb$MedX[i] = median(seqdep)
      buscodb$DensX[i] = densModeZoom(seqdep,adjust=settings$adjust)
    }else{
      logWrite(paste("#WARN BUSCOID",buscodb$BuscoID[i],"has zero-length depth vector"))
    }
    cat(".", file = stderr())
  }
  cat("\n", file = stderr())
  # Return updated data
  logWrite("#BUSCO BUSCO depth calculations complete.")
  # Add Homology / Kmer data if present
  buscodb$DensK = NA
  buscodb$SelfK = NA
  buscodb$HomPC = NA
  writeLines("Calculating BUSCO kmers/homology...")
  # Calculate region stats
  prevname = ""
  for(i in 1:nrow(buscodb)){
    seqname = buscodb$Contig[i]
    if(seqname != prevname){
      prevname = seqname
      if(settings$debug){
        cat(seqname)
      }
    }
    cat(".", file = stderr())
    seqi = buscodb$Start[i]
    seqj = buscodb$End[i]
    if(file.exists(katfile)){
      #i# katlist = raw read kmer frequencies
      seqdep = depVector(seqname,seqi,seqj,katlist)
      buscodb$DensK[i] = densModeZoom(seqdep,adjust=settings$adjust)
    }
    if(file.exists(katself)){
      #i# katself = assembly kmer frequencies
      seqdep = depVector(seqname,seqi,seqj,selflist)
      buscodb$SelfK[i] = densModeZoom(seqdep,adjust=settings$adjust)
    }
    if(file.exists(homfile)){
      #i# homlist = Assembly vs self homology
      seqdep = depVector(seqname,seqi,seqj,homlist)
      #buscodb$HomX[i] = densModeZoom(seqdep,adjust=settings$adjust)
      buscodb$HomPC[i] = 100.0 * sum(seqdep > 0) / length(seqdep)
    }
  }
  cat("\n", file = stderr())
  logWrite("#BUSCO BUSCO kmer/homology calculations complete.")
  return( list(depvec=depvec, buscodb=buscodb, buscovec=buscovec, scdepth=scdepth))
}

#!# Cut down fields and reformat long to calculate CN easily.
buscoCN <- function(buscodat){
  scdepth = buscodat$scdepth
  buscodb = buscodat$buscodb
  # # Reshape the table
  # buscodb = pivot_longer(buscodb,colnames(buscodb)[8:ncol(buscodb)],
  #                        names_to = "Method",
  #                        values_to = "XDepth"
  # ) %>% select(BuscoID, Contig, Start, End, Length, Method, XDepth)
  # Calculate CN
  buscodb$CN = buscodb$DensX / scdepth
  # Calculate CI
  buscoMean = mean(buscodb$DensX)
  buscoSD = sd(buscodb$DensX)
  buscodb$CIsyst = 1.96 * buscoSD * sqrt(buscodb$DensX / (buscoMean ** 3))
  buscodb$CIrand = 1.96 * buscoSD * buscodb$DensX / (buscoMean ** 2)
  return(buscodb)
}

### ~ Region XDepth and CN calclations ~~~~~~~~~~~~~~~ ###
depVector <- function(seqname,posi,posj,deplist){
  return(deplist[[seqname]][posi:posj])
}
regCN <- function(regdb,buscoMean=0,buscoSD=0){
  regdb$MeanX = 0
  regdb$MedX = 0  
  regdb$ModeX = 0  
  regdb$DensX = 0  
  regdb$CN = 0
  regdb$CIsyst = 0.0
  regdb$CIrand = 0.0
  logWrite("Calculating region CN...")
  # Calculate gene stats
  for(i in 1:nrow(regdb)){
    seqname = regdb[i,reghead[1]]
    seqi = regdb[i,reghead[2]]
    seqj = regdb[i,reghead[3]]
    seqdep = depVector(seqname,seqi,seqj,deplist)
    if(length(seqdep) >= 1){
      regdb$MeanX[i] = mean(seqdep)
      regdb$ModeX[i] = getmode(seqdep)
      regdb$MedX[i] = median(seqdep)
      regdb$DensX[i] = densModeZoom(seqdep,adjust=settings$adjust)
      #i# Performing per row calculations for reporting
      regdb$CN[i] = regdb$DensX[i] / scdepth
      regdb$CIsyst[i] = 1.96 * buscoSD * sqrt(regdb$DensX[i] / (buscoMean ** 3))
      regdb$CIrand[i] = 1.96 * buscoSD * regdb$DensX[i] / (buscoMean ** 2)
      writeLines(paste0(seqname,":",seqi,"-",seqj," = ",round(regdb$DensX[i],2),"X = ",round(regdb$CN[i],2),"N +/- ",round(regdb$CIrand[i],2)))
    }else{
      logWrite(paste0("#WARN Region ",seqname,":",seqi,"-",seqj," has zero-length depth vector"))
    }
    
  }
  logWrite("#CNCALC Region CN calculations complete.")
  # Add Homology / Kmer data if present
  regdb$DensK = NA
  regdb$SelfK = NA
  #regdb$HomX = NA
  regdb$HomPC = NA
  logWrite("Calculating region kmers/homology...")
  # Calculate region stats
  prevname = ""
  for(i in 1:nrow(regdb)){
    seqname = regdb[i,reghead[1]]
    if(seqname != prevname){
      prevname = seqname
      if(settings$debug){
        cat(seqname)
      }
    }
    cat(".", file = stderr())
    seqi = regdb[i,reghead[2]]
    seqj = regdb[i,reghead[3]]
    if(file.exists(katfile)){
      #i# katlist = raw read kmer frequencies
      seqdep = depVector(seqname,seqi,seqj,katlist)
      regdb$DensK[i] = densModeZoom(seqdep,adjust=settings$adjust)
    }
    if(file.exists(katself)){
      #i# katself = assembly kmer frequencies
      seqdep = depVector(seqname,seqi,seqj,selflist)
      regdb$SelfK[i] = densModeZoom(seqdep,adjust=settings$adjust)
    }
    if(file.exists(homfile)){
      #i# homlist = Assembly vs self homology
      seqdep = depVector(seqname,seqi,seqj,homlist)
      #regdb$HomX[i] = densModeZoom(seqdep,adjust=settings$adjust)
      if(length(seqdep) > 0){
        regdb$HomPC[i] = 100.0 * sum(seqdep > 0) / length(seqdep)
      }
    }
  }
  cat("\n", file = stderr())
  logWrite("#KCALC Region kmer/homology calculations complete.")

  return(regdb)
}


################# ::: RUN CODE ::: ######################
### ~ Key inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Depth file of read depths per base
depfile = settings$depfile
logWrite(paste("Depth File:",depfile))
if(! file.exists(depfile)){
  quit("no",2)  
}
#i# BUSCO table or scdepth
buscofile = settings$buscofile
scdepth = as.numeric(settings$scdepth)
if(file.exists(buscofile)){
  logWrite(paste("BUSCO Full File:",buscofile))
}else{
  logWrite(paste("SC Depth:",round(scdepth,2)))
}
#i# Region file
regfile = settings$regfile
reghead = c("SeqName","Start","End")
if(file.exists(regfile)){
  logWrite(paste("Region File:",regfile))
}
#i# Window-based analysis
if(settings$winsize>0){
  if(settings$winstep<=1){
    settings$winstep = settings$winstep * settings$winsize
  }
  logWrite(paste("Window size:",settings$winsize,"bp"))
  logWrite(paste("Window step:",settings$winstep,"bp"))
}
logWrite('#RCODE Setup complete.')

### ~ Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# 1. Load the depth data list
deplist = depthList(depfile)

## ~ Additional vectors for cross-plotting ~ ##
#i# katfile=FILE, katself=FILE, homfile=FILE
katfile = settings$katfile
logWrite(paste("KAT File:",katfile))
if(file.exists(katfile)){
  katlist = depthList(katfile)
  # Pad kmer counts
  for(seqname in names(katlist)){
    seqlen = length(deplist[[seqname]])
    seqkat = katlist[[seqname]]
    k = seqlen - length(seqkat)
    if(settings$debug){
      writeLines(paste(seqname,"Pad",k,"kmers"))
    }
    i = as.integer(k/2)
    if(k > 0){
      katlist[[seqname]] = c( rep(seqkat[1],i), seqkat, rep(seqkat[length(seqkat)],k-i) )
    }
  }
}
#i# Assembly versus self KAT
katself = settings$katself
logWrite(paste("Self-KAT File:",katself))
if(file.exists(katself)){
  selflist = depthList(katself)
  # Pad kmer counts
  for(seqname in names(selflist)){
    seqlen = length(deplist[[seqname]])
    seqkat = selflist[[seqname]]
    k = seqlen - length(seqkat)
    if(settings$debug){
      writeLines(paste(seqname,"Pad",k,"kmers"))
    }
    i = as.integer(k/2)
    if(k > 0){
      selflist[[seqname]] = c( rep(seqkat[1],i), seqkat, rep(seqkat[length(seqkat)],k-i) )
    }
  }
}
#i# Assembly vs self homology
homfile = settings$homfile
logWrite(paste("Self-Homology File:",homfile))
if(file.exists(homfile)){
  homlist = depthList(homfile)
}

### ~ Calculate SC Depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# 2. Calculate SC read depth if mot provided
buscoMean = 0
buscoSD = 0
if(file.exists(buscofile)){
  #!# Add ability to re-use loaded BUSCO data and given scdepth=NUM value
  buscodat = buscoDepData(buscofile,deplist)
  scdepth = buscodat$scdepth
  logWrite("Calculating BUSCO gene depth stats...")
  buscodat = buscoDepCalc(buscodat)
  buscodb = buscoCN(buscodat)
  buscoMean = mean(buscodb$DensX)
  buscoSD = sd(buscodb$DensX)
  
  ### ~ Save BUSCO table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  #i# Output table to verify details
  outfile = basename(buscofile)
  outvec = strsplit(outfile,".",TRUE)[[1]]
  outfile = paste(c(outvec[1:length(outvec)-1],"regcnv.tsv"),sep=".",collapse=".")
  logWrite(paste("#SAVE",nrow(buscodb),"BUSCO genes output to",outfile))
  write.table(tidyTable(buscodb),outfile,sep="\t",quote=FALSE,row.names=FALSE)
  
  ### ~ Duplicated BUSCO genes ~~~~~~~~~~~~~~~~~~~~~~~~ ###
  dupdb = buscoDupTable(buscofile)
  #i# Calculate CN per gene
  oldhead = reghead
  reghead = c("Contig","Start","End")
  dupdb <- regCN(as.data.frame(dupdb),buscoMean,buscoSD)
  reghead = oldhead
  ### ~ Save CN table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  #i# Output table to with CN data
  outfile = paste(c(outvec[1:length(outvec)-1],"dupcnv.tsv"),sep=".",collapse=".")
  logWrite(paste("#SAVE",nrow(dupdb),"Duplicated BUSCO genes output to",outfile))
  write.table(tidyTable(dupdb),outfile,sep="\t",quote=FALSE,row.names=FALSE)
}

### ~ Calculate Region CN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# 3. Load region file
if(file.exists(regfile)){
  #i# GFF file
  if(sum(endsWith(regfile,c("gff","gff3")))){
    gtype = settings$gfftype
    regdb = gffTable(regfile,gtype)
  }else{
    reghead = strsplit(settings$reghead,",",TRUE)[[1]]
    regdb = regionTable(regfile)
  }
#i# 4. Calculate CN per region
  regdb <- regCN(regdb,buscoMean,buscoSD)

### ~ Save CN table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  #i# Output table to with CN data
  outfile = basename(regfile)
  outvec = strsplit(outfile,".",TRUE)[[1]]
  outfile = paste(c(outvec[1:length(outvec)-1],"regcnv.tsv"),sep=".",collapse=".")
  logWrite(paste("#SAVE",nrow(regdb),"CN predictions output to",outfile))
  write.table(tidyTable(regdb),outfile,sep="\t",quote=FALSE,row.names=FALSE)
  
}

### ~ Calculate Window CN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Make window regions table
if(settings$winsize>0){
  windb = tibble(SeqName=c(),Start=c(),End=c(),SeqLen=c())
  for(seqname in names(deplist)){
    seqlen = length(deplist[[seqname]])
    if(seqlen <= settings$winsize){
      windb <- bind_rows(windb,
                         tibble(SeqName=seqname,Start=1,End=seqlen,SeqLen=seqlen))
      logWrite(paste0(seqname,"(",seqlen," bp) => 1 window"))
    }else{
      winx = as.integer(seqlen / settings$winstep) + 1
      starts = 0:winx * settings$winstep + 1
      ends = starts + settings$winsize - 1
      seqdb <- tibble(Start=starts,End=ends) %>% 
        filter(End <= seqlen) 
      seqdb <- seqdb %>%
        add_row(Start=max(seqdb$Start)+settings$winstep,
                End=min(seqlen,max(seqdb$End)+settings$winstep)) %>%
        mutate(SeqName=seqname,SeqLen=seqlen)
      windb <- bind_rows(windb,seqdb) %>% select(SeqName,Start,End,SeqLen)
      writeLines(paste0(seqname,"(",seqlen," bp) => ",nrow(seqdb)," windows"))
    }
  }
  #i# Calculate CN per region
  oldhead = reghead
  reghead = c("SeqName","Start","End")
  windb <- regCN(as.data.frame(windb),buscoMean,buscoSD)
  reghead = oldhead
  ### ~ Save CN table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  #i# Output table to with CN data
  outfile = paste0(settings$basefile,".windows.regcnv.tsv")
  logWrite(paste("#SAVE",nrow(windb),"CN predictions output to",outfile))
  write.table(tidyTable(windb),outfile,sep="\t",quote=FALSE,row.names=FALSE)
}


### ~ Calculate SeqStats CN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Make whole sequences regions table
if(settings$seqstats){
  scaffdb = tibble(SeqName=c(),Start=c(),End=c())
  for(seqname in names(deplist)){
    seqlen = length(deplist[[seqname]])
    scaffdb <- bind_rows(scaffdb,
                         tibble(SeqName=seqname,Start=1,End=seqlen))
  }
  #i# Calculate CN per sequence
  oldhead = reghead
  reghead = c("SeqName","Start","End")
  scaffdb <- regCN(as.data.frame(scaffdb),buscoMean,buscoSD)
  reghead = oldhead
  ### ~ Save CN table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  #i# Output table to with CN data
  outfile = paste0(settings$basefile,".sequences.regcnv.tsv")
  logWrite(paste("#SAVE CN statistics for",nrow(scaffdb),"sequences output to",outfile))
  write.table(tidyTable(scaffdb),outfile,sep="\t",quote=FALSE,row.names=FALSE)
}

### ~ Generate Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Merge buscodb and regdb with Dataset field
plotdb = tibble(Dataset=c(),SeqName=c(),Start=c(),End=c(),MeanX=c(),MedX=c(),ModeX=c(),DensX=c(),CN=c(),DensK=c(),SelfK=c(),HomPC=c())
if(file.exists(regfile)){
  plotdb <- bind_rows(plotdb, 
                      regdb %>% rename(SeqName=reghead[1],Start=reghead[2],End=reghead[3]) %>%
                        mutate(Dataset="Regions") %>%
                        select(Dataset,SeqName,Start,End,MeanX,MedX,ModeX,DensX,CN,DensK,SelfK,HomPC))
}
if(file.exists(buscofile)){
  plotdb <- bind_rows(plotdb, 
                      buscodb %>% rename(SeqName=Contig) %>%
                        mutate(Dataset="BUSCO") %>%
                        select(Dataset,SeqName,Start,End,MeanX,MedX,ModeX,DensX,CN,DensK,SelfK,HomPC))
  plotdb <- bind_rows(plotdb, 
                      dupdb %>% rename(SeqName=Contig) %>%
                        mutate(Dataset="Duplicated") %>%
                        select(Dataset,SeqName,Start,End,MeanX,MedX,ModeX,DensX,CN,DensK,SelfK,HomPC))
}
if(settings$seqstats){
  plotdb <- bind_rows(plotdb, 
                      scaffdb %>% 
                        mutate(Dataset="Sequences") %>%
                        select(Dataset,SeqName,Start,End,MeanX,MedX,ModeX,DensX,CN,DensK,SelfK,HomPC))
}
#i# Window analysis, incorporating chromcheck
if(settings$winsize>0){
  #i# Setup chromcheck
  chromcheck = settings$chromcheck
  chromsize = as.integer(chromcheck[1])
  if(! chromcheck[1] %in% windb$SeqName & ! is.na(chromsize) & chromsize > 0 & max(windb$End) >= chromsize){
    #i# Set chromosome check to those sequences exceeding length threshold
    logWrite(paste("Chromosome check:",chromsize,"bp+"))
    chromcheck = unique(windb[windb$End >= chromsize,]$SeqName)
  }
  #i# Chromcheck subsets
  if(length(chromcheck) > 0 & chromcheck[1] %in% windb$SeqName){
    logWrite(paste("Chromosome check:",length(chromcheck),"sequences"))
    for(chrom in chromcheck){
      plotdb <- bind_rows(plotdb, 
                          windb %>% filter(SeqName==chrom) %>%
                            rename(SeqName=reghead[1],Start=reghead[2],End=reghead[3]) %>%
                            mutate(Dataset=chrom) %>%
                            select(Dataset,SeqName,Start,End,MeanX,MedX,ModeX,DensX,CN,DensK,SelfK,HomPC))
    }
    plotdb <- bind_rows(plotdb, 
                        windb %>% filter(SeqName!=chromcheck) %>%
                          rename(SeqName=reghead[1],Start=reghead[2],End=reghead[3]) %>%
                          mutate(Dataset="Other") %>%
                          select(Dataset,SeqName,Start,End,MeanX,MedX,ModeX,DensX,CN,DensK,SelfK,HomPC))
  }else{
    #i# Full Windows data
    plotdb <- bind_rows(plotdb, 
                        windb %>% rename(SeqName=reghead[1],Start=reghead[2],End=reghead[3]) %>%
                          mutate(Dataset="Windows") %>%
                          select(Dataset,SeqName,Start,End,MeanX,MedX,ModeX,DensX,CN,DensK,SelfK,HomPC))
  }
}

#i# Violin plot function
vioPlot <- function(plotdb,plotfield){
  #i# Violin plots by Method
  #!# Add labels of n=X for each Dataset
  labels = levels(factor(plotdb$Dataset))
  labels = paste("n =",table(factor(plotdb$Dataset)))
  nx = 1:length(labels)
  ny = rep(0,length(labels))
  p <- ggplot(plotdb, aes(factor(Dataset), .data[[plotfield]])) + 
    aes(fill = factor(Dataset)) +
    geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + 
    xlab("Dataset") + ylab(plotfield) +
    labs(fill="Dataset") +
    #geom_text(x=nx, y=ny, label=labels, vjust=1) +
    theme_bw(base_size = settings$pointsize)
  
  if(settings$ggstatsplot){
    splotdb = plotdb
    splotdb$plotfield = splotdb[[plotfield]]
    p <- ggbetweenstats(
      data = splotdb,
      results.subtitle = length(labels) <= 4,
      x = Dataset,
      y = plotfield
    ) + ylim(0,median(splotdb$plotfield)*4) +
      #geom_hline(yintercept=1, color="black") +
      labs(
        #title = "DepthSizer accuracy versus reference assembly size",
        x = "Dataset",
        y = plotfield
      ) +
      theme(text = element_text(size=settings$pointsize))

  }

  return(p)
}

#i# Generate PNG of violin plots
plotfields = c("MeanX","MedX","ModeX","DensX","CN")
if(file.exists(katfile)){
  plotfields = c(plotfields, "DensK")
}
if(file.exists(katself)){
  plotfields = c(plotfields, "SelfK")
}
if(file.exists(homfile)){
  plotfields = c(plotfields, "HomPC")
}

if(nrow(plotdb) > 0){
  for(plotfield in plotfields){
    pngfile = paste(settings$basefile,plotfield,"png",sep=".")
    pngfile = file.path(settings$pngdir,pngfile)
    png(pngfile,width=settings$pngwidth,height=settings$pngheight,pointsize=settings$pointsize*2)
    print(vioPlot(plotdb,plotfield))
    dev.off()
    logWrite(paste(plotfield,"violin plot output to:",pngfile))
  }
}

### ~ Depth / Kmer plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
hexPlot <- function(plotdb,plotfield){
  p = ggplot(plotdb, aes(x=CN, y=.data[[plotfield]]) ) +
  geom_hex(bins = 100) + #scale_x_log10() + scale_y_log10() +
  xlab("CN") + ylab(plotfield) +
  #coord_trans(x = "log2", y = "log2") +
  scale_x_continuous(trans = "log2") + 
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept=1, color="black") +
  # geom_hline(yintercept=0.5, linetype=2, color="#d29f23") +
  # geom_vline(xintercept=0.5, linetype=2, color="#d29f23") +
  # geom_hline(yintercept=2, linetype=2, color="#d29f23") +
  # geom_vline(xintercept=2, linetype=2, color="#d29f23") +
  theme_bw(base_size = settings$pointsize)
  
  if(! plotfield %in% c("HomPC")){
    p = p +  scale_y_continuous(trans = "log2") 
  }
  
  return(p)
}
#i# Generate PNG of hex plots
for(plotfield in plotfields[plotfields != "CN"]){
  for(dset in unique(plotdb$Dataset)){
    datdb = plotdb %>% filter(Dataset==dset)
    if(nrow(datdb) > 0){
      pngfile = paste(settings$basefile,dset,plotfield,"png",sep=".")
      pngfile = file.path(settings$pngdir,pngfile)
      png(pngfile,width=settings$pngwidth,height=settings$pngheight,pointsize=settings$pointsize*2)
      print(hexPlot(datdb,plotfield))
      dev.off()
      logWrite(paste(dset,plotfield,"hex plot output to:",pngfile))
    }
  }
}
logWrite(paste0('#PLOTS PNG graphics output to ',settings$pngdir,"/",settings$basefile,".*"))
if(settings$ggstatsplot){
  logWrite('#CITE If using violin plots, please cite:  Patil, I. (2021). J. Open Source Software, 6(61), 3167, doi:10.21105/joss.03167')
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE DepthCopy.R finished.")
quit("no",0)

