## ----pre,echo=FALSE,results='hide'-----------------------------------------
library(knitr)
opts_chunk$set(warning=FALSE,message=FALSE,cache=TRUE)

## ----style,echo=FALSE,results='asis'---------------------------------------
BiocStyle::markdown()

## ----loadLibrary-----------------------------------------------------------
library(dupRadar)

## ----eval=FALSE------------------------------------------------------------
#  # call the duplicate marker and analyze the reads
#  bamDuprm <- markDuplicates(dupremover="bamutil",
#  			   bam="test.bam",
#  	    	   path="/opt/bamUtil-master/bin",
#  			   rminput=FALSE)

## --------------------------------------------------------------------------
if(suppressWarnings(require(AnnotationHub))) {
    ah = AnnotationHub()
    query(ah, c("ensembl", "80", "Takifugu", "gtf"))    # discovery
    cache(ah["AH47101"])                                # retrieve / file path
}

## --------------------------------------------------------------------------
attach(dupRadar_examples)

## ----eval=FALSE------------------------------------------------------------
#  # The call parameters:
#  bamDuprm <- "test_duprm.bam"	# the duplicate marked bam file
#  gtf <- "genes.gtf"	# the gene model
#  stranded <- 2		# '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
#  paired   <- FALSE	# is the library paired end?
#  threads  <- 4		# number of threads to be used
#  
#  # Duplication rate analysis
#  dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)

## ----fig.width=14,fig.height=7---------------------------------------------
## make a duprate plot (blue cloud)
par(mfrow=c(1,2))
duprateExpDensPlot(DupMat=dm)		# a good looking plot
title("good experiment")
duprateExpDensPlot(DupMat=dm.bad)	# a dataset with duplication problems
title("duplication problems")

## duprate boxplot
duprateExpBoxplot(DupMat=dm)		# a good looking plot
duprateExpBoxplot(DupMat=dm.bad)	# a dataset with duplication problems

## ----fig.width=7,fig.height=7----------------------------------------------
duprateExpDensPlot(DupMat=dm)

# or, just to get the fitted model without plot
fit <- duprateExpFit(DupMat=dm)
cat("duprate at low read counts: ",fit$intercept,"\n",
	"progression of the duplication rate: ",fit$slope,fill=TRUE)

## ----eval=FALSE------------------------------------------------------------
#  ## INTERACTIVE: identify single points on screen (name="ID" column of dm)
#  duprateExpPlot(DupMat=dm)		# a good looking plot
#  duprateExpIdentify(DupMat=dm)

## ----fig.width=14,fig.height=7---------------------------------------------
par(mfrow=c(1,2))
cols <- colorRampPalette(c("black","blue","green","yellow","red"))
duprateExpPlot(DupMat=dm,addLegend=FALSE)
duprateExpPlot(DupMat=dm.bad,addLegend=FALSE,nrpoints=10,nbin=500,colramp=cols)

## ----fig.width=7,fig.height=7----------------------------------------------
readcountExpBoxplot(DupMat=dm)

## ----fig.width=7,fig.height=7----------------------------------------------
expressionHist(DupMat=dm)

## --------------------------------------------------------------------------
# calculate the fraction of multimappers per gene
dm$mhRate <- (dm$allCountsMulti - dm$allCounts) / dm$allCountsMulti
# how many genes are exclusively covered by multimappers
sum(dm$mhRate == 1, na.rm=TRUE)

# and how many have an RPKM (including multimappers) > 5
sum(dm$mhRate==1 & dm$RPKMMulti > 5, na.rm=TRUE)

# and which are they?
dm[dm$mhRate==1 & dm$RPKMMulti > 5, "ID"]

## ----fig.width=7,fig.height=7----------------------------------------------
hist(dm$mhRate, 
     breaks=50, 
     col="red", 
     main="Frequency of multimapping rates", 
     xlab="Multimapping rate per gene", 
     ylab="Frequency")


## ----fig.width=7,fig.height=7----------------------------------------------
# comparison of multi-mapping RPK and uniquely-mapping RPK
plot(log2(dm$RPK), 
     log2(dm$RPKMulti), 
     xlab="Reads per kb (uniquely mapping reads only)",
     ylab="Reads per kb (all including multimappers, non-weighted)"
)

## ---- eval=F---------------------------------------------------------------
#  identify(log2(dm$RPK),
#           log2(dm$RPKMulti),
#  	 labels=dm$ID)

## ----eval=F----------------------------------------------------------------
#  library(dupRadar)
#  library(biomaRt)
#  
#  ## for detailed explanations on biomaRt, please see the respective
#  ## vignette
#  
#  ## set up biomart connection for mouse (needs internet connection)
#  ensm <- useMart("ensembl")
#  ensm <- useDataset("mmusculus_gene_ensembl", mart=ensm)
#  
#  ## get a table which has the gene GC content for the IDs that have been
#  ## used to generate the table (depends on the GTF annotation that you
#  ## use)
#  tr <- getBM(attributes=c("mgi_symbol", "percentage_gc_content"),
#              values=TRUE, mart=ensm)
#  
#  ## create a GC vector with IDs as element names
#  mgi.gc <- tr$percentage_gc_content
#  names(mgi.gc) <- tr$mgi_symbol
#  
#  ## using dm demo duplication matrix that comes with the package
#  ## add GC content to our demo data and keep only subset for which we can
#  ## retrieve data
#  keep <- dm$ID %in% tr$mgi_symbol
#  dm.gc <- dm[keep,]
#  dm.gc$gc <- mgi.gc[dm.gc$ID]
#  
#  ## check distribution of annotated gene GC content (in %)
#  boxplot(dm.gc$gc, main="Gene GC content", ylab="% GC")

## ----eval=F----------------------------------------------------------------
#  par(mfrow=c(1,2))
#  
#  ## below median GC genes
#  duprateExpDensPlot(dm.gc[dm.gc$gc<=45,], main="below median GC genes")
#  
#  ## above median GC genes
#  duprateExpDensPlot(dm.gc[dm.gc$gc>=45,], main="above median GC genes")

## ----eval=F----------------------------------------------------------------
#  #!/usr/bin/env Rscript
#  
#  ########################################
#  ##
#  ## dupRadar shell script
#  ## call dupRadar R package from the shell for
#  ## easy integration into pipelines
#  ##
#  ## Holger Klein & Sergi Sayols
#  ##
#  ## https://sourceforge.net/projects/dupradar/
#  ##
#  ## input:
#  ## - _duplicate marked_ bam file
#  ## - gtf file
#  ## - parameters for duplication counting routine:
#  ##   stranded, paired, outdir, threads.
#  ##
#  ########################################
#  
#  library(dupRadar)
#  
#  ####################
#  ##
#  ## get name patterns from command line
#  ##
#  args   <- commandArgs(TRUE)
#  
#  ## the bam file to analyse
#  bam <- args[1]
#  ## usually, same GTF file as used in htseq-count
#  gtf <- gsub("gtf=","",args[2])
#  ## no|yes|reverse
#  stranded <- gsub("stranded=","",args[3])
#  ## is a paired end experiment
#  paired   <- gsub("paired=","",args[4])
#  ## output directory
#  outdir   <- gsub("outdir=","",args[5])
#  ## number of threads to be used
#  threads  <- as.integer(gsub("threads=","",args[6]))
#  
#  if(length(args) != 6) {
#    stop (paste0("Usage: ./dupRadar.sh <file.bam> <genes.gtf> ",
#  			   "<stranded=[no|yes|reverse]> paired=[yes|no] ",
#  			   "outdir=./ threads=1"))
#  }
#  
#  if(!file.exists(bam)) {
#    stop(paste("File",bam,"does NOT exist"))
#  }
#  
#  if(!file.exists(gtf)) {
#    stop(paste("File",gtf,"does NOT exist"))
#  }
#  
#  if(!file.exists(outdir)) {
#    stop(paste("Dir",outdir,"does NOT exist"))
#  }
#  
#  if(is.na(stranded) | !(grepl("no|yes|reverse",stranded))) {
#    stop("Stranded has to be no|yes|reverse")
#  }
#  
#  if(is.na(paired) | !(grepl("no|yes",paired))) {
#    stop("Paired has to be no|yes")
#  }
#  
#  if(is.na(threads)) {
#    stop("Threads has to be an integer number")
#  }
#  
#  stranded <- if(stranded == "no") 0 else if(stranded == "yes") 1 else 2
#  
#  ## end command line parsing
#  ##
#  ########################################
#  
#  ########################################
#  ##
#  ## analyze duprates and create plots
#  ##
#  cat("Processing file ", bam, " with GTF ", gtf, "\n")
#  
#  ## calculate duplication rate matrix
#  dm <- analyzeDuprates(bam,
#                        gtf,
#                        stranded,
#                        (paired == "yes"),
#                        threads)
#  
#  ## produce plots
#  
#  ## duprate vs. expression smooth scatter
#  png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_drescatter.png"),
#      width=1000, height=1000)
#  duprateExpDensPlot(dm, main=basename(bam))
#  dev.off()
#  
#  ## expression histogram
#  png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_ehist.png"),
#      width=1000, height=1000)
#  expressionHist(dm)
#  dev.off()
#  
#  ## duprate vs. expression boxplot
#  png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_drebp.png"),
#      width=1000, height=1000)
#  par(mar=c(10,4,4,2)+.1)
#  duprateExpBoxplot(dm, main=basename(bam))
#  dev.off()

## ----citation--------------------------------------------------------------
citation("dupRadar")

## ----echo=FALSE------------------------------------------------------------
sessionInfo()

