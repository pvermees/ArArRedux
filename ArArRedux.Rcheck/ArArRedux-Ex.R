pkgname <- "ArArRedux"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "ArArRedux-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('ArArRedux')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Melbourne")
### * Melbourne

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Melbourne
### Title: An example dataset
### Aliases: Melbourne

### ** Examples

data(Melbourne)
plotcorr(Melbourne$X)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Melbourne", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("average")
### * average

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: average
### Title: Calculate the arithmetic mean
### Aliases: average

### ** Examples

data(Melbourne)
K <- average(Melbourne$X,grep("K:",Melbourne$X$labels),newlabel="K-glass")
plotcorr(K)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("average", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("averagebyday")
### * averagebyday

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: averagebyday
### Title: Average all the data collected on the same day.
### Aliases: averagebyday

### ** Examples

dfile <- system.file("Calibration.csv",package="ArArRedux")
md <- loaddata(dfile)
ld <- fitlogratios(blankcorr(md))
d <- averagebyday(ld,"DCAL")
plotcorr(d)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("averagebyday", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("blankcorr")
### * blankcorr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blankcorr
### Title: Apply a blank correction
### Aliases: blankcorr blankcorr.default blankcorr.timeresolved
###   blankcorr.PHdata blankcorr.WiscAr

### ** Examples

samplefile <- system.file("Samples.csv",package="ArArRedux")
m <- loaddata(samplefile) # samples and J-standards
blanklabel <- "EXB#"
l <- fitlogratios(blankcorr(m,blanklabel),"Ar40")
plotcorr(l)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blankcorr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calibration")
### * calibration

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calibration
### Title: Detector calibration
### Aliases: calibration calibration.default calibration.WiscAr

### ** Examples

data(Melbourne)
C <- calibration(Melbourne$X,"DCAL")
plotcorr(C)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calibration", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("clcorrection")
### * clcorrection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: clcorrection
### Title: Cl-interference correction
### Aliases: clcorrection

### ** Examples

data(Melbourne)
Cl <- clcorrection(Melbourne$X,Melbourne$irr)
plotcorr(Cl)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("clcorrection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("concat")
### * concat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: concat
### Title: Merge a list of logratio data
### Aliases: concat concat.default concat.WiscAr

### ** Examples

samplefile <-  system.file("Samples.csv",package="ArArRedux")
kfile <- system.file("K-glass.csv",package="ArArRedux")
cafile <- system.file("Ca-salt.csv",package="ArArRedux")
dfile <- system.file("Calibration.csv",package="ArArRedux")
blanklabel <- "EXB#"
Jpos <- c(3,15)
 
m <- loaddata(samplefile) # samples and J-standards
mk <- loaddata(kfile) # K-interference data
mca <- loaddata(cafile) # Ca interference data
md <- loaddata(dfile) # detector intercalibrations
 
# form and fit logratios
l <- fitlogratios(blankcorr(m,blanklabel),"Ar40")
lk <- fitlogratios(blankcorr(mk,blanklabel,"K:"),"Ar40")
k <- getmasses(lk,"Ar39","Ar40") # subset on the relevant isotopes
lca <- fitlogratios(blankcorr(mca,blanklabel,"Ca:"),"Ar37")
ca <- getmasses(lca,c("Ar36","Ar39"),c("Ar37","Ar37")) # subset
ld <- fitlogratios(blankcorr(md))
d <- averagebyday(ld,"DCAL")

# merge all data (except air shots) into one big logratio structure
X <- newredux(concat(list(l,k,ca,d)),Jpos)
data(Melbourne)
if (isTRUE(all.equal(Melbourne$X,X))) {
   print("We just reconstructed the built-in dataset Melbourne$X")}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("concat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("decaycorrection")
### * decaycorrection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: decaycorrection
### Title: Correct for radioactive decay occurred since irradiation
### Aliases: decaycorrection

### ** Examples

data(Melbourne)
C <- calibration(Melbourne$X,"DCAL")
A <- massfractionation(C,Melbourne$fract)
D9 <- decaycorrection(A,Melbourne$irr,"Ar39")
plotcorr(D9)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("decaycorrection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitlogratios")
### * fitlogratios

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitlogratios
### Title: Extrapolation to 'time zero'
### Aliases: fitlogratios fitlogratios.default fitlogratios.timeresolved
###   fitlogratios.PHdata fitlogratios.WiscAr

### ** Examples

samplefile <- system.file("Samples.csv",package="ArArRedux")
m <- loaddata(samplefile) # samples and J-standards
blanklabel <- "EXB#"
l <- fitlogratios(blankcorr(m,blanklabel),"Ar40")
plotcorr(l)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitlogratios", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fractionation")
### * fractionation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fractionation
### Title: Compute the mass fractionation correction
### Aliases: fractionation

### ** Examples

data(Melbourne)
fd37file <- system.file("AirL2.csv",package="ArArRedux")
fd40file <- system.file("AirH1.csv",package="ArArRedux")
fract <- list(fractionation(fd37file,"L2"),
              fractionation(fd40file,"H1"))
if (isTRUE(all.equal(Melbourne$fract,fract))){
  print("We just re-created the fractionation correction for the Melbourne dataset")
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fractionation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get4039")
### * get4039

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get4039
### Title: Calculate the 40Ar*/39ArK-ratios
### Aliases: get4039

### ** Examples

data(Melbourne)
R <- get4039(Melbourne$X,Melbourne$irr)
plotcorr(R)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get4039", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getJfactors")
### * getJfactors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getJfactors
### Title: Calculate the irradiation parameter ('J factor')
### Aliases: getJfactors

### ** Examples

data(Melbourne)
R <- get4039(Melbourne$X,Melbourne$irr)
J <- getJfactors(R)
plotcorr(J)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getJfactors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getages")
### * getages

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getages
### Title: Calculate 40Ar/39Ar ages
### Aliases: getages

### ** Examples

data(Melbourne)
R <- get4039(Melbourne$X,Melbourne$irr)
J <- getJfactors(R)
ages <- getages(J)
plotcorr(ages)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getages", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getmasses")
### * getmasses

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getmasses
### Title: Select a subset of isotopes from a dataset
### Aliases: getmasses getmasses.default getmasses.timeresolved
###   getmasses.WiscAr getmasses.logratios getmasses.redux

### ** Examples

kfile <- system.file("K-glass.csv",package="ArArRedux")
mk <- loaddata(kfile)
lk <- fitlogratios(blankcorr(mk,"EXB#","K:"),"Ar40")
k <- getmasses(lk,"Ar39","Ar40") # subset of the relevant isotopes
plotcorr(k)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getmasses", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("interference")
### * interference

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: interference
### Title: define the interference corrections
### Aliases: interference

### ** Examples

samplefile <- system.file("Samples.csv",package="ArArRedux")
irrfile <- system.file("irradiations.csv",package="ArArRedux")
X <- read(samplefile,blabel="EXB#",Jpos=c(3,15))
irr <- loadirradiations(irrfile)
# assume log(36Ar/37Ar) = log(39Ar/37Ar) = 1 in co-irradiate Ca-salt
# with variances of 0.0001 and zero covariances
ca <- interference(intercepts=c(1,1),
                   covmat=matrix(c(0.001,0,0,0.001),nrow=2),
                   num=c("Ar39","Ar36"),den=c("Ar37","Ar37"),
                   irr=X$irr[1],label="Ca-salt")
# assume log(39Ar/40Ar) = 4.637788 in co-irradiate K-glass
# with variance 7.9817e-4
k <- interference(intercepts=4.637788,covmat=7.9817e-4,
                  num="Ar39",den="Ar40",irr=X$irr[1],
                  label="K-glass")
ages <- process(X,irr,ca=ca,k=k)
summary(ages)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("interference", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("isoratios")
### * isoratios

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: isoratios
### Title: Isochron ratios
### Aliases: isoratios

### ** Examples

data(Melbourne)
IR <- isoratios(Melbourne$X,irr=Melbourne$irr,fract=Melbourne$fract)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("isoratios", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("loaddata")
### * loaddata

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: loaddata
### Title: Load mass spectrometer data
### Aliases: loaddata

### ** Examples

samplefile <- system.file("Samples.csv",package="ArArRedux")
masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
m <- loaddata(samplefile,MS='ARGUS-VI') # samples and J-standards
plot(m,"MD2-1a","Ar40")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("loaddata", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("loadirradiations")
### * loadirradiations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: loadirradiations
### Title: Load the irradiation schedule
### Aliases: loadirradiations

### ** Examples

irrfile <- system.file("irradiations.csv",package="ArArRedux")
irr <- loadirradiations(irrfile)
str(irr)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("loadirradiations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("massfractionation")
### * massfractionation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: massfractionation
### Title: Apply the mass fractionation correction
### Aliases: massfractionation

### ** Examples

graphics.off()
data(Melbourne)
C <- calibration(Melbourne$X,"DCAL")
A <- massfractionation(C,Melbourne$fract)
 plotcorr(A)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("massfractionation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("param")
### * param

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: param
### Title: Set or get Ar-Ar_Redux parameters
### Aliases: param

### ** Examples

data(Melbourne)
param(Melbourne$X)$air
Y <- param(Melbourne$X,air=295.5)
param(Y)$air



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("param", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot")
### * plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.timeresolved
### Title: Plot a time resolved mass spectrometry signal
### Aliases: plot.timeresolved plot.WiscAr plot.PHdata

### ** Examples

samplefile <- system.file("Samples.csv",package="ArArRedux")
mMC <- loaddata(samplefile)
plot(mMC,"MD2-1a")
mPH <- loaddata(samplefile)
plot(mPH,"MD2-1a","Ar40")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotcorr")
### * plotcorr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotcorr
### Title: Plot a matrix with correlation coefficients
### Aliases: plotcorr

### ** Examples

graphics.off()
data(Melbourne)
plotcorr(Melbourne$X)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotcorr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("process")
### * process

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: process
### Title: Process logratio data and calculate 40Ar/39Ar ages
### Aliases: process

### ** Examples

data(Melbourne)
ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
summary(ages)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("process", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read")
### * read

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read
### Title: Read mass spectrometer data
### Aliases: read

### ** Examples

samplefile <- system.file("Samples.csv",package="ArArRedux")
kfile <- system.file("K-glass.csv",package="ArArRedux")
cafile <- system.file("Ca-salt.csv",package="ArArRedux")
dfile <- system.file("Calibration.csv",package="ArArRedux")
X <- read(samplefile,blabel="EXB#",Jpos=c(3,15),
          kdat=kfile,cadat=cafile,ddat=dfile,MS='ARGUS-VI')
plotcorr(X)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("subset")
### * subset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: subset.timeresolved
### Title: Select a subset of some data
### Aliases: subset.timeresolved subset.logratios subset.redux
###   subset.results

### ** Examples

data(Melbourne)
ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
MD <- subset(ages,labels=c("MD2-1","MD2-2","MD2-3","MD2-4","MD2-5"))
plotcorr(MD)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("subset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary")
### * summary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.results
### Title: Summary table
### Aliases: summary.results

### ** Examples

data(Melbourne)
ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
summary(ages)[1:5,]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("weightedmean")
### * weightedmean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: weightedmean
### Title: Calculate the weighted mean age
### Aliases: weightedmean

### ** Examples

data(Melbourne)
ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
weightedmean(ages,"MD2-")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("weightedmean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
