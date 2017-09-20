#!/usr/bin/env Rscript

args  <- commandArgs(TRUE)
a<-args[1]

source("parameter.R")
source("Routines.R")

######################################################################
## Read parameter fille and check mandatory commands
######################################################################

if (!exists("ntrait")){
	cat('ntrait should be specified\n')
	quit(save='no') 
} 
if (!exists("PHENO")){
	cat('The PHENO parameter is missing. Phenotype data should be provided\n')
	quit(save='no')
} 
if (!exists("GENO")){
	cat('The GENO parameter is missing. Genotype data should be provided\n')
	quit(save='no')
} 
if (!exists("RAN")){
	cat('The RAN parameter is missing. The column for random effects shoud be provided\n')
	quit(save='no')
} 

if (!exists("VAR")){
	cat('The VAR parameter is missing. Variance components should be provided\n')
	quit(save='no')
}
if (!exists("REL")){
	warning('REL is not specfied for any of the random effects. Independence will be assumed\n')
}  

nfix=0
if (exists("FIX")){
	nfix <- length(FIX)
}

ncov=0
if (exists("COV")){
	ncov <- length(COV)
}

nran <- length(RAN)
mvar<-((ntrait*(ntrait+1))/2)
nvar<-(nran+1)*mvar

if (length(VAR)< nvar) {
	cat('The variance components provided in VAR is not as it should be\n')
	quit(save='no')
} 

if (length(TRLOC)!= ntrait) {
	cat('The number of traits in TRLOC should be equal to ntrait\n')
	quit(save='no')
} 

################################################################
###		Calling other routines                       ###
################################################################

ReadFiles(PHENO[1], GENO[1], a)


######################################################################
### 	Prepare input file for wombat				   ###
######################################################################

ped=read.table("pheno.dat",header=FALSE, sep=" ")
nfam<-length(unique(ped[,2]))

rec=ped[,-1]

N=nrow(rec)
k= ntrait     # number of traits                     
K=c(1:k)
l=TRLOC  # start of wave number

trait_no=rep(1:k, each=N)     # trait number
obs1=matrix(0, nrow=N*k, ncol=1)
raneff<- array(0,dim=c(N*k,nran))

for( i in 1:nran){
	raneff[,i]=rep(rec[,RAN[i]],k)     
}


if(ncov > 0){
	coveff<- array(0,dim=c(N*k,ncov))
	for( i in 1:ncov){
		coveff[,i]=rep(rec[,COV[i]],k)    
	}
}

if(nfix > 0){
	fixeff<- array(0,dim=c(N*k,nfix))
	for( i in 1:nfix){
		fixeff[,i]=rep(rec[,FIX[i]],k)    
	}
}

if(ncov > 1 && nfix > 1) {
	FFeff <- data.frame(raneff, coveff, fixeff)
} else if (ncov > 1 && nfix == 0) {
	FFeff <- data.frame(raneff, coveff)
} else if (ncov == 0 && nfix > 0) {
	FFeff <- data.frame(raneff, fixeff)
} else if (ncov == 0 && nfix == 0) {
	FFeff <- raneff
}

S=rep(0,k)     # initial vector for start
E=S             # initial vector to end

# loop to creat strat location  while extracting wave numbers
  for (i in 1:k) {
    S[i]=(i*N)-N +1
    }

# loop to creat end location while extracting wave numbers
  for ( i in 1:k)  {
    E[i]=i*N
  }

# loop to reorder bservatoins in  WOMBAT format
  for (i in K) {
    j=l[i]
    obs1[S[i]:E[i]]=rec[,j]
      }

# combine all information needed in one date frame
obs=data.frame(trait_no, FFeff, obs1)
obs <- obs[!(obs$obs1 > 900.0), ]

# order the data frame first by animal number and then by trait number within animals

k0=order(obs[,2])
obs2=obs[k0,]                      

write.table(obs2,file="pop_100.d", col.names=FALSE, row.names=FALSE, sep="\t")

######################################################################
### Write out wombat parameter file				   ###
######################################################################

if (exists("COM")){
	write (paste("COM", COM, sep= " "),file="wf.par", append = F)
} else {
	write ("COM Blup Analysis",file="wf.par", append = F)
}
	
if ( ntrait == 1){
	write ("ANA UNI",file="wf.par", append = T)
} else {
	write (paste("ANA MUV", ntrait, sep=" "),file="wf.par", append = T)
}

if (exists("PED")){
 write ("PED pedigree.ped",file="wf.par", append = T)
}
 
write (paste("DATA", "pop_100.d", "GRP", sep=" "), file="wf.par", append = T)
cat (c(" TRN", 1:ntrait), file="wf.par", append = T,"\n") 

if ( ntrait > 1){
	write (paste(" traitno", ntrait,sep=" "), file="wf.par", append = T)
}

for( i in 1:nran){
	write (paste("",PHENO[RAN[i]+1], length(unique(ped[,RAN[i]]))+10, sep=" "), file="wf.par", append = T)   
}

if(ncov > 0){
	for( i in 1:ncov){
		coveff[,i]=rep(rec[,COV[i]],k)
	write (paste("",PHENO[COV[i]+1], length(unique(ped[,COV[i]]))+10, sep=" "), file="wf.par", append = T)     
	}
}
if(nfix > 0){
	for( i in 1:nfix){
		fixeff[,i]=rep(rec[,FIX[i]],k)
	write (paste("",PHENO[FIX[i]+1], length(unique(ped[,FIX[i]]))+10, sep=" "), file="wf.par", append = T)     
	}
}
cat (c(" NAM", PHENO[TRLOC+1]), file="wf.par", append = T,"\n") 
write ("END", file="wf.par", append = T)
write (" ",file="wf.par", append = T)

write ("MODEL", file="wf.par", append = T)
if (exists("FIX")){
	for( i in 1:length(FIX)){
		write (paste(" FIX", PHENO[FIX[i]+1],sep=" "), file="wf.par", append = T)
	}
}
if (exists("COV")){
	for( i in 1:length(COV)){
		write (paste(" COV", PHENO[COV[i]+1],"(1)",sep=" "), file="wf.par", append = T)
	}
}
for( i in 1:length(RAN)){
     if (exists("REL")){
	      write (paste(" RAN", PHENO[RAN[i]+1], REL[i], sep=" "), file="wf.par", append = T)
     } else  {
       write (paste(" RAN", PHENO[RAN[i]+1], "IDE", sep=" "), file="wf.par", append = T)
     }
}
for( i in 1:ntrait){
	write (paste(" TR", PHENO[TRLOC[i]+1], i, sep=" "), file="wf.par", append = T)
}

write ("END MODEL", file="wf.par", append = T)
write (" ",file="wf.par", append = T)

for( i in 0:length(RAN)){
  if (i == 0){
		write (paste("VAR Residual", ntrait ,sep=" "),file="wf.par", append = T)
    for (j in i+1:mvar) {
            write (VAR[j],file="wf.par", append = T)
    }
	} else {
		write (paste("VAR", PHENO[RAN[i]+1], ntrait ,sep=" "),file="wf.par", append = T)
    for (j in ((mvar*i)+1):(mvar*(i+1))) {
            write (VAR[j],file="wf.par", append = T)
    }
	}
}

###############################################################################
##   Write out the parameter file:for the Gmatrix                          ####
###############################################################################

write ("$MAPFILE",file="par.dat", append = F)
write ('"map.dat"',file="par.dat", append = T)
write ("$MARKERFILE",file="par.dat", append = T)
write ('"rep.mrk"',file="par.dat", append = T)

write ("$PEDIGREEFILE",file="par.dat", append = T)
write ('" "',file="par.dat", append = T)
write ("$IDBASEFILE",file="par.dat", append = T)
write ('" "',file="par.dat", append = T)

write ("$GMATRIXFILE",file="par.dat", append = T)
if (exists("GMATOUT")){    
   write ("gmat.out",file="par.dat", append = T)
}  else {
   write (" ",file="par.dat", append = T)
}
write ("$IGMATRIXFILE",file="par.dat", append = T)
write ('"animal.gin"',file="par.dat", append = T)

write ("$AMATRIXFILE",file="par.dat", append = T)
write ('" "',file="par.dat", append = T)
write (" ",file="par.dat", append = T)

write ("%PARSTART",file="par.dat", append = T)
write (" ",file="par.dat", append = T)

write ("$MINMAF",file="par.dat", append = T) 
if (exists("MAF")){    
   write (MAF,file="par.dat", append = T)
}  else {
   write (0.05,file="par.dat", append = T)
} 

write ("$FREQMETHOD",file="par.dat", append = T) 
if (exists("FREQMETHOD")){    
   write (FREQMETHOD,file="par.dat", append = T)
}  else {
   write (1,file="par.dat", append = T)
}

write ("$SCALEMETHOD",file="par.dat", append = T) 
if (exists("SCALEMETHOD")){    
   write (SCALEMETHOD,file="par.dat", append = T)
}  else {
   write (1,file="par.dat", append = T)
} 

write ("$DIAG_ADD",file="par.dat", append = T) 
if (exists("DIAG")){    
   write (DIAG,file="par.dat", append = T)
}  else {
   write (0.01,file="par.dat", append = T)
}

write ("$G_ADD",file="par.dat", append = T) 
write (0,file="par.dat", append = T)

write ("$DIAG_ONE",file="par.dat", append = T) 
write (0,file="par.dat", append = T)

write ("$AGSCALE",file="par.dat", append = T) 
write (0,file="par.dat", append = T) 

write ("$PROP_A_to_G",file="par.dat", append = T) 
write (0,file="par.dat", append = T) 

write ("$CAL_DET",file="par.dat", append = T) 
write (1,file="par.dat", append = T) 

write ("$OUT_GMATRIX",file="par.dat", append = T)
if (exists("GMATOUT")){    
   write (GMATOUT,file="par.dat", append = T)
}  else {
   write (0,file="par.dat", append = T)
} 

write ("$OUT_IGMATRIX ",file="par.dat", append = T) 
if (exists("IGMATOUT")){    
   write (IGMATOUT,file="par.dat", append = T)
}  else {
   write (1,file="par.dat", append = T)
} 

write ("$OUT_AMATRIX",file="par.dat", append = T) 
write (0,file="par.dat", append = T) 
