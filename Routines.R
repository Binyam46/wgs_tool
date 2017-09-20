#!/usr/bin/env Rscript

ReadFiles = function (phenoFile, genoFile, nfam) {

##########################################
## Read genotype and phenotype files #####
##########################################

ped=read.table(phenoFile,header=FALSE, sep=" ")
gen=read.table(genoFile,header=FALSE)

fam<-nfam
if (fam > 0) {
	ped=ped[ped[,2]==fam,]	
	gen<-gen[c(ped[,3],ped[,4],ped[,1]),]
}


gen[,1]<-c(1,2, gen[3:nrow(gen),1])

all_id<-gen[,1]
nsnp=(ncol(gen)-1)/2

ped<-data.frame(ped[,1], ped)

ped[,2]<-match(ped[,2], unique(sort(all_id)))
gen[,1]<-match(gen[,1], unique(sort(all_id)))

ko = order(gen[,1])
gen1<-gen[ko,]

aa<-data.frame(c(1:nsnp), c(1:nsnp))
write.table(aa, file="map.dat", quote=FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

write.table(gen1,"rep.mrk", quote=F,row.names=F,col.names=F)
write.table(ped,"pheno.dat", quote=F,row.names=F,col.names=F)

}


htd = function (n) {
i=n
ped=read.table("pheno.dat",header=FALSE, sep=" ")
testIndexes <- which(ped[,14]==i,arr.ind=TRUE)
candi<-ped[testIndexes, c(1,2,3,10:14)]
write.table(candi,file="candidates.dat", col.names=FALSE, row.names=FALSE, sep="\t")

rec=ped[-testIndexes,-1]

N=nrow(rec)
k= 3     # number of traits                     
K=c(1:k)
l=c(9,11,12)  # start of wave number

trait_no=rep(1:k, each=N)     # trait number
obs1=matrix(0, nrow=N*k, ncol=1)

animal_ID=rep(rec[,1],k)     # multiple animal ID
fam=rep(rec[,2],k)     # multiple family effect
don=rep(rec[,8],k)     # multiple doner family effect
pen=rep(rec[,7],k)     # multiple pin effect
trial=rep(rec[,6],k)     # multiple trial effect
sex=rep(rec[,5],k)     # multiple sex effect

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
obs=data.frame(animal_ID,fam, don, pen, trial, sex,trait_no,obs1)
obs <- obs[!(obs$obs1 > 900.0), ]

# order the data frame first by animal number and then by trait number within animals

obs2=obs[do.call(order,obs),]                        
dat_womb=data.frame(obs2$trait_no,obs2$animal_ID, obs2$fam, obs2$don, obs2$pen, obs2$trial, obs2$sex,obs2[,8])

write.table(dat_womb,file="pop_100.d", col.names=FALSE, row.names=FALSE, sep="\t")
}









