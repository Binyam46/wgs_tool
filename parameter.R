GENO <-"fam_60.gen"
COM <- "Univariate"
PHENO <- c("pheno.dat","ID", "FamID","sire","dam","trial","Sex","pen","don",
	"STW","Wegight","Res","Tol")
ntrait<-2
VAR <-c(1000,3,1999,1000,2,1800) # list all variances here accordingly. 
  #the first variace is for the residual and the order must be according 
  #to the random effects given in the "RAN" and the order of the traits
RAN<-c(1) 	# location of random effects
FIX<-c(6,7,8) 		# location of fixed effects
#COV<- c(1,2)		 # location of covariable effects
TRLOC <-c(9,11)	# location of the traits
REL<-c("GIN")
MAF <-0.05
FREQMETHOD <- 1
DIAG <-0.02
GMATOUT<-0
IGMATOUT<-1
SCALEMETHOD<-1
