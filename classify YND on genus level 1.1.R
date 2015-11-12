#classify YND using AAMM on genus level
#first load YND data
library(Biostrings)
library(flexmix)
library(graph)
library(seqinr)
library(QuasiAlign)
db <- openGenDB("F:/green.db")
ge <- sort(table(getRank(db, rank="genus",all=TRUE)),decreasing=TRUE) #get all genus name and sort by genus size
genus.name <- rownames(ge)

ff <- function(){
x<-readDNAStringSet("F:/Gene data/gene16s YND.wri")
name.x <- names(x)
loc.superkingdom <- vector()
loc.phylum <- vector()
loc.class <- vector()
loc.order<- vector()
loc.family <- vector()
loc.genus <- vector()
loc.species <- vector()

for(i in 1:6026){
  a <- which(strsplit(name.x[i],":")[[1]]==" superkingdom")
  loc.superkingdom[i] <- strsplit(name.x[i],": ")[[1]][a+1]
  loc.phylum[i] <- strsplit(name.x[i],": ")[[1]][a+3]
  loc.class[i] <- strsplit(name.x[i],": ")[[1]][a+5]
  loc.order[i] <- strsplit(name.x[i],": ")[[1]][a+7]
  loc.family[i] <- strsplit(name.x[i],": ")[[1]][a+9]
  loc.genus[i] <- strsplit(name.x[i],": ")[[1]][a+11]
  loc.species[i] <- strsplit(name.x[i],": ")[[1]][a+13]}

return(loc.genus)}
loc.genus <- ff()
x<-readDNAStringSet("F:/Gene data/gene16s YND.wri")
#====================================================================
id <-match(names(sort(table(loc.genus),decreasing=TRUE)[1:10]),genus.name)
#decide test genus
test.genus.nam <- names(ge[id])
test.genus <- table(loc.genus[which(loc.genus %in% test.genus.nam)])
test.genus.seq <- tolower(x[which(loc.genus %in% test.genus.nam)])
correct.class <- loc.genus[which(loc.genus %in% test.genus.nam)]
#load target AAMMs for supervised learning
#decide the target AAMMs 
target.genus.nam <- names(ge[id])
M <- 2:16
#============================================================================================================================
set.seed(12345)
seed <- c(83512, 49648, 66009, 80846, 34734,sample(100000,195))

score.matrix.ynd <- matrix(nrow=length(test.genus.seq),ncol=length(id))
mvalue <- c()
count <- 1
for(i in id){
  cat("Model creation for",genus.name[i],"\n")
  nam.AH    <- paste(paste0(strsplit(genus.name[i],"")[[1]],collapse=""),".AH",sep="")
  assign("a", getSequences(db, rank = "genus", name = genus.name[i]))
  set.seed(seed[i-2])
  position <- sample(ge[i],2*ceiling(0.1*ge[i]),replace=FALSE)
  assign("b",ab.construct(k=2,a[position]))
  
  s <- c()
  cat("M-cut selection for model:",genus.name[i],"\n")
  for(m in 1:length(M)){
    cat("Score calculation against model using M-cut:",M[m],"\n")
    s <- cbind(s,Score.seq(tolower(a[position]),a[position],b,mcut=Mcut(m,b)))
    }
  mvalue[count] <- mm <- M[which.max(apply(s,2,mean)[-1])]
  score.matrix.ynd[,count] <- Score.seq(tolower(test.genus.seq),a[position],b,mcut=Mcut(mm,b))
  count <- count +1
  assign(nam.AH,b)
}

class.rate <-c()
for(i in 1:length(id)){
rw <- match(which(loc.genus==test.genus.nam[i]),which(loc.genus %in% test.genus.nam))
class.rate[i] <- sum(target.genus.nam[apply(score.matrix.ynd[rw,],1,which.max)]==correct.class[rw])/length(correct.class[rw])
}

#=========================================================================
for(i in id){
  s <- c()
  cat("M-cut selection for model:",genus.name[i],"\n")
  assign("a", getSequences(db, rank = "genus", name = genus.name[i]))
  set.seed(seed[i-2])
  position <- sample(ge[i],2*ceiling(0.1*ge[i]),replace=FALSE)
  assign("b",get(paste0(genus.name[i],".AH")))
  for(m in 1:length(M)){
    s <- cbind(s,Score.seq(tolower(a[position]),a[position],b,mcut=Mcut(m,b)))
  }
  mm <- M[which.max(apply(s,2,mean)[-1])]
 print(mm)}
#=========================================================================
#Classification
#build score matrix
#using which target genus

plot(1:16,class.rate3,type='o',xlab="M-cut",ylab="Classification Rate",main="Classification Rate using AAMM"
     ,pch=16,lwd=2,col='red',cex=0.8,axes=FALSE,lty=3)
axis(1,seq(1,16,3),seq(30,45,3))
axis(2,seq(0,1.0,0.2))
points(class.rate2,col='blue',type='o',pch=16,lwd=2,cex=0.8,lty=3)
points(class.rate1,col='green',type='o',pch=16,lwd=2,cex=0.8,lty=3)
points(class.rate,col='black',type='o',pch=16,lwd=2,cex=0.8,lty=1)
legend(10,0.8,legend=c("Corynebacterium","Staphylococcus","Pseudomonas","Average"),col=c('red','blue','green','black'),
       lty=c(3,3,3,1),box.col='white',cex=1.2,lwd=2)
#===================================================================

