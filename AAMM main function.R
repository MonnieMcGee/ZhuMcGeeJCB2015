# revise the function of seaching parent for each node
library(Biostrings)
library(flexmix)
library(graph)
library(seqinr)
#tolower
# revise function oligonucleotideFrequency1
# to make it calculate the k-grams frequecy of multiple seq.datauences

oligonucleotideFrequency1 <- function(x,w){
freq <- oligonucleotideFrequency(DNAString(paste0(x,collapse="")),width=w)
  return(freq)}

#Calculate the context of all k-grams with respect to D
ContestDS <- function(seq.data,k){ #seq.data
  freq <- oligonucleotideFrequency1(seq.data,w=k+1)
  contest.matrix <- matrix(rep(NA,4^k*4),nrow=4)
  for(i in 1:4){
    for(j in 1:(4^k)){
      contest.matrix[i,j] <- (1+sum(freq[4*j-(4-i)]))/(4+sum(freq[(4*j-3):(4*j)]))
    }
  }
  rownames(contest.matrix) <- c("A","C","G","T")
  colnames(contest.matrix) <- names(oligonucleotideFrequency(DNAString(""),width=k))
  return(contest.matrix)}

#----------------------------------------
#Calculate the context of an abstraction aj
# abstraction = kgrams
ContestDA <- function(abstraction,k,seq.data){ #input abstraction = mcut
  abstraction1 <- unlist(strsplit(abstraction,c(",")))
  contest.ds <- ContestDS(seq.data,k)
  freq1 <- oligonucleotideFrequency1(seq.data,w=k)
  weight <- freq1/sum(freq1)
  kgrams <- names(oligonucleotideFrequency(DNAString(""),width=k))
  id <- match(abstraction1,kgrams)
  new.weight <- weight[id]/sum(weight[id])
  contest.da <- contest.ds[,id]%*%t(t(new.weight))
  rownames(contest.da) <- c("A","C","G","T")
  result <- list(contest.da,weight);names(result) <- c("Contest of Abstractions","Prior Probability")
  return(result)
}
#-------------------------------
#pi1,pi2 are prior probability of ai and aj
#a1 and a2 
#weighted Jensen-Shannon divergence between two probability distributions
WeightJSdis <- function(pi1,pi2,p1,p2){
  pbar <- (pi1/(pi1+pi2))*p1+(pi2/(pi1+pi2))*p2
  js <- (pi1/(pi1+pi2))*KLdiv(cbind(p1,pbar))+(pi2/(pi1+pi2))*KLdiv(cbind(p2,pbar))
  return(js[1,2])
}
#----------------------------------------
#using R create stack
#A stack implementation consists of three main components:
#a container variable --- a.k.a. the stack
#a push method to add elements
#a pop method to remove elements
#Last in first out stack
push <- function(x, value, ...) UseMethod("push")
pop  <- function(x, ...) UseMethod("pop")
push.stack <- function(x, value, ...) x$push(value)
pop.stack  <- function(x) x$pop()

new_stack <- function() { 
  stack <- new.env()
  stack$.Data <- vector()
  stack$push <- function(x) .Data <<- c(.Data,x)
  stack$pop  <- function() {
    tmp <- .Data[1]
    .Data <<- .Data[-1]
    return(tmp)
  }
  environment(stack$push) <- as.environment(stack)
  environment(stack$pop) <- as.environment(stack)
  class(stack) <- "stack"
  stack
}

#---------------------------------------
#Abstraction Construction
ab.construct <- function(k,seq.data)
{
  kgrams <- A <- names(oligonucleotideFrequency(DNAString(""),width=k))
  Tree <- new_stack()
  push(Tree,kgrams)
  
  w <- oligonucleotideFrequency1(seq.data,w=k)
  weight <- w/sum(w)
  
  for (ite in 1:(length(kgrams)-1)){
    cat("Abstraction Construction Iteration #",ite,"\n")
    n <- length(A)
    d.matrix <- matrix(rep(NA,n^2),nrow=n)
    for( num.au in 1:(n-1)){
      for(num.av in n:(num.au+1)){
        d.matrix[num.au,num.av] <- sum(weight[c(num.au,num.av)])*WeightJSdis(weight[num.au],weight[num.av],ContestDA(A[num.au],k,seq.data)[[1]],ContestDA(A[num.av],k,seq.data)[[1]])
      }
    }
    rownames(d.matrix) <- A;colnames(d.matrix) <- A
    argmin<- which(d.matrix==min(d.matrix,na.rm=T),arr.ind=TRUE)
    au <- argmin[1];av <- argmin[2]
    aw <- paste0(c(A[au],",",A[av]),collapse="")
    A <- c(A[-c(au,av)],aw);add.weight <- sum(weight[c(au,av)])
    weight <- c(weight[-c(au,av)],add.weight);names(weight)[length(weight)] <- aw
    push(Tree,aw)
  }
  Tree
}
#===========================
#find parent of each node in the Tree
  find.parent <- function(Tree){
  #find k by the size of tree  
  k <- length(s2c(Tree$.Data[1]))
  N <- 4^k
  spl <- function(x){unlist(strsplit(x,c(",")))}
  #find the parent of each node
  parent <- vector("list", length=2*N-1);names(parent) <- as.character(1:(2*N-1))
  for(i in c(1:(2*N-2))){
    if(i <= 16) j=17 else j=i+1
    while(is.null(parent[[i]])){
      if(sum(spl(Tree$.Data[i]) %in% spl(Tree$.Data[j]))!=0) parent[[i]] <- as.character(j)
      j <- j+1}}
  return(parent)}
#==============================
#find m-cut
Mcut <- function(m=3,Tree,ab.name=TRUE){
  tree.parent <- find.parent(Tree)
  l <- length(tree.parent)+1
  mcut <- c()
  for(i in 1:(l-2)){ if(tree.parent[[i]] %in% as.character((l-m+1):(l-1))) mcut <- c(mcut,i)}
  if (ab.name == TRUE) return(Tree$.Data[mcut[1:m]]) else return(mcut[1:m])}
#======================================
# plot Abstraction Hierarchy T
plot.AH <- function(Tree){
  k <- length(s2c(Tree$.Data[1]))
  N <- 4^k
  V <- as.character(1:(2*N-1))
  edL <- find.parent(Tree)
  gR <- graphNEL(nodes=V, edgeL=edL,edgemode='directed')
  return(plot(gR))
}
#=====================================
#Score a sequence
Score.seq <- function(new.seq,data,Tree,mcut,log=TRUE){
  k <- length(s2c(Tree$.Data[1]))
  kgrams <- names(oligonucleotideFrequency(DNAString(""),width=k))
  
  m <- match(mcut,Tree$.Data)        
  mcut.group <- c()
  contestda.matrix <- matrix(nrow=4,ncol=length(m))
  for(i in 1:length(m)) {mcut.group[which(kgrams %in% strsplit(mcut[i], ",")[[1]])] <- i
                         a <- ContestDA(mcut[i],k=k,data)
                         contestda.matrix[,i] <- a[[1]];weight <-a[[2]]}
  
  score.v <- matrix(ncol=1,nrow=length(new.seq))
  for(l in 1:length(new.seq)){
    seq1 <- s2c(new.seq[l]);seq.num1 <- match(seq1,c("a","c","g","t"))
    na.id <- which(is.na(seq.num1));seq1[na.id] <- sample(c("a","c","g","t"),length(na.id),replace=TRUE)
    seq.num1 <- match(seq1,c("a","c","g","t"))
    seq.kgram <- c()
    for(i in 1:(length(seq1)-(k-1))) {seq.kgram[i] <-paste0(seq1[i:(i+k-1)],collapse="")}
    seq.num <- match(toupper(seq.kgram),kgrams)
    
    ln.prob <- log(weight[seq.num[1]])
    for(j in (k+1):length(seq1)){
      ln.prob <- ln.prob+log(contestda.matrix[seq.num1[j],mcut.group[seq.num[j-2]]])
    }
    score.v[l,1] <- ln.prob
  }
  
  
  if(log==TRUE) return(score.v) else return(exp(score))
}