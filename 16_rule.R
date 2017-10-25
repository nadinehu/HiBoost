set.seed(48)
library(Matrix)
library(MASS)
source("../boosting_functions.R")
original=TRUE
STOP=FALSE
FLIP=TRUE
text="rand.txt" if STOP = true, give it the input file as text, and if STOP=FALSE give it the output file  to write the iterations in as text 
# probabilty with which the regs are turned on
regruleon=0.5
regruleon2=0.2
regruleon3=0.02
# probabilty with which the motifs are turned on
motruleon=0.5
motruleon2=0.05

niter_MAX=50
niter=50
ndatasets=100 

# model params
nreg=45
nmotif=22
ngenes=600
ncelltypes=16
intnodes=ncelltypes-1

errors=matrix(0,2,niter)

# generate cell types hierarchy 
# f leaves
final_leaves=matrix(0,intnodes,ncelltypes)
final_leaves[1,]=c(1:16); final_leaves[2,c(1:9)]=c(1:9); final_leaves[3,1:7]=c(10:16); final_leaves[4,c(1:8)]=c(2:9)
final_leaves[5,1:3]=c(10:12); final_leaves[6,c(1:4)]=c(13:16); final_leaves[7,1:4]=c(2:5); final_leaves[8,1:4]=c(6:9)
final_leaves[9,1:2]=c(11,12); final_leaves[10,1:2]=c(13,14); final_leaves[11,1:2]=c(15,16); final_leaves[12,1:2]=c(2,3)
final_leaves[13,1:2]=c(4,5); final_leaves[14,1:2]=c(6,7); final_leaves[15,1:2]=c(8,9)

# d children
alf=-1000 # for all leaves
direct_children=matrix(0,intnodes,3) # bin tree here - last indicates the leaves
direct_children[1,1:2]=c(2,3); direct_children[2,]=c(4,1,-1); direct_children[3,1:2]=c(5,6); direct_children[4,1:2]=c(7,8);
direct_children[5,]=c(9,10,-1); direct_children[6,1:2]=c(10,11); direct_children[7,1:2]=c(12,13); direct_children[8,1:2]=c(14,15);
direct_children[9,]=c(11,12,alf); direct_children[10,]=c(13,14,alf); direct_children[11,]=c(15,16,alf); 
direct_children[12,]=c(2,3,alf);direct_children[13,]=c(4,5,alf); direct_children[14,]=c(6,7,alf);direct_children[15,]=c(8,9,alf);



nrules=16 
rules=matrix(0,nrules,5) # motif, reg, reg sign, internal node, parent in ADT
rules[,3]=sample(0:1,nrules,replace=T) # random sign assignment. 
rules[which(rules[,3]==0),3]=-1

# # associate rules to internal nodes in hierarchy 
rules[1,-3]=c(1,1,1,0); rules[2,-3]=c(2,2,1,0); rules[3,-3]=c(6,2,1,0); rules[4,-3]=c(7,6,1,3);
rules[5,-3]=c(8,7,1,4); rules[6,-3]=c(3,3,2,0); rules[7,-3]=c(5,4,2,1); rules[8,-3]=c(4,4,2,1);
rules[9,-3]=c(4,5,2,7); rules[10,-3]=c(10,10,3,0);
rules[11,-3]=c(11,11,6,10); rules[12,-3]=c(11,12,6,10);
rules[13,-3]=c(11,13,6,10); rules[14,-3]=c(12,11,6,10); rules[15,-3]=c(12,12,6,10); rules[16,-3]=c(12,13,6,10);

# generate/set random alpha for each rule. 
alphas=runif(nrules, -10, 10)
alphas[1]=-3.5; alphas[2]=4; alphas[3]=-1;
alphas[3]=3 # maybe try 3 later on. 
alphas[4]=-3; alphas[5]=-5; alphas[6]=-2; #alphas[6]=-3
alphas[7:14]=c(-20,17,-10,-13,-13,20,18,-15)

recov=matrix(0,ndatasets,7)
iteration_stop=rep(0,ndatasets)
if(STOP==TRUE)
  iteration_stop=scan(file=text)
iteration_stop=as.vector(iteration_stop)

for(l in 1:ndatasets){
  cat("\n dataset ", l)
  mot_gen=matrix(0, nmotif,ngenes)
  reg_ctype=matrix(0, nreg, ncelltypes)
  
  # creating the regulator matrix 
  for(i in 1:ncelltypes){
    # satisfy most of the relevant rules to this cell type 
    rulestocheck=rep(0,nrules); 
    for(k in 1:nrules){
      if(i %in% final_leaves[rules[k,4],])
        rulestocheck[k]=1
    }
    check=which(rulestocheck==1)
    for(k in 1:length(check)){
      reg_ctype[rules[check[k],2],i]=rules[check[k],3]
      # turning on that relevant regulator for that cell type with proba (0.1,0.9)
      zero=sample(0:1,1,replace=T,prob=c(1-regruleon,regruleon)) 
      if(zero==0) # setting it to zero
        reg_ctype[rules[check[k],2],i]=0
    }
    
    # satisfy other rule regulators with lower proba
    all=c(1:nrules); rest=all[-check]; 
    for(k in 1:length(rest)){
      reg_ctype[rules[rest[k],2],i]=rules[rest[k],3]
      zero=sample(0:1,1,replace=T,prob=c(1-regruleon2,regruleon2)) # turning on a regulator that is relevant in some rules but not that cell type 
      if(zero==0) # setting it to zero
        reg_ctype[rules[rest[k],2],i]=0
    }
    
    # add some random regulators with even lower proba 
    for(k in 1:nreg){
      isreg_inrule= (k %in% rules[,2])
      if(isreg_inrule==FALSE ) # then it is a regulator that is not present in any rule 
        reg_ctype[k,i]= sample(c(-1,0,1),1,replace=T,prob=c(regruleon3,1-regruleon3*2,regruleon3))
    }
    
  }
  
  # creating the gene by motif matrix 
  for(j in 1:ngenes){
    # make it satisfy motifs in all rules with a lower probablity - so that using the hierarchy gives 
    # an advantage
    all=c(1:nrules); 
    for(k in 1:length(all)){
      mot_gen[rules[all[k],1],j]=sample(0:1,1,replace=T,prob=c(1-motruleon,motruleon))# prob=c(0.1,0.9)
    }
    
    # also turn on some that arent in any rules at all 
    for(k in 1:nmotif){
     ismot_inrule=(k %in% rules[,1])
     if(ismot_inrule==FALSE ) # then it is a regulator that is not present in any rule 
        mot_gen[k,j]= sample(0:1,1,replace=T,prob=c(1-motruleon2,motruleon2))
    }
  }
  
  # dists
  cat("\n motifs per gene")
  print(table(colSums(mot_gen)))
  
  cat("\n genes per motif")
  print(table(rowSums(mot_gen)))
  
  cat("\n cell types per reg")
  print(table(rowSums(abs(reg_ctype))) )
  
  cat("\n regs per cell types")
  print(table(colSums(abs(reg_ctype))) )
  
  # generate labels by following the ADT structure 
  labels=matrix(0,ngenes,ncelltypes)
  for(i in 1:ncelltypes){
    rulestocheck=rep(0,nrules); 
    for(k in 1:nrules){
      if(i %in% final_leaves[rules[k,4],])
        rulestocheck[k]=1
    }
    check=which(rulestocheck==1)
    # cat("\n print relevant rules numbers for obsevations in col ", i); print(check )
    
    for(j in 1:ngenes){
      count=0
      val=0 # alpha to start for root.
      rpassed=c(0)
      for(k in 1:length(check)){ # only want to check rules that are relevant to this cell type. + check that parent also passes rule
        #  cat("\n checking that rule ", check[k]); cat(" 's parent ",  rules[check[k],5]) ;cat(" has been reached ")
        if( mot_gen[rules[check[k],1],j]== 1 && reg_ctype[rules[check[k],2],i]== rules[check[k],3] && ( rules[check[k],5] %in% rpassed )  )
          # and need to check passes the rules before it in the hierarchy otherwise wrong application of ADT. 
        { val=val+alphas[check[k]]; count=count+1 ; rpassed=c(rpassed,check[k]);  }
      }
      labels[j,i]=sign(val)
      if(labels[j,i]!=0)
      {
        # cat("\n for gene ", j); cat(" and cell type ", i); cat(" passes rules ", count)
        # cat("\n val ", val)
        # cat("\n passed rules number "); print(rpassed); #cat("\n rules are ") ; print(rules[rpassed[-1],]);
        
      }
    }
  }
  
  cat("\n labels distribution \n")
  print(table(labels))
  
  # or add noise to the observations. 
  
  # add noise to test - flip some labels -1 & 1. 
  if (FLIP){
     noise=matrix(runif(ngenes*ncelltypes), ncol=ncelltypes) 
     noise=which(noise>0.9,arr.ind=TRUE)
     noise=sparseMatrix(i=noise[,1], j=noise[,2],x=1, dims = c(ngenes,ncelltypes))
     pos=which(noise>0&labels>0,arr.ind=TRUE)
     neg=which(noise>0&labels<0,arr.ind=TRUE)
     labels_noisy=labels 
     labels_noisy[pos]=-1; labels_noisy[neg]=1 # 10perc of noise obs flipped. 
  }
  
  #write.table(final_leaves, file="~/Dropbox/reg_networks/code_boosting/sim/fleaves_sim.txt", row.names=FALSE, col.names=FALSE)
  #write.table(direct_children, file="~/Dropbox/reg_networks/code_boosting/sim/dchildren_sim.txt", row.names=FALSE, col.names=FALSE)
  #write.table(labels, file="~/Dropbox/reg_networks/code_boosting/small_5rule_sim/labels_sim.txt", row.names=FALSE, col.names=FALSE)
  #write.table(mot_gen, file="~/Dropbox/reg_networks/code_boosting/small_5rule_sim/gen_mot_sim.txt", row.names=FALSE, col.names=FALSE)
  #write.table(reg_ctype, file="~/Dropbox/reg_networks/code_boosting/small_5rule_sim/reg_ctype_sim.txt", row.names=FALSE, col.names=FALSE)
  
  if(STOP==TRUE)
    niter=iteration_stop[l]
  
  boosted=boost(as.matrix(final_leaves), as.matrix(direct_children), intnodes,
                as.matrix(mot_gen), as.matrix(t(reg_ctype)), as.matrix(labels), niter , rules, original)
  cat("\n boosting training loss ", boosted$trl[niter]); cat(" test loss ", boosted$tsl[niter])
  recov[l,]=boosted$recovery_Eval
  ll=length(c( boosted$trl))
  errors[1,1:ll]=c( boosted$trl)+ errors[1,1:ll]; errors[2,1:ll]=c( boosted$tsl)+errors[2,1:ll]
  
  iteration_stop[l]=boosted$iteration_stop
}

if(STOP==FALSE)
  write.matrix(iteration_stop,text)

cat("\n recovery over datasets \n")
cat("\n rules, correct rules, wrong rules, correct with int node, correct with parent \n")
print(recov)
cat("\n summing the recovery values \n ")
print(colSums(recov)/ndatasets)
plot(c(1:niter_MAX),errors[1,]/ndatasets)
lines(c(1:niter_MAX),errors[1,]/ndatasets)
lines(c(1:niter_MAX),errors[2,]/ndatasets,col="blue")
title(main = list("Average Training and Test errors", cex = 1.5, col = "red", font = 3))

