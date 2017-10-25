set.seed(13)
library(Matrix)
library(MASS)
source("boosting_allfunctions.R")
original=FALSE
FLIP=TRUE
pos=c(1:500) # for the gene number 
motif_names="motifs.txt"
reg_names="regs.txt"

maxchildren=5
niter=100
ndatasets=100  

# data dims 
nreg=32
nmotif=22
ngenes=500
ncelltypes=7
intnodes=3
errors=matrix(0,2,niter)
recov=matrix(0,ndatasets,10) # was 7, added total chosen with an int node, ones with a correct motif there, and a correct reg there separately, without checking that the matching is correct. 
rules=matrix(0,32,3) # just for the code but no recovery to compare for now. 
rules[,3]=1 # for now just assume that all rules require a positively expressed regulator 
# now in the fourth entry we specify whether or not the rule is activating or repressing 
rules[1:10,1]=c(1:10); rules[1:10,2]=c(1:10) # these are activating rules with a direct motif
rules[11:20,1]=c(11:20); rules[11:20,2]=c(11:20) # these are activating rules with a direct motif
rules[21:32,1]=c(1,1,1,2,2,2,3,3,3,4,4,4); rules[21:32,2]=c(21:32) # these are new regulators with previous motifs, pathways. 
# have some repetitions because some motifs have multiple regulators binding to them. 

rulesf=matrix(0,40,5) # in 4 put act/deact sign and in last the cell internal node 
regruleon3=0.3 
#flipreg=0.3
motruleon2=0.02
 
# implementing the hierarchy 
direct_children=matrix(0,intnodes, maxchildren)
subtree_nodes=matrix(0,ncelltypes, ncelltypes)
parent_nodes=c(0,1,1,2,2,3,3)

direct_children[1,1:2]=c(2,3)
direct_children[2,1:2]=c(4,5)
direct_children[3,1:2]=c(6,7)

subtree_nodes[1,1:7]=c(1,2,3,4,5,6,7)
subtree_nodes[2,1:3]=c(2,4,5)
subtree_nodes[3,1:3]=c(3,6,7)
subtree_nodes[4,1]=c(4)
subtree_nodes[5,1]=c(5)
subtree_nodes[6,1]=c(6)
subtree_nodes[7,1]=c(7)

# transitions between the cell types in the hierarchy 
transitions=matrix(0,6,2)
transitions[1,]=c(1,2);transitions[3,]=c(2,4);
transitions[2,]=c(1,3);transitions[4,]=c(2,5);
transitions[5,]=c(3,6);transitions[6,]=c(3,7);

rules_per_cell=matrix(0,4*ncelltypes,20) #  have 4 rules per cell type. 
# and each of these rules can at most be the AND of 5 rules (each requiring 4 )
rules_per_cell[1,1:4]=c(1,1,1,1) # one activating rule for initial cell type. 
rules_per_cell[3,1:4]=c(11,11,1,-1) # the deactivating rule in the initial cell type
rulesf[1,]=c(rules_per_cell[1,1:4],1)
rulesf[2,]=c(rules_per_cell[3,1:4],1)
h=2

# sample some rules for the cell types: to add the dependencies as well with and rules 
for(i in 2:ncelltypes)  { # ncelltypes
  cat("\n cell type ", i )
  # sample the direct activating rules which are in 1:10
  rules_per_cell[4*(i)-3,1:4]=c(rules[sample(c(2:10), 1, replace=T, prob=rep(1/9,9) ),],1) # first 4 entries
  h=h+1
  rulesf[h,]=c(rules_per_cell[4*(i)-3,1:4],i)
  
  # sample the activating rules with ANDs from the parent node rules 
  rules_per_cell[4*(i)-2,1:4]=c(rules[sample(c(21:26), 1, replace=T, prob=rep(1/6,6) ),],1) # first 4 entries 
  h=h+1
  rulesf[h,]=c(rules_per_cell[4*(i)-2,1:4],i)
  # now for the and rule first select which rule of the parent we will add 
  if(i==2 || i==3)
    num=2 # do not sample, only have one activating rule for now in the first cell . 
  else
   num=sample(c(1:2), 1,replace=T, prob=rep(1/2,2) ) # sample which one of the two parent activating rules to add 
  l=length(which(rules_per_cell[4*parent_nodes[i]-(1+num),]!=0) )
  rules_per_cell[4*(i)-2,5:(5+l-1)]=rules_per_cell[4*parent_nodes[i]-(1+num),1:l] # add the and rule 
  vec=rules_per_cell[4*(i)-2,5:(5+l-1)]
  ll=1;
  while(ll <= (l/4)){
   h=h+1
   inds=c((ll*4-3):(ll*4))
   rulesf[h,]=c(vec[inds],i)
   ll=ll+1
  }

  # sample the deactivating rule which are in 11:20
  rules_per_cell[4*(i)-1,1:4]= c(rules[sample(c(12:20), 1, replace=T, prob=rep(1/9,9) ),], -1)
  h=h+1
  rulesf[h,]=c(rules_per_cell[4*(i)-1,1:4],i)
  
  # sample the deactivating rule with an AND from the parent node rules
  rules_per_cell[4*(i),1:4]=c(rules[sample(c(27:32), 1, replace=T, prob=rep(1/6,6) ),], -1) # first 4 entries from pathway rule 
  # now for the 'and' rule first select which rule of the parent we will add  
  if(i==2 || i==3)
    num=2 # do not sample, only have one activating rule for now in the first cell . 
  else
   num=sample(c(1:2), 1,replace=T, prob=rep(1/2,2) ) # sample which one of the two parent activating rules to add 
  l=length(which(rules_per_cell[4*parent_nodes[i]-(num-1),]!=0) )
  rules_per_cell[4*(i),5:(5+l-1)]=rules_per_cell[4*parent_nodes[i]-(num-1),1:l] # add the and rule 
  vec=rules_per_cell[4*(i),5:(5+l-1)]
  ll=1; 
  while(ll <= (l/4)){
    h=h+1
    inds=c((ll*4-3):(ll*4))
    rulesf[h,]=c(vec[inds],i)
    ll=ll+1
  }
}

cat("\n check printing the rules per cell here \n ")
print(rules_per_cell)

# # the regulator matrix is fixed in this setting - minus the noise.
# # get it from the rules associated with the cells. 
reg_ctype=matrix(0, nreg, ncelltypes) 
for(i in 1:ncelltypes){
  for(j in 1:4)
  {
    done=FALSE;k=2
    while(!done){
      if(rules_per_cell[(i-1)*4+j,k]!=0){
        reg=rules_per_cell[(i-1)*4+j,k]
        if(rules_per_cell[(i-1)*4+j,k+1]==1)
         reg_ctype[reg,i]=1
        else
         reg_ctype[,i]=-1      
      }
      k=k+4
      if( rules_per_cell[(i-1)*4+j,k]==0)
        done=TRUE
    }
  }
}
cat("\n check printing the regs ")
print(reg_ctype)


maxrulespercell=4
for(l in 1:ndatasets){
  cat("\n dataset ", l)
 
 # # generating the random labels in c1 and their transitions. want roughly 30percent of the genes active in every cell type
 cat("\n generating labels ")
 labels=matrix(0,ngenes,ncelltypes)
 labels[,1]=sample(c(-1,0,1), ngenes, replace=T, prob=c(0.15,0.7,0.15) )
 for(i in 2:3)
   labels[,i]=sample(c(-1,0,1), ngenes, replace=T, prob=c(0.25,0.5,0.25) )
 for(i in 4:7)
   labels[,i]=sample(c(-1,0,1), ngenes, replace=T, prob=c(0.3,0.4,0.3) )


 # creating the motif by gene matrix from the labels
  mot_gen=matrix(0, nmotif,ngenes)
  
  # code up the rules in a different way. and generate the labels matrix first and reg by cell type 
  # then generate the motif by gene matrix afterwards 
  for(i in 1:ngenes){ # activate some motif for each gene
     #cat("\n gene ", i)
    # activate in cell type 1
    if(labels[i,1]==1)  { # if the gene is on in cell type 1.
      count=length(which(rules_per_cell[1:(maxrulespercell/2),1]>0)) # check how many activating rules there are 
      if(count==1) # if only one motif is responsible, then activate it 
        mot_gen[rules_per_cell[1,1],i]=1 # activate it since its the only possible one 
      else{ # if have multiple rules (motifs) that can be responsible, activate one of them 
        rule=sample(c(1:count), 1,replace=F,prob=rep(1/count,count) )  # selects which rule to satisfy 
        len=length(which(rules_per_cell[rule,]!=0)) # check how many and rules in that one 
        ind= seq(from = 1, to = len, by = 4) # get indices of corresponding motifs 
        motifs=rules_per_cell[rule,ind] # get all motifs corresponding to that rule 
        for(k in 1:length(motifs)){ # go over motifs to activate 
          mot_gen[motifs[k],i]=1 # activate it for that gene 
         # if(i==1)
        #    cat("\n for gene 1 activating motif ", motifs[k])
        }
      }
    }
    else if (labels[i,1]==-1){ # if gene is inactive in cell type 1 
      count=length(which(rules_per_cell[(maxrulespercell/2+1):maxrulespercell,1]>0)) # check how activating rules there are  
      if(count==1) # if only one motif is responsible, then activate it  
        mot_gen[rules_per_cell[(maxrulespercell/2)+1,1],i]=1 # activate it since its the only possible one 
      else{ # if have multiple rules (motifs) that can be responsible, activate one of them 
        # get the motifs from those rules, each just containing one mot-tf interaction
        rule=sample(c((maxrulespercell/2+1):(maxrulespercell)), 1,replace=F,prob=rep(1/count,count) )  # selects which rule to satisfy 
        len=length(which(rules_per_cell[rule,]!=0))
        ind= seq(from = 1, to = len, by = 4)
        motifs=rules_per_cell[rule,ind]
        for(k in 1:length(motifs)){ # go over motifs to activate 
          mot_gen[motifs[k],i]=1 # activate it 
          #if(i==1)
          #  cat("\n for gene 1 activating motif ", motifs[k])
        }
      }
    }
    
    # activate for transitions
   # cat("\n activating in transitions")
    for(j in 1:nrow(transitions)){ # depending on which transitions change its state 
      if((labels[i,transitions[j,1]]==-1 && labels[i,transitions[j,2]]==+1)
         ||  (labels[i,transitions[j,1]]==0 && labels[i,transitions[j,2]]==+1)
         ||  (labels[i,transitions[j,1]]==-1 && labels[i,transitions[j,2]]==0)
      ){
        # cat("\n activated when leaving ", j)
        tocell=transitions[j,2] # transitioning to this cell 
        startind=(tocell)*4-3 # index of rules for that cell 
        endind=startind+(maxrulespercell/2)-1 # 
        len=endind-startind+1
        x=sample(startind:endind,1,replace=F,prob=rep(1/len,len) )  # sample one rule to activate  
        # get its motifs
        len=length(which(rules_per_cell[x,]!=0)); # get how many ANDs in this rule 
        ind= seq(from = 1, to = len, by = 4) # get indices of all relevant motifs 
        motifs=rules_per_cell[x,ind] # 
        mot_gen[motifs,i]=1 # activate them all 
       }
      else if( ( labels[i,transitions[j,1]]==+1 && labels[i,transitions[j,2]]==-1)
        ||  (labels[i,transitions[j,1]]==0 && labels[i,transitions[j,2]]==-1)
      ||  (labels[i,transitions[j,1]]==1 && labels[i,transitions[j,2]]==0)
      ){
        tocell=transitions[j,2] # transitioning to this cell 
        startind=(tocell)*4-maxrulespercell/2+1 # index of rules for that cell 
        endind=(tocell)*4
        len=endind-startind+1
        x=sample(startind:endind,1,replace=F,prob=rep(1/len,len) ) 
        len=length(which(rules_per_cell[x,]!=0)) # get how many ANDs in this rule 
        ind= seq(from = 1, to = len, by = 4)  # get indices of all motifs in that rule 
        motifs=rules_per_cell[x,ind]
        mot_gen[motifs,i]=1
      }
    } # end loop over transition
  }# end loop over genes
   

  for(k in 1:nreg)
  {
    rel=FALSE
    if(length(which(reg_ctype[k,]!=0))>0) # check if that regulator is ever relevant
     rel=TRUE
    for(u in 1:ncelltypes){
      if(rel==TRUE && reg_ctype[k,u]==0)
       reg_ctype[k,u]=sample( c(-1,0,1),1,prob=c(regruleon3,1-2*regruleon3,regruleon3) )
    }

  }
  
  for(u in 1:nmotif){
    for(k in 1:ngenes){
      if(mot_gen[u,k]==0)
        mot_gen[u,k]= sample(c(0,1),1,replace=T,prob=c(1-motruleon2,motruleon2))
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

#   # add noise to test - flip some labels -1 & 1. 
   if(FLIP){  
    noise=matrix(runif(ngenes*ncelltypes), ncol=ncelltypes) 
    noise=which(noise>0.9,arr.ind=TRUE)
    noise=sparseMatrix(i=noise[,1], j=noise[,2],x=1, dims = c(ngenes,ncelltypes))
    pos=which(noise>0&labels>0,arr.ind=TRUE)
    neg=which(noise>0&labels<0,arr.ind=TRUE)
    labels_noisy=labels 
    labels_noisy[pos]=-1; labels_noisy[neg]=1 # 10perc of noise obs flipped.
   }
   
  boosted=boost_complete(as.matrix(subtree_nodes), as.matrix(direct_children), intnodes,
                         as.matrix(mot_gen), as.matrix(t(reg_ctype)), as.matrix(labels), niter , rulesf, original,
                         pos, reg_names, motif_names)
  cat("\n boosting training loss ", boosted$trl[niter]); cat(" test loss ", boosted$tsl[niter])
  recov[l,]=boosted$recovery_Eval
  errors[1,]=c( boosted$trl)+ errors[1,]; errors[2,]=c( boosted$tsl)+errors[2,]
}

cat("\n recovery over datasets \n")
cat("\n rules, correct rules, wrong rules, correct with int node, correct with parent \n")
print(recov)
cat("\n summing the recovery values \n ")
print(colSums(recov)/ndatasets)
plot(c(1:niter),errors[1,]/ndatasets, ylim=c(0.0,0.45) )
lines(c(1:niter),errors[1,]/ndatasets)
lines(c(1:niter),errors[2,]/ndatasets,col="blue")
title(main = list("Average Training and Test errors", cex = 1.5, col = "red", font = 3))
print(errors/ndatasets)
