update_weights<-function(predtr,labels,idxtrall){
  prd = -labels*predtr # entrywise mult
  w = exp(prd)*idxtrall; # mult by the mat indicating indices of included training examples the corresponding scores
  w = w/sum(sum(w))
  return(w)
}

check<-function(score,splitx2sign, rules) {
  checks=FALSE
  if( (rules[1]==1 && rules[2]==1 && score > 0 && splitx2sign>0 )
     || (rules[1]==1 && rules[2]==1 && score< 0 && splitx2sign<0 )  )  # reg sign is positive and activating 
    checks=TRUE
  else if ((rules[1]==1 && rules[2]==-1 && score < 0 && splitx2sign>0 )
    || (rules[1]==1 && rules[2]==-1 && score> 0 && splitx2sign<0 )  )
     checks=TRUE
  else if ((rules[1]==-1 && rules[2]==-1 && score < 0 && splitx2sign<0 )
  || (rules[1]==-1 && rules[2]==-1 && score > 0 && splitx2sign>0 )  )
     checks=TRUE
  else if ( (rules[1]==-1 && rules[2]==1 && score > 0 && splitx2sign<0 )
    || (rules[1]==-1 && rules[2]==1 && score < 0 && splitx2sign>0 )  )
      checks=TRUE
  return(checks)
}

allequal<-function(v1,v2){
  eq=TRUE
  if(length(v1)!=length(v2)){
    cat("\n diff lengths")
    return(FALSE)
  }
  for(i in 1:length(v1)){
    if(v1[i]!=v2[i])
      eq=FALSE
  }
  return(eq)
}

get_candidate<-function(w,idxpredtmp,idxtrup,idxtrdown,wtmp,x1,x2,nx1,nx2,epsilon){  
  
  wtmp = w*idxpredtmp; # ws of training examples that pass parent node
  wpos = wtmp*idxtrup; # ws of +ve training examples that pass parent node - dim of the labels matrix 
  wneg = wtmp*idxtrdown; # ws of -ve training examples that pass parent node
  # separate x2 into positive and negative examples. get wposup and wposdown - get wnegup and wnegdown
  # get 2 wother matrices and 2 loss matrices and then get the minimum loss out of all the losses 
  
  # cat("\n check wpos \n")
  # print(wpos[1:2,1:2])
  x2pos=x2;ind1=which(x2==-1,arr.ind=TRUE);x2pos[ind1]=0 # set negative ones to zero
  
  #  cat("\n check dims ")
  #  print(dim(x1)); print(dim(wpos)); print(dim(x2pos)) 
  wposRegup = x1%*%wpos%*%x2pos; # W+ for feat pairs. mult by x2 to select regs, and by x1 to select the weights of motifs. 
  # cat("\n checks ")
  #print(wposRegup[1:2,1:2])
  # so that each entry in that matrix contains the sum of the weights of all the observations that satisfy the rule (= feature pair)
  # corresponding to that entry
  wnegRegup = x1%*%wneg%*%x2pos; # W- for all feature pairs
  #print(wnegRegup[1:2,1:2])
  wotherRegup = 1 - wposRegup - wnegRegup; # W0 for all feature pairs
  
  x2neg=x2;ind2=which(x2==1,arr.ind=TRUE);x2neg[ind2]=0; # set positive ones to zero
  ind2=which(x2==-1,arr.ind=TRUE);x2neg[ind2]=1;
  wposRegdown = x1%*%wpos%*%x2neg; # W+ for all feature pairs. 
  wnegRegdown = x1%*%wneg%*%x2neg; # W- for all feature pairs
  wotherRegdown = 1 - wposRegdown - wnegRegdown; # W0 for all feature pairs
  
  # loss= matrix(0,nx1,nx2); #loss matrix (store for all obs) for each feat pair (x1,x2) in new node  
  lossup = 2*sqrt(wposRegup*wnegRegup) + wotherRegup; # Z = W0 + 2*sqrt(W+.W-)o
  # loss evaluated at optimized alphas 
  lossdown = 2*sqrt(wposRegdown*wnegRegdown) + wotherRegdown;

  
  lsbestup = min(min(lossup)); lsbestdown =  min(min(lossdown));
  # cat("\n printing loss tables "); print(table(as.vector(lossup)) ); print(table(as.vector(lossdown)) )
  # cat("\n w pos and w neg " ); print(table(as.vector(wposRegdown)));  print(table(as.vector(wnegRegdown)));
  ind=which.min(c(lsbestup,lsbestdown));lsbest=min(lsbestup,lsbestdown); 
  if(ind==1) inds = which(lossup==lsbestup,arr.ind=TRUE); # just take a random one satisfying being the min. 
  if(ind==2) inds = which(lossdown==lsbestdown,arr.ind=TRUE);
  # take = ceil(rand*length[i]);
  mbest = inds[1,1];
  pbest = inds[1,2];
  if(ind==2) ind=-1;
  return(list(mbest=mbest,pbest=pbest,lsbest=lsbest,ind=ind)) # returns ind=1 if reg2 is up in added rule
}

cormat=function(m){
  vars=ncol(m)
  mat=matrix(0,vars,vars)
  cov=(t(m)%*%m)/vars
  means=colMeans(m, na.rm = FALSE, dims = 1) # mean by col
  m1=matrix(rep(means,each=vars),nrow=vars); m2=t(m1) # by col then by row. 
  sds=apply(m, 2, sd)
  sd1=matrix(rep(sds,each=vars),nrow=vars); sd2=t(sd1)
  cormat=(cov-m1*m2)/(sd1*sd2)
  
  cormat=cov(m,m)/(sd1*sd2) # equiv to R version of it : fct cor
  return(cormat)
}

# returns the final training and test error at last iteration, with the ADT recovery evaluation
boost=function(final_leaves, direct_children, intnodes, motif_mat, tf_mat, labels, niter, rules, original ){
  all=c(1:(intnodes+1) )
  nx1 = dim(motif_mat)[1];
  nx2 = dim(tf_mat)[2];
  nrY = dim(labels)[1];
  ncY = dim(labels)[2];
  
  # test set to holdout. read in as a matrix of experiments to keep. 
  holdout=matrix(runif(ncY*nrY), ncol=ncY) 
  holdout=which(holdout>0.6,arr.ind=TRUE)
  holdout=sparseMatrix(i=holdout[,1], j=holdout[,2],x=1)
  idxtstup = (holdout==1 & labels==1) # get matrix of FALSE / TRUE. 
  idxtstdown = (holdout==1 & labels==-1)
  idxtstall = idxtstup+idxtstdown
  
  idxtrup=(holdout==0 & labels==1) # matrix with only upreg entries
  idxtrdown=(holdout==0 & labels==-1) # only updown
  idxtrall=idxtrup+idxtrdown # all up or down.
  
  ntr=sum(sum(idxtrall))
  ntest=sum(sum(idxtstall))
  npredtot=niter+1 #adding a node per iteration. This is the number of nodes in the ADT. 
  
  idxpredtr=vector(mode="list",npredtot); # get vector of matrices. 
  idxpredtst=vector(mode="list",npredtot);
  
  # prediction scores F(x). Updated after adding each node to the ADT. 
  predtr =  as(matrix(0,nrY,ncY) , "sparseMatrix") 
  predtst = as(matrix(0,nrY,ncY) , "sparseMatrix") #sparseMatrix(dims=c(nrY,ncY) )
  
  splitx1 = rep(0,niter); splitx2 = rep(0,niter); splitx2sign = rep(0,niter);
  intnode = rep(0,niter); # to keep track of what internal node in the hierarchy we are at. 
  abovepred = rep(0,npredtot); # save splitter index for every prediction node (to be able to reconstruct the tree)
  abovesplit = rep(0,niter); # save prediction node index for every splitter node
  score = rep(0,npredtot); # save prediction score 'alpha' at every prediction node
  
  trloss = rep(0,niter);
  tstloss = rep(0,niter);
  
  w=(1/ntr)*idxtrall; # just sets them directly in the matrix entries.
  eps=1/ntr; # smoothing
  
  scroot=(1/2)*log( sum(sum(w*idxtrup))/ sum(sum(w*idxtrdown)) )
  
  # update training and test score matrices, just same score for all now since constant root that all pass. 
  predtr = predtr + scroot*idxtrall
  predtst = predtst + scroot*idxtstall
  
  # update the weights for next optim of loss
  npred=1; nsplit=0  
  score[1]=scroot
  w=update_weights(predtr,labels,idxtrall) # 6099*161 matrix of weights for obs
  
  # first root node is constant 1, all training&test pass it.
  idxpredtr[[1]]= idxtrall; # just indicating all observations pass the first node. 
  idxpredtst[[1]]= idxtstall;
  intnode[1]=1 # set the first internal node to just 1 for the root. 
  
  for(i in 1:niter){
    if(i>1){ cat("\n iteration " , i); cat(" test loss " , tstloss[i-1]*100);
             cat(" tr loss ", trloss[i-1]*100) 
             if(trloss[i-1]==0)
               break
    }
    
    tmpls = Inf*rep(1,npred); # min loss if optimally split each current prediction node in ADT. 
    tmpidxm = rep(0,npred); # index of best x1 motif feature at current pred nodes
    tmpidxp = rep(0,npred); # index of best x2 regulator feature at current pred nodes
    tmpidxpsign = rep(0,npred); # sign of best x2 regulator feature at current pred nodes in added rule
    tmpidxc = rep(0,npred); # new node's internal node number in the cell hierarchy if splitting node i 
    
    # get best position and best feature pair for each node in ADT now.
    for(j in 1:npred) # go over all current prediction nodes
    {
      idxpredtemp=idxpredtr[[j]]; # indices of the observations that pass until node
      childnum=1; # so when original equals true it only checks the current node.  
      if(original==FALSE && intnode[j]>0 ) # so not reached a leaf nodes yet which would have a negative sign flag. 
        childnum = childnum+length( which(direct_children[intnode[j],]>0) ) # get children num in hierarchy
      
      loss=rep(0,childnum); motif=rep(0,childnum); reg=rep(0,childnum); regsign=rep(0,childnum);
      
      for(k in 1:childnum){
        #cat("\n checking child  ", k)
        idxpredtempS= idxpredtemp; # nodes that pass until the parent's node only.
        if(k>1) # filter out columns not in leaf subset for that child if not at a leaf node. 
        { 
          child=direct_children[intnode[j],k-1] #one direct child.(k-1 because first option is to stay at same int node) 
          if( length(which(direct_children[intnode[j],] <= -1)) > 0 ) # if some direct children are leaves, keep that node.
          { 
            if( length(which(direct_children[intnode[j],] == -1000))==1 ) # then all are leaf nodes - just keep that leaf
              keep_leaves=child  
            else # Some are leaves (order will be first internal nodes, then leaves. )
            {
              val=direct_children[intnode[j],which(direct_children[intnode[j],] <= -1)] # get number of leaves in those children. 
              ind=length(which(direct_children[intnode[j],]>0) ) # get number of children. 
              if((k-1)<=(ind-abs(val) ) )
              {  keep_leaves=final_leaves[child,] # get still another internal node. 
              }
              else{
                keep_leaves=child # it is a leaf   
              }
            }   
          }
          else # all are internal nodes
            keep_leaves=final_leaves[child,(final_leaves[child,]>0)] # get its leaves in subtree of that child.
          idxpredtempS[,-keep_leaves]=0; # reducenodes that pass for new loss computation
        }
        fit=get_candidate(w,idxpredtempS,idxtrup,idxtrdown,wtmp,motif_mat,tf_mat,nx1,nx2,eps); #optimal motif-reg pair 
        loss[k]=fit$lsbest # optimal loss  for best candidate for node in hierarchy
        motif[k]=fit$mbest; # index of best motif
        reg[k]=fit$pbest; # index of best regulator
        regsign[k]=fit$ind; # sign of best regulator in rule
        
      } # END OF GOING OVER CHILDREN FOR THAT NODE. 
      
      idx=which.min(loss) #  which child for that nodes's split is best  
      # cat("\n updating tmpls ", j); cat("\n from the losses "); print(loss) # get losses of only 1 when nothing passes that rule. 
      tmpls[j]=loss[idx]; tmpidxm[j]=motif[idx]; tmpidxp[j]=reg[idx]; tmpidxpsign[j]=regsign[idx]; 
      children=direct_children[intnode[j],]
      
      if(idx==1 || childnum==1) # no additional split wrt parent. stick to same internal node in cell hierarchy
        tmpidxc[j]=intnode[j] 
      else if ( (length(which(children <=-1)) > 0 ) && # some leaves 
                  ( (idx -1)> ( length( which(children>0)  ) - abs(children[which(children <= -1)]) ) ) ) # optimal is leaf
      {
        if(intnode[j]>0) # first encounter with this leaf node. 
          tmpidxc[j]=-direct_children[intnode[j],idx-1] # give leaf node negative index in int node list 
      }      
      else #  move to new internal node in the tree hierarchy (none are leaves, or some are but optimal is internal)
      { tmpidxc[j]= direct_children[intnode[j],idx-1] ; #cat("\n move to new internal node ");
      }
    } # end of going over current prediction nodes. 
    
    #cat("\n all losses for ADT nodes possible splits ") ;print(tmpls) ; 
    minindex=which.min(tmpls) #index of pred node whose split gives best loss, set the feature splitters
    m= tmpidxm[minindex] # optimal motif of optimal new rule
    p= tmpidxp[minindex] # optimal reg
    psign= tmpidxpsign[minindex] # opt sign for new reg rule
    child= tmpidxc[minindex] # new internal node index (negative if its a leaf node)
    
    # add new splitter node to current ADT structure.
    nsplit=nsplit+1 # num of nodes
    splitx1[nsplit]=m; splitx2[nsplit]=p; splitx2sign[nsplit]=psign; # store new motif and regulator for that node. 
    abovesplit[nsplit]=minindex # store the parent of the new node, or which node we split to add the new rule
    intnode[npred+1]=child # corresponding index for the new node's internal node index in the output hierarchy. 
    cat("\n set the int node of the newly added splitter node at ", (npred+1)); cat("\n with value ", child )
    
    
    ##### GET EXAMPLES THAT PASS THROUGH NEW SPLITTER NODE THAT WAS ADDED. ##### 
    x = which(idxpredtr[[minindex]]>0, arr.ind=TRUE); # array indices passing prev nodes of pred node we split. 
    validr = which(motif_mat[m,]!=0); validc = which(tf_mat[,p]==psign); # genes that have that motif (gives rows), celltypes that have that regulator.
    tf = match(x[,1],validr,nomatch = FALSE) & match(x[,2],validc,nomatch = FALSE) # get out of obs that satisfy the rule which passed the parent. 
    r = x[tf,1]; c = x[tf,2]; # only select ones that pss the parent and the new rule 
    yestr = sparseMatrix(i=r,j=c,x=1,dims=c(nrY,ncY) );  # create matrix with 1 for (i,j) if row/col vals of obs pass new splitter node) 
    if(child>0 && child!=1) # if not at a leaf node, need to filter more. Or if first time reach that leaf node 
    {  all=c(1:ncY); remove=all[-final_leaves[child,which(final_leaves[child,]>0)]]
       yestr[,remove] = 0 # remove final leaves, unless at root of hierarchy in which case keep all. 
       # cat("\n removed leaves ") ; print(remove)
    }
    else if (child < 0 )# if at a leaf node, set all the other columns to zero. 
    {  yestr[,child] = 0 # child negative so directly set all else to 0
       #cat("\n only keeping ", -child)
    }
    idxpredtr[[npred+1]] = yestr;
    
    #### GET THE ONES FROM THE TEST SET THAT PASS THE NEW SPLITTER NODE. ####
    x = which(idxpredtst[[minindex]]>0, arr.ind=TRUE); 
    tf = match(x[,1],validr,nomatch = FALSE) & match(x[,2],validc,nomatch = FALSE);
    r = x[tf,1]; c = x[tf,2];
    yestst = sparseMatrix(i=r,j=c,x=1,dims=c(nrY,ncY) ); 
    if(child>0 & child!=1) # case with 1 just keep everything. 
      yestst[,-final_leaves[child,which(final_leaves[child,]>0)]] = 0
    else if(child<0) # cause -child in that case contains the column number corresponding to that leaf/cell type. 
      yestst[,child] = 0 
    idxpredtst[[npred+1]] = yestst;
    
    # add new prediction node 
    npred=npred+1
    abovepred[npred]=nsplit
    wpos=w*idxtrup; wneg=w*idxtrdown;
    score[npred]= 1/2*log((sum(sum(wpos*yestr)) + eps)/(sum(sum(wneg*yestr)) + eps));
    
    # update PREDICTION SCORE computed over the ADT as F_T
    predtr= predtr + yestr*score[npred]; # update matrix of prediction for train and test sets. 
    # scores for each observation that passes the rules until that newly added node. 
    predtst= predtst + yestst*score[npred];
    
    # update weights of all observations. # no modification. already contains correct prediction scores 
    w=update_weights(predtr,labels,idxtrall);
    
    # calculate training and test error by summing over all observations. 
    trloss[i] = sum(sum(labels*predtr<0))/ntr;
    tstloss[i] = sum(sum(labels*predtst<0))/ntest;
    # or diff with other measure: just the same but with loss computed sep on pos and neg ex
  }
  
  # stats for ADT recovery - number of rules, number of correct rules, number of rules with wrong parent, number of rules with wrong internal node
  # size of training data for dataset, final training error, 
  nrules=niter; 
  if(i <= niter)
    nrules=i-1;
  
  got=rep(0,nrow(rules)); wrongrules=0;gotwrule=rep(0,nrow(rules)); gotparent=rep(0,nrow(rules))
  for(i in 1:nrules){# go over all added rules
    found=FALSE
    for(j in 1:nrow(rules) ){ # go over rules from ADT
      vec=c(splitx1[i],splitx2[i],splitx2sign[i])
      if(allequal(vec,rules[j,1:3]) )
      {  found=TRUE
         got[j]=1
         if(intnode[i+1]==rules[j,4]) # so if once gets it with the correct child -- then counting it as correct for now ..
           gotwrule[j]=1
         if( (abovesplit[i]==1 && rules[j,5]==0 ) )
         { gotparent[j]=1 ;}
         else if( abovesplit[i]>1 && rules[j,5]>0 ){
           if( splitx1[abovesplit[i]-1] == (rules[(rules[j,5]), 1]) 
               && splitx2[abovesplit[i]-1] == (rules[(rules[j,5]), 2]) )
             gotparent[j]=1 ;
         }
         break
      }
      
    }
    if(!found)
      wrongrules=wrongrules+1
  }
  
  correct=sum(got); rule_int_nodes=sum(gotwrule); rule_parents=sum(gotparent);
  recovery_Eval=c( nrules,  correct, wrongrules, rule_int_nodes, rule_parents, trloss[nrules], tstloss[nrules])
  
  par(mfrow=c(1, 1))
  plot(tstloss, ylim=c(min(trloss), max(trloss) ) )
  lines(tstloss)
  lines(trloss,col="blue")
  title("Training and test error versus iteration number ")  
  cat("\n test losses ")
  print(tstloss)
  cat("\n split x1 ")
  print(splitx1); cat("\n split x2 and sign "); print(splitx2)
  print(splitx2sign)
  cat("\n above pred node")
  print(abovepred)  
  cat("\n above split ");print(abovesplit) 
  cat("\n intnodes "); print(intnode)
  cat("\n split motif "); print(splitx1)
  cat("\n split reg "); print(splitx2)
  cat("\n reg sign "); print(splitx2sign)
  
  ####
  genenames=c(1:500)#readLines("~/Dropbox/reg_networks/BoostingTestData/original/target_gene_names_G.txt") # 6099
  regnames=c(1:45)#readLines("~/Dropbox/reg_networks/BoostingTestData/original/reg_names_R.txt") # 472
  motifnames=c(1:22)#readLines("~/Dropbox/reg_networks/BoostingTestData/original/motif_names_M.txt") # 203
  
  cat("\n root: score ", score[1]) ; cat(" int node: 1 ")
  ###### output for latex table: ######
  for(i in 1:nrules){
    cat("\n ", (i+1) )
    cat(" & " , motifnames[splitx1[i]])
    cat(" & " , regnames[splitx2[i]])
    cat(" & " , splitx2sign[i])
    cat(" & " , score[i+1] ) # first is root score. 
    cat(" & ", abovesplit[i])
    cat(" & " , intnode[i+1])
    cat(" \\\\ "); cat('\\hline')
  }
  
  return(list(trl=trloss,tsl=tstloss,recovery_Eval=recovery_Eval))
}

# returns the final training and test error at last iteration, with the ADT recovery evaluation
# done - requirement : numbering of cell types is from 1...ncelltypes
# direct children are going to be the direct node numbers (unique values)
# final leaves: contains all the internal nodes and leaves in a subtree. If looking at final leaves from
# internal node 2: the format should be [2,etc ]. The node itself needs to be included as well because it belongs 
# the cell types that need to be checked. 
boost_complete=function(final_leaves, direct_children, intnodes, motif_mat, tf_mat, labels, niter, rules, original, file1peakgenes,
                        file2regulators, file3motifs){
  cat("\n boosting with ADT full node tree ")
  all=c(1:(intnodes+1) )
  nx1 = dim(motif_mat)[1];
  nx2 = dim(tf_mat)[2];
  nrY = dim(labels)[1];
  ncY = dim(labels)[2];
  print(nrY)
  print(ncY)
  
  # test set to holdout. read in as a matrix of experiments to keep. 
  cat("\n holding out ")
  holdout=matrix(runif(ncY*nrY), ncol=ncY) 
  #print(dim(holdout))
  holdout=which(holdout>0.6,arr.ind=TRUE)
  #print(dim(holdout))
  #print(holdout)
  holdout=sparseMatrix(i=holdout[,1], j=holdout[,2],x=1,dims=c(nrY,ncY))
  #print(dim(holdout))
  #print(dim(labels))
  #write.matrix(holdout, file = "~/Dropbox/reg_networks/holdout_var020.txt", sep = " ")
  idxtstup = (holdout==1 & labels==1) # get matrix of FALSE / TRUE. 
  idxtstdown = (holdout==1 & labels==-1) # WAS -1
  idxtstall = idxtstup+idxtstdown
  
  idxtrup=(holdout==0 & labels==1) # matrix with only upreg entries
  idxtrdown=(holdout==0 & labels==-1) # only updown
  idxtrall=idxtrup+idxtrdown # all up or down.
  
  ntr=sum(sum(idxtrall))
  ntest=sum(sum(idxtstall))
  npredtot=niter+1 #adding a node per iteration. This is the number of nodes in the ADT. 
  
  idxpredtr=vector(mode="list",npredtot); # get vector of matrices. 
  idxpredtst=vector(mode="list",npredtot);
  
  # prediction scores F(x). Updated after adding each node to the ADT. 
  predtr =  as(matrix(0,nrY,ncY) , "sparseMatrix") 
  predtst = as(matrix(0,nrY,ncY) , "sparseMatrix") #sparseMatrix(dims=c(nrY,ncY) )
  
  splitx1 = rep(0,niter); splitx2 = rep(0,niter); splitx2sign = rep(0,niter);
  intnode = rep(0,niter); # to keep track of what internal node in the hierarchy we are at. 
  abovepred = rep(0,npredtot); # save splitter index for every prediction node (to be able to reconstruct the tree)
  abovesplit = rep(0,niter); # save prediction node index for every splitter node
  score = rep(0,npredtot); # save prediction score 'alpha' at every prediction node
  
  trloss = rep(0,niter);
  tstloss = rep(0,niter);
  trba = rep(0,niter);
  tstba = rep(0,niter);
  
  w=(1/ntr)*idxtrall; # just sets them directly in the matrix entries.
  eps=1/ntr; # smoothing
  scroot=(1/2)*log( sum(sum(w*idxtrup))/ sum(sum(w*idxtrdown)) )
  
  # update training and test score matrices, just same score for all now since constant root that all pass. 
  predtr = predtr + scroot*idxtrall
  predtst = predtst + scroot*idxtstall
  
  # update the weights for next optim of loss
  npred=1; nsplit=0  
  score[1]=scroot
  
  w=update_weights(predtr,labels,idxtrall) # 6099*161 matrix of weights for obs

  
  # first root node is constant 1, all training&test pass it.
  idxpredtr[[1]]= idxtrall; # just indicating all observations pass the first node. 
  idxpredtst[[1]]= idxtstall;
  intnode[1]=1 # set the first internal node to just 1 for the root. 
  
  for(i in 1:niter){
    if(i>1){ cat("\n iteration " , i); cat(" test loss " , tstloss[i-1]*100);
             cat(" tr loss ", trloss[i-1]*100) 
             cat("\n iteration " , i); cat(" test Balanced Accuracy " , tstba[i-1]);
             cat(" tr Balanced Accuracy ", trba[i-1]) 
             if(trloss[i-1]==0)
               break
    }
    
    tmpls = Inf*rep(1,npred); # min loss if optimally split each current prediction node in ADT. 
    tmpidxm = rep(0,npred); # index of best x1 motif feature at current pred nodes
    tmpidxp = rep(0,npred); # index of best x2 regulator feature at current pred nodes
    tmpidxpsign = rep(0,npred); # sign of best x2 regulator feature at current pred nodes in added rule
    tmpidxc = rep(0,npred); # new node's internal node number in the cell hierarchy if splitting node i 
    
    # get best position and best feature pair for each node in ADT now.
    for(j in 1:npred) # go over all current prediction nodes
    {
      idxpredtemp=idxpredtr[[j]]; # indices of the observations that pass until node
      childnum=1; # so when original equals true it only checks the current node.  
      subtree=which(final_leaves[intnode[j],]>0)
      if(original==FALSE && length(subtree)>1 ) # so not reached a leaf node yet which would have a negative sign flag. 
        childnum = childnum+length( which(direct_children[intnode[j],]>0) ) # get children num in hierarchy
      
      #cat("\n child num to check for that node ")
      loss=rep(0,childnum); motif=rep(0,childnum); reg=rep(0,childnum); regsign=rep(0,childnum);
      
      for(k in 1:childnum){
        #cat("\n checking child  ", k)
        
        idxpredtempS= idxpredtemp; # nodes that pass until the parent's node only.
        if(k>1) # filter out columns not in leaf subset for that child if not at a leaf node. 
        { 
          child=direct_children[intnode[j],k-1] #one direct child.(k-1 because first option is to stay at same int node) 
          keep_leaves=final_leaves[child,(final_leaves[child,]>0)] # get its leaves in subtree of that child.
          idxpredtempS[,-keep_leaves]=0; # reduce nodes that pass for new loss computation
          #cat("\n for child ", child)
          #cat("\n keep leaves ")
          #print(keep_leaves)
        }

        
        #cat("\n in get candidtate")
        fit=get_candidate(w,idxpredtempS,idxtrup,idxtrdown,wtmp,motif_mat,tf_mat,nx1,nx2,eps); #optimal motif-reg pair 
        #cat("\n out of get candidate")
        loss[k]=fit$lsbest # optimal loss for best candidate for node in hierarchy
        motif[k]=fit$mbest; # index of best motif
        reg[k]=fit$pbest; # index of best regulator
        regsign[k]=fit$ind; # sign of best regulator in rule
        
      } # END OF GOING OVER CHILDREN FOR THAT NODE. 
      
      idx=which.min(loss) #  which child for that nodes's split is best  
      #cat("\n updating tmpls ", j); cat("\n from the losses "); print(loss) # get losses if only 1 when nothing passes that rule. 
      tmpls[j]=loss[idx]; tmpidxm[j]=motif[idx]; tmpidxp[j]=reg[idx]; tmpidxpsign[j]=regsign[idx]; 
      #children=direct_children[intnode[j],]
      
      if(idx==1 || childnum==1) # no additional split wrt parent. stick to same internal node in cell hierarchy
        tmpidxc[j]=intnode[j]       
      else #  move to new internal node in the tree hierarchy (none are leaves, or some are but optimal is internal)
        tmpidxc[j]= direct_children[intnode[j],idx-1] ; #cat("\n move to new internal node ");
      
    } # end of going over current prediction nodes. 
    
    #cat("\n all losses for ADT nodes possible splits ") ; print(tmpls) ; 
    minindex=which.min(tmpls) #index of pred node whose split gives best loss, set the feature splitters
    m= tmpidxm[minindex] # optimal motif of optimal new rule
    p= tmpidxp[minindex] # optimal reg
    psign= tmpidxpsign[minindex] # opt sign for new reg rule
    child= tmpidxc[minindex] # new internal node index (negative if its a leaf node)
    
    # add new splitter node to current ADT structure.
    nsplit=nsplit+1 # num of nodes
    splitx1[nsplit]=m; splitx2[nsplit]=p; splitx2sign[nsplit]=psign; # store new motif and regulator for that node. 
    abovesplit[nsplit]=minindex # store the parent of the new node, or which node we split to add the new rule
    intnode[npred+1]=child # corresponding index for the new node's internal node index in the output hierarchy. 
    #cat("\n set the int node of the newly added splitter node at ", (npred+1)); cat("\n with value ", child )
    
    
    ##### GET EXAMPLES THAT PASS THROUGH NEW SPLITTER NODE THAT WAS ADDED. ##### 
    x = which(idxpredtr[[minindex]]>0, arr.ind=TRUE); # array indices passing prev nodes of pred node we split. 
    validr = which(motif_mat[m,]!=0); validc = which(tf_mat[,p]==psign); # genes that have that motif (gives rows), celltypes that have that regulator.
    tf = match(x[,1],validr,nomatch = FALSE) & match(x[,2],validc,nomatch = FALSE) # get out of obs that satisfy the rule which passed the parent. 
    r = x[tf,1]; c = x[tf,2]; # only select ones that pss the parent and the new rule 
    yestr = sparseMatrix(i=r,j=c,x=1,dims=c(nrY,ncY) );  # create matrix with 1 for (i,j) if row/col vals of obs pass new splitter node) 
    if( child!=1) # if NOT at a leaf node, need to filter more. Or if first time reach that leaf node .Just root node needs no filter. 
    {  all=c(1:ncY); remove=all[-final_leaves[child,which(final_leaves[child,]>0)]]
       yestr[,remove] = 0 # remove final leaves, unless at root of hierarchy in which case keep all. 
       # cat("\n removed leaves ") ; print(remove)
    }
    idxpredtr[[npred+1]] = yestr;
    
    #### GET THE ONES FROM THE TEST SET THAT PASS THE NEW SPLITTER NODE. ####
    x = which(idxpredtst[[minindex]]>0, arr.ind=TRUE); 
    tf = match(x[,1],validr,nomatch = FALSE) & match(x[,2],validc,nomatch = FALSE);
    r = x[tf,1]; c = x[tf,2];
    yestst = sparseMatrix(i=r,j=c,x=1,dims=c(nrY,ncY) ); 
    if( child!=1) # case with 1 just keep everything. 
      yestst[,-final_leaves[child,which(final_leaves[child,]>0)]] = 0
    idxpredtst[[npred+1]] = yestst;
    
    # add new prediction node 
    npred=npred+1
    abovepred[npred]=nsplit
    wpos=w*idxtrup; wneg=w*idxtrdown;
    score[npred]= 1/2*log((sum(sum(wpos*yestr)) + eps)/(sum(sum(wneg*yestr)) + eps));
    
    # update PREDICTION SCORE computed over the ADT as F_T
    predtr= predtr + yestr*score[npred]; # update matrix of prediction for train and test sets. 
    # scores for each observation that passes the rules until that newly added node. 
    predtst= predtst + yestst*score[npred];
    
    # update weights of all observations. # no modification. already contains correct prediction scores 
    w=update_weights(predtr,labels,idxtrall);
    
    # calculate training and test error by summing over all observations. 
    trloss[i] = sum(sum(labels*predtr<0))/ntr;
    tstloss[i] = sum(sum(labels*predtst<0))/ntest;  

    
    posl=which(holdout*labels>0); negl=which(holdout*labels<0);
    posp=which(holdout*predtst>0); negp=which(holdout*predtst<0); 
    TP=length(which(match(posl,posp)!='NA') ); TN=length(which(match(negp,negl)!='NA'))
    FP=length(which(match(posp,negl)!='NA')); FN=length(which(match(posl,negp)!='NA'))
    tstba[i] = 0.5*( TP/(TP+FN) + TN/(FP+TN) ) ;
  
    
    posl=which((1-holdout)*labels>0); negl=which((1-holdout)*labels<0);
    posp=which((1-holdout)*predtr>0); negp=which((1-holdout)*predtr<0); 
    TP=length(which(match(posl,posp)!='NA') ); TN=length(which(match(negp,negl)!='NA'))
    FP=length(which(match(posp,negl)!='NA')); FN=length(which(match(posl,negp)!='NA'))
    trba[i] = 0.5*( TP/(TP+FN) + TN/(FP+TN) ) ;
    
    # or diff with other measure: just the same but with loss computed sep on pos and neg ex
  }
  
  # stats for ADT recovery - number of rules, number of correct rules, number of rules with wrong parent, number of rules with wrong internal node
  # size of training data for dataset, final training error, 
  nrules=niter; 
  if(i <= niter)
    nrules=i-1;
  
  got=rep(0,nrow(rules)); wrongrules=0;gotwrule=rep(0,nrow(rules)); gotparent=rep(0,nrow(rules))
  countintnodetotal=0; countintnodecorrectmotif=0; countintnodecorrectreg=0; # just counts if the motif or the reg placed there
    # is actually placed correctly. 
  for(i in 1:nrules){# go over all added rules
    found=FALSE
    for(j in 1:nrow(rules) ){ # go over rules from ADT
      vec=c(splitx1[i],splitx2[i],splitx2sign[i])
      if( #allequal(vec,rules[j,1:3]) || 
        # this is just checking the motif regulator matching . 
         ( allequal(vec[1:2],rules[j,1:2]) && check(score[i+1],splitx2sign[i], rules[j,3:4])  )  )
      {  found=TRUE
         got[j]=1
         if(intnode[i+1]==rules[j,4]) # so if once gets it with the correct child -- then counting it as correct for now ..
           gotwrule[j]=1       
         if( (abovesplit[i]==1 && rules[j,5]==0 ) )
         { gotparent[j]=1 ;}
         else if( abovesplit[i]>1 && rules[j,5]>0 ){
           if( splitx1[abovesplit[i]-1] == (rules[(rules[j,5]), 1]) 
               && splitx2[abovesplit[i]-1] == (rules[(rules[j,5]), 2]) )
             gotparent[j]=1 ;
         }
         break
      }   
    }
    if(!found)
      wrongrules=wrongrules+1
    if( intnode[i+1]>1) # if it has a pointer. need to check that the motif and regulator (separately) are 
       # in fact in that subtree 
    {
      #cat("\n for the rule "); print(intnode[i+1]); print(splitx1[i]); print(splitx2[i])
      countintnodetotal=countintnodetotal+1
      relintnodes=final_leaves[intnode[i+1],]
      ind=which(rules[,5]  %in% relintnodes ) #==intnode[i+1]) # get the rules that are placed for that cell type
      #cat("\n relevant cell types ", relintnodes)
        if( splitx1[i] %in% rules[ind,1]) # actually just want to check if correct side of the tree 
                                          # so for that need to use the subtree nodes. 
        {  countintnodecorrectmotif=countintnodecorrectmotif+1; 
           #cat("\n motif corrrect ")
        }
        if( splitx2[i] %in% rules[ind,2])
        {  countintnodecorrectreg=countintnodecorrectreg+1;
           #cat("\n reg correct ")
        }
    }
  }
  
  correct=sum(got); rule_int_nodes=sum(gotwrule); rule_parents=sum(gotparent);
  recovery_Eval=c( nrules,  correct, wrongrules, rule_int_nodes, rule_parents, trloss[nrules], tstloss[nrules],
                   countintnodetotal,countintnodecorrectmotif,countintnodecorrectreg)
  
  par(mfrow=c(1, 1))
  plot(tstloss, ylim=c(0, 0.5 ) )
  lines(tstloss)
  lines(trloss,col="blue")
  title("Training and test error versus iteration number ")  
  cat("\n test losses ")
  print(tstloss)
  plot(tstba, ylim=c(0,1 ) )
  lines(tstba)
  lines(trba,col="blue")
  title("Training and test balanced accuracy versus iteration number ")  
  cat("\n test ba ")
  print(tstba)
  
  cat("\n split x1 ")
  print(splitx1); cat("\n split x2 and sign "); print(splitx2)
  print(splitx2sign)
  cat("\n above pred node")
  print(abovepred)  
  cat("\n above split ");print(abovesplit) 
  cat("\n intnodes "); print(intnode)
  cat("\n split motif "); print(splitx1)
  cat("\n split reg "); print(splitx2)
  cat("\n reg sign "); print(splitx2sign)
  
  
  genenames=file1peakgenes #"~/Dropbox/reg_networks/BoostingTestData/20150415_forPeyton/peak.headers.txt") # 6099
  regnames=readLines(file2regulators)#"~/Dropbox/reg_networks/BoostingTestData/original/reg_names_R.txt") # 472
  motifnames=readLines(file3motifs)#"~/Dropbox/reg_networks/BoostingTestData/20150415_forPeyton/annotationMatrix.headers.txt") # 203
  
  cat("\n root: score ", score[1]) ; cat(" int node: 1 ")
  ###### output for latex table: ######
  for(i in 1:nrules){
    cat("\n ", (i+1) )
    cat(" & " , motifnames[splitx1[i]])
    cat(" & " , regnames[splitx2[i]])
    cat(" & " , splitx2sign[i])
    cat(" & " , score[i+1] ) # first is root score. 
    cat(" & ", abovesplit[i])
    cat(" & " , intnode[i+1])
    cat(" \\\\ "); cat('\\hline')
  }
  
  return(list(trl=trloss,tsl=tstloss,trba=trba,tstba=tstba , recovery_Eval=recovery_Eval))
}

