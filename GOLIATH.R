


library("RColorBrewer")
library("ggplot2")
library("beepr")
library("mclust")
library("ROSE")
library("kernelboot")
library("dbscan")
library("greybox")
library("pracma")
library("geometry")
library("TDA")
library("reticulate")
library("uniformly")
library("mvtnorm")
library("dplyr")
library("data.table")
library("FNN")
library("rgl")
library("philentropy")
library("mixtools")
library("MASS")
#library("lcmix")
library("Hmisc")
library("simukde")
library("np")
library("ks")
library("Ake")
library("rtemis")
library("randomForest")


# Algorithmes

## GOLIATH Generator



GOLIATH = function(X,mode="mix",type,Y=NULL,method_Y=5,clustering=F,seed=NULL,synth=NULL,m=nrow(X),w=rep(1,nrow(X)),N=nrow(X), CSB_param=NULL, NCSB_param=NULL, NNSB_param=NULL, BSB_param=NULL){
  n = nrow(X)
  if (is.null(Y)){train0=X} else {train0=cbind(X,Y)}
  if (clustering=="train") {
    clust = cluster(train0)
  } else if (clustering=="Y"){
    test_clust_w = endo                                                         # cluster sur le Y et w
    test_clust_w$w = w
    clust = cluster(test_clust_w)
  }else {
    clust = data.frame("clust" = rep(1,n))
    row.names(clust) = row.names(X)
  }
  
  X_clust = cbind(X,clust)
  nb_clust = max(clust)
  
  # tirage des graines
  if (is.null(seed)) {
    seed = sample(seq(nrow(X)),N, replace=T, prob=w)
    seed = row.names(X[seed,])
    seed = as.character(floor(as.numeric(seed))) 
  }
  seed_clust = data.frame("seed"= seed,"clust" = X_clust[seed,"clust"])
  
  # generation données synthétique 
  if (m<n){                                                                     # Mixture Model approaches (w is estimated)
    if (type == "GaussMM"){                                                     ## Gaussian Mixture Model : GMM de MClust
      print("Gaussian Mixture Model Generator")
      synth = MM_GaussMM(train0,m,N)
    }
  } else {
    for (cl in seq(nb_clust)){
      X_cl = X[clust==cl,] # OK CAR CLUSTER en fonction de train0 en fonction de X
      seed_cl = seed_clust[seed_clust$clust==cl,"seed"]
      N_cl = length(seed_cl)
      w_cl = w[clust==cl] # OK CAR w en fonction de X
      
      # Kernel approaches (m=n)
      if (type == 'CSB'){                                                         ## (Weighted) Classical Smoothed Boostrap
        if (CSB_param$kern  %in% c("multivariate", "gaussian", "epanechnikov", "rectangular", "triangular","biweight", "cosine", "optcosine", "none")){
          print("Classical Smoothed Boostrap Generator")                          # Classical Smoothed Bootstrap (H : Silverman)
          synth_cl = GOLIATH_CSB(X=X_cl,N=N_cl, w=w_cl,seed=seed_cl, kern=CSB_param$kern, hmult=CSB_param$hmult)
        }
        if (CSB_param$kern == 'ROSE'){                                            ## Gaussian Smoothed Bootstrap (H : ROSE (Bowman-Azzalini))
          print("ROSE Generator")
          synth_cl = GOLIATH_ROSE(X=X_cl,N=N_cl, w=w_cl,seed=seed_cl,hmult=CSB_param$hmult, w_sd=CSB_param$w_sd[clust==cl])
        }
        if (CSB_param$kern == 'GN'){                                              ## Gaussian Noise (H : pert x Sigma)
          print("Gaussian Noise Generator")
          synth_cl = GOLIATH_GN(X=X_cl,N=N_cl, w=w_cl,seed=seed_cl,pert=CSB_param$pert, w_sd = CSB_param$w_sd[clust==cl])
        }
        if (CSB_param$kern == 'FGSB'){                                            ## Full bandwidth matrices (Non-diagonal) Gaussian SB 
          print("Full bandwidth matrice Gaussian Smoothed Bootstrap Generator")
          synth_cl = GOLIATH_FGSB(X=X_cl,N=N_cl,w=w_cl,seed=seed_cl)
        }
      }
      if (type == 'NCSB'){                                                        ## (weighted) Non-Classical Smoothed Bootstrap
        print("Non-Classical Smoothed Boostrap Generator")
        synth_cl = GOLIATH_NCSB(X=X_cl,N=N_cl,w= w_cl,seed=seed_cl,hmult=NCSB_param$hmult,NCSB_param=NCSB_param)
      }
      if (type == 'NNSB'){                                                        ## Nearest-Neighbors Smoothed Bootstrap
        if (is.null(NNSB_param$Ex_sig)){
          print("SMOTE-Family Generator")
          synth_cl = GOLIATH_SMOTEf(X=X_cl,N=N_cl,w=w_cl,seed=seed_cl, k=NNSB_param$k, pi_ij=NNSB_param$pi_ij, fun_interp=NNSB_param$fun_interp, dist=NNSB_param$dist, Beta_alpha=NNSB_param$Beta_alpha, Beta_beta=NNSB_param$Beta_beta, LogNorm_mu=NNSB_param$LogNorm_mu, LogNorm_sig=NNSB_param$LogNorm_sig)
        }
        else {
          print('Extended-SMOTE')
          synth_cl = GOLIATH_ExSMOTEf(X=X_cl,N=N_cl,w=w_cl,seed=seed_cl, k=NNSB_param$k, pi_ij=NNSB_param$pi_ij, fun_interp=NNSB_param$fun_interp, dist=NNSB_param$dist, Beta_alpha=NNSB_param$Beta_alpha, Beta_beta=NNSB_param$Beta_beta, LogNorm_mu=NNSB_param$LogNorm_mu, LogNorm_sig=NNSB_param$LogNorm_sig, Ex_sig=NNSB_param$Ex_sig)
        }
      }  
      if (type %in% c('BSB')){                                                    ## Balloon Smoothed Bootstrap
        if (min(w)==max(w)){  
          print("kNN Density  Estimate Generator")
          synth_cl = GOLIATH_kNN_DE(X_cl,N_cl,k = BSB_param$k, n_grid=BSB_param$n_grid, p_sd=BSB_param$p_sd, d=BSB_param$dist, p_mink=BSB_param$p_mink, exp_dist=BSB_param$exp_dist,BSB_param$d_min)}
        else {
          print("Weighted kNN Density  Estimate Generator")
          synth_cl = GOLIATH_CkNN_DE(X_cl,N_cl,w_cl,seed_cl,k = BSB_param$k, n_grid=BSB_param$n_grid, p_sd=BSB_param$p_sd, d=BSB_param$dist, p_mink=BSB_param$p_mink, exp_dist=BSB_param$exp_dist,BSB_param$d_min,cp_knn=BSB_param$cp_knn)
          
        }
      }
      
      # Generation of Y
      if (is.null(Y)==F){                                                           
        synth_cl[,colnames(Y)] = GOLIATH_Y(X_cl,subset(Y,subset= (clust==cl)),synth_cl,seed_cl,method=method_Y, pert=CSB_param$pert,w_sd=CSB_param$w_sd[clust==cl])
      }
      
      if (cl==1){
        synth = synth_cl
      } else {synth = rbind(synth,synth_cl)}
    }
  }
  
  
  # Restitution
  if (type == "GaussMM"){
    train_synth = synth                                                         # si generateur = GMM, renvoi du jeu synth
  } else {
    if (mode == "synth"){
      train_synth = synth} # full synth : renvoie synth
    else if (mode == "augment"){
      train_synth = rbind(train0,synth)} # augmented dataset : renvoie (X,Y)  + synth
    else if (mode == "mix"){
      train0=X[seed,]
      if (!is.null(Y)){train0[,colnames(Y)] =Y[seed,]} # (X,Y)-seed original : over/under sampling
      train_synth = synth # synth on (X,Y)-seed
      occ = (duplicated(seed))*1   # = 1 à partir seconde occurence, 0 pour la première
      for (lig in seq(length(seed))){
        if (occ[lig]==0){train_synth[lig,]=train0[lig,]} # si première occurence alors on remplace la synth par seed
      }
    }
  }
  
  list("synth"= train_synth, "seed"=seed)
}



  ## Generation of $Y$ : Nadaraya-Watson Estimator (kernel regression) 
  
  
GOLIATH_Y = function(X,Y,synth,seed, method=5, pert=0.1,w_sd=rep(1,nrow(X))){
  print(paste("method_Y : ", method))
  x_seed = X[seed,]
  y_seed = Y[seed,]
  train = X
  Y_lab = colnames(Y)
  train[,Y_lab] =Y
  Y_w = Y[w_sd,]
  
  # RF predicton 
  model = randomForest(formula(paste0(Y_lab ,"~ .")), data = train)
  pred = predict(model,synth, predict.all = T)
  pred_seed = predict(model,x_seed, predict.all = T)
  
  ###  perturbation de la graine y
  synth[,Y_lab] = 0 # à modifier dans goliath
  
  # method 1 : ajout du delta de prediction (mean) sur la graine
  if (method == 0){ 
    synth[,Y_lab] =  y_seed
  }
  # method 1 : ajout du delta de prediction (mean) sur la graine
  else if (method == 1){ 
    for (i in seq(nrow(synth))){
      eps = pred$aggregate[i] - pred_seed$aggregate[i]  
      synth[i,Y_lab] =  y_seed[i]  + eps
    }
  }
  # method 2 : tirgae selon distrib erreur et génération d'un y : wild-smoothed boostrap / essayer une smoothboot
  else if (method == 2){  
    for (i in seq(nrow(synth))){
      eps = y_seed[i] - pred_seed$individual[i,]  
      h = sqrt(as.numeric(bw.silv(as.data.frame(eps))))
      synth[i,Y_lab] =  y_seed[i]  +  as.numeric(abs(rnorm(1, 0 ,h)) * sign(pred$aggregate[i]-pred_seed$aggregate[i]))
    }
  }
  # method 3 : bruitage noyau de la distribution de prediction
  else if (method == 3){  
    for (i in seq(nrow(synth))){
      d_pred = pred_seed$individual[i,]  
      h = sqrt(as.numeric(bw.silv(as.data.frame(d_pred))))
      synth[i,Y_lab] = y_seed[i] + as.numeric(abs(rnorm(1, 0 ,h)) * sign(pred_seed$aggregate[i] - pred$aggregate[i])) 
    }
  }
  # method 4 : bruitage GN sur Y (pour comparer les méthodes au-dessus)
  else if (method == 4){  
    for (i in seq(nrow(synth))){
      sig = sd(Y_w, na.rm = T)  
      synth[i,Y_lab] =  y_seed[i] + as.numeric(abs(rnorm(1,0,sig*pert)) * sign(pred_seed$aggregate[i] - pred$aggregate[i]))
    }
  }
  # method 5 : classic bootstrap dans la distribution des erreurs de prediction
  else if (method == 5){  
    for (i in seq(nrow(synth))){
      eps = y_seed[i] - pred_seed$individual[i,]  
      b_eps = sample(eps,1,replace=T)
      synth[i,Y_lab] = y_seed[i]  +  as.numeric(abs(b_eps) * sign(pred_seed$aggregate[i] - pred$aggregate[i]))
    }
  }
  # method 6 : bruitage par noyau de la distribution de prediction signé
  else if (method == 6){  
    for (i in seq(nrow(synth))){
      d_pred = pred_seed$individual[i,]  
      h = sqrt(as.numeric(bw.silv(as.data.frame(d_pred))))
      synth[i,Y_lab] = y_seed[i]  + as.numeric(abs(rnorm(1,0,h)) *sign(pred_seed$aggregate[i] - pred$aggregate[i]))
    }
  }
  # method 7 : perturbation du Y par noyau calibré sur le Y 
  else if (method == 7){  
    for (i in seq(nrow(synth))){
      h = sqrt(as.numeric(bw.silv(as.data.frame(Y_w))))
      synth[i,Y_lab] =  y_seed[i]  + as.numeric(abs(rnorm(1,0,h)) * sign(pred_seed$aggregate[i] - pred$aggregate[i]))
    }
  }
  # method 4 : bruitage GN sur Y (pour comparer les méthodes au-dessus)
  else if (method == 8){  
    for (i in seq(nrow(synth))){
      sig = sd(Y_w, na.rm = T)  
      synth[i,Y_lab] =  y_seed[i] + rnorm(1,0,sig*pert)
    }
  }
  # method 5 : classic bootstrap dans la distribution des erreurs de prediction
  else if (method == 9){  
    for (i in seq(nrow(synth))){
      eps = y_seed[i] - pred_seed$individual[i,]  
      b_eps = sample(eps,1,replace=T)
      synth[i,Y_lab] = y_seed[i]  +  as.numeric(abs(b_eps*pert) * sign(pred_seed$aggregate[i] - pred$aggregate[i]))
    }
  }
  synth[, Y_lab]
}





## Clustering


cluster=function(train){
  
  DMC = densityMclust(train, plot=F,G=1:10)
  #cbind(train,"cluster" = DMC$classification)
  DMC$classification
}




## Mixture Model

### Gaussian Mixture Model


MM_GaussMM = function(X,m,N){                                                      # m peut être un vecteur, l'algo va alors chercher à optimiser le nb de composantes  
  DMC = densityMclust(X, plot=F,G=1:m)
  ech_MC = cbind(X,"cluster" = DMC$classification)
  DS_GMM = as.data.frame(mclust::sim(modelName=DMC$modelName,parameters=DMC$parameters,n=N))
  DS_GMM$group=NULL
  colnames(DS_GMM) = colnames(X)
  DS_GMM
} 



## KDE


### Classical Smoothed Bootstrap



# Attention, la matrice de lissage ne tient pas compte des poids.
GOLIATH_CSB = function(X,N, w,seed=sample(seq(nrow(X)),N, replace=T, prob=w), kern,hmult){
  # if (kern=="gaussian") {
  #   data.frame(rmvg(N, X, weights =  w))
  # }else {
  
  synth = kernelboot(data=X,function(dat) dat,kernel=kern,weights = w,R=5,adjust = hmult)
  synth$boot.samples[seed,]
  # }
}



### Full bandwidth matrices (Non-diagonal) Smoothed Bootstrap


GOLIATH_FGSB = function(X,N, w,seed=sample(seq(nrow(X)),N, replace=T, prob=w)){
  H = ks::kde(x=X,w=w*nrow(X))$H       # bw.silv(X)
  synth = data.frame(rmvg(N, X, bw=H, weights =  w))          # X = seed
  synth[seed,]
}



### ROSE


GOLIATH_ROSE = function(X,N, w,seed=sample(seq(nrow(X)),N, replace=T, prob=w),hmult=1, w_sd=rep(1,nrow(X))){
  if (is.null(hmult)){hmult=1}
  if (is.null(w_sd)){w_sd=rep(1,nrow(X))}
  n = nrow(X)
  q = ncol(X)
  kern <- seed
  cons.kernel <- (4/((q+2)*n))^(1/(q+4))
  if(q!=1){
    H <- hmult*cons.kernel*diag(apply(X, 2, function(z) sqrt(wtd.var(z,w_sd))), q)
  }else {
    H <- hmult*cons.kernel*sqrt(wtd.var(X,w))}
  Xnew.num <- matrix(rnorm(N*q), N, q)%*%H
  data.frame(Xnew.num + X[kern,])
}


### Gaussian Noise


GOLIATH_GN = function(X,N, w,seed=NULL,pert=0.5,w_sd=rep(1,nrow(X))){
  if (is.null(pert)){pert=0.1}
  n = nrow(X)
  q = ncol(X)
  kern <- seed
  if(q!=1){
    H <- pert*diag(apply(X, 2, function(z){sqrt(wtd.var(z,w_sd))}), q)
  }else {
    H <- pert*sqrt(wtd.var(X,w_sd))}
  Xnew.num <- matrix(rnorm(N*q), N, q)%*%H
  data.frame(Xnew.num + X[kern,])  
}



### Non-classical kernels



GOLIATH_NCSB = function(X,N, w,seed=sample(seq(nrow(X)),N, replace=T, prob=w), hmult=1,NCSB_param=list('min'=rep(NA,ncol(X)),'max'=rep(NA,ncol(X)))){
  if (is.null(NCSB_param$min)){NCSB_param$min=rep(NA,ncol(X))}
  if (is.null(NCSB_param$max)){NCSB_param$max=rep(NA,ncol(X))}
  Synth = X[seed,]
  
  n = nrow(X)
  q = ncol(X)
  cons.kernel <- (4/((q+2)*n))^(1/(q+4)) / (4/((1+2)*n))^(1/(1+4)) # constante de silverman
  hmult = hmult * cons.kernel
  
  for (j in (seq(q))){
    if (class(X[,j]) == "numeric"){
      mini = NCSB_param$min[j]
      maxi = NCSB_param$max[j]
      if (is.na(mini) & is.na(maxi)){
        print(paste0("Gaussian kernel for the variable ", as.character(colnames(X)[j])))
        #h = npudens(X[,j])$bw
        h = hmult * sqrt(bw.silv(subset(X,select=j)))
        Synth[,j] = sapply(Synth[,j],function(x) {x + rnorm(1)*h})
      }
      else if(is.na(mini) & !is.na(maxi)){
        print(paste0("Negative Gamma kernel for the variable ", as.character(colnames(X)[j])))
        h = hmult * npuniden.boundary(X=-(X[,j]-maxi),kertype = "gamma",a=0)$h
        Synth[,j] = sapply(Synth[,j],function(x) {-rgamma(1,shape = -(x-maxi)/h+1,scale=h)+maxi })
      }
      else if(!is.na(mini) & is.na(maxi)){
        print(paste0("Gamma kernel for the variable ", as.character(colnames(X)[j])))
        h = hmult  * npuniden.boundary(X=X[,j]-mini,kertype = "gamma",a=0,bwmethod = c("cv.ml"))$h
        Synth[,j] = sapply(Synth[,j],function(x) {rgamma(1,shape = (x-mini)/h+1,scale=h)+mini})
        # h = npuniden.boundary(X=X[,j],kertype = "gamma",a=mini,bwmethod = c("cv.ml"))$h
        # Synth[,j] = sapply(Synth[,j],function(x) {rgamma(1,shape = (x/h+1) * hmult,scale=h * hmult)})
      }
      else if(!is.na(mini) & !is.na(maxi)){
        if (mini==0 & maxi==0){
          print(paste0("Dirac kernel (classic bootstrap) for the variable ", as.character(colnames(X)[j])))
        } 
        else if (mini==0 & maxi==1){
          print(paste0("Beta kernel for the variable ", as.character(colnames(X)[j])))
          h = hmult * npuniden.boundary(X=X[,j],kertype = "beta1",a=0,b=1)$h
          Synth[,j] = sapply(Synth[,j],function(x) {rbeta(1,x/h+1,(1-x)/h+1)})
        } 
        else {
          print(paste0("Truncated Gaussian kernel for the variable ", as.character(colnames(X)[j])))
          h = hmult * npuniden.boundary(X=X[,j],kertype = "gaussian1",a=mini,b=maxi)$h
          u = data.frame('bb'=rep(0,nrow(Synth)),'bh'=rep(0,nrow(Synth)))
          u$bb = sapply(Synth[,j],function(x) {pnorm(mini,x,h)})
          u$bh = sapply(Synth[,j],function(x) {pnorm(maxi,x,h)})
          u$tir = apply(u,1,function(x){runif(1,x[1],x[2])})
          Synth[,j] = qnorm(u$tir,Synth[,j],h)
        }
      }
    } 
    else if (class(X[,j]) == "integer"){
      print(paste0("Binomial kernel for the variable ", as.character(colnames(X)[j])))
      h = hcvd.fun(X[,j],ker="bino")$hcv
      Synth[,j] = sapply(Synth[,j],function(x) {rbinom(1,x+1,(x+h*hmult)/(x+1))})
    }
  }
  Synth
} 


# X = exo
# j = 1
# mini=min(data$lread)
# h = npuniden.boundary(X=X[,j]-mini,kertype = "gamma",a=0,bwmethod = c("cv.ml"))$h
# hmult=1
# x=exo$lread[1];x
# gam = rgamma(1000,shape = (x-mini)/(h*hmult)+1,scale=h*hmult) + mini
# hist(gam,breaks=100)

## Nearest-Neighbors Bootstrap

### SMOTE-Family 


standardiz = function(v,m,s){
  (v-m)/s
}




GOLIATH_SMOTEf = function(X,N,w,seed=sample(seq(nrow(X)),N, replace=T, prob=w),dist="euclidean", k, pi_ij='Unif',fun_interp='Unif',Beta_alpha=1,Beta_beta=1,LogNorm_mu=0,LogNorm_sig=0.5){
  if (is.null(dist)){dist='euclidean'}
  if (is.null(pi_ij)){pi_ij='Unif'}
  if (is.null(fun_interp)){fun_interp='Unif'}
  if (is.null(Beta_alpha)){Beta_alpha=1}
  if (is.null(Beta_beta)){Beta_beta=1}
  if (is.null(LogNorm_mu)){LogNorm_mu=0}
  if (is.null(LogNorm_sig)){LogNorm_sig=0.5}
  if (is.null(seed)){sample(nrow(X),N,prob = w, replace=TRUE)}
  # handling of confused neighbors (to avoid zero distance)
  XX = cbind(X,w)
  XX = aggregate(w ~ .,data=XX,FUN=sum)
  w = XX$w
  XX$w = NULL
  # definition k-nn
  Xst = apply(X,2,function(z) {standardiz(z, mean(z),sd(z))})
  nn_array = kNN_array(Xst,dist,k)#kNN(XX,k)                                                           
  d = nn_array$dist                                                             
  d = 1/d
  d = d / rowSums(d)
  nn_array = nn_array$id
  row.names(d) = seed
  row.names(nn_array) = seed
  # drawing
  kern = seed
  kern_nn = nn_array[kern,]
  # interpolation
  if (pi_ij=='Unif'){
    nn = apply(kern_nn, 1, function(x) sample(x,1,replace=TRUE,prob=rep(1/k,k)))
  } else if (pi_ij=='Dist') {                                                   # inverse distance weighting
    nn = as.numeric(sample(kern_nn[1,],1, replace=TRUE,prob = d[kern[1],]))
    for (i in seq(2,N)){nn = c(nn, as.numeric(sample(kern_nn[i,],1, replace=TRUE,prob = d[kern[i],])))} # à optimiser avec un apply
  } # ajouter poids gaussien : https://towardsdatascience.com/make-your-knn-smooth-with-gaussian-kernel-7673fceb26b9
  
  synth = X[kern,]*0
  if (fun_interp=='Unif'){lambda = runif(N)}
  else if (fun_interp == 'Beta') {lambda = rbeta(N,Beta_alpha,Beta_beta)}
  else if (fun_interp == 'LogNorm') {lambda = rlogitnorm(N,LogNorm_mu,LogNorm_sig)}
  
  synth=X[kern,] + lambda * (X[nn,]-X[kern,])
  synth
}


### (Space-)Extended Nearest-Neighbors Bootstrap : E-SMOTE 




GOLIATH_ExSMOTEf = function(X,N,w,seed=sample(seq(nrow(X)),N, replace=T, prob=w), k,dist="euclidean", pi_ij='Unif',fun_interp='Unif',Beta_alpha=1,Beta_beta=1,LogNorm_mu=0,LogNorm_sig=0.5,Ex_sig=0.5){
  if (is.null(dist)){dist='euclidean'}
  if (is.null(pi_ij)){pi_ij='Unif'}
  if (is.null(fun_interp)){fun_interp='Unif'}
  if (is.null(Beta_alpha)){Beta_alpha=1}
  if (is.null(Beta_beta)){Beta_beta=1}
  if (is.null(LogNorm_mu)){LogNorm_mu=0}
  if (is.null(LogNorm_sig)){LogNorm_sig=0.5}
  if (is.null(Ex_sig)){Ex_sig=0.5}
  if (is.null(seed)){sample(nrow(X),N,prob = w, replace=TRUE)}
  # handling of confused neighbors (to avoid zero distance)
  XX = cbind(X,w)
  XX = aggregate(w ~ .,data=XX,FUN=sum)
  w = XX$w
  XX$w = NULL
  # definition k-nn
  Xst = apply(X,2,function(z) {standardiz(z, mean(z),sd(z))})
  nn_array = kNN_array(Xst,dist,k)#kNN(XX,k)                                                           
  d = nn_array$dist                                                             
  d = 1/d
  d = d / rowSums(d)
  nn_array = nn_array$id
  row.names(d) = seed
  row.names(nn_array) = seed
  # drawing
  kern = seed
  kern_nn = nn_array[kern,]
  # interpolation
  if (pi_ij=='Unif'){
    nn = apply(kern_nn, 1, function(x) sample(x,1,replace=TRUE,prob=rep(1/k,k)))
  } else if (pi_ij=='Dist') {                                                   # inverse distance weighting
    nn = as.numeric(sample(kern_nn[1,],1, replace=TRUE,prob = d[kern[1],]))
    for (i in seq(2,N)){nn = c(nn, as.numeric(sample(kern_nn[i,],1, replace=TRUE,prob = d[kern[i],])))} # à optimiser avec un apply
  } # ajouter poids gaussien : https://towardsdatascience.com/make-your-knn-smooth-with-gaussian-kernel-7673fceb26b9
  
  synth = X[kern,]*0
  if (fun_interp=='Unif'){lambda = runif(N)}
  else if (fun_interp == 'Beta') {lambda = rbeta(N,Beta_alpha,Beta_beta)}
  else if (fun_interp == 'LogNorm') {lambda = rlogitnorm(N,LogNorm_mu,LogNorm_sig)}
  # tirage d'une graine selon smotefamily
  synth_kern = X[kern,] + lambda * (X[nn,]-X[kern,])
  # perturbation de la graine selon GN 
  n = nrow(X)
  q = ncol(X)
  H = Ex_sig * sqrt(bw.silv(X))
  H[is.na(H)]=0
  #   if(q!=1){
  #   H <- Ex_sig*diag(apply(X, 2, sd), q)
  # }else {
  #   H <- Ex_sig*sd(X)}
  Xnew.num <- matrix(rnorm(N*q), N, q)%*%H
  data.frame(Xnew.num + synth_kern)
}



kNN_array = function(X,d="euclidean",k){
  n = nrow(X)
  nn_array = data.frame(matrix(rep(seq(n-1),n),n,n-1,byrow = T))
  distX = philentropy::distance(X,method=d)
  nn_array_id=t(apply(distX,1,order))
  nn_array_d=t(apply(distX,1,sort))
  nn_array_id=nn_array_id[,-1][,1:k]
  nn_array_d=nn_array_d[,-1][,1:k]
  list('id'=nn_array_id,'dist'=nn_array_d)
}
