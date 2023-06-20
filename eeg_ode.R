

ssa.svd<-function(df,L,chans,eigs){
  ## chose based on stationarity - e.g. nutrlan or not not - eigen, additionally töplitz only for stationary
  
  if(ncol(df)>1){
    df.ssa<-ssa(df[,chans],L=L,neig=eigs,kind="mssa",svd.method="propack")
  }else{
    df.ssa<-ssa(df[,chans],L=L,neig=eigs,kind="1d-ssa",svd.method="propack")
  }
  
  return(df.ssa)
}

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

rsvd<-function(df,L,chans,eigs){
  H = matrix(nrow=(length(df)+L), ncol=L)
  for(k in seq(1:L)){
    H[-k, k] = df
  }
  #   H = H((L):end,:);
  #   H=transpose(H);
  # ny = size(X,2)
  # P = randn(ny,r+p)
  # Z = df*P
  #   for k=1:q
  #   Z = X*(X'*Z);
  #     end
  #     % Gram-Schmidt or QR-Algorithm
  # %     for j=1:n 
  # %         v=A(:,j); % restricted to rows I v=A(I,j);
  # %         for i=1:j-1
  # %             R(i,j)=Q(:,i); % Q unitary at start -
  # %             v=v-R(i,j)*Q(:,i);
  # %         end
  # %         R(j,j)=norm(v); % extension formulae v|I = sqrt(v*v)
  # %         Q(:,j)=v/R(j,j);
  # %     end
  #     [Q,R] = qr(Z,0);
  #     % Step 2: Compute SVD on projected Y=Qâ€™*X;
  #      Y = Q'*X;
  #          %     eig(Y*Y');
  # %     UY=real(Y);
  # %     eig(Y'*Y);
  #          %     V=real(Y);
  #          %     S=Y*UY/V;
  #          tic;
  #          [UY,S,V] = svd(Y,'econ');
  #          U = Q*UY; 
  
}


grouping.ssa<-function(df,g,nc){
  sig.comps  <- grouping.auto.wcor(df,nclust=nc,groups=g,method="complete")
  return(sig.comps)
  ## Reconstructin by Components,
  ## e.g. seasonality, sine-cosine frequency components (use period) and noise (by correlation matrix)
  #residuals(eeg_klt_re)
  # klt_dat<-x
  # for(i in 1:as.numeric(input$klt_clust)){
  #   klt_dat<-cbind(klt_dat,as.numeric(eeg_klt_re[[as.character(i)]]))
  #   colnames(klt_dat)[i+1]<-paste("PC-Signal",i)
  # }
  # klt_dat<-cbind(klt_dat,"PC-uncovered"=as.numeric(residuals(eeg_klt_re)))
  # klt_dat<-as.ts(klt_dat,frequency=200)
  # klt_dat[,-1][,1:(ncol(klt_dat)-(2+input$final_component))]
}

recons<-function(issa,gclust,stacked=F){
  if(!is.list(gclust)){
    print("Clusters must be list")
    return(issa)
  }
  # residuals first
  tmp=residuals(issa)
  colnames(tmp)[1]<-"Unexplained"
  ## test for full recons
  ## remember sum eigenvalues
  gcls=c()
  eig_rank=c()
  for (i in gclust) {
    gcls=c(gcls,i)
    # print(i)
    eig_rank=c(eig_rank,sum(issa$sigma[i]))
    
  }
  
  if(!all(1:max(gcls) %in% gcls)){
    sings=(1:max(gcls))[!(1:max(gcls) %in% gcls)]
    clustsigs=list()
    for(i in 1:length(gclust)){
      clustsigs[[i]]=gclust[[i]]}
    singsigs=as.list(sings)
    clustsigs=append(clustsigs,singsigs)
    sig.split <- reconstruct(issa,clustsigs,drop.attributes=F)
  }else{
    sig.split <- reconstruct(issa,gclust,drop.attributes=F)
  }
  
  # if incomplete clustering, reconstruct so
  # if(){
  # sig.split <- reconstruct(df,g,drop.attributes=T)
  # }else{
  #
  # }
  
  ## add signals according to eigvalues and return as cumsum timeseries
  if(stacked){
    ##wrong!
    for(i in 1:length(gclust)){
      # tmp2=tmp[,i]
      # for(j in 1:i){
      tmp2<- reconstruct(issa,list(unlist(gclust[1:i])),drop.attributes=F)$F1
      # tmp2=tmp2+sig.split[[order(eig_rank,decreasing=T)[j]]]
      # }
      colnames(tmp2)<-paste0("Clust",i)
      tmp=merge(tmp,tmp2)
    }
  }else{
    for(i in 1:length(sig.split)){
      tmp2=sig.split[[order(eig_rank,decreasing=T)[i]]]
      colnames(tmp2)<-paste0("Clust",order(eig_rank,decreasing=T)[i])
      tmp=merge(tmp,tmp2)
    }
  }
  return(tmp)
  
}






# kalman_filter<-function(df,autoarma=T){
#   if(autoarma){
#     buildFun <- function(x) {
#       res<-auto.arima(df,allowdrift = F)
#       ar_sum<-sum(grepl("ar",names(res$coef)))
#       ma_sum<-sum(grepl("ma",names(res$coef)))
#       if(ar_sum+ma_sum!=length(res$coef))print("ERROR!!!")
#       m <- dlmModARMA(ar=res$coef[1:ar_sum], ma=res$coef[(1+ar_sum):(length(res$coef))])
#       return(m)
#     }
#   }else{
#     buildFun <- function(x,p,q) {
#       res<-arima(df,order=c(p,0,q))
#       ar_sum<-sum(grepl("ar",names(res$coef)))
#       ma_sum<-sum(grepl("ma",names(res$coef)))
#       if(ar_sum+ma_sum!=length(res$coef))print("ERROR!!!")
#       m <- dlmModARMA(ar=res$coef[1:ar_sum], ma=res$coef[(1+ar_sum):(length(res$coef))])
#       return(m)
#     }
#   }
#   fit <- dlmMLE(df,parm = rep(0,8), build = buildFun)
#   
#   dlmExG <- buildFun(fit$par)
#   
#   filt_dat<-dlmFilter(df, dlmExG)
#   
#   # p_dat<-dlmSmooth(filt_dat)
#   
#   return(fil_dat$s[-1])
# }

# 
rssa <- function(x,
                L = (N + 1) %/% 2,
                neig = NULL,
                fmask = NULL, wmask = NULL,
                column.projector = "none", row.projector = "none",
                column.oblique = "identity", row.oblique = "identity",
                ...,
                kind = c("1d-ssa", "2d-ssa", "nd-ssa", "toeplitz-ssa", "mssa", "cssa"),
                circular = FALSE,
                svd.method = c("auto", "nutrlan", "propack", "svd", "eigen", "rspectra", "primme"),
                force.decompose = TRUE) {
  svd.method <- match.arg(svd.method)
  N <- dim(x)[1]
  # Squeeze the attributes
  xattr <- attributes(x)
  iattr <- NULL
  # Grab class separately. This way we will capture the inherit class as well
  xclass <- class(x)

  call <- match.call(); cargs <- as.list(call)[-1]
  ## wmask is special and will be treated separately later
  cargs$wmask <- NULL
  ecall <- do.call("call", c("ssa", lapply(cargs, eval, parent.frame())))

  ## Provide some sane defaults, e.g. complex inputs should default to cssa
  if (missing(kind)) {
    if (is.complex(x))
      kind <- "cssa"
    else if (inherits(x, "mts") || inherits(x, "data.frame") || inherits(x, "list") || inherits(x, "series.list"))
      kind <- "mssa"
    else if (is.matrix(x))
      kind <- "2d-ssa"
    else if (is.array(x))
      kind <- "nd-ssa"
    else
      kind <- "1d-ssa"
  }
  kind <- match.arg(kind)

  # Do the fixups depending on the kind of SSA.
  weights <- NULL
  if (identical(kind, "1d-ssa") || identical(kind, "toeplitz-ssa")) {
    ## Nothing special here (yet!)
  } else if (identical(kind, "2d-ssa") || identical(kind, "nd-ssa")) {
    # 2d-SSA is just a special case of nd-ssa
    if (length(dim(x)) == 2)
      kind <- c("2d-ssa", "nd-ssa")
    else
      kind <- "nd-ssa"
  } else if (identical(kind, "mssa")) {
    ## Nothing special here (yet!)
  } else if (identical(kind, "cssa")) {
    ## Nothing special here (yet!)
  } else {
    N <- -1;
    fmask <- NULL
    stop("invalid SSA kind")
  }

  if (!identical(column.projector, "none") || !identical(row.projector, "none")) {
    # Add `pssa` class if appropriate implementation exists

    if (!any(kind %in% c("1d-ssa", "2d-ssa", "nd-ssa"))) {
      stop("SSA with projection is not implemented for such SSA kind yet")
    }

    kind <- c("pssa", paste("pssa", kind, sep = "-"), kind)
  }

  if (!identical(column.oblique, "identity") || !identical(row.oblique, "identity")) {
    # Add `wossa` class if appropriate implementation exists

    if (!any(kind %in% c("1d-ssa", "2d-ssa", "nd-ssa", "pssa"))) {
      stop("SSA with weights is not implemented for such SSA kind yet")
    }

    # TODO: Accept only row.oblique or column.oblique
    if (identical(column.oblique, "identical") || identical(row.oblique, "identical")) {
      stop("Both column.oblique and row.oblique must be numeric")
    }

    kind <- c("wossa", paste("wossa", kind, sep = "-"), kind)
  }

  # Normalize the kind to be used
  kind <- gsub("-", ".", kind, fixed = TRUE)

  # Create information body
  this <- list(call = call, ecall = ecall,
               kind = kind,
               svd.method = svd.method)

  # Create data storage
  this <- Rssa:::.create.storage(this)

  # Save the names of the essential fields
  this$fields <- c("F",
                   "wmask", "fmask", "weights", "circular",
                   "Fattr", "Fclass", "Iattr",
                   "column.projector", "row.projector",
                   "column.oblique", "row.oblique")

  # Make this S3 object
  class(this) <- c(kind, "ssa")

  ## Perform additional init steps, if necessary. We cannot simply eval .init in
  ## the current environment because we're using S3 dispatch at the same
  ## time... UseMethod uses NSE.
  ## NOTE: This will modify the *current* environment (local vars of the function)
  parent.env <- parent.frame()
  eval(Rssa:::.init.fragment(this))

  # Save attributes
  Rssa:::.set(this, "Fattr", xattr)
  Rssa:::.set(this, "Fclass", xclass)
  Rssa:::.set(this, "Iattr", iattr)

  # Deprecated stuff
  Rssa:::.deprecate(this, "lambda", "sigma")

  ## Window and series length should be ready by this moment
  this$length <- N
  this$window <- L

  ## Save series
  Rssa:::.set(this, "F", x)

  ## Save masks, weights and topology
  Rssa:::.set(this, "wmask", wmask)
  Rssa:::.set(this, "fmask", fmask)
  Rssa:::.set(this, "weights", weights)
  Rssa:::.set(this, "circular", circular)

  ## Store projectors
  Rssa:::.set(this, "column.projector", column.projector)
  Rssa:::.set(this, "row.projector", row.projector)

  ## Store oblique matrices
  Rssa:::.set(this, "column.oblique", column.oblique)
  Rssa:::.set(this, "row.oblique", row.oblique)

  # If 'neig' is specified, then we need to decompose
  if (!is.null(neig) && !force.decompose) {
    warning("`force.decompose = FALSE` is ignored because number of eigentriples is specified")
    force.decompose <- TRUE
  }

  # Determine the desired number of eigentriples, if necessary
  if (is.null(neig))
    neig <- Rssa:::.default.neig(this, ...)

  # Fix SVD method
  if (identical(svd.method, "auto"))
    svd.method <- Rssa:::.determine.svd.method(this, kind = kind, neig = neig, ...)

  this$svd.method <- svd.method

  # Decompose, if necessary
  if (force.decompose) {
    if (!is.null(weights) && all(weights == 0))
      stop("Nothing to decompose: the given field shape is empty")

    this <- decompose.ssa(this, neig = neig, ...);
  }

  this;
}

decompose.ssa <- function(x,
                          neig = NULL,
                          ...,
                          force.continue = FALSE) {
  ## Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0 &&
      !capable(x, "decompose.continue"))
    stop("Continuation of decomposition is not yet implemented for this method.")
  
  if (is.null(neig))
    neig <- Rssa:::.default.neig(x, ...)
  
  if (identical(x$svd.method, "svd")) {
    S <- svd(as.matrix(Rssa:::.get.or.create.trajmat.1d.ssa(x)), nu = neig, nv = neig)
    Rssa:::.set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "eigen")) {
    S <- eigen(tcrossprod(Rssa:::.get.or.create.trajmat.1d.ssa(x)), symmetric = TRUE)
    
    ## Fix small negative values
    S$values[S$values < 0] <- 0
    
    Rssa:::.set.decomposition(x,
                       sigma = sqrt(S$values[1:neig]),
                       U = S$vectors[, 1:neig, drop = FALSE])
  } else if (identical(x$svd.method, "nutrlan")) {
    S <- trlan.svd( Rssa:::.get.or.create.trajmat.1d.ssa(x), neig = neig, ...,
                    lambda = Rssa:::.sigma(x), U = Rssa:::.U(x))
    Rssa:::.set.decomposition(x, sigma = S$d, U = S$u)
  } else if (identical(x$svd.method, "propack")) {
    S <- propack.svd(Rssa:::.get.or.create.trajmat.1d.ssa(x), neig = neig, ...)
    Rssa:::.set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "rspectra")) {
    if (!requireNamespace("RSpectra", quietly = TRUE))
      stop("RSpectra package is required for SVD method `rspectra'")
    h <- Rssa:::.get.or.create.trajmat.1d.ssa(x)
    A <- function(x, args) ematmul(args, x)
    Atrans <- function(x, args) ematmul(args, x, transposed = TRUE)
    S <- RSpectra::svds(A, k = neig, Atrans = Atrans, dim = dim(h), args = h, ...)
    Rssa:::.set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "primme")) {
    if (!requireNamespace("PRIMME", quietly = TRUE))
      stop("PRIMME package is required for SVD method `primme'")
    h <- Rssa:::.get.or.create.trajmat.1d.ssa(x)
    pA <-function(x, trans) if (identical(trans, "c")) crossprod(h, x) else h %*% x
    S <- PRIMME::svds(pA, NSvals = neig, m = nrow(h), n = ncol(h), isreal = TRUE, ...)
    Rssa:::.set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else
    stop("unsupported SVD method")
  
  x
}


.get.or.create.hmat <- function(x) {
  .get.or.create(x, "hmat",
                 new.hmat(.F(x), L = x$window, circular = x$circular,
                          wmask = x$wmask, fmask = x$fmask, weights = x$weights,
                          fft.plan = .get.or.create.fft.plan(x)))
}

.get.or.create <- function(x, name, default) {
  (if (.exists.non.null(x, name))
    get(name, envir = .storage(x))
   else
     assign(name, default, envir = .storage(x), inherits = FALSE))
}

.set <- function(x, name, value)
  assign(name, value, envir = .storage(x), inherits = FALSE);

.exists <- function(x, name)
  exists(name, envir = .storage(x), inherits = FALSE);

.exists.non.null <- function(x, name) {
  ret <- FALSE
  
  if (exists(name, envir = .storage(x), inherits = FALSE)) {
    val <- get(name, envir = .storage(x))
    ret <- !is.null(val) && (typeof(val) != "externalptr" || !.is.extptrnull(val)) &&
      (!isS4(val) || !inherits(val, "externalptr") || !.is.extptrnull(val@.xData))
  }
  ret
}

fft.plan.1d <- function(N, L, circular = FALSE,
                        wmask = NULL, fmask = NULL, weights = NULL) {
  storage.mode(N) <- "integer"
  
  if (!is.null(wmask)) {
    storage.mode(wmask) <- "logical"
  }
  
  if (!is.null(fmask)) {
    storage.mode(fmask) <- "logical"
  }
  
  if (is.null(weights)) {
    if (!circular) {
      weights <- .hweights.default(N, L)
    } else {
      weights <- rep(L, N)
    }
  }
  storage.mode(weights) <- "integer"
  
  .Call("initialize_fft_plan", N, wmask, fmask, weights)
}
