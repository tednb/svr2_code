### epidish.R

### Summary: a reference-based function to infer the proportions of a priori known cell subtypes present in a sample representing a mixture of such cell-types. Inference proceeds via one of 3 methods (Robust Partial Correlations-RPC, Cibersort (CBS), Constrained Projection (CP)), as determined by user.
### Author: Andrew E Teschendorff (a.teschendorff@ucl.ac.uk)
### Date: 19th May 2016

#### INPUT arguments:
### 参考：ref.m: a matrix of reference "centroids", i.e. representative molecular profiles, for a number of cell subtypes.
### rows label molecular features (e.g. CpGs,...) and columns label the cell-type. IDs need to be provided as rownames and colnames, respectively.
### No missing values are allowed, and all values in this matrix should be positive or zero. For DNAm data, values should be beta-values.

### 研究对象：avdata.m: a data matrix with rows labeling the molecular features (should use same ID as in cent.m) 
### and columns labeling samples (e.g. primary tumour specimens). 
### No missing values are allowed and all values should be positive or zero. In the case of DNA methylation, these are beta-values.


### OUTPUT arguments

### (CP-mode):
### a list with the following entries:
### estF: the estimated cell fraction matrix
### ref: the reference centroid matrix used
### dataREF: the input data matrix over the probes defined in the reference matrix

### (CBS-mode):
### a list with the following entries:
### estF: the estimated cell fraction matrix
### nu: a vector of "best" nu-parameter for each sample
### ref: the reference centroid matrix used
### dataREF: the input data matrix over the probes defined in the reference matrix

### (RPC-mode):
### a list with the following entries:
### estF: the estimated cell fraction matrix
### ref: the reference centroid matrix used
### dataREF: the input data matrix over the probes defined in the reference matrix


epidish <- function(avdata.m,ref.m,method=c("RPC","CBS","CP"),maxit=50,nu.v=c(0.25,0.5,0.75)){

    if(method=="RPC"){
       out.o <- DoRPC(avdata.m,ref.m,maxit);
    }
    else if (method=="CBS"){
       out.o <- DoCBS(avdata.m,ref.m,nu.v);
    }
    else if (method=="CP"){
       out.o <- DoCP(avdata.m,ref.m);
    }
    else {
        print("Input a valid method!");
    }

    return(out.o);
}
    
### Reference-based methods

### RPC
DoRPC <- function(avdata.m,ref.m,maxit){
    require(MASS);
    # Map the row names of ref.m to the corresponding rows in avdata.m
    map.idx <- match(rownames(ref.m),rownames(avdata.m)); 
    # Extract the rows that have a match and store the results in data2.m and ref2.m
    rep.idx <- which(is.na(map.idx)==FALSE); # isn't NA indexes in map.idx
    data2.m <- avdata.m[map.idx[rep.idx],];
    ref2.m <- ref.m[rep.idx,];
    # Create a matrix "est.m" with dimensions equal to the number of columns of data2.m by the number of columns of ref2.m
    est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m));
    # Set the column and row names of est.m, col is cell type and row is sample name
    colnames(est.m) <- colnames(ref2.m);
    rownames(est.m) <- colnames(data2.m);
    # Loop through each column in data2.m
    for(s in 1:ncol(data2.m)){
      # Run a robust regression analysis (rlm) for each column in data2.m against all columns in ref2.m
      rlm.o <- rlm( data2.m[,s] ~ ref2.m ,maxit=maxit) #maxit: maximum number of iterations for the regression algorithm (rlm)
      # Store the regression coefficients in coef.v
      coef.v <- summary(rlm.o)$coef[2:(ncol(ref2.m)+1),1];
      # Ensure the coefficients are non-negative
      coef.v[which(coef.v<0)] <- 0;
      # Normalize the coefficients so they sum to 1
      total <- sum(coef.v);
      coef.v <- coef.v/total;
      # Store the normalized coefficients in est.m
      est.m[s,] <- coef.v;
    }
    # Return a list containing the estimated regression coefficients (est.m), the reference data (ref2.m), and the corresponding data from the avdata.m matrix (dataREF)
    return(list(estF=est.m,ref=ref2.m,dataREF=data2.m));
}



###  CIBERSORT

require(e1071); # Load the e1071 library, used for SVM analysis

map.idx <- match(rownames(ref.m),rownames(avdata.m)); # Match the gene names in ref.m and avdata.m
rep.idx <- which(is.na(map.idx)==FALSE); # Get the index of non-NA matches

data2.m <- avdata.m[map.idx[rep.idx],]; # Get the matched data from avdata.m
ref2.m <- ref.m[rep.idx,]; # Get the matched data from ref.m

est.lm <- list(); # Initialize a list to store the estimated coefficients
nui <- 1; # Initialize a variable to index the list of coefficients

for(nu in nu.v){ # Loop over the values in nu.v
  est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m)); # Initialize a matrix to store the estimated coefficients for each value of nu
  colnames(est.m) <- colnames(ref2.m); # Set the column names of est.m to the column names of ref2.m
  rownames(est.m) <- colnames(data2.m); # Set the row names of est.m to the column names of data2.m
  
  for(s in 1:ncol(data2.m)){ # Loop over the columns of data2.m
    svm.o <- svm(x=ref2.m, y=data2.m[,s], scale = TRUE, type="nu-regression", kernel ="linear", nu = nu); # Fit a SVM regression model with nu = nu
    coef.v <- t(svm.o$coefs) %*% svm.o$SV; # Compute the coefficients
    coef.v[which(coef.v<0)] <- 0; # Set negative coefficients to 0
    total <- sum(coef.v); # Compute the total of the coefficients
    coef.v <- coef.v/total; # Normalize the coefficients
    est.m[s,] <- coef.v; # Store the normalized coefficients in est.m
  }
  est.lm[[nui]] <- est.m; # Store the estimated coefficients for this value of nu in est.lm
  nui <- nui+1; # Increment the index of est.lm
}
### select best nu
  # Initialize a matrix to store the root mean square error (RMSE) for each column in avdata.m and each value of nu in nu.v
  rmse.m <- matrix(NA,nrow=ncol(avdata.m),ncol=length(nu.v));
  # Loop through each value of nu in nu.v
  for(nui in 1:length(nu.v)){
    # Calculate the reconstruction of data2.m using the estimated coefficients from the SVM model trained with the current value of nu
    reconst.m <- ref2.m %*% t(est.lm[[nui]]);
    # Loop through each column in avdata.m
    for(s in 1:ncol(avdata.m)){
      # Calculate the RMSE between the actual data and the reconstructed data using the current value of nu
      rmse.m[s,nui] <- sqrt(mean((data2.m[,s] - reconst.m[,s])^2));
      }
    # Print the index of the current value of nu in nu.v
    print(nui);
    }
  # Label the columns of rmse.m with the values of nu in nu.v
  colnames(rmse.m) <- nu.v;
  # Find the index of the minimum RMSE value for each column in avdata.m
  nu.idx <- apply(rmse.m,1,which.min);
  # Initialize a matrix to store the final estimated coefficients
  estF.m <- est.m;    
  # Loop through each row of est.m
  for(s in 1:nrow(estF.m)){
    # Replace each row of estF.m with the row of the est.lm matrix corresponding to the minimum RMSE value for the current column in avdata.m
    estF.m[s,] <- est.lm[[nu.idx[s]]][s,];
  }
  # Return the final estimated coefficients, the selected value of nu, the reference matrix (ref2.m), and the reference data (data2.m)
  return(list(estF=estF.m,nu=nu.v[nu.idx],ref=ref2.m,dataREF=data2.m));
    

### Houseman CP
  DoCP <- function(avdata.m, ref.m){
    
    # Load the quadprog library
    require(quadprog);
    
    ### Define D matrix
    # Number of columns in the reference matrix
    nCT <- ncol(ref.m);
    # Initialize the D matrix with size nCT x nCT and filled with NAs
    D <- matrix(NA, nrow=nCT, ncol=nCT);
    # Fill the D matrix with the dot product of each pair of columns in the reference matrix
    for(j in 1:nCT){
      for(k in 1:nCT){
        D[j,k] <- 2 * sum(ref.m[,j] * ref.m[,k])
      }
    }
    
    ### Define constraints
    # Initialize the A matrix with size (nCT+1) x nCT and filled with zeros
    A.m <- matrix(0, nrow=nCT+1, ncol=nCT);
    # Set the first row of the A matrix to -1
    A.m[1,] <- -1;
    # Set the diagonal elements of the A matrix to 1
    for(i in 1:nCT){
      A.m[1+i,i] <- 1;
    }
    # Transpose the A matrix
    A.m <- t(A.m);
    # Initialize the b0 vector with length nCT + 1
    b0.v <- c(-1, rep(0,nCT));
    
    ### Define d-vector and solve for each sample
    # Number of columns in the avdata matrix
    nS <- ncol(avdata.m);
    # Initialize the westQP matrix with size nS x nCT and filled with NAs
    westQP.m <- matrix(NA, ncol=ncol(ref.m), nrow=nS);
    # Set the column names of the westQP matrix to the column names of the reference matrix
    colnames(westQP.m) <- colnames(ref.m);
    # Set the row names of the westQP matrix to the column names of the avdata matrix
    rownames(westQP.m) <- colnames(avdata.m);
    # Get the indices of the rows in the reference matrix that match the rows in the avdata matrix
    match(rownames(ref.m), rownames(avdata.m)) -> map.idx;
    # Get the indices of the non-NA elements in the map.idx vector
    rep.idx <- which(is.na(map.idx)==FALSE);
    ### loop through each sample in the available data
    for(s in 1:nS){
      # extract the data for the current sample
      tmp.v <- avdata.m[,s];
      # calculate the d-vector for the current sample
      d.v <- as.vector(2*matrix(tmp.v[map.idx[rep.idx]],nrow=1) %*% ref.m[rep.idx,]);
      # solve the quadratic program for the current sample
      qp.o <- solve.QP(D,d.v,A.m,b0.v,meq=0);
      # store the solution for the current sample
      westQP.m[s,] <- qp.o$sol;
      # print the current sample number
      print(s);
    }
    
    ### return the estimated functional data and the reference data used
    return(list(estF=westQP.m,ref=ref.m[rep.idx,],dataREF=avdata.m[map.idx[rep.idx],]));
}
