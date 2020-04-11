standardize<- function(X)
{
  X = as.matrix(X)
  n = dim(X)[1]
  p = dim(X)[2]
  
  X = scale(X, center = TRUE, scale=FALSE)
  X = scale(X,center=FALSE, scale=sqrt(apply(X^2,2,sum)))

  # m = apply(X,2,mean)
  # st = sqrt(apply(X^2,2,sum));
  # st_mat = matrix(st, nrow = length(st), ncol = dim(X)[2], byrow=FALSE)
  # X2 = X / st_mat
  
}


 hidden_convariate_linear <- function(F,Y,k = 10 ,lambda = 5,lambda2 = 1, lambda3 =1 , iter = 1000) {
   #
   #
   # function [Z,B,U,o,error,error1,error2,dz,db,du] = hidden_covariate_linear(F,Y,k,lambda,lambda2,lambda3,iter);
   # input:
   #      F: a matrix nxd of known covariates, where n is the number of
   #      subjects and d is the number of known covariates. *must be standardize (columns have 0 mean and constant SS).
   #      Y: a matrix of nxg of expression data (must be standardized (columns
   #      scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
   #      k: number of inferred hidden components (k is an integer)
   #      lambda, lambda2, lambda3 are model parameters
   #      (optional) iter: number of iterations (default = 100);
   #
   #      note: k>0, lambda>0, lambda2>0, lambda3>0 must be set by the user based on the data at hand. one can set these values
   #      using cross-validation, by evaluating the "performance" of the  resulting residual data on a desired task.
   #      typically, if lambda>5, then hidden factors match the known covariates closely. 
   #
   # objective:
   #
   # this function solves the following problem:
   # argmin_{Z,B,U}   ||Y-Z*B||_2 + \lambda*||Z-F*U||_2 + \lambda2*||B||_2 + \lambda_3||U||_2
   #
   # output:
   #      Z: matrix of hidden components, dimensionality: nxk
   #      B: matrix of effects of hidden components, dimensionality: kxg
   #      o: value of objective function on consecutive iterations.
   #
   # to use the residual data: Residual = Y - Z*B
   library(MASS)
   library(pracma)
   
   tol = 1e-6;
   
   U = matrix(0, nrow=dim(F)[2],k)
   Z = matrix(0, nrow=dim(F)[1],k)
   B= matrix(runif(dim(Z)[2]*dim(Y)[2]), nrow=dim(Z)[2], ncol=dim(Y)[2])
   F = as.matrix(F)
   
   n1 = dim(F)[1];
   d1 = dim(F)[2]
 
   n2 = dim(Y)[1]
   d2 = dim(Y)[2]
   
   if(n1!=n2)    stop("number of rows in F and Y must agree")
   
   
   if (k<1 | lambda<1e-6 | lambda2<1e-6 | lambda3<1e-6 ) {
     stop("lambda, lambda2, lambda3 must be positive and/or k must be an integer");
   }
   
   o = vector(length=iter)
   
   for (ii in 1:iter) {
     o[ii] = sum((Y - Z%*%B)^2) + sum((Z -  F%*%U)^2)*lambda + (sum(B^2))*lambda2 + lambda3*(sum(U^2));
     Z = (Y %*% t(B) + lambda * F %*%U) %*% ginv(B %*% t(B) + lambda * diag(dim(B)[1]))
     B = mldivide(t(Z) %*% Z + lambda2 * diag(dim(Z)[2]), (t(Z) %*% Y))
     U = mldivide(t(F) %*% F * lambda + lambda3 * diag(dim(U)[1]), lambda * t(F) %*% Z)
     
     if(ii > 1 &&  (abs(o[ii]-o[ii-1])/o[ii]) < tol)  break
   }
         
   error =  sum((Y - Z%*%B)^2) / sum(Y^2)  + sum((Z - F%*%U)^2)/sum((F%*%U)^2)
   error1 = sum((Y - Z%*%B)^2) / sum(Y^2);
   error2 = sum((Z - F%*%U)^2) / sum((F%*%U)^2);
   
   dz = Z%*%(B%*%t(B) + lambda*diag(dim(B)[1]))-(Y%*%t(B) + lambda*F%*%U);
   db = (t(Z)%*%Z + lambda2*diag(dim(Z)[2]))%*%B - t(Z)%*%Y;
   du = (t(F)%*%F*lambda + lambda3*diag(dim(U)[1]))%*%U-lambda*t(F)%*%Z;
                    
           
   dataout = list(Z = Z, B = B, U = U)
   return(dataout)
                    
  
   
 }
