#' Implementing the Power Method
#'
#' Given a square matrix A, implement the Power Method to approximate the dominant or minimal eigenvector of A up to a specified error
#'
#' @param A square numeric matrix
#' @param dominant TRUE/FALSE choose to find the dominant or minimal eigenvector
#' @param v numeric vector with length equal to the number of rows of A
#' @param epsilon desired maximum of approximation
#' @param max_iter maximum number of iterations before outputting
#' @param plot TRUE/FALSE choose to plot the error versus iterations
#' @return power.method() returns a list which includes a vector of errors at each iteration, vectors of the approximate eigenvectors/eigenvalues at each iteration, and the dominant/minimal eigenvector and eigenvalue of A
#' @export
#' @examples A = rbind(c(1,3),c(4,5))
#' power.method(A)
#' eigen(A)
#'
power.method = function(A, dominant = TRUE, v = NULL, epsilon= 1e-06, max_iter = 100, plot = FALSE){

  if(!is.square.matrix(A)){
    stop("power.method() requires a square numeric matrix")
  }
  if(!is.null(v)){
    if(!is.vector(v) || !is.numeric(v)){
      stop("power.method() requires 'v' to be a numeric vector")
    }
    if(length(v) != nrow(A)){
      stop(" The vector `v` is not conformable to the matrix A in power.method()")
    }
  }

  # Finds the dominant or minimal eigenvalue and eigenvector using eigen
  if(dominant == TRUE){
    eigen_A = eigen(A)
    dom_eigenvalue = which.max(abs(eigen_A$values))
    dom_eigenvector = eigen_A$vectors[,dom_eigenvalue]
  }else{
    A = solve(A)
    eigen_A = eigen(A)
    dom_eigenvalue = which.max(abs(eigen_A$values))
    dom_eigenvector = eigen_A$vectors[,dom_eigenvalue]
  }

  # Initialise random eigenvector and correspondong eigenvalue
  b_old = c(runif(nrow(A),0,1))
  lambda_old = (A %*% b_old)[1]/b_old[1]

  # Calculate error of initial guess
  error = norm(abs(dom_eigenvector) - abs(b_old),type="2")


  # Initialise output vectors
  error_vec = c(error)
  vector_approx = list(b_old)
  value_approx = c(lambda_old)
  iter = 0

  # Iterate the eigenvalues/vectors using the Power Method
  while(error_vec[iter+1] > epsilon){
    iter = iter + 1
    if(iter == max_iter){
      warning("The power method has failed to converge to a vector within the desired error")
      break
    }

    b_new = A %*% b_old
    b_new = b_new/sqrt(sum(b_new^2))

    if(dominant){
      lambda_new = (A %*% b_new)[1]/b_new[1]
    }
    else{lambda_new = 1/((A %*% b_new)[1]/b_new[1])}

    error_vec[iter + 1] = norm(abs(dom_eigenvector) - abs(b_new),type="2")
    vector_approx[[iter + 1]] = b_new

    value_approx[iter+1] = lambda_new

    b_old = b_new
  }

  # Output contextually for minimal versus dominate
  if(dominant){
    dom_eigenvalue = (A %*% b_new)[1]/b_new[1]
    vector_approx = do.call(cbind, vector_approx)
    output = list(error_vec, vector_approx, value_approx, dom_eigenvector,dom_eigenvalue)
    names(output) = c("errors", "eigenvectors", "eigenvalues", "dom_eigenvector","dom_eigenvalue")
  }else{
    dom_eigenvalue = 1/((A %*% b_new)[1]/b_new[1])
    vector_approx = do.call(cbind, vector_approx)
    output = list(error_vec, vector_approx, value_approx, dom_eigenvector,dom_eigenvalue)
    names(output) = c("errors", "eigenvectors", "eigenvalues", "min_eigenvector","min_eigenvalue")
  }


  # Produce plot of the error versus iteration
  if(plot){
    plot(seq(0,iter),error_vec,type="l",xlab="Iteration Number", ylab="Error")
  }

  return(output)

}
