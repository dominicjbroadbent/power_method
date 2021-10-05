#' Implementing the Power Method
#'
#' Given a square matrix A, implement the Power Method to approximate the dominant eigenvector of A up to a specified error
#'
#' @param A square numeric matrix
#' @param v numeric vector with length equal to the number of rows of A
#' @param epsilon desired maximum of approximation
#' @param max_iter maximum number of iterations before outputting
#' @param plot TRUE/FALSE plot or not
#' @return power.method() returns a list which includes a vector of errors at each iteration, a vector of the approximate eigenvector at each iteration and the dominant eigenvector of A
#' @export
#' @examples A = rbind(c(1,3),4,5)
#' power.method(A)
#' eigen(A)
#'
power.method = function(A,v = NULL, epsilon= 1e-06, max_iter = 100, plot = FALSE){

  if(!is_square_matrix(A)){
    stop("power.method() requires a square numeric matrix")
  }
  if(!is.null(v)){
    if(!is.vector(v) || !is.numeric(v)){
      stop("power.method() requires 'v' to be a numeric vector")
    }
    if(len(v) != nrow(A)){
      stop(" The vector `v` is not conformable to the matrix A in power.method()")
    }
  }


  eigen_A = eigen(A)
  dominant = eigen_A$vectors[,which.max(abs(eigen_A$values))]

  b_old = c(runif(nrow(A),0,1))
  error = len(abs(dominant) - abs(b_old))

  error_vec = c(error)
  vector_approx = list(b_old)
  iter = 0

  while(error_vec[iter+1] > epsilon){
    iter = iter + 1
    if(iter == max_iter){
      warning("The power method has failed to converge to a vector within the desired error")
      break
    }

    b_new = A %*% b_old
    b_new = b_new/len(b_new)

    error_vec[iter + 1] = len(abs(dominant) - abs(b_new))
    vector_approx[[iter + 1]] = b_new

    b_old = b_new
  }

  vector_approx = do.call(cbind, vector_approx)
  output = list(error_vec, vector_approx,dominant)
  names(output) = c("errors", "approximations","dominant")

  if(plot){
    plot(seq(0,iter),error_vec,type="l",xlab="Iteration Number", ylab="Error")
  }

  return(output)

}
