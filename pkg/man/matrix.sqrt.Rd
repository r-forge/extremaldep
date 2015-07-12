\name{matrix.sqrt}
\alias{matrix.sqrt}
\alias{print.matrix.sqrt}

\title{Matrix square root.}

\description{Matrix square root using singular value decomposition (svd).}
  
\usage{
matrix.sqrt(A)
}

\arguments{
	\item{A}{A symmetric, positive definite matrix of size \eqn{n \times n}{n x n}.}
}

\details{
	Using svd, the matrix A can be written as 
	$$ A = U S V^T, $$
	where 
	\itemize{
		\item $U$ is called the left-singular vectors of $A$, which corresponds to the eigenvectors of $AA^T$, 
		\item $S$ has entries $S_{i,j}$ which are the singular values of $A$, 
		\item $V$ is called the right-singular vectors of $A$, which corresponds to the eigenvectors of $A^T A$.
	}
	The non-zero singular values of $A$ are the square roots of the non-zero eignevalues of both $AA^T$ and $A^T A$.	
	Due to the symmetry and positive definiteness of $A$, without loss of generality we can set $U=V$.
	Thus $A = V S V^T$, and the columns of $V$ are eignevectors of $A$ with $j-th$ eignevalue begin $s_j$.
	Define
	$$ S^{1/2} = diag((s_1^{1/2},\ldots,s_p^{1/2})) $$
	and hence we have 
	$$ B = V S^{1/2} V^T $$, a symmetric \eqn{n \times n}{n x n} matrix such that $A = BB$. 
}

\value{
	Return a \eqn{n \times n}{n x n} matrix.
}

\author{
	Simone Padoan, \email{simone.padoan@unibocconi.it},
	\url{http://faculty.bocconi.it/simonepadoan};
	Boris Beranger, \email{borisberanger@gmail.com}
}
\examples{
	A <- matrix(c(1,2,2,3),ncol=2,byrow=TRUE)
	print(A)
	B <- matrix.sqrt(A)
}
