\name{RVinePDF}        
\alias{RVinePDF}                

\title{PDF of an R-Vine Copula Model}

\description{
This function calculates the probability density function of a d-dimensional R-vine copula.
}

\usage{
RVinePDF(newdata, RVM)
}

\arguments{
  \item{newdata}{An N x d data matrix that specifies where the density shall be evaluated.}
  \item{RVM}{An \code{\link{RVineMatrix}} object including the structure and the pair-copula families and parameters.}
}


\details{
The density of a \eqn{d}-dimensional R-vine copula with \eqn{d-1} trees and corresponding edge sets \eqn{E_1,...,E_{d-1}} is given by
\deqn{
\prod_{\ell=1}^{d-1} \prod_{e\in E_\ell}
c_{j(e),k(e)|D(e)}\left(F(u_{j(e)}|\boldsymbol{u}_{D(e)}),F(u_{k(e)}|\boldsymbol{u}_{D(e)})|\boldsymbol{\theta}_{j(e),k(e)|D(e)}\right),
}{
=\prod_{k=1}^{d-1} \prod_{e\in E_k}
c_{j(e),k(e)|D(e)}(F(u_{j(e)}|u_{D(e)}),F(u_{k(e)}|u_{D(e)})|\theta_{j(e),k(e)|D(e)}),
}
where \eqn{\boldsymbol{u}=(u_{1},...,u_{d})^\prime\in[0,1]^d}{u=(u_{1},...,u_{d})'\in[0,1]^d}.
Further \eqn{c_{j(e),k(e)|D(e)}} denotes a bivariate copula density associated to an edge \eqn{e} and with parameter(s) \eqn{\boldsymbol{\theta}_{j(e),k(e)|D(e)}}{\theta_{j(e),k(e)|D(e)}}.
Conditional distribution functions such as \eqn{F(u_{j(e)}|\boldsymbol{u}_{D(e)})}{F(u_{j(e)}|u_{D(e)})} are obtained recursively using the relationship
\deqn{
h(u|\boldsymbol{v},\boldsymbol{\theta}) := F(u|\boldsymbol{v}) =
\frac{\partial C_{uv_j|\boldsymbol{v}_{-j}}(F(u|\boldsymbol{v}_{-j}),F(v_j|\boldsymbol{v}_{-j}))}{\partial F(v_j|\boldsymbol{v}_{-j})},
}{
h(u|v,\theta) := F(u|v) =
d C_{uv_j|v_{-j}}(F(u|v_{-j}),F(v_j|v_{-j}))/d F(v_j|v_{-j}),
}
where \eqn{C_{uv_j|\boldsymbol{v}_{-j}}}{C_{uv_j|v_{-j}}} is a bivariate copula distribution function with parameter(s) \eqn{\boldsymbol{\theta}}{\theta}
and \eqn{\boldsymbol{v}_{-j}}{v_{-j}} denotes a vector with the \eqn{j}-th component \eqn{v_j} removed.
The notation of h-functions is introduced for convenience. For more details see Dissmann et al. (2013).

The function is actually just a wrapper to \code{\link{RVineLogLik}}.
}

\references{
Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.
}

\author{Thomas Nagler}

\seealso{\code{\link{BiCopHfunc}}, \code{\link{RVineMatrix}}, \code{\link{RVineMLE}}, \code{\link{RVineAIC}}, \code{\link{RVineBIC}}}

\examples{
# define 5-dimensional R-vine tree structure matrix
Matrix <- c(5, 2, 3, 1, 4,
            0, 2, 3, 4, 1,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 1)
Matrix <- matrix(Matrix, 5, 5)

# define R-vine pair-copula family matrix
family <- c(0, 1, 3, 4, 4,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 3,
            0, 0, 0, 0, 0)
family <- matrix(family, 5, 5)

# define R-vine pair-copula parameter matrix
par <- c(0, 0.2, 0.9, 1.5, 3.9,
         0, 0, 1.1, 1.6, 0.9,
         0, 0, 0, 1.9, 0.5,
         0, 0, 0, 0, 4.8,
         0, 0, 0, 0, 0)
par <- matrix(par, 5, 5)

# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 5, 5)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3", "V4", "V5"))
                  
# compute the density at (0.1, 0.2, 0.3, 0.4, 0.5)
RVinePDF(c(0.1, 0.2, 0.3, 0.4, 0.5), RVM)
}
