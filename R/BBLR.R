## Bayesian Bridge
#' @title Bayesian Bridge Linear Regression
#' @description This function implements the Bayesian bridge linear regression, a statistical method for linear regression analysis in high-dimensional data sets.
#'
#' @param y A numeric vector of the response variable with length n.
#' @param X A numeric matrix with n rows (observations) and p columns (predictors).
#' @param step_sizes_tuning_iterations The number of iterations dedicated to tuning the MCMC sampling algorithm for optimal step sizes.
#' @param Tb The tuning block size, specifying how often the algorithm adjusts the step sizes.
#' @param burn The number of initial MCMC iterations to discard, allowing the algorithm to reach stationary distribution.
#' @param nmc The total number of MCMC iterations to perform post-tuning and burn-in, for final analysis.
#' @param thin Thinning interval for the MCMC sampling to reduce autocorrelation in the chain.
#' @param method.alpha Specifies the method to update the alpha parameter. It can be "fixed" for a constant alpha throughout the sampling, or "beta" to sample alpha from a Beta distribution.
#' @param alpha_sig Significance level used for computing credible intervals of the estimated coefficients.
#'
#' @return A list containing the mean estimated coefficients (BetaHat), their credible intervals (LeftCI, RightCI), the median of the coefficient estimates (BetaMedian), the estimated variance of the error term (Sigma2Hat), the estimated parameter of the bridge prior (LambdaHat), the estimated shape parameter of the bridge prior (AlphaHat), and samples from the posterior distributions of the coefficients (BetaSamples), error variance (Sigma2Samples), lambda (LambdaSamples), and alpha (AlphaSamples).
#' @export
#'
#' @examples
#' \donttest{
#' # Assuming y is your response variable and X is your matrix of predictors
#' result <- BBLR(y, X, tuning=1000, Tb=100, burn=500, nmc=5000, thin=5, method.alpha="beta")
#' }
BBLR <- function(y,X,step_sizes_tuning_iterations,Tb,burn,nmc,thin,method.alpha,alpha_sig = 0.05){
    N = step_sizes_tuning_iterations + burn + nmc # the length of sampling
    effsamp = nmc/thin # ultimate chain length we get
    n = nrow(X) # sample size
    p = ncol(X) # dimension

    ### initial values of the chain
    Beta = rep(0.5, p)
    sigmasq = 1
    lambda = 1
    alpha = 0.5

    ### storing the chain
    betaout = matrix(0, p, effsamp)
    sigmasqout = rep(0, effsamp)
    lambdaout = rep(0, effsamp)
    alphaout = rep(0, effsamp)

    ### defining variables for tuning
    step.beta = rep(0.75,p)
    step.alpha = 0.25

    acc.beta = rep(0,p)
    acc.alpha = 0

    ### some pre calc
    Xy = t(X)%*%y
    XX = t(X) %*% X

    ### sampling
    for (i in 1:N) {

        ## sampling beta
        Beta.new = rep(0,p)
        for (j in 1:p) {
            muj = (Xy[j] - sum(XX[j,-j]*Beta[-j]))/XX[j,j]
            sigmaj = sigmasq/XX[j,j]
            beta.star = Beta[j] + rnorm(1,0,step.beta[j])
            T1 = lambda*((abs(Beta[j]))^alpha - (abs(beta.star))^alpha) # prior
            T2 = (0.5/sigmaj)*( (Beta[j] - muj)^2 - (beta.star - muj)^2 ) # likelihood
            MH = min(1,exp(T1+T2))
            U = runif(1,0,1)
            if(U<=MH){Beta.new[j] = beta.star;acc.beta[j]=acc.beta[j]+1}
            else{Beta.new[j] = Beta[j];acc.beta[j]=acc.beta[j]}
        }
        Beta = Beta.new

        ## sample alpha
        if(method.alpha=="beta"){ #Beta(0.5,0.5)
            Z0 = tan(pi*(alpha-0.5))
            Z.new = Z0 + rnorm(1,0,step.alpha)
            alpha.new = 0.5+atan(Z.new)/pi
            T1 = (1+Z0^2)/(1+Z.new^2)
            T2 = dbeta(alpha.new,0.5,0.5,log = T) - dbeta(alpha,0.5,0.5,log = T)
            T3 = log(lambda)*p*(1/alpha.new - 1/alpha) + p*(lgamma(1+1/alpha)-lgamma(1+1/alpha.new)) + lambda*( sum((abs(Beta))^alpha) - sum((abs(Beta))^alpha.new))
            MH = min(1,T1*exp(T2+T3))
            Ua = runif(1,0,1)
            if(Ua<=MH){alpha = alpha.new;acc.alpha = acc.alpha + 1}
        }
        if(method.alpha=="fixed"){
            alpha = alpha
        }

        ## sampling lambda
        lambda = rgamma(1,0.1+p/alpha,0.1+sum((abs(Beta))^alpha))

        ## sampling sigmasq
        sigmasq = 1/rgamma(1,0.5*n,0.5*sum((y-as.vector(X%*%Beta))^2))

        ## print the progress
        if (i%%1000 == 0) {print(i)}

        ## tuning
        if (i<=step_sizes_tuning_iterations & i%%Tb == 0){
            ar.beta = acc.beta/Tb
            ar.alpha = acc.alpha/Tb
            for (il in 1:p) {
                if(ar.beta[il]>0.5){step.beta[il]=step.beta[il]*1.1}
                if(ar.beta[il]<0.3){step.beta[il]=step.beta[il]*0.91}
            }
            if(ar.alpha>0.5){step.alpha = step.alpha*1.1}
            if(ar.alpha<0.3){step.alpha = step.alpha*0.91}
            acc.beta = rep(0,p)
            acc.alpha = 0
        }

        ## store
        if (i > (burn+step_sizes_tuning_iterations) && (i-burn-step_sizes_tuning_iterations)%%thin == 0) {
            betaout[, (i - burn - step_sizes_tuning_iterations)/thin] = Beta
            lambdaout[(i - burn - step_sizes_tuning_iterations)/thin] = lambda
            alphaout[(i - burn - step_sizes_tuning_iterations)/thin] = alpha
            sigmasqout[(i - burn - step_sizes_tuning_iterations)/thin] = sigmasq
        }
    }

    ## inference
    pMean = apply(betaout, 1, mean)
    pMedian = apply(betaout, 1, stats::median)
    pSigma = mean(sigmasqout)
    pLambda = mean(lambdaout)
    pAlpha = mean(alphaout)

    left <- floor(alpha_sig * effsamp/2)
    right <- ceiling((1 - alpha_sig/2) * effsamp)
    BetaSort <- apply(betaout, 1, sort, decreasing = F)
    left.points <- BetaSort[left, ]
    right.points <- BetaSort[right, ]

    ## return results
    result = list(BetaHat = pMean, LeftCI = left.points, RightCI = right.points,
                  BetaMedian = pMedian, Sigma2Hat = pSigma, LambdaHat = pLambda, AlphaHat = pAlpha,
                  BetaSamples = betaout,  Sigma2Samples = sigmasqout,
                  LambdaSamples = lambdaout, AlphaSamples = alphaout)
    return(result)
}
