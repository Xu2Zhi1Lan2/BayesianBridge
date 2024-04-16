## Bayesian Bridge Normal Mean Model
#' @title Bayesian Bridge Normal Mean Model
#' @description This function implements the Bayesian bridge normal mean model, which is an adaptation of the Bayesian bridge linear regression for scenarios where the response variable y follows a normal distribution with mean values influenced by the predictor variables. It is designed to handle high-dimensional datasets and employs a Markov Chain Monte Carlo (MCMC) sampling method for inference.
#'
#' @param y A numeric vector of the response variable, assumed to be normally distributed, with length p.
#' @param step_sizes_tuning_iterations The number of iterations dedicated to tuning the MCMC sampling algorithm for optimal step sizes.
#' @param Tb The frequency of tuning adjustments within the initial tuning phase.
#' @param burn The number of initial MCMC iterations to discard, allowing the algorithm to reach stationary distribution.
#' @param nmc The total number of MCMC iterations to perform post-tuning and burn-in, for final analysis.
#' @param thin Thinning interval for the MCMC sampling to reduce autocorrelation in the chain.
#' @param method.alpha Specifies the method to update the alpha parameter. It can be "fixed" for a constant alpha throughout the sampling, or "beta" to sample alpha from a Beta distribution.
#' @param alpha_sig The significance level for calculating credible intervals for the estimated parameters.
#'
#' @return Returns a list of outputs including the mean estimates of coefficients (BetaHat), credible intervals for these estimates (LeftCI, RightCI), the median of the coefficients (BetaMedian), the estimated error variance (Sigma2Hat), the estimated lambda parameter of the bridge prior (LambdaHat), the estimated alpha parameter of the bridge prior (AlphaHat), and samples from their posterior distributions (BetaSamples, Sigma2Samples, LambdaSamples, AlphaSamples).
#' @export
#'
#' @examples
#' \donttest{
#' # Assuming y is your normally distributed response variable
#' result <- BBNM(y, tuning=1000, Tb=100, burn=500, nmc=5000, thin=5, method.alpha="beta")
#' }

## Bayesian Bridge normal mean model
BBNM <- function(y,step_sizes_tuning_iterations,Tb,burn,nmc,thin,method.alpha,alpha_sig = 0.05){
    N = step_sizes_tuning_iterations + burn + nmc # the length of sampling
    effsamp = nmc/thin # ultimate chain length we get
    p = length(y) # dimension

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

    ### sampling
    for (i in 1:N) {

        ## sampling beta
        Beta.new = rep(0,p)
        for (j in 1:p) {
            muj = y[j]
            sigmaj = sigmasq
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
        if(method.alpha=="fixed"){ # fix alpha @ 0.5
            alpha = alpha
        }

        ## sampling lambda
        lambda = rgamma(1,0.1+p/alpha,0.1+sum((abs(Beta))^alpha))

        ## sampling sigmasq
        sigmasq = 1/rgamma(1,0.5*p,0.5*sum((y-Beta)^2))

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
