#' Empirical effect sizes based on latent trait estimates
#'
#' Description FIXME!
#'
#' @param mod a multipleGroup object which estimated only 2 groups
#' @param focal_items a numeric vector indicating which items to include the tests. The
#'   default uses all of the items. Selecting fewer items will result in tests of
#'   'differential bundle functioning' when \code{DIF = FALSE}
#' @param npts number of points to use in the integration. Default is 50
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is
#'   used in conjunction with \code{npts}
#' @param Theta.focal an optional matrix of Theta values from the focal group to be evaluated. If not supplied
#'   the default values to \code{\link{fscores}} will be used in conjunction with the \code{...}
#'   arguments passed
#' @param DIF logical; return a data.frame of item-level imputation properties? If \code{FALSE},
#'   DBF and DTF statistics will be computed
#' @param ref.group either 1 or 2 to indicate which group is considered the 'reference' group. Default
#'   is 1
#' @param plot logical; plot the effects? FIXME!
#' @param ... additional arguments to be passed to \code{\link{fscores}}
#'
#' @author Adam Meade and Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Meade, A. (2010)...... #FIXME
#'
#' @export empirical_ES
#' @examples
#' \dontrun{
#'
#' #no DIF
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('dich', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('Ref', N), rep('Focal', N))
#'
#' mod <- multipleGroup(dat, 1, group = group,
#'    invariance = c(colnames(dat)[1:5], 'free_means', 'free_var'))
#' coef(mod, simplify=TRUE)
#'
#' empirical_ES(mod)
#' empirical_ES(mod, DIF=FALSE)
#'
#' empirical_ES(mod, plot=TRUE)
#' empirical_ES(mod, plot=TRUE, DIF=FALSE)
#'
#' ###---------------------------------------------
#' #DIF
#' set.seed(12345)
#' a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
#' a2[10:15,] <- a2[10:15,] + rnorm(6, 0, .3)
#' d2[10:15,] <- d2[10:15,] + rnorm(6, 0, .3)
#' itemtype <- rep('dich', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a1, d1, N, itemtype)
#' dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('Ref', N), rep('Focal', N))
#'
#' mod <- multipleGroup(dat, 1, group = group,
#'    invariance = c(colnames(dat)[1:5], 'free_means', 'free_var'))
#' coef(mod, simplify=TRUE)
#'
#' empirical_ES(mod)
#' empirical_ES(mod, DIF=FALSE)
#'
#' empirical_ES(mod, plot=TRUE)
#' empirical_ES(mod, plot=TRUE, DIF=FALSE)
#'
#' }
empirical_ES <- function(mod, Theta.focal = NULL, focal_items = 1L:extract.mirt(mod, 'nitems'),
                 DIF = TRUE, npts = 61, theta_lim=c(-6,6), ref.group = 1, ...){
    stopifnot(extract.mirt(mod, 'nfact') == 1L)
    stopifnot(extract.mirt(mod, 'ngroups') == 2L)
    ref <- extract.group(mod, ref.group)
    focal <- extract.group(mod, ifelse(ref.group == 1, 2, 1))
    focal_select <- extract.mirt(mod, 'group') != extract.mirt(mod, 'groupNames')[ref.group]
    if(is.null(Theta.focal)){
        Theta <- fscores(mod, full.scores = TRUE, full.scores.SE = FALSE, ...)
        Theta.focal <- Theta[focal_select, , drop = FALSE]
    } else Theta.focal <- as.matrix(Theta.focal)
    if(sum(focal_select) == nrow(Theta.focal))
        stop('Theta elements do not match the number of individuals in the focal group')

    # add the rest here
    if(DIF){

        ret <- data.frame()
        if(plot){
            ## return plot object. mirt imports lattice so those functions can be used if you are comfortable
            plot <- lattice::xyplot()
            return(plot)
        }
    } else { #DTF/DBF

        ret <- data.frame()
        if(plot){
            ## return plot object. mirt imports lattice, so this can be used if you are comfortable
            plot <- lattice::xyplot()
            return(plot)
        }
    }
    ret
}
