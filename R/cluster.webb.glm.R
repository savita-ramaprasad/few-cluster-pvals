#' @export
#' @importFrom stats as.formula coefficients glm predict quantile residuals update
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom sandwich sandwich estfun
#' @importFrom lmtest coeftest

cluster.webb.glm <- function (mod, dat, cluster, vars.boot = NULL, ci.level = 0.95, impose.null = TRUE,  boot.reps = 1000, report = TRUE, prog.bar = TRUE, output.replicates = FALSE) 
{
    if (mod$family[1] != "gaussian" | mod$family[2] != "identity") {
        stop("Use only with gaussian family models with a linear link")
    }
    if (output.replicates == TRUE & impose.null == TRUE) {
        stop("Recovering bootstrap replicates requires setting impose.null = FALSE")
    }
    form <- mod$formula
    variables <- all.vars(form)
    clust.name <- all.vars(cluster)
    used.idx <- which(rownames(dat) %in% rownames(mod$model))
    dat <- dat[used.idx, ]
    clust <- as.vector(unlist(dat[[clust.name]]))
    G <- length(unique(clust))
    ind.variables <- attr(mod$terms, "term.labels")
    "%w/o%" <- function(x, y) x[!x %in% y]
    dv <- variables %w/o% all.vars(update(form, 1 ~ .))
    ind.variables.data <- all.vars(update(form, 1 ~ .))
    ind.variables.names <- names(coefficients(mod))
    dat$dv.new <- mod$y
    form.new <- update(form, dv.new ~ .)
    fac <- c()
    for (i in 1:length(ind.variables.data)) {
        fac[i] <- is.factor(dat[, ind.variables.data[i]])
    }
    fac <- max(fac)
    if (fac == 1 & impose.null == TRUE) {
        cat("\n", "\n", "Note: null not imposed (factor variables are present).", 
            "\n")
        impose.null <- FALSE
    }
    cl <- function(dat, fm, cluster) {
        M <- length(unique(cluster))
        N <- length(cluster)
        K <- fm$rank
        dfc <- (M/(M - 1))
        uj <- apply(estfun(fm), 2, function(x) tapply(x, cluster, 
            sum))
        vcovCL <- dfc * sandwich(fm, meat. = crossprod(uj)/N)
        coeftest(fm, vcovCL)
    }
    se.clust <- cl(dat, mod, clust)[ind.variables.names, 2]
    beta.mod <- coefficients(mod)[ind.variables.names]
    w <- beta.mod/se.clust
    
    if (!is.null(vars.boot)){
    	ind.variables.boot <- match(vars.boot, ind.variables)
    	} else {
		ind.variables.boot <- 1:length(ind.variables)		
	}    
    if (impose.null == TRUE) {
        p.store <- c()
        w.store <- matrix(data = NA, nrow = boot.reps, ncol = length(ind.variables))
        if (attr(mod$terms, "intercept") == 1) {
            offset <- 1
        } else {
            offset <- 0
        }
        if (prog.bar == TRUE) {
            cat("\n")
        }        
        for (j in ind.variables.boot) {
            if (prog.bar == TRUE) {
                cat("Independent variable being bootstrapped: ", 
                  ind.variables[j], "\n")
            }
            form.null <- as.formula(paste(variables[1], "~", 
                paste(ind.variables[1:length(ind.variables) %w/o% 
                  j], collapse = " + ")))
            mod.null <- glm(form.null, data = dat, family = mod$family)
            null.resid <- residuals(mod.null)
            boot.dat <- dat
            wald.store <- c()
            if (prog.bar == TRUE) {
                pb <- txtProgressBar(min = 0, max = boot.reps, 
                  initial = 0, style = 3)
            }
            for (i in 1:boot.reps) {
                if (prog.bar == TRUE) {
                  setTxtProgressBar(pb, value = i)
                }
                weight <- sample(c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5)), 
                  size = G, replace = T, prob = rep(1/6, 6))[match(clust, unique(clust))]
                pseudo.resid <- null.resid * weight
                pseudo.dv <- predict(mod.null) + pseudo.resid
                boot.dat[, "dv.new"] <- pseudo.dv
                boot.mod <- glm(form.new, data = boot.dat, family = mod$family)
                se.boot <- cl(boot.dat, boot.mod, clust)[offset + 
                  j, 2]
                beta.boot <- coefficients(boot.mod)[offset + 
                  j]
                wald.store[i] <- beta.boot/se.boot
            }
            if (prog.bar == TRUE) {
                close(pb)
            }
            p.store[j] <- min(sum(abs(w[offset + j]) > abs(wald.store)), sum(abs(w[offset + j]) < abs(wald.store)))/boot.reps
            w.store[, j] <- wald.store
        } 
        p.store <- p.store[ind.variables.boot]
        w.store <- w.store[ ,ind.variables.boot]
        out <- matrix(p.store, ncol = 1)
        colnames(out) <- c("wild cluster BS p-value")
        rownames(out) <- vars.boot
        ci.lo = NULL
        ci.hi = NULL
        print.ci = NULL
        out.ci = NULL
        } else {
        if (prog.bar == TRUE) {
            cat("Wild Cluster bootstrapping w/o imposing null...", 
                "\n")
        }
        boot.dat <- dat
        w.store <- matrix(data = NA, nrow = boot.reps, ncol = length(ind.variables.names))
        rep.store <- matrix(data = NA, nrow = boot.reps, ncol = length(beta.mod))
        colnames(rep.store) <- ind.variables.names
        resid <- residuals(mod)
        if (prog.bar == TRUE) {
            pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, 
                style = 3)
        }
        for (i in 1:boot.reps) {
            if (prog.bar == TRUE) {
                setTxtProgressBar(pb, value = i)
            }
            weight <- sample(c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5)),  
              size = G, replace = T, prob = rep(1/6, 6))[match(clust, unique(clust))]
            pseudo.resid <- resid * weight
            pseudo.dv <- predict(mod) + pseudo.resid
            boot.dat[, "dv.new"] <- pseudo.dv
            boot.mod <- glm(form.new, data = boot.dat, family = mod$family)
            se.boot <- cl(boot.dat, boot.mod, clust)[, 2]
            beta.boot <- coefficients(boot.mod)
            w.store[i, ] <- (beta.boot - beta.mod)/se.boot
            rep.store[i, ] <- beta.boot
        }
        if (prog.bar == TRUE) {
            close(pb)
        }
        comp.fun <- function(vec2, vec1) {
            as.numeric(vec1 > vec2)
        }
        p.store.s <- t(apply(X = abs(w.store), FUN = comp.fun, 
            MARGIN = 1, vec1 = abs(w)))
        p.store <- 1 - (colSums(p.store.s)/dim(w.store)[1])
        crit.t <- apply(X = abs(w.store), MARGIN = 2, FUN = quantile, 
            probs = ci.level)
        ci.lo <- beta.mod - crit.t * se.clust
        ci.hi <- beta.mod + crit.t * se.clust
        print.ci <- cbind(ind.variables.names, ci.lo, ci.hi)
        print.ci <- rbind(c("variable name", "CI lower", "CI higher"), 
            print.ci)
        out.ci <- cbind(ci.lo, ci.hi)
        rownames(out.ci) <- ind.variables.names
        colnames(out.ci) <- c("CI lower", "CI higher")
        
        out <- matrix(p.store, ncol = 1)
        colnames(out) <- c("wild cluster BS p-value")
        rownames(out) <- ind.variables.names

    }
     out.list <- list()
     out.list[["pvalues"]] <- out
     out.list[["ci"]] <- out.ci
     if (output.replicates == TRUE) {
        out.list[["replicates"]] <- rep.store
    }
    if(report == TRUE){
      	print(out.list)
      }
      return(invisible(out.list))
}
   
