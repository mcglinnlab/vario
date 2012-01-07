'census' = function(sim, snap = length(sim$snaps), type = "census") 
{
    ## Function from package neutral.vp written by Tyler Smith
    ## Citation: 
    ## Smith, T. W. and Lundholm, J. T. 2010. Variation partitioning as a tool
    ##   to distinguish between niche and neutral processes. – Ecography 33:
    ##   648-655. 
    pop = sim$snaps[[snap]]
    dim(pop) = c(sim$p$S, sim$p$M, sim$p$M)
    if (type == "census") {
        output = t(pop[, , 1])
        for (i in 2:sim$p$M) output <- rbind(output, t(pop[, 
            , i]))
    }
    else if (type == "richness") 
        output = apply(pop, MARGIN = 2:3, FUN = function(x) sum(as.logical(x)))
    else if (type == "abundance") 
        output = apply(pop, MARGIN = 2:3, FUN = sum)
    else stop("Invalid type")
    return(output)
}

