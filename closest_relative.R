library(ape)
library(parallel)
library(getopt)

get_options <- function() {

    spec <- matrix(c(
        'input', 'i', 1, 'charcter', 'Input cgMLST calls'
    ), byrow = T, ncol = 5)

    opt <- getopt(spec)

    opt
}

get_dm <- function(calls) {

    d <- dist.gene(calls, method = 'percentage')
    m <- as.matrix(d)

    m
}

nonself_min <- function(n, dm) {

    others <- colnames(dm)[which(colnames(dm) != n)]

    other_calls <- dm[n, others]
    names(other_calls) <- others

    d <- min(other_calls)
    # names(d) <- names(which(other_calls == d))

    d
}

main <- function() {

    opt <- get_options()

    calls <- read.csv(opt$input, row.names = 1)
    dm <- get_dm(calls)

    minima <- sapply(rownames(dm), nonself_min, dm = dm)

    mean(as.integer(minima))
}

main()
