library(ape)
library(parallel)
library(getopt)

get_options <- function() {

    spec <- matrix(c(
        'help',  'h', 0, 'logical',   'Print this help and exit',
        'input', 'i', 1, 'character', 'Input cgMLST calls',
        'cores', 'c', 2, 'integer',   'Number of cores to use [1]'
    ), byrow = T, ncol = 5)

    opt <- getopt(spec)

    if (!is.null(opt$help)) {
        cat(getopt(spec, usage = TRUE))
        q(status = 1)
    }

    if (is.null(opt$cores)) {
        opt$cores <- 1
    }
    opt
}

get_dm <- function(calls) {

    f <- compiler::cmpfun(dist.gene)
    d <- f(calls, method = 'percentage')
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

    minima <- mclapply(rownames(dm), nonself_min, dm = dm, mc.cores = opt$cores)

    mean(as.numeric(minima))
}

main()
