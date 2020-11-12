suppressMessages(library(metafor))
suppressMessages(library(weights))

# x$N <- ifelse(x$N == 4, 5, x$N)
# x$vi <- escalc(measure= "COR", ri  = x$Effect.r, ni = x$N)$vi
# x$id <- 1:nrow(x)

Table1 <- data.frame(   Construct = character(),
                        k.ES      = numeric(),
                        k.indSamp = numeric(),
                        k.papers  = numeric(),
                        N.observa = numeric(),
                        N.indpend = numeric(),
                        r         = numeric(),
                        lb        = numeric(),
                        ub        = numeric(),
                        pMR       = numeric(),
                        Q         = numeric(),
                        tau       = numeric(),
                        I2        = numeric(),
                        rankTestR = numeric(),
                        rankTestP = numeric(),
                        stringsAsFactors = FALSE)
# For testing #
data <- x
construct <- "Effect.of.exposure"
line <- 1
  
buildTable <- function(construct, line = 1, data = x) {
  
    if (length(table(data[, "Sample"])) > 1) {
    res.mv0 <- rma.mv(yi, vi,
                    random =  ~  1 | factor(Author)/factor(Sample)/factor(id),
                    data = data,
                    control=list(optimizer="optim", optmethod="Nelder-Mead"))
    print(res.mv0, digits = 2)
    Table1[line, "k.ES"]      <- res.mv0$s.nlevels.f [3]
    Table1[line, "k.indSamp"] <- res.mv0$s.nlevels.f [2]
    Table1[line, "k.papers"]  <- res.mv0$s.nlevels.f [1]
    } else {
    res.mv0 <- rma.mv(Effect.r, vi,
                        random =  ~  1 | FullArticle / factor(id),
                        data = data)
    print(res.mv0, digits = 2)  
    Table1[line, "k.ES"]      <- res.mv0$s.nlevels.f [2]
    Table1[line, "k.indSamp"] <- res.mv0$s.nlevels.f [1]
    Table1[line, "k.papers"] <- res.mv0$s.nlevels.f [1]  
      
    }
  
    # Calculate I-squared
    # http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate#multilevel_models
    W  <- diag(1/data$vi)
    X  <- model.matrix(res.mv0)
    P  <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2 <- 100 * sum(res.mv0$sigma2) / 
        (sum(res.mv0$sigma2) + (res.mv0$k-res.mv0$p)/sum(diag(P)))

    Table1[line, "Construct"] <- construct

    Table1[line, "N.observa"] <- format(sum(data$N), big.mark = ",")
    Ns                        <- data[, c("Author", "Sample", "N")]
    Ns[, c("Author", "Sample")] <- apply(Ns[, c("Author", "Sample")], 2, as.character)
    mins <- aggregate(Ns, list(data[, "Author"], data[, "Sample"]), min)
    Table1[line, "N.indpend"] <- format(sum(mins$N), big.mark = ",")
    Table1[line, "r"]         <- rd(res.mv0$b, 2)
    Table1[line, "lb"]        <- rd(res.mv0$ci.lb, 2)
    Table1[line, "ub"]        <- rd(res.mv0$ci.ub, 2)
    Table1[line, "Q"]         <- format(round(res.mv0$QE, 1), big.mark = ",")
    Table1[line, "tau"]       <- ifelse(sqrt(sum(res.mv0$sigma2) < .001), ".00",
                                             rd(sqrt(sum(res.mv0$sigma2)), 2))
    Table1[line, "I2" ]       <- round(I2, 1)
    Table1[line, "rankTestR"] <- suppressWarnings(ifelse
                                 (ranktest(res.mv0)$tau == 0, ".00",
                                 rd(ranktest(res.mv0)$tau, 2)))
    Table1[line, "rankTestP"] <- suppressWarnings(ifelse
                                 (ranktest(res.mv0)$pval < .01, ".01", 
                                 rd(ranktest(res.mv0)$pval, 2)))
    return(Table1)
}
# buildTable ("Effect.of.exposure", line = 1, data = x)
