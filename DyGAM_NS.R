library(magrittr)
library(foreach)
library(doParallel)
library(mgcv)

make_DESeq2 = function(geneCounts, colData, design){
    dds = DESeq2::DESeqDataSetFromMatrix(
        countData=geneCounts,
        colData=colData,
        design=design
    ) %>% DESeq2::DESeq()
    keep = rowSums(DESeq2::counts(dds) >= 10) >= 3
    return(dds[keep, ])
}

make_data_list = function(dds){
    cts = DESeq2::counts(dds, normalized=TRUE) %>% round()
    purrr::map2(
        t(cts) %>% tibble::as_tibble(), rownames(cts),
        function(x, y){
            tibble::tibble(counts=x, time=dds$Time, geneId=y)
        }
    )
}

make_theta_list = function(dds){
    1 / DESeq2::dispersions(dds)
}

make_knots_list = function(timePoints){
    timePointsCount = length(timePoints)
    firstTimePoint = timePoints[1]
    lastTimePoint = timePoints[timePointsCount]
    knotsList = purrr::map(
        1:(timePointsCount-2),
        ~ combn(
            timePoints[2:(timePointsCount-1)],
            .x, simplify=FALSE
        )
    ) %>%
        purrr::reduce(c) %>%
        purrr::map(~ c(firstTimePoint, .x, lastTimePoint))
    return(knotsList)
}

DyGAM_NS_optimize_knots = function(data, knotsList, theta){
    fitList = purrr::map(
        knotsList,
        ~ mgcv::gam(
            counts ~ s(time, bs="cr", k=length(.x)),
            family=mgcv::negbin(theta),
            data=data, method="REML",
            knots=list(time=.x)
        )
    )
    AICcList = purrr::map_dbl(
        fitList,
        AICcmodavg::AICc
    )
    bestFit = fitList[[which.min(AICcList)]]
    return(bestFit)
}

DyGAM_NS = function(dataList, knotsList, thetaList){
    n = length(dataList)
    bestFit = foreach (i=1:n, .verbose=TRUE) %dopar% {
        DyGAM_NS_optimize_knots(data=dataList[[i]], knotsList=knotsList, theta=thetaList[[i]])
    }
    tibble::tibble(
        geneId = names(dataList),
        data = dataList,
        bestFit = bestFit,
        pvalue = purrr::map_dbl(bestFit, ~ broom::tidy(.x)$p.value),
        padj = p.adjust(pvalue, method="BH"),
        rsquare = purrr::map_dbl(bestFit, ~ summary(.x)$r.sq),
        devExpl = purrr::map_dbl(bestFit, ~ summary(.x)$dev.expl),
        theta = purrr::map_dbl(bestFit, ~ .x$family$getTheta()),
        knots = purrr::map_chr(
            bestFit,
            ~ stringr::str_c(.x$smooth[[1]]$xp, collapse=",")
        ),
        converged = purrr::map_lgl(bestFit, ~ .x$converged)
    )
}

filter_DTIGs = function(x){
    dplyr::filter(x, pvalue < 0.05, padj < 0.05, devExpl > 0.7)$geneId
}



registerDoParallel(80)

knotsList = make_knots_list(c(0, 1, 3, 6, 9, 12, 24, 36, 48, 72))

dds = make_DESeq2(cts, colData, ~ time)

dygamns_fit = DyGAM_NS(
    make_data_list(dds),
    knotsList,
    make_theta_list(dds)
)

DTIGs = filter_DTIGs(dygamns_fit)

stopImplicitCluster()
