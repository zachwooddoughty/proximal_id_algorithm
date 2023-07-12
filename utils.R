# suppressMessages(library(R.utils, warn.conflicts=FALSE))  # for arg parsing
suppressMessages(library(jsonlite))
library(stringr)

df_make = function(mat, names){
    df <- as.data.frame(mat)
    colnames(df) <- names
    return(df)
}

lst_to_df <- function(data) {
    mat <- c()
    for (name in names(data)) {
        mat <- cbind(mat, unlist(data[[name]]))
    }
    df <- data.frame(mat)
    names(df) <- names(data)

    return(df)
}

df_resample = function(df){
    # Resample dataframe (for bootstrapping)
    nrow <- dim(df)[1]
    ncol <- dim(df)[2]
    sampled_df <- df_make(matrix(nrow=nrow, ncol=ncol), names(df))

    sampled_rows <- sample(c(1:nrow), nrow, replace=TRUE)
    for (i in c(1:nrow)) {
        sampled_df[i, ] <- df[sampled_rows[i],]
    }
    
    return(sampled_df)
}

bootstrap_ci = function(df, k, f, ...){
    # Use nonparametric bootstrap to call function f
    # on k different bootstrap samples from dataset df.
    # Each time, sample a new df of the same length.
    stats <- c()
    failures <- 0
    for (i in c(1:k)){
        tryCatch(
        {
            df1 <- df_resample(df)
            result <- f(df1, ...)
            stats <- rbind(stats, result)
        }, error=function(cond){
            failures <- failures + 1
        })
    }

    if (failures > 0){
        warning(paste(failures, "bootstrap failures"))
    }
    return(stats)
}

interval_stats <- function(arr, truth) {
    # Compute bootstrap interval coverage and width
    quantiles <- quantile(arr, c(0.025, 0.975))
    coverage <- quantiles[1] < truth & truth < quantiles[2]
    width <- quantiles[2] - quantiles[1]
    return(c(as.integer(coverage), width))
}

clip_weights <- function(weights) {
    # Clip weights to the middle 95% interval
    percentiles <- quantile(weights, c(0.025, 0.975))
    low <- percentiles[1]
    high <- percentiles[2]
    return(low + (weights - low > 0) * (weights - low)
           - (weights - high > 0) * (weights - high))
}

printout <- function(name, arr, cols=NULL) {
    # Helper function to print things nicely
    tmp <- as.data.frame(t(arr))
    colnames(tmp) <- cols
    rownames(tmp) <- name
    print(format(tmp, digits=3, nsmall=2))
}
