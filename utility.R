##
binnit <- function (model) {
    predictions <- round(predict(model, type = "fitted.ind"),2)
    predictions.cut <- cut2(predictions, cuts = seq(0,1,0.05))
    binsdat <- deri2[,.(predictions.cut, Reference = factor(model$y))]
    ggplot(binsdat,aes(predictions.cut)) + geom_bar(aes(fill=Reference), position="dodge")
}

n_percent_format <- function(tabl, event="1") {
    tcols <- colnames(tabl)
    n_events <- tabl[which(rownames(tabl)==event),] # a named vector
    n_all <- apply(tabl,2,sum) # a named vector
    v <- unlist(lapply(tcols,function(x) paste0(n_events[x],"/",n_all[x], " (",round(n_events[x]/n_all[x],2)*100,"%)")))
    names(v) <- tcols
    return(v)
}

#' @param x vector of character names 
#' @return converted namesl 
#' @export
convert_names <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub('\\s|-+','_',x)
    return(x)
}
