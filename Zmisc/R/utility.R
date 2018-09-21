#' create a 3-group class. receives a data frame with obs prob columns
create_pred_table <- function(DF,thl,thh, levs = c("low","mid","high")) {
    DF <- as.data.table(DF)
    DF[prob<thl, cl := levs[1]]
    DF[prob >= thl & prob < thh, cl := levs[2]]
    DF[prob >= thh, cl := levs[3]]
    DF$cl <- factor(DF$cl)
    DF
}

#' just summary after factorization
#' @export
summaryf <- function(vec) {
    summary(factor(vec))
}

#' summary for prediction table
#' @export
summtab <- function(pred_table) {
    ta <- addmargins(table(pred_table[,c(1,3)]))
    ta <- ta[,c(2,3,1,4)]
    p <- round(prop.table(ta,1)[,1],2)
    ta <- rbind(ta, "Percent Yes" = round(ta[1,]/ta[3,],2))
    rbind(ta, "Percent Cat" = round(ta[3,]/ta[3,4],2))
}

####
#' defined for a two level factor, will omit the rest.
#' @export
summ <- function (v, pos = "yes", neg = "no") {
    ta <- as.data.frame(t(as.matrix(table(as.factor(v)))))  
    ta <- ta[,c(pos,neg)]
    ta["Total"] <- ta[,pos]+ta[,neg] 
    ta["Percent"] <- round(ta[,pos]/ta$Total*100,2)
    ta
}


get_col_perc <- function (v) {
    l <- do.call(rbind.data.frame,strsplit(v," "))
    names(l) <- c("n1","prop1")
    l[,2] <- gsub("[()]", "", l[,2])
    l <- apply(l,2,as.character)
    l
}

#' compare two proportions
compare_prop <- function(n1,prop1,n2,prop2) {
    prop1 = prop1/100
    prop2 = prop2/100
    tot1 = round(n1/prop1)
    tot2 = round(n2/prop2)
    ret <- tryCatch(
         chisq.test(as.table(matrix(c(n1,tot1-n1,n2,tot2-n2),ncol=2)))$p.value,error=function(cond) return("NA"))
    return(format.pval(as.numeric(ret),digits=3,eps=0.001))
}

#'
comp_col <- function (t1,col1,col2) {
    l <- data.frame(get_col_perc(t1[[col1]]),get_col_perc(t1[[col2]]))
    l <- as.data.frame(apply(l,2,as.numeric))
    names(l) <- c("n1","prop1","n2","prop2")
    ch <- apply(l, 1, function(x)do.call(compare_prop, as.list(x)))
    ch
}

#' get fit from caret and extract ppv/npv
get_ppv_npv_res <- function (fit) {
        res <- fit$resampledCM
        pr <- as.data.table(fit$pred)
        if (fit$method == "svmRadial") {
            bestTune <- fit$bestTune
            res <- res[res$sigma==bestTune[,"sigma"] & res$C==bestTune[,"C"],]
            pr <- pr[pr$sigma==bestTune[,"sigma"] & pr$C==bestTune[,"C"],]
        }
        res <- colMeans(res[1:4])
        ppv <- unname(res[1]/(res[1]+res[3]))
        npv <- unname(res[4]/(res[2]+res[4]))
        pr <- pr[,.(prob = mean(Yes)), by = .(rowIndex,obs)][,.(obs,prob)]
        list("NPV" = npv, "PPV" = ppv,"probs" = pr)
}


dtoc <- function (d) {
    (1+d/100)/2
}

getppv <- function (sen,spe) {
    sen/(sen*(1-spe))
}

switch_table <- function (conftab) {
    confv <- as.vector(conftab)
    matrix(c(confv[4],confv[3],confv[2],confv[1]),ncol=2)
}

#' bin a model
#' @param model lm object type
#' @return cut values 
#' @export
binnit_model <- function (model) {
    predictions <- round(predict(model, type = "fitted.ind"),2)
    predictions.cut <- cut2(predictions, cuts = seq(0,1,0.05))
}

#' create a bin
#' @param x predictions  
#' @param y reference 
#' @return a ggplot 
#' @export
binnit <- function (x,y,cuts = seq(0,1,0.05),g = 10, pos = "dodge") {
    if (is.numeric(cuts)) {
        x.cut <- cut2(x, cuts = cuts)
    }
    else {
        x.cut <- cut2(x, g=g)
    }
    binsdat <- data.frame(x.cut, Reference = factor(y))
    ggplot(binsdat,aes(x.cut)) + geom_bar(aes(fill=Reference), position=pos)
}

#' format text to be used in Rmd's mostly
#' @param a matrix with observations and events
#' @param event the way events are coded
#' @return converted namesl 
#' @export
n_percent_format <- function(tabl, event="1") {
    tcols <- colnames(tabl)
    n_events <- tabl[which(rownames(tabl)==event),] # a named vector
    n_all <- apply(tabl,2,sum) # a named vector
    v <- unlist(lapply(tcols,function(x) paste0(n_events[x],"/",n_all[x], " (",round(n_events[x]/n_all[x],2)*100,"%)")))
    names(v) <- tcols
    return(v)
}

#' Take care of variabe names and try to make them consistent
#' @param x vector of character names 
#' @return converted namesl 
#' @export
convert_names <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub('\\s|-+','_',x)
    return(x)
}

#' Format regression results using a table
#' @param obj object create by glm() or coxph()
retablef <- function(obj) {
    require(broom)
    retable <- tidy(summary(pool(glmobj)))[,c(1,2,6,7,8)]
    retable <- cbind(retable[,1],round(exp(retable[,c(2,4,5)]),2),format.pval(retable[,3],digits=3,eps=0.001))
    retable["95% CI"] <- paste0(retable$lo.95,"-",retable$hi.95)
    retable <- retable[,c(1,2,6,5)]
    colnames(retable) <- c("Names","Odds Ratio", "95% CI", "P-value")
    names(retable) <- capitalize(names(retable))
    return(retable)
}

#' Print only the cox table from cph
#' @param x object create by glm() or coxph()
#' @return pretty printed talbe
cph.table <- function(x, conf.int = 0.95, digits=4, eps=0.001) { 
    coefs = x$coef
    se   = sqrt(diag(x$var))
    Z    <- coefs/se
    P    <- 1 -pchisq(Z ^ 2, 1)
    pdigits <- nchar(sub('\\.','',as.character(0.001)))-1
    ta <- cbind(round(cbind(exp(coefs),se,Z),digits),format.pval(P,digits=pdigits,eps=eps))
    colnames(ta)  <- c("Hazard Ratio","Standard Error", "Z", "p-value")
    if (is.numeric(conf.int)) {
        zcrit <- qnorm((1 + conf.int)/2)
        lo <- round(exp(coefs - zcrit * se),digits)
        hi <- round(exp(coefs + zcrit * se),digits)
        ta <- cbind(ta,paste(lo," - ",hi))
        colnames(ta)[dim(ta)[2]] <- "95% CI"
    }
    print(ta)
}

#' pretty print regression table
#' @export
lrm_modfit <- function(x, conf.int=0.95, digits=2) { 
    lang = 'plain'
    obj <- list(coef=x$coef,se=sqrt(diag(x$var)))
    beta <- obj$coef
    se   <- obj$se
    Z    <- beta/se
    P    <- 1 - pchisq(Z ^ 2, 1)
    DF <- data.frame("Odds ratio" = round(exp(beta),2), "S.E" = round(se,2), "Wald Z" = round(Z,2), "p-value" = format.pval(P,eps=0.001,digits=digits))
    if (is.numeric(conf.int)) {
        zcrit <- qnorm((1 + conf.int)/2)
        lo <- round(exp(beta - zcrit * se),digits)
        hi <- round(exp(beta + zcrit * se),digits)
        DF <- cbind(DF,paste(lo," - ",hi))
        colnames(DF)[dim(DF)[2]] <- "95% CI"
        DF[,c(1,5,2,3,4)]
    }
    DF
}

#'
#' @export
central_disp <- function(param, cent_type = mean, disp_type = sd, digits = 2) {
    disp <- disp_type(param,na.rm=T)
    if (is.numeric(disp)) {
        disp <- round(disp,digits)
    }
    paste0(round(cent_type(param,na.rm=T),digits)," (",disp,")")
}

#' like t() but makes sure a data.table is a dataframe and that row names are inserted
#' @export
transp <- function(DF) {
    DF <- as.data.frame(DF)
    rownames(DF) <- DF[[1]]
    DF[,1] <- NULL
    t(DF)
}

#' like IQR but formatted as range
#' @export
iqr <- function(x,na.rm=F) {
    paste0(quantile(x,1/4,na.rm=na.rm)," - ",quantile(x,3/4,na.rm=na.rm))
}

#' take two factors and arrange in a proportion table
#' @export
table_by <- function(fac1,fac2,dat,yes="Yes") {
     dat <- as.data.frame(dat)
     ta <- table(dat[[fac1]], dat[[fac2]])
     ta <- ta[, !(colnames(ta) %in% "")]
     cbind(ta, total = margin.table(ta, 1), percent_yes = round(prop.table(ta, 
         1)[, yes], 3))
}

#' similar to table_by but accomodates multilevel factors.
#' @export
table_prop <- function(fac1, fac2,dat) {
    tab <- table(dat[[fac1]],dat[[fac2]])
    pr <- prop.table(tab,1)[,2]
    tab <- as.data.frame.matrix(tab)
    tab$percent <- round(pr,2)
    tab
}

#' take a 2-level vector and code to numeric
#' @export
yesno_to_num <- function (v,yes = "Yes", no ="No") {
    v[v==yes & !is.na(v)] <- 1
    v[v==no & !is.na(v)] <- 0
    as.integer(v)
}

#' turn a two level num vector to yes/no
#' @export
num_to_yesno <- function (v, yes = "Yes", no = "No") {
    v[v==1 & !is.na(v)] <- yes
    v[v==0 & !is.na(v)] <- no
    factor(as.character(v))

}

#' turn a factor back ti numeric
#' @export
asnf <- function (v) {
    as.numeric(levels(v))[v]
}

