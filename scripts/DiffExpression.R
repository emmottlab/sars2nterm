
# Peptide differential expression analysis


library(limma)

load('ProtFiltered.RData')

# Differential expression of peptides

dea1 = run_dea_24h(dat, "mlab")
dea_report(dea1)

dea2 = run_dea_24h(dat, "imp")
dea_report(dea2)

dea3 = run_dea_24h(datF, "mlab")
dea_report(dea3)

dea4 = run_dea_24h(datF, "imp")
dea_report(dea4)


## The NeoNT peptides have been re-assigned in a more restrictive way (motif_analysis.R)
## and this analysis uses that assignment
datF_ff = datF
datF_ff$A549$annot$NeoNT = 'Original'
datF_ff$A549$annot$NeoNT[datF_ff$A549$pep$Sequence %in% pep_set_A549$peptide_seq] = 'NeoNT'

datF_ff$Vero$annot$NeoNT = 'Original'
datF_ff$Vero$annot$NeoNT[datF_ff$Vero$pep$Sequence %in% pep_set_vero$peptide_seq] = 'NeoNT'

dea5 = run_dea_24h(datF_ff, 'imp')
dea_report(dea5)

res_A549 = dea5$A549$neopeps_24h$res
sig_pep_A549 = datF_ff$A549$pep[datF_ff$A549$annot$NeoNT == 'NeoNT',]
sig_pep_A549 = sig_pep_A549[res_A549$q.values <= 0.05,]

res_vero = dea5$Vero$neopeps_24h$res
sig_pep_vero = datF_ff$Vero$pep[datF_ff$Vero$annot$NeoNT == 'NeoNT',]
sig_pep_vero = sig_pep_vero[res_vero$q.value <= 0.05,]

save.image('DiffExprAnalysis.RData')

## -------------------------------------------------- 
## Functions
## -------------------------------------------------- 

# Limma 2 groups test (24h)
limma_2_treatments_24h <- function(xpm, tgt, varlabs, ph, tp) {
    require(limma)
    tgt = tgt[ tgt$timep == tp ,]
    grp = factor(tgt$treatment, levels = ph)
    dm = model.matrix(~grp)
    mat = xpm[,tgt$tmt]
    fit = eBayes(lmFit(mat , dm))
    res = topTable(fit, coef=ncol(dm), number = Inf, sort.by = "none")
    res = cbind(res, peptides = varlabs)
    coefs = colnames(fit$coefficients)
    
    resl = list(tgt = tgt,
                grp = grp,
                dm = dm,
                mat = mat,
                fit = fit,
                res = res,
                coefs = coefs)
    
    return(resl)
    
}

# Check rows with constant values
check_constant_values <- function(mat, tgt, ph, tp) {
    tgt = tgt[ tgt$timep == tp ,]
    grp = factor(tgt$treatment, levels = ph)
    mat = mat[,tgt$tmt]
    i = 0
    for(r in 1:nrow(mat)) {
        if(var(mat[r,]) == 0) {
            # print(paste0("Row ", r, ": ", paste(mat[r,], collapse = ", ")))
            i = i+1
        }
    }
    return(i)
}

# t-test: replicate MATLAB results
ttest_as_matlab_24h <- function(xpm, tgt, ph, tp) {
    tgt = tgt[ tgt$timep == tp ,]
    grp = factor(tgt$treatment, levels = ph)
    mat = xpm[,tgt$tmt]
    pvs = c()
    for(r in 1:nrow(mat)) {
        g1 = mat[r,tgt$treatment == ph[1]]
        g2 = mat[r,tgt$treatment == ph[2]]
        res = try(t.test(g1, g2), TRUE)
        if(length(res)>1) {
            pvs[r] = res$p.value
        } else {
            pvs[r] = NA
        }
    }
    return(pvs)
}

# Run DEA
run_dea_24h <- function(datl, mtype) {
    
    dea = list()
    
    ph = c("Mock", "Exp")
    tp = 24
    
    # Mock 24h vs Infected 24h
    for(s in samples) {
        dea[[s]][["allpeps_24h"]] = limma_2_treatments_24h(datl[[s]][[mtype]],
                                                           datl[[s]]$des,
                                                           datl[[s]]$annot$Leading.razor.protein,
                                                           ph, tp)
    }
    
    # Mock 24h vs Infected 24h (Neo-NT only, as done in the manuscript)
    for(s in samples) {
        ri = rownames(datl[[s]]$annot[ datl[[s]]$annot$NeoNT == "NeoNT" , ])
        dea[[s]][["neopeps_24h"]] = limma_2_treatments_24h(datl[[s]][[mtype]][ri,],
                                                           datl[[s]]$des,
                                                           datl[[s]]$annot[ri,]$Leading.razor.protein,
                                                           ph, tp)
    }
    
    ### Note about the warning from Limma:
    # This just means that for at least one gene the log ratio is identical for all samples
    # Since this will give a zero variance (which will end up in the denominator of your statistic
    # and could possibly result in an infinite value for your test statistic)
    # it has been offset to a small value to prevent that possibility.
    
    # This has happened because some values are constant
    # In MATLAB or t.test in R this is not permitted
    # In MATLAB results is NaN and and error in R
    
    for(s in samples) {
        nc = check_constant_values(datl[[s]][[mtype]], datl[[s]]$des, ph, tp)
        print(paste(s, "has", nc, "rows with constant values."))
    }
    
    # Compare with MATLAB p-values (neoNT only)
    for(s in samples) {
        ri = rownames(datl[[s]]$annot[ datl[[s]]$annot$NeoNT == "NeoNT" , ])
        plot(dat[[s]]$mlab_p[ri], dea[[s]]$neopeps_24h$res$P.Value)    
    }
    
    
    # T-tests
    ttest = list()
    
    for(s in samples) {
        ri = rownames(datl[[s]]$annot[ datl[[s]]$annot$NeoNT == "NeoNT" , ])
        ttest[[s]][["neopeps_24h"]] = ttest_as_matlab_24h(datl[[s]][[mtype]][ri,], datl[[s]]$des, ph, tp)
    }
    
    # Compare with MATLAB p-values (neoNT only)
    for(s in samples) {
        ri = rownames(datl[[s]]$annot[ datl[[s]]$annot$NeoNT == "NeoNT" , ])
        plot(dat[[s]]$mlab_p[ri], ttest[[s]]$neopeps_24h)    
    }
    
    # Apply Storey's FDR (as MATLAB code)
    
    library(qvalue)
    
    for(s in samples) {
        for(e in names(dea[[s]])) {
            res = dea[[s]][[e]][["res"]]
            qobj = qvalue(res$P.Value)
            res = cbind(res, q.values = qobj$qvalues)
            dea[[s]][[e]][["res"]] = res
        }
    }
    
    return(dea)
    
}

# Volcano Plots for DEA
MyVolcanoPlot <- function(res, pep, main) {
    
    q = 0.05
    fc = 0.5
    
    require(EnhancedVolcano)
    
    keyvals <- ifelse(pep[rownames(res),"db"] == "SARS2", "red", "blue")
    names(keyvals)[keyvals == 'red'] <- 'Viral'
    names(keyvals)[keyvals == 'blue'] <- 'Cellular'
    
    p = EnhancedVolcano(res,
                        x = 'logFC',
                        y = 'q.values',
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        lab = rep("", nrow(res)),
                        ylab = bquote(~-Log[10] ~ italic(Q)),
                        colCustom = keyvals,
                        title = main,
                        subtitle = "q cutoff: 0.05 / LogFC: 0.5")
    return(p)
    
}

# Report DEA
dea_report <- function(dea) {
    
    # How many DE proteins for each analysis (q<0.05)?
    for(s in samples){
        for(e in names(dea[[s]])) {
            res = dea[[s]][[e]][["res"]]
            nsig = nrow(res[ res$q.values < 0.05 , ])
            print(paste(s, e, ":", nsig, "hits from", nrow(res), "proteins."))
            
        }
    }
    
    for(s in samples) {
        for(e in names(dea[[s]])) {
            print(MyVolcanoPlot(dea[[s]][[e]][["res"]],
                                dat[[s]][["pep"]],
                                paste(s, e)))    
        }  
    }
    
}
