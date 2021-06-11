

library(clusterProfiler)

# Enrichment analysis using results of the subset enrichment only

load('MotifAnalysis.RData')

path = "../data/genesets/"
gmt_files = list.files(path, "\\.gmt$")

genesets = list() # List containing all gene sets

for(f in gmt_files) {
    genesets[[f]] = read.gmt(paste0(path, f))
}

# List gene set names
for(gs in names(genesets)) {
    print(gs)
}

enr1 = do_enrich(dea5, datF_ff)
report_enrich(enr1)


dotplot(enr1$A549$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt, showCategory=20)
dotplot(enr1$A549$neopeps_24h$c5.bp.v7.1.entrez.gmt, showCategory=20)
dotplot(enr1$A549$neopeps_24h$c5.cc.v7.1.entrez.gmt, showCategory=20)
dotplot(enr1$A549$neopeps_24h$c5.mf.v7.1.entrez.gmt, showCategory=20)

cnetplot(enr1$A549$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt)
cnetplot(enr1$A549$neopeps_24h$c5.bp.v7.1.entrez.gmt, showCategory=20)
cnetplot(enr1$A549$neopeps_24h$c5.cc.v7.1.entrez.gmt, showCategory=20)
cnetplot(enr1$A549$neopeps_24h$c5.mf.v7.1.entrez.gmt, showCategory=20)

emapplot(enr1$A549$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt)
emapplot(enr1$A549$neopeps_24h$c5.bp.v7.1.entrez.gmt, showCategory=20)
emapplot(enr1$A549$neopeps_24h$c5.cc.v7.1.entrez.gmt, showCategory=20)
emapplot(enr1$A549$neopeps_24h$c5.mf.v7.1.entrez.gmt, showCategory=20)


dotplot(enr1$Vero$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt, showCategory=20)
dotplot(enr1$Vero$neopeps_24h$c5.cc.v7.1.entrez.gmt, showCategory=20)
dotplot(enr1$Vero$neopeps_24h$c5.mf.v7.1.entrez.gmt, showCategory=20)

cnetplot(enr1$Vero$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt)
cnetplot(enr1$Vero$neopeps_24h$c5.cc.v7.1.entrez.gmt, showCategory=20)
cnetplot(enr1$Vero$neopeps_24h$c5.mf.v7.1.entrez.gmt, showCategory=20)

save.image('Final.RData')

##-------------------------------------------------- 
# Functions
##-------------------------------------------------- 

get_entrez_gls <- function(res, annot, eidcol) {
    sig = res[res$q.values < 0.05, ]
    hits = annot[rownames(sig),eidcol]
    bg = annot[,eidcol] # All peptides in the enriched dataset
    return(list(hits = hits, bg = bg))
}

# Enrichment in all collections of gene sets
batch_enrich <- function(gl, gsets, cutoff) {
    
    if(any(is.na(gl$hits))) {
        message("Warning: NA gene names.")
    }
    if(any(duplicated(gl$hits))) {
        message("Warning: Duplicated gene names.")
    }
    
    resl = list()
    
    for(gs in names(gsets)) {
        res = enricher(gene = gl$hits,
                       TERM2GENE = gsets[[gs]],
                       minGSSize = 10,
                       maxGSSize = 500,
                       universe = gl$bg,
                       pvalueCutoff = cutoff,
                       pAdjustMethod = "BH")
        if(nrow(res@result) > 0) {
            resl[[gs]] = res 
        }
    }
    return(resl)
}

# Do Functional Enrichment
do_enrich <- function(dea, datl) {
    
    gl = list()
    
    for(s in samples) {
        for(e in names(dea[[s]])) {
            gl[[s]][[e]] = get_entrez_gls(dea[[s]][[e]][["res"]],
                                          datl[[s]]$annot,
                                          ifelse(s == "Vero", "entrezid.hs", "entrezid"))
        }
    }
    
    enr = list()
    
    for(s in samples) {
        for(e in names(gl[[s]])) {
            enr[[s]][[e]] = batch_enrich(gl[[s]][[e]], genesets, 0.05)
        }
    }
    
    return(enr)
}

report_enrich <- function(enr) {
    
    # Significant gene sets (q<0.05)
    for(s in samples){
        for(e in names(enr[[s]][2])) {
            for(gs in names(enr[[s]][[e]])) {
                res = enr[[s]][[e]][[gs]]@result
                res = res[ res$qvalue < 0.05 , ]
                if(nrow(res)>0) {
                    cat(sprintf('\n\n\n'))
                    print(paste(s, e, gs, ":", nrow(res), "."))
                    print(res[,c("GeneRatio", "BgRatio", "pvalue", "qvalue", "Count")])
                }
            }
        }
    }
    
}
