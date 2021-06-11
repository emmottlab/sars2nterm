

## Data reprocessing and annotation

##-------------------------------------------------- 
## Loading datasets
##-------------------------------------------------- 

dat = list()

samples = c('A549','Vero')

# Load the peptide files
dat$A549$pep <- read.table('../data/tA549_Enrich/peptides.txt', header = T, sep = '\t')
dat$Vero$pep <- read.table('../data/tVero_Enrich/peptides.txt', header = T, sep = '\t')

# Load in the TMTpro design matrices containing the randomised layouts
# A549-ACE2
dat$A549$tmt <- read.csv('../data/SARS2a549tmtlabelling_20200507.csv', header = F)
# VeroE6
dat$Vero$tmt <- read.csv('../data/Verotmtlabelling_20200511.csv', header = F)

# Reorder the corrected RI intensity channels according to randomised layouts
reorder_randomcols <- function(pep, tmt) {
    ric_cols = colnames(pep)
    ric_cols = ric_cols[grepl("Reporter.intensity.corrected.", ric_cols)]
    ric_n = as.numeric(sapply(strsplit(ric_cols, "\\."), function (x) x[4]))
    ric_n = tmt[ 2 , match(ric_n, tmt[1,])]
    ric_ord = paste0("Reporter.intensity.corrected.", ric_n)
    pep[,ric_cols] = pep[,ric_ord]
    return(pep)
}

for(s in samples) {
    dat[[s]][['pep']] = reorder_randomcols(dat[[s]][['pep']], dat[[s]][['tmt']])  
}

# Design Matrix
for(s in samples){
    dat[[s]][['des']] = data.frame(tmt = paste0("Reporter.intensity.corrected.", 1:16),
                                   timep = c(rep(0, 4), rep(6, 3), rep(12, 3), rep(24, 6)),
                                   treatment = c("Mock", rep("Exp", 12), rep("Mock", 3)),
                                   replicate = c(1, rep(1:3, 5)))
    dat[[s]][['des']]$expid = paste0(dat[[s]][['des']]$treatment, "_",
                                     "t", dat[[s]][['des']]$timep, "h", "_", dat[[s]][['des']]$replicate)
}

# Previously normalised data as in the manuscript
# I have re-run MATLAB github code and expoorted normalised matrix before neo-N-termini filtering on line 75
dat$A549$mlab = as.matrix(read.csv("emmott_git/sars2nterm-matlab/exports/dat.mat.A549.norm.txt", header = F))
dat$Vero$mlab = as.matrix(read.csv("emmott_git/sars2nterm-matlab/exports/dat.mat.Vero.norm.txt", header = F))
# I have re-run MATLAB github code and expoorted normalised matrix before knnimpute
dat$A549$mlab2 = as.matrix(read.csv("emmott_git/sars2nterm-matlab/exports/dat.mat.A549.beforeimp.txt", header = F))
dat$Vero$mlab2 = as.matrix(read.csv("emmott_git/sars2nterm-matlab/exports/dat.mat.Vero.beforeimp.txt", header = F))

# Add colnames to MATLAB matrices
for(s in samples) {
    for(m in c("mlab", "mlab2")) {
        colnames(dat[[s]][[m]]) = paste0("Reporter.intensity.corrected.", 1:16)
    }
}

# Previously generated P-values
dat$A549$mlab_p = as.numeric(readLines("emmott_git/sars2nterm-matlab/exports/dat.p.A549.txt"))
dat$Vero$mlab_p = as.numeric(readLines("emmott_git/sars2nterm-matlab/exports/dat.p.Vero.txt"))

rm(s, m)

##-------------------------------------------------- 
## Normalisation as in the manuscript
##-------------------------------------------------- 

normTMT <- function(datl, nacutoff) {
    
    # Peptides Table
    pep <- datl$pep
    
    ### QC
    # Remove Reverse hits, Contaminants and low PEP identifications
    
    # Remove PEP <= 0.02
    pep <- pep[ pep$PEP <= 0.02 , ]
    # Remove Reverse hits
    pep <- pep[ pep$Reverse != '+' , ]
    # Remove potential contaminants
    pep <- pep[ pep$Potential.contaminant != '+' , ]
    
    ### Expression Matrix
    
    # Columns with RI data
    ric_cols = colnames(pep)[grepl("Reporter.intensity.corrected.", colnames(pep))]
    # Convert table columns to a matrix
    mat <- as.matrix(pep[,ric_cols])
    # Convert 0 to NaN
    mat[mat == 0] <- NA
    
    ### Filter out rows above missing data threshold
    
    # comment from github code: identify rows with >75% NaN 
    # however, code was actually identifying rows with > 0.8125% NaN
    # I am reproducing what the original code does in here
    allNA <- (apply(is.na(mat), 1, sum)/ncol(mat)) > nacutoff
    mat <- mat[!allNA , ]
    pep <- pep[!allNA , ]
    
    ### Normalisation and imputation
    
    # Median normalise to control for protein loading
    mat <- med_norm_by_col(mat)
    
    ######## Until here I am getting same results as MATLAB code!  
    
    # KNNimpute missing data (k calculated as sqrt of column n)
    # Default: >50% missing values it performs mean imputation instead
    library(impute)
    imp <- impute.knn(mat, k=4)
    imp <- imp$data
    
    # Row normalise by mean (whole dataset)
    imp <- imp / rowMeans(imp)
    
    # Adjust rownames of MATLAB imported matrices
    rownames(datl$mlab) = rownames(mat)
    rownames(datl$mlab2) = rownames(mat)
    
    # Return datasets
    datl$pep = pep
    datl$mat = mat
    datl$imp = imp
    return(datl)
    
}

# Normalise by the column median
med_norm_by_col <- function(m) {
    medcols <- apply(m, 2, median, na.rm = T)
    medmat <- t(t(m)/medcols)
    return(medmat)
}

# Do preprocessing as in the manuscript
# Missing data threshold: 81.25%
for (s in samples){
    dat[[s]] = normTMT(dat[[s]], 0.8125)
}

##-------------------------------------------------- 
### Annotation of Neo-N-Termini
##-------------------------------------------------- 

label_neoNT <- function(datl, start, exp) {
    # Labelling peptides starting at a given position (start)
    datl$pep$NeoNT = ifelse(datl$pep$Start.position >= start, "NeoNT", "Original")
    # Comparison with MATLAB exported table with Neo-N-termini
    path = paste0("emmott_git/sars2nterm-matlab/exports/dat.pep.", exp ,".neoNTs.txt")
    if(sum(datl$pep$NeoNT == "NeoNT") == length(readLines(path)) - 1) {
        print(paste0("Neo-N-termini labelling same as MATLAB output for ", exp, "."))
    } else {
        print("Warning! Different number of identified Neo-N-termini peptides.")
    }
    return(datl)
}

# Neo-N-termini: position in protein is greater than or equal 2 (comment from MATLAB code)
# Warning! Code is actually selecting greater than 2
# If Neo-N-peptides also start at position 2, these are being ignored
# Labelling Neo-N-termini here as manuscript (starting position = 3)
for(s in samples) {
    dat[[s]] = label_neoNT(dat[[s]], 3, s)
}

### Annotate p-values from MATLAB with rownames (neo-N-termini only)
for(s in samples) {
    names(dat[[s]]$mlab_p) = rownames( dat[[s]]$pep[ dat[[s]]$pep$NeoNT == "NeoNT" , ] )
}

rm(s)

##-------------------------------------------------- 
### Protein names, gene names and IDs annotation
##-------------------------------------------------- 

gn = list()

# Import gene name/signal peptide data used in MATLAB github code (Uniprot)
# reference file: EnrichedCellularNterm_20200923.m

# Human
gn$A549$up <- read.csv('emmott_git/sars2nterm-matlab/data/human_acc_to_gn_signaltransit.csv',
                       header = T, na.strings = "")
colnames(gn$A549$up) = c("up_id", "symbols", "signalpep", "transitpep")
# Chlorocebus sp
gn$Vero$up <- read.csv('emmott_git/sars2nterm-matlab/data/vero_acc_to_gn_signaltransit.csv',
                       header = T, na.strings = "")
colnames(gn$Vero$up) = c("up_id", "symbols", "signalpep", "transitpep")

# Import Uniprot to Entrez IDs mapping (From Uniprot website) based on IDs in tables imported above
# Note: Not all Uniprot IDs had a matching EntrezID

# Human
gn$A549$eid <- read.table('data/Human_upID_2_EntrezID_fromUniprot.txt', sep = "\t", header = T)
colnames(gn$A549$eid) = c("up_id", "entrezid")
# Chlorocebus sp
gn$Vero$eid <- read.table('data/Vero_upID_2_EntrezID_fromUniprot.txt', sep = "\t", header = T)
colnames(gn$Vero$eid) = c("up_id", "entrezid")

# Fix genes "date" names ruined by Excel (such as MARCH7 as 7-Mar)
for(s in samples) {
    pat = "[0-9]{1,2}-Mar"
    fix = grep(pat, gn[[s]]$up$symbols)
    if(length(fix) > 0) {
        gn[[s]]$up$symbols[fix] = gsub("-Mar", "", gn[[s]]$up$symbols[fix])
        gn[[s]]$up$symbols[fix] = paste0("MARCH", gn[[s]]$up$symbols[fix])
    }
}

# Split IDs from peptides table
split_up_id <- function(pep) {
    uids = pep$Leading.razor.protein
    annot = as.data.frame(do.call(rbind, strsplit(uids, split = "\\|")))
    colnames(annot) = c("db", "up_id", "up_acc")
    pep = cbind(pep, annot)
    return(pep)
}

for(s in samples) {
    dat[[s]]$pep = split_up_id(dat[[s]]$pep)
}

# Warning in case not all ids from peptides tables in Uniprot annot tables
peptide_has_annot <- function(datl, gnl) {
    peps = datl$pep[ datl$pep$db != "SARS2" , ] # Exclude viral prots
    if(!all(peps$up_id %in% gnl$up$up_id)) {
        warning("IDs in peptides table not in Uniprot signal peptide table.")
    }
}

for(s in samples) {
    peptide_has_annot(dat[[s]], gn[[s]])  
}


### Map Vero (Chlorocebus sp) gene names to human gene names (Orthologs)

# Filter out Vero Uniprot IDs not in peptides tables
gn$Vero$up = gn$Vero$up[ gn$Vero$up$up_id %in% dat$Vero$pep$up_id , ]
gn$Vero$eid = gn$Vero$eid[ gn$Vero$eid$up_id %in% dat$Vero$pep$up_id , ]

# Split gene names in multiple rows
library(tidyr)
for(s in samples) {
    gn[[s]]$up2 = separate_rows(gn[[s]]$up, symbols, sep = " ")
}

# Any Gene Name not in Human Uniprot table?
gn$Vero$up2$symbols[!gn$Vero$up2$symbols %in% gn$A549$up2$symbols]

# Chsa-A is HLA-A orthologue. Checking name in human Uniprot table
gn$A549$up$symbols[grep("HLA-A", gn$A549$up$symbols)]




# Create another a human orthologue name column, replacing Chsa-A to HLA-A
gn$Vero$ort = gn$Vero$up2
gn$Vero$ort$symbols.hs = replace(gn$Vero$ort$symbols,
                                 which(gn$Vero$ort$symbols == "Chsa-A"), "HLA-A")

# Relationship table: map Vero Uniprot IDs to Human Uniprot IDs
gn$Vero$ort = merge(gn$Vero$ort[,c("symbols.hs", "up_id", "symbols")],
                    gn$A549$up2[,c("up_id", "symbols")],
                    all.x = T,
                    by.x = "symbols.hs",
                    by.y = "symbols",
                    suffixes = c(".v", ".hs"),
                    incomparables = NA)


### Final annotation tables

# Merge
gn$A549$annot = merge(gn$A549$up, gn$A549$eid, all.x = T, by = "up_id")
gn$Vero$annot = merge(gn$Vero$up, gn$Vero$ort[,c("up_id.v", "up_id.hs")],
                      all.x = T, by.x = "up_id", by.y = "up_id.v")
gn$Vero$annot = merge(gn$Vero$annot, gn$Vero$eid, all.x = T, by = "up_id")
gn$Vero$annot = merge(gn$Vero$annot, gn$A549$annot[,c("up_id", "symbols", "entrezid")],
                      all.x = T, by.x = "up_id.hs", by.y = "up_id", suffixes = c("", ".hs"))

# Filter out Human Uniprot IDs not in peptides tables
gn$A549$annot = gn$A549$annot[ gn$A549$annot$up_id %in% dat$A549$pep$up_id , ]

# Clean-up and try fill fill-out missing Entrez IDs
search_Entrez <- function(tbl) {
    # Load db
    require(org.Hs.eg.db)
    hs = org.Hs.eg.db
    # Entrez ID column
    if("entrezid.hs" %in% colnames(tbl)) {
        sf = ".hs"
    } else {
        sf = ""
    }
    # Check for Uniprot IDs with not matching Entrez ID
    misstbl = tbl[ is.na(tbl[,paste0("entrezid",sf)]) &
                       !is.na(tbl[,paste0("symbols",sf)]) , ]
    if(nrow(misstbl) > 0) {
        # Get first gene name
        syms = strsplit(misstbl[,paste0("symbols", sf)], split = " ")
        syms = sapply(syms, function (x) x[1])
        # Search hs db for matching EntrezID
        eids = AnnotationDbi::select(hs, syms, "ENTREZID", "SYMBOL")
        # Replace annotation values
        if(nrow(eids) == nrow(misstbl)) {
            tbl[rownames(misstbl),paste0("entrezid",sf)] <- eids$ENTREZID
        } else {
            message("'select()' returned 1:many relationships. Verify manually.")
        }
    }
    # Missing Entrez IDs to be matched?
    miss = tbl[ is.na(tbl[,paste0("entrezid",sf)]) &
                    !is.na(tbl[,paste0("symbols",sf)]) , paste0("symbols",sf)]
    if(length(miss) > 0) {
        message(paste("No Entrez IDs found for", paste(miss, collapse = ", ")))
    }
    return(tbl)
}

for(s in samples) {
    gn[[s]]$annot = search_Entrez(gn[[s]]$annot)
}

# ANKRD20A3 in Human is a pseudogene


### Full annotation table matching peptides table
library(limma)

for(s in samples) {
    dat[[s]]$annot = gn[[s]]$annot[ match(dat[[s]]$pep$up_id, gn[[s]]$annot$up_id) , ]
    dat[[s]]$annot = cbind(dat[[s]]$pep[,c("Leading.razor.protein", "NeoNT")], dat[[s]]$annot)
    # Gene to be used in plots
    dat[[s]]$annot$rep.gene = strsplit2(dat[[s]]$annot$symbols, " ")[,1]
}

rm(fix, pat, s)

##-------------------------------------------------- 
## QC
##-------------------------------------------------- 

library(finalfit)
library(ggplot2)

# Dims
dim(dat$A549$mat)
dim(dat$Vero$mat)

# Missing Plot after filtering
missplot <- function(datl, title) {
    mat = datl$mat
    des = datl$des
    df = as.data.frame(mat)
    colnames(df) = des$expid[ match(colnames(df), des$tmt)]
    rownames(df) = NULL
    print(missing_plot(df, title = title, plot_opts = list(xlab("Peptides"))))
}

for(s in samples) {
    missplot(dat[[s]], paste0(s, " Missing Plot"))  
}


### Missing Data Percentage

missperc <- function(datl, exp) {
    mat = datl$mat
    des = datl$des
    # Total missing data %
    mpct = round(sum(is.na(mat))/length(mat), 2) * 100
    print(paste0(exp, " missing data percentage: ", mpct, "%"))
    # Missing data % by rows
    mpctr = (apply(is.na(mat), 1, sum)/ncol(mat)) * 100
    mpctr = mpctr[order(mpctr, decreasing = T)]
    mpctr = mpctr[mpctr != 0]
    #barplot(mpctr, names.arg = NA, space = 0, border = NA, col = "#339966")
    plot(mpctr, type = 'l', ylab = "missing data %",
         main = paste0(exp, " missing data % by peptides"),
         sub = "Complete rows not shown")
    # Missing data % by columns
    mpctc = (apply(is.na(mat), 2, sum)/nrow(mat)) * 100
    mpctc = mpctc[order(mpctc, decreasing = T)]
    names(mpctc) = des$expid[match(names(mpctc), des$tmt)]
    barplot(mpctc, las = 2, cex.names = 0.8,
            main = paste0(exp, " missing data % by sample"))
}

for(s in samples) {
    missperc(dat[[s]], s)
}

### Comparison between current and MATLAB outputs

compare_dats <- function(datl, exp) {
    mat1 = datl$mat
    mat2 = datl$mlab2
    imp1 = datl$imp
    imp2 = datl$mlab # From MATLAB code, after imputation
    namap = which(is.na(mat1))
    
    # Before imputation
    plot(mat1, mat2, pch=".", cex = 2, xlab = "current", ylab = "MATLAB",
         main = paste(exp, "comparison of RI values before KNN impute"))
    
    # After imputation
    plot(imp1[-namap], imp2[-namap], pch=".", xlab = "current", ylab = "MATLAB",
         main = paste(exp, "comparison of RI values after KNN impute"), col = "blue")
    points(imp1[namap], imp2[namap], pch=".", col = "red")
    legend("topright", legend = c("not imputed", "imputed"), col = c("blue", "red"), pch = 16)
}

for(s in samples) {
    compare_dats(dat[[s]], s)
}

rm(s)

##-------------------------------------------------- 
## Checking normalisation
##-------------------------------------------------- 


# Log-transform datasets
for(s in samples) {
    dat[[s]]$mlab = log2(dat[[s]]$mlab)
    dat[[s]]$imp = log2(dat[[s]]$imp)
}

# Boxplot samples in imputed dataset
for(s in samples) {
    xlab = dat[[s]]$des$expid[ match(colnames(dat[[s]]$imp), dat[[s]]$des$tmt) ]
    boxplot(dat[[s]]$imp, las = 2, names = xlab, cex.axis = 0.8, pch = ".",
            main = paste(s, "normalised and imputed R"))
}

# Boxplot samples in MATLAB dataset
for(s in samples) {
    xlab = dat[[s]]$des$expid[ match(colnames(dat[[s]]$mlab), dat[[s]]$des$tmt) ]
    boxplot(dat[[s]]$mlab, las = 2, names = xlab, cex.axis = 0.8, pch = ".",
            main = paste(s, "normalised from MATLAB"))
}


##-------------------------------------------------- 
### PCA
##-------------------------------------------------- 

library(factoextra)

plot_pca <- function(mat, tgt, main) {
    groups = paste(tgt$treatment, tgt$timep, sep = "_")
    pca = prcomp(t(mat), scale = F) # no scaling: same measurement type
    p = fviz_pca_ind(pca, habillage = groups, axes = c(1,2),
                     repel = T, pointsize = 3, mean.point = F,
                     geom.ind = "point",
                     title = main)
    print(p)
}

# R imputed datasets
for(s in samples) {
    plot_pca(dat[[s]]$imp, dat[[s]]$des, paste(s, "R imputed"))
}

# MATLAB datasets
for(s in samples) {
    plot_pca(dat[[s]]$mlab, dat[[s]]$des, paste(s, "MATLAB imputed"))
}

# Mapping rows with more than one imputed values for the groups of interest
filter_TooImputed <- function(ph, tp, datl) {
    # Select samples based on time and treatment
    tgt = datl$des
    tgt = tgt[tgt$timep %in% tp & tgt$treatment %in% ph , ]
    tgt$grs = paste(tgt$treatment, tgt$timep, sep = "_")
    # Select matrix columns based on target tbl
    mat = datl$mat
    filt = list()
    for(g in unique(tgt$grs)) {
        msub = mat[ , tgt$tmt[tgt$grs == g] ]
        filt[[g]] = rowSums(is.na(msub))
    }
    map = do.call(cbind, filt)
    lrows = apply(map, 1, function(x) any(x>1))
    datl$pep = datl$pep[!lrows,]
    datl$mat = datl$mat[!lrows,]
    datl$imp = datl$imp[!lrows,]
    datl$mlab = datl$mlab[!lrows,]
    datl$mlab2 = datl$mlab2[!lrows,]
    datl$annot = datl$annot[!lrows,]
    return(datl)
}

datF = list()

# 24h Mock vs 24h Exp
for(s in samples) {
    datF[[s]] = filter_TooImputed(c("Mock", "Exp"), 24, dat[[s]])  
}

save.image('ProcAnnot.RData')
