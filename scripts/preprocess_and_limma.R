library(affy)
library(limma)
library(AnnotationDbi)
library(drosophila2.db)

file_dir <- "raw_data/E-MTAB-1620/"

raw <- ReadAffy(celfile.path =file_dir)
sampleNames(raw)
length(sampleNames(raw))

#ExpressionSet
eset <-rma(raw) #background correction, quantile normalization, summarization, log2 transformation
expr <- exprs(eset) #access numeric expression values

#link probes to genes
probes <- rownames(expr)
ann <- AnnotationDbi::select(
  drosophila2.db,
  keys     = probes,
  keytype  = "PROBEID",
  columns  = c("SYMBOL","FLYBASE","GENENAME","ENTREZID")
)

sym_map <- setNames(ann$SYMBOL, ann$PROBEID)
SYMBOL  <- unname(sym_map[probes])
fData(eset) <- data.frame(PROBEID = probes, SYMBOL = SYMBOL, row.names = probes)
keep <- !grepl("_x_at$", probes)       # drop cross-hyb
expr  <- expr[keep, ]
SYMBOL <- SYMBOL[keep]
valid <- !is.na(SYMBOL) & SYMBOL != ""
expr_gene <- avereps(expr[valid, ], ID = SYMBOL[valid], fun = median) #average replicates - for probes that mapped to a single gene

#information about each .CEL file is here: https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-1620/sdrf
region <- factor(c("Crop", "Crop", "Crop", "R1", 
                   "R1", "R1", "R2", "R2", "R2", 
                   "R3", "R3", "R3", "R4", "R4", "R4", "R5", "R5", "R5", 
                   "Hindgut", "Hindgut", "Hindgut", "WholeGut", "WholeGut", "WholeGut"))
region <- factor(region, levels = c("Crop","R1","R2","R3","R4","R5","Hindgut", "WholeGut"))
design <- model.matrix(~0 + region)
colnames(design) <- levels(region)

fit <- lmFit(expr_gene, design)

lvls <- levels(region)
others_mean <- function(target, all_lvls) {
  others <- setdiff(all_lvls, target)
  sprintf("%s - (%s)/%d", target, paste(others, collapse = "+"), length(others))
}
cons_str <- sapply(lvls, others_mean, all_lvls = lvls)
names(cons_str) <- paste0(lvls, "VsOthers")   # <-- names become column names

cont <- makeContrasts(contrasts = cons_str, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
colnames(fit2$coefficients)
r3 <- topTable(fit2, coef = "R3 - (Crop+R1+R2+R4+R5+Hindgut+WholeGut)/7", number = Inf, adjust.method = "BH")
