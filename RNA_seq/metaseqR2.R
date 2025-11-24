# ------------------------------------------------------------
# Set seed for reproducibility
# ------------------------------------------------------------
set.seed(1234)

# ------------------------------------------------------------
# Define input/output paths and output name
# the.path: directory containing input target file and data
# save: directory to store output results
# out_name: name for output folder / report prefix
# ------------------------------------------------------------
the.path <- "/media/samba/hatzis_lab/iliopoulos2023/analysis"
save <- "/media/samba/hatzis_lab/iliopoulos2023/analysis"
out_name="delidakis_all_comp"

# ------------------------------------------------------------
# Define contrasts to be analyzed with metaseqR2
# Each contrast corresponds to a pairwise differential expression comparison
# ------------------------------------------------------------
the.contrasts <- c("trxRNAi-Suc_vs_cnt-Suc",
                   "trxRNAi-DSS_vs_cnt-DSS",
                   "cnt-Suc_vs_cnt-DSS",
                   "trxRNAi-Suc_vs_trxRNAi-DSS")

# Load metaseqR2 package
library(metaseqR2)

# ------------------------------------------------------------
# Run metaseqR2 pipeline
#   - sampleList: target sample annotation file
#   - org: organism reference (dm6 for Drosophila melanogaster)
#   - countType: using UTR count files
#   - normalization: DESeq normalization
#   - statistics: list of statistical methods for differential expression
#   - metaP: P-value integration method ("pandora")
#   - qcPlots: quality control plots to generate
#   - export options: which data tables and values to export
#   - reportDb: reporting system ("dexie")
#   - createTracks: UCSC genome browser tracks
#   - trackInfo: browser track settings
# ------------------------------------------------------------
metaseqr2(sampleList = file.path(the.path, "targets_delidakis.txt"), 
          contrast = the.contrasts, 
          org = "dm6",
          countType = "utr", 
          #annotation="embedded",
          normalization = "deseq", 
          statistics = c("deseq", "deseq2", "edger", "nbpseq", "noiseq"), 
          metaP = "pandora", 
          figFormat = c("png", "pdf"), 
          exportWhere = file.path(save, print(out_name)), 
          restrictCores = 0.25,
          qcPlots = c("mds", "biodetection", "countsbio", "readnoise", 
                      "filtered", "correl", "boxplot", "gcbias", "lengthbias", 
                      "meandiff", "volcano", "biodist", "mastat"), 
          exonFilters = NULL, 
          pcut = 0.05, 
          exportStats= "cv",
          exportWhat=c("annotation", "p_value", "adj_p_value", "meta_p_value" ,"fold_change", "stats", "counts", "flags"),
          exportScale = c("natural","log2", "rpgm"), 
          exportValues = "normalized", 
          exportCountsTable = TRUE,
          saveGeneMode = TRUE,
          reportTop = 0.1,
          localDb="/media/raid/users/agalaras/annotation.sqlite",
          reportDb= "dexie",
          createTracks = TRUE, 
          overwrite = TRUE, 
          trackInfo = list(stranded = TRUE, 
                           normTo = 1e+08, 
                           urlBase = "http://epigenomics.fleming.gr/~alexandros/delidakis", 
                           hubInfo = list(name = "delidakis Quant-seq", 
                                          shortLabel = "delidakis Quant-seq", 
                                          longLabel = "delidakis Quant-seq", 
                                          email = "galaras@fleming.gr")))
