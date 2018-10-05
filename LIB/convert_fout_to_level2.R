# v0.02 MM
# This script generates the second level output with some more information into the table
#
# Author: Margherita Mutarelli

#R CMD BATCH --no-save --no-restore '--args <aname> <dbhost> <dbname> <dbuser> <dbpw> ' convert_fout_to_level2.R

args =  commandArgs(trailingOnly = TRUE)
# print(args)
aname = args[1]#FRANK
aname
dbhost <- args[2]#FRANK
dbhost
dbname <- args[3]#FRANK
dbname
dbuser <- args[4]#FRANK
dbuser
dbpw <- args[5]#FRANK
dbpw
workdir <- args[6]#FRANK
workdir
vargenius_out <- args[7]#FRANK
vargenius_out
out_2nd_lev <- args[8]#FRANK
out_2nd_lev
rprep_file <- args[9]#FRANK
rprep_file


####################################################
### chunk number 1: general directories and files
###################################################
#rprep_file <- '/pico/work/TELET_UDP/vargenius_analyses/DATA/prepare_UDPpico_final.RData'#FRANK
#base_analysis_dir <- '/pico/work/TELET_UDP/vargenius_analyses'
# base_website_dir <- '/home/users/ngs/exomeAnalysis/exoMeDiData/results/UDP'

###################################################
### chunk number 2: load libraries, functions and options
###################################################
options(java.parameters = "-Xmx16G")
options(stringsAsFactors = FALSE)

library('RPostgreSQL')

###################################################
### chunk number 3: set general variables
###################################################
#dbhost <- 'postgresql.pico.cineca.it'#FRANK
#dbname <- 'telethon'#FRANK
#dbuser <- 'ttudp001'#FRANK
#dbpw <- '!nj5DRMi'#FRANK

# finalout_dir <- 'finalout_out'
#txtsuffix <- '.hg19.final_out_resgr14.txt'#FRANK
#outsuffix <- '.hg19.final_out_resgr14_level2'#FRANK

# frequency column names
frecols <- c('exac_all', 'exac_max', 'esp6500siv2_all','db1000g2015aug_all', 'db1000g2015aug_afr', 'db1000g2015aug_amr', 'db1000g2015aug_eas', 'db1000g2015aug_eur', 'db1000g2015aug_sas', 'db1000g_max', '_1000g2015aug_all', '_1000g2015aug_afr', '_1000g2015aug_amr', '_1000g2015aug_eas', '_1000g2015aug_eur', '_1000g2015aug_sas', '_1000g_max', 'gnomad_exome_all', 'gnomad_genome_all', 'gnomad_ex_max', 'gnomad_gen_max')

# chrcol <- 'chr' # alcuni 'chrom'
chrcol <- 'chrom'
poscol <- 'pos'
refcol <- 'ref'
altcol <- 'alt'
keycols <- c(chrcol, poscol, refcol, altcol)

clinvarcols <- c('clinsig', 'clndbn', 'clnacc', 'clndsdb', 'clndsdbid')
clinsigcol <- clinvarcols[1]
clinphencol <- 'ClinVar_PHEN'
hgmdcols2remove <- c('CLASS_hgmd', 'GENE_hgmd', 'MUT_hgmd', 'PHEN_hgmd', 'PMID_hgmd', 'hgmd_info')

funcol <- 'func_refgene'
# spdcol <- 'splvar_dist_from_exon'
spdcol <- c('splvar_dist_from_exon', 'refgene_splvar_dist_from_exon') 

gencol <- 'gene_refgene'
splth <- 2
fth <- 0.01
filtcol = 'filter'
 
genefuncionclasscol <- 'GeneFunctionClass'
lastgenefuncionclasslab <- 'Other'
customsplicecol <- 'canonSplSites'

impactstring <- list(
	ncRNA = c('ncRNA')
	, LossOfFunction = c('stop', 'frameshift', 'canonical')
	, CodingChange = c('nonframeshift', 'nonsynonymous')
	)#,
	# Other = c())
refseqol <- 'refgene_transcript'
gendetcol <- 'genedetail_refgene'

###################################################
### chunk number 4: get analysis specific variables
###################################################

# aname = 'UD_NA011_QUARTET'
#an_out_dir <- paste(base_analysis_dir, aname, 'finalout_out', sep='/')#FRANK
#outxtf <- paste(aname, txtsuffix, sep='')#FRANK
outxtf <- vargenius_out#FRANK
# pedf <- paste(aname, '.ped', sep='')

#setwd(an_out_dir)#FRANK
setwd(workdir)#FRANK
#parsedtxtf <- paste(aname, outsuffix, '.txt', sep='')#FRANK
#parsedRf <- paste(aname, outsuffix, '.RData', sep='')#FRANK
parsedtxtf <- paste(out_2nd_lev, '.txt', sep='')#FRANK
parsedRf <- paste(out_2nd_lev, '.RData', sep='')#FRANK

### load annotation tables and functions
load(rprep_file)


### read variations table
outxt <- read.delim(outxtf, na='-', header=FALSE, skip=1)
colnames(outxt) <- scan(outxtf, what='ch', nlines=1)
cat(aname, 'has', nrow(outxt), 'variants\n')
outxt <- removeEmptyCols(outxt, aname=aname)
# dim(outxt)

### add key
cat('- in', aname, 'adding key\n')
# chrcol <- colnames(outxt)[grep('chr', colnames(outxt))]
# keycols <- c(chrcol, poscol, refcol, altcol)
outxt$key <- apply(outxt[, keycols], 1, paste, collapse='_')
outxt$key <- gsub(' ', '', outxt$key)


### summarize frequency 
cat('- in', aname, 'adding Freq_Max\n')
frecols <- intersect(frecols, colnames(outxt))
outxt$Freq_Max <- apply(outxt[, frecols], 1, function(x) {max(x, na.rm=TRUE)} )

### add ExAC genotype counts 
cat('- in', aname, 'adding exac_n_hom\n')
# disp(hhcnts)
# # 10134070        4
# # ExAC_AC_Het ExAC_AC_Hom ExAC_AC_Hemi            key
hhcnts <- unique(hhcnts[hhcnts$key %in% outxt$key, ])
cat('- in', aname, 'found', nrow(hhcnts), 'variants with ExAC N homozygous count\n')
outxt <- merge(outxt, hhcnts, by='key', all=TRUE)
# dim(outxt)
# # 85050   101

### add denovo-db
cat('- in', aname, 'adding de novo annotation from denovo-db\n')
# disp(dndb2ann)
 # # 75472 4 
  # # key denovodb_PubmedID denovodb_PrimaryPhenotype  denovodb_Validation
dndb2ann <- unique(dndb2ann[dndb2ann$key %in% outxt$key, ])
# disp(dndb2ann)
if(nrow(dndb2ann) > 0) {
	cat('- in', aname, 'found', nrow(dndb2ann), 'variants in denovo-db\n')
	outxt <- merge(outxt, dndb2ann, by='key', all=TRUE)
	# dim(outxt)
} else {
	cat('- in', aname, 'no variant found in denovo-db\n')
}

### reorder freq columns
# summary(outxt[, c(frecols, 'Freq_Max')])
if(nrow(dndb2ann) > 0) {
	frecol2keep <- c('exac_n_hom', frecols[grep('_all$', frecols)], 'Freq_Max', colnames(outxt)[grep('^denovodb', colnames(outxt))])
} else {
	frecol2keep <- c('exac_n_hom', frecols[grep('_all$', frecols)], 'Freq_Max')
}
# if(nrow(dndb2ann) > 0) {
# 	frecol2keep <- c(colnames(hhcnts)[-match('key', colnames(hhcnts))], frecols[grep('_all$', frecols)], 'Freq_Max', tdbcols2keep[-match('key', tdbcols2keep)], c('IVF', 'freq_factors'), colnames(outxt)[grep('^denovodb', colnames(outxt))])
# } else {
# 	frecol2keep <- c(colnames(hhcnts)[-match('key', colnames(hhcnts))], frecols[grep('_all$', frecols)], 'Freq_Max', tdbcols2keep[-match('key', tdbcols2keep)], c('IVF', 'freq_factors'))
# }
if(nrow(dndb2ann) > 0) {
	frecol2keep <- c(colnames(hhcnts)[-match('key', colnames(hhcnts))], frecols[grep('_all$', frecols)], 'Freq_Max', c('IVF', 'freq_factors'), colnames(outxt)[grep('^denovodb', colnames(outxt))])
} else {
	frecol2keep <- c(colnames(hhcnts)[-match('key', colnames(hhcnts))], frecols[grep('_all$', frecols)], 'Freq_Max', c('IVF', 'freq_factors'))
}
# disp(outxt[,frecol2keep])
ixf <- match(frecols, colnames(outxt))
ixf2keep <- match(frecol2keep, colnames(outxt))
ix2keepfinal <- c(1:(ixf[1]-1), ixf2keep, setdiff((ixf[length(ixf)]+1):ncol(outxt), ixf2keep))
# setdiff(colnames(outxt), colnames(outxt)[ix2keepfinal])
# disp(outxt[, ix2keepfinal])
outxt <- outxt[, ix2keepfinal]
for (fqcl in frecol2keep[grep('all', frecol2keep)]) {
	outxt[is.na(outxt[, fqcl]), fqcl] <- -1
}
# summary(outxt[, frecol2keep])
rm(ixf, ixf2keep, ix2keepfinal, fqcl)


### summarize ClinVar
cat('- in', aname, 'summarizing ClinVar\n')
outxt$ClinVar_CLASS <- ifelse(grepl('Pathogenic', outxt[, clinsigcol]), 'Pathogenic', ifelse(grepl('pathogenic', outxt[, clinsigcol]), 'Likely pathogenic', ifelse(grepl('drug response', outxt[, clinsigcol]), 'drug_response', ifelse(grepl('Likely benign', outxt[, clinsigcol]), 'Likely benign', ifelse(grepl('Benign', outxt[, clinsigcol]), 'Benign', ifelse(grepl('Uncertain significance', outxt[, clinsigcol]) | grepl('not provided', outxt[, clinsigcol]) | grepl('other', outxt[, clinsigcol]), 'Not_informative', ''))))))
outxt$ClinVar_Details <- apply(outxt[, clinvarcols[-2]], 1, paste, collapse='|')
  outxt$ClinVar_Details[outxt$ClinVar_Details == 'NA|NA|NA|NA'] <- NA
  colnames(outxt)[match(clinvarcols[2], colnames(outxt))] <- clinphencol
cat('- in', aname, 'found', sum(!is.na(outxt$ClinVar_Details)), 'variants present in ClinVar\n')

### add HGMD
cat('- in', aname, 'adding HGMD\n')
hgmd <- hgmd[hgmd$key %in% outxt$key,]
cat('- in', aname, 'found', nrow(hgmd), 'variants present in HGMD\n')
outxt <- merge(outxt, hgmd, by='key', all=TRUE)
dim(outxt)

### reorder known clinical variants columns
ixf <- match(clinphencol, colnames(outxt))
ixf2keep <- c(grep('ClinVar', colnames(outxt)), grep('HGMD', colnames(outxt)))
ix2keepfinal <- c(1:(ixf[1]-1), ixf2keep, setdiff((ixf[length(ixf)]+1):ncol(outxt), ixf2keep))
# colnames(outxt)[ix2keepfinal]
outxt <- outxt[, ix2keepfinal]
col2remove <- c(clinvarcols, hgmdcols2remove)
# colnames(outxt)[setdiff(match(col2remove, colnames(outxt)), NA)]
outxt <- outxt[,-setdiff(match(col2remove, colnames(outxt)), NA)]
rm(ixf, ixf2keep, ix2keepfinal, col2remove)


### reclassify canonical splicing
cat('- in', aname, 'reclassify canonical splicing\n')
spdcol <- intersect(spdcol, colnames(outxt))
cat('using column', spdcol, 'for exon edge distance\n')
canonical_splicing <- abs(outxt[, spdcol]) <= splth & !is.na(outxt[, spdcol] <= splth)
outxt[canonical_splicing, funcol] <- paste(outxt[canonical_splicing, funcol], '_canonical', sep='')
# mytable(outxt[,funcol])


### cleanup databases after use to free memory
rm(hgmd, hhcnts, dndb2ann)
gc()


### create GeneFunctionClass
outxt[, genefuncionclasscol] <- rep(lastgenefuncionclasslab, nrow(outxt))
# outxt$GeneFunctionClass = rep('', nrow(outxt))
# mytable(outxt[, genefuncionclasscol])
x <- outxt[, funcol]

# ncRNA Class
impactclass <- names(impactstring)[1]
currentstring <- impactstring[[impactclass]]
ispos <- rep(FALSE, length(x))
for (i in seq(length(currentstring))) {
	cat(i, '-', currentstring[i], '\n')
	isposel <- grepl(currentstring[i], x)
	print(mytable(isposel))
	print(mytable(x[isposel]))
	ispos <- ispos | isposel
}
# mytable(ispos)
# print(mytable(x[ispos]))
# print(mytable(x[!ispos]))
# mytable(outxt[ispos, genefuncionclasscol])
outxt[ispos, genefuncionclasscol] <- impactclass
# mytable(outxt[, genefuncionclasscol])
rm(ispos, i, isposel)

# LossOfFunction Class
# impactclass <- 'LossOfFunction'
impactclass <- names(impactstring)[2]
currentstring <- impactstring[[impactclass]]
ispos <- rep(FALSE, length(x))
for (i in seq(length(currentstring))) {
	cat(i, '-', currentstring[i], '\n')
	isposel <- grepl(currentstring[i], x)
	print(mytable(isposel))
	print(mytable(x[isposel]))
	ispos <- ispos | isposel
}
# mytable(ispos)
# print(mytable(x[ispos]))
# print(mytable(x[!ispos]))
# mytable(outxt[ispos, genefuncionclasscol])
outxt[ispos, genefuncionclasscol] <- impactclass
# mytable(outxt[, genefuncionclasscol])
rm(ispos, i, isposel)

# CodingChange Class
impactclass <- names(impactstring)[3]
currentstring <- impactstring[[impactclass]]
ispos <- rep(FALSE, length(x))
for (i in seq(length(currentstring))) {
	cat(i, '-', currentstring[i], '\n')
	isposel <- grepl(currentstring[i], x)
	print(mytable(isposel))
	print(mytable(x[isposel]))
	ispos <- ispos | isposel
}
# mytable(ispos)
# print(mytable(x[ispos]))
# print(mytable(x[!ispos]))
# mytable(outxt[ispos, genefuncionclasscol])
outxt[ispos, genefuncionclasscol] <- impactclass
# mytable(outxt[, genefuncionclasscol])
rm(ispos, i, isposel)

# custom_splice_column-based class
if (!is.na(match(customsplicecol, colnames(outxt)))) {
	custom_canonical_splicing <- !is.na(outxt[, 'canonSplSites'])
	# mytable(custom_canonical_splicing)
	outxt[custom_canonical_splicing, genefuncionclasscol] <- names(impactstring)[grep('[L,l]oss', names(impactstring))][1]
	# mytable(outxt[, genefuncionclasscol])
}


### add coordinates to visualize in IGV
outxt <- addCoords(outxt, chrcol='chrom', stcol = 'pos', encol = NULL, appEnd = FALSE)


### add variant type
outxt$vartype <- ifelse(nchar(outxt$alt) != 1, 'INSDEL', ifelse(nchar(outxt$ref) != 1, 'INSDEL', 'SNP'))


### create Quality Class
gqcols <- colnames(outxt)[grep('GQ', colnames(outxt))]
ns <- length(gqcols)
samps <- sub('_GQ', '', gqcols)

cat('Quality Class\n\n')
chkperc <- data.frame(matrix(NA, nrow=nrow(outxt), ncol=ns))
colnames(chkperc) <- samps                     
chkqual <- chkperc 
colnames(chkqual) <- paste(samps, 'QualityClass', sep='_')
for (j in seq(ns)) {
	# j = 1
	sampname <- samps[j]
	gqcol <- gqcols[grep(sampname, gqcols)]
	tmp <- checkPerc(outxt, samp=sampname, hethlo=0.25, pasfx='PercALT_reads', nrefx='REF_reads', naltx='ALT_reads', sep='_', plot=TRUE)
	# tmp[outxt[, zygcol] %in% c('HOMREF', '') | is.na(outxt[, zygcol])] <- 'NO_CALL'
	chkperc[, j] <- tmp
	chkqual[, j] <- chkperc[, j] == 'OK' & outxt[, filtcol] %in% c('PASS', '.') & (outxt[, gqcol] >= 50 & !is.na(outxt[, gqcol]))
	rm(sampname, tmp)
}
outxt <- cbind(chkqual, outxt)


### recover empty RefSeq for some categories
narefseq2replace <- is.na(outxt[, refseqol]) & !is.na(outxt[, gendetcol]) & outxt[, funcol] != 'intergenic'
tmp <- outxt[narefseq2replace, gendetcol]
tmp <- sub(':', '|', tmp)
# disp(tmp)
# 43947
tmp2 <- text2columns(tmp, str='|', nc=2)
outxt$RefSeq <- NA
outxt$RefSeq[narefseq2replace] <- tmp2[,1]
outxt$RefSeq[is.na(outxt[, 'RefSeq']) & !is.na(outxt[, refseqol])] <- outxt[is.na(outxt[, 'RefSeq']) & !is.na(outxt[, refseqol]), refseqol]


### add gene annotation
#disp(anno_table)
# 202091 11 
cat('- in', aname, 'adding gene annotation\n')
anno2table <- anno_table[anno_table$RefSeq %in% outxt$RefSeq, ]
anno2table <- anno2table[!is.na(anno2table$RefSeq), ]
cat('- in', aname, 'found', nrow(anno2table), 'transcripts and', length(unique(anno2table$Symbol)), 'genes with annotation\n')
#dim(outxt)
outxt <- merge(outxt, anno2table, by = 'RefSeq', all = TRUE)
#dim(outxt)
#length(unique(anno_table, c('RefSeq','GeneID','Symbol')

### write output 
cat('- saving', parsedRf, 'RData file\n')
save(outxt, file = parsedRf)	
cat('- saving', parsedtxtf, 'txt file\n')
write.table(outxt, file = parsedtxtf, row.names=FALSE, na = '', sep = '\t', quote=F)

cat('-', aname, 'done\n')


