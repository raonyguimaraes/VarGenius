
#VarGenius - Variant Discovery and Annotation Tool
#Copyright (C) <2017>  <Francesco Musacchia>

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
#R CMD BATCH --no-save --no-restore '--args <table1> <table2> <outname>' r_utils.R


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

job_type = args[1] 
job_type

#FUNCTIONS

#This function binds the columns of the two input tables.
cbind_tables = function (table1,table2,out_name){
	t1 <- read.table(table1,stringsAsFactors=F,head=T,sep="\t",quote="")
	t2 <- read.table(table2,stringsAsFactors=F,head=T,sep="\t",quote="")
	
	cbind_table <- cbind(t1,t2)
	#Write output
	write.table(cbind_table,file=out_name,sep="\t",row.names=F,col.names=T,quote=F)
}

#This function merges two tables using a given column
merge_tables = function (table1,table2,selcol,out_name){
	t1 <- read.table(table1,stringsAsFactors=F,head=T,sep="\t",quote="")
	t2 <- read.table(table2,stringsAsFactors=F,head=T,sep="\t",quote="")
	
	merged_tables <- merge(t1,t2,by=selcol)
	#Write output
	write.table(merged_tables,file=out_name,sep="\t",row.names=F,col.names=T,quote=F)
}

#This function merges two tables using a given column
merge_tables2 = function (table1,table2,selcol,out_name){
	t1 <- read.table(table1,stringsAsFactors=F,head=T,sep="\t",quote="")
	t2 <- read.table(table2,stringsAsFactors=F,head=T,sep="\t",quote="")
	
	merged_tables <- merge(t1,t2,by=selcol,all=TRUE)
	#Write output
	write.table(merged_tables,file=out_name,sep="\t",row.names=F,col.names=T,quote=F)
}

#This function filters those rows where the fields containing the pattern are
#all contained into list
filter_no_variant = function (table1,pattern,values,out_name){
	t1 <- read.table(table1,stringsAsFactors=F,head=T,sep="\t",quote="")
	#Get the GT columns indexes
	gt_cols <- grep(pattern,names(t1))
	#Get a new table by filtering
	fout <- t1[apply(t1[,gt_cols],1,function(row) !(all(row %in% values)) ),]	
	#Write output
	write.table(fout,file=out_name,sep="\t",row.names=F,col.names=T,quote=F)
}

#This function Given panel genes and candidate genes to belong to the panel, will add
#a column with the panel name and the field contains TRUE if the gene
#is one from the panel and CAND if belong to the candidates list
#Panel and candidate list must have column name called "Gene"
tag_panel_genes = function (table1,panel_name,panel_path,cand_path,empty_val,gene_field,out_name) {
	trueflag = "TRUE"
	candflag = "CAND"
	panel <- read.table(panel_path,stringsAsFactors=F,quote="",head=T,sep="\t")
	cand <- read.table(cand_path,stringsAsFactors=F,quote="",head=T,sep="\t")
	
	a<-read.delim(table1,stringsAsFactors=F,quote="",head=T,sep="\t")
	#Add the flag-column with the panel name
	a[[panel_name]] <-  ifelse(a[[gene_field]] %in% panel$Gene, trueflag, ifelse(a[[gene_field]] %in% cand$Gene, candflag, empty_val))
	#Print out a table
	write.table(a,file=out_name,sep="\t",quote=F,row.names=F)
}

#JOB SELECTION
if (job_type == "CBIND"){
	
	length(args)
	if (length(args) < 4){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <table1> <table2> <outname>' r_utils.R")
	}
	table1 = args[2] #First table
	table1
	table2 = args[3] #Second table
	table2
	out_name = args[4] #Output name
	out_name
	cbind_tables(table1,table2,out_name)
}

#JOB SELECTION
if (job_type == "MERGE"){
	
	length(args)
	if (length(args) < 5){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <table1> <table2> <selcol> <outname>' r_utils.R")
	}
	table1 = args[2] #First table
	table1
	table2 = args[3] #Second table
	table2
	selcol = args[4] #Selected column to use for merge
	selcol
	out_name = args[5] #Output name
	out_name
	merge_tables(table1,table2,selcol,out_name)
}

#JOB SELECTION
if (job_type == "MERGE2"){
	
	length(args)
	if (length(args) < 5){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <table1> <table2> <selcol> <outname>' r_utils.R")
	}
	table1 = args[2] #First table
	table1
	table2 = args[3] #Second table
	table2
	selcol = args[4] #Selected column to use for merge
	selcol
	out_name = args[5] #Output name
	out_name
	merge_tables2(table1,table2,selcol,out_name)
}

#When the Joint output has been filtered per set of samples, the table still contains
#all the variants obtained from the joint analysis.
#Here I write a function to filter this table removing all those cases
#where must not happen that all the samples of the family have GT=0/0 or 0|0
#pattern = "_GT" and values = c("0/0","0|0")
#hence there is no variant
if (job_type == "FILTER_NO_VARIANTS"){
	length(args)
	if (length(args) < 3){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <table1> <pattern> <out_name>' r_utils.R")
	}
	table1 = args[2] #input table
	table1
	out_name = args[3] #Output name
	out_name

	pattern = "_GT" #pattern to find
	pattern	
	values = c("0/0","0|0") #Values that the pattern should not be
	values	
	
	filter_no_variant(table1,pattern,values,out_name)

}


#Given panel genes and candidate genes to belong to the panel, will add
#a column with the panel name and the field contains TRUE if the gene
#is one from the panel and CAND if belong to the candidates list
if (job_type == "TAG_PANEL_GENES"){
	length(args)
	if (length(args) < 8){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <table1> <panel_name> <panel_path> <cand_path> <empty_val> <gene_field> <out_name>' r_utils.R")
	}
	table1 = args[2] #input table
	table1
	panel_name = args[3] #genes list
	panel_name 
	panel_path = args[4] #genes list
	panel_path
	cand_path = args[5] #candidates list
	cand_path
	empty_val = args[6] # empty value
	empty_val
	gene_field = args[7] # gene field
	gene_field
	out_name = args[8] #Output name
	out_name

	tag_panel_genes(table1,panel_name,panel_path,cand_path,empty_val,gene_field,out_name)

}
