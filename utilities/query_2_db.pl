#!/usr/bin/perl

#Given a list of compid returns the information about what are the samples continaing it and the sequencing type

use strict;
use warnings;
use Getopt::Long;


#vargenius.pl path
use lib '/pico/work/TELET_UDP/VarGeniusBeta/';

#Using a library for database management
use LIB::db_management qw(get_id_if_exists_from_db do_query_select_all get_id_if_exists_from_db_woconn
												fetch_all_rows get_count_from_selected do_fetch_row_array_woconn
												do_query_select_all_woconn);
#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash initialize_folders print_and_log separate_input_ids);

#Using a library to manage files
use LIB::files_management qw( list_to_array list_to_hash);


my $foldersFile = "folders.txt";
my $config_folder = "CONFIGURATION";
my $program_name = "query_2_db.pl";
my $program_version = "0.1";

						
											
my $func = "";#function to execute			
my $user_config = "";
my $input_f = "";
my $suppl_file = "";
my $out_file = "";
my $log_file = "";
my $parameters = "";

#Get the working folder and vargenius	
my ($working_folder,$program_folder) = initialize_folders($foldersFile);
my $program_config = $program_folder."/".$config_folder."/program_config.txt";

#die $usage."too less parameters\n" unless scalar(@ARGV) >= 2;#If arguments are less than 2, die


parse_command_line_args();

my $cfg_hash;

configFile2Hash($program_config,\$cfg_hash);
configFile2Hash($user_config,\$cfg_hash);


#Given a file with compids returns, for each compid, the list of samples which have those variants
#The compid is: CHR_POS_ALT_REF
if ( $func eq 'VARIANTS'){
	die "ERROR: input file  $input_f does not exist\n" unless (-e $input_f);#If input file does not exist
	print "Getting info for all variants of $input_f\n";

 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
		
	open (FILE, "<$input_f") or die "Cannot open $input_f\n";
	print "Variant\tAnalysis\tSample\tGT\tDP\tFilter\tSequencingType\n";
	while (my $row = <FILE>){
			chop($row);
			print "Analysing $row...\n";
			#Get the variant id
			my $varid = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_variants_table'},$cfg_hash->{'db_var_id'},
																$cfg_hash->{'db_var_compid'},"'".$row."'");	

			#Get all samples and analyses with the given variant id
			my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_genotype_sample_table'}.
			" WHERE ".$cfg_hash->{'db_var_ids'}." LIKE ANY ( values('$varid'),('$varid,%'),('%,$varid'),('%,$varid,%') ) ;";	
			#print "Executing: $query\n";#DEBUGCODE
			my $sample_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

			foreach my $sample_id (keys %{$sample_ids}){
					my $analysisid = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_analysis_id'},
																$cfg_hash->{'db_sample_id'},$sample_id);	
					my $sample_name = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
																$cfg_hash->{'db_sample_id'},$sample_id);
					my $seq_type = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
																$cfg_hash->{'db_analysis_id'},$analysisid);	
					my $analysis_name = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysisid);	
					my $qual = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_var_statistics_table'},$cfg_hash->{'db_qual'},
																$cfg_hash->{'db_analysis_id'},$analysisid);
					my $filter = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_var_statistics_table'},$cfg_hash->{'db_filter'},
																$cfg_hash->{'db_analysis_id'},$analysisid);
																
																
																 			
					#Get the genotype info
					my $var_query = "SELECT ".lc($cfg_hash->{'db_genotype_gt'})." FROM ".$cfg_hash->{'db_genotype_sample_table'}." WHERE ".$cfg_hash->{'db_var_ids'}." LIKE ".
				"ANY ( values('$varid'),('$varid,%'),('%,$varid'),('%,$varid,%') ) AND sampleid=$sample_id ;";	

					my $genotypes = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$var_query,lc($cfg_hash->{'db_genotype_gt'}));
					my $genotype = "";	
					
					foreach my $gen (keys %{$genotypes}){
							$genotype = $gen;
							#print $gen."\n";#DEBUGCODE
					}
					print "$row\t$analysis_name\t$sample_name\t$genotype\t$qual\t$filter\t$seq_type\n"									
			}	
	}
	close(FILE);
	
	#Disconnect db
  $dbh->disconnect()	
  
}



#Given a gene name, extracts a file with all the variants present on that gene
if ( $func eq 'VARS_ON_GENE'){
	

	if ( $parameters eq ''){ die "ERROR: the gene name  is needed. Please specify one using -p XXX. Exiting..\n";}
	#if ( $out_file eq '' ){ die "ERROR: out file is needed. Please specify one using -o XXX. Exiting..\n";}

	my $gene = $parameters;
	
	#open(OUT,">$out_file") or die "Cannot open file $out_file\n";

	#Get all the variant ids into an hash
	my $variants;
	#print "Get all the variant ids\n";
	my $query = "SELECT ".$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_var_compid'}." FROM ".$cfg_hash->{'db_variants_table'}.";";							
	#print "Executing: $query\n";
	#fetch_all_rows
	my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
								$cfg_hash->{'db_pass'},$query);
	foreach my $var_field (@$res) {
		my @arr = @$var_field;
		#Get the database id of the variant
		my $varid = $arr[0];
		my $compid = $arr[1];
		$variants->{$varid} = $compid;
	}				
  
	#Get all the variants for the gene
	$query = "SELECT ".$cfg_hash->{'db_var_id'}." FROM  ".$cfg_hash->{'db_var_annotations_table'}." WHERE ".
	            " ".$cfg_hash->{'gene_field'}." = '$gene';";	
	#print "Executing: $query\n";
	my $variant_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_var_id'});

	if (scalar (keys %{$variant_ids}) == 0){die "No gene associated with: $gene\n";}
	
	 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
		
		
	#Print the file header
	print  "variant\t".$cfg_hash->{"db_exome_allele_freq"}."\t".$cfg_hash->{"db_exome_freq_factors"}.
										"\t".$cfg_hash->{"db_targeted_allele_freq"}."\t".$cfg_hash->{"db_targeted_freq_factors"}."\t";
	print join("\t",split(",",$cfg_hash->{'fields_to_import'}	))."\n";
	foreach my $varid (keys %{$variant_ids}){
		print  $variants->{$varid}."\t";

		my $var_query = "SELECT ".$cfg_hash->{"db_exome_allele_freq"}.",".$cfg_hash->{"db_exome_freq_factors"}.
										",".$cfg_hash->{"db_targeted_allele_freq"}.",".$cfg_hash->{"db_targeted_freq_factors"}.
										" FROM ".$cfg_hash->{'db_variants_table'}." WHERE ".$cfg_hash->{'db_var_id'}." = $varid;";	
		#print_and_log( "Executing: $var_query\n",$log_file);		
		my @res = ();
		do_fetch_row_array_woconn($dbh,$var_query,\@res);
		if ( scalar(@res) > 0 ){
			print $res[0]."\t".$res[1]."\t".$res[2]."\t".$res[3];
		}else{
			print "-\t-\t-\t-\n";
		}

		#Get information from the annotation table
		$var_query = "SELECT ".$cfg_hash->{'fields_to_import'}.
										" FROM ".$cfg_hash->{'db_var_annotations_table'}." WHERE ".$cfg_hash->{'db_var_id'}." = $varid;";	
		#print_and_log( "Executing: $var_query\n",$log_file);		
		@res = ();
		do_fetch_row_array_woconn($dbh,$var_query,\@res);
		my @fields = split(",",$cfg_hash->{'fields_to_import'});
		my $num_fields = scalar(@fields);
		if ( scalar(@res) > 0 ){
			
			for (my $i = 0; $i < $num_fields; $i++){
				print "\t".$res[$i];
		  }
		}else{
			for (my $i = 0; $i < $num_fields; $i++){
				print "\t-";
		  }
		}
		print "\n";		
	}			
	
	#Disconnect db
  $dbh->disconnect()		
}




#Given a list of sample identifiers and a gene name, returns for each sample the coverage of that gene
if ( $func eq 'SAMPLES_GENE_COVERAGE'){
	
	if ( $input_f eq ''){ die "ERROR: the list of sample ids is needed. Please specify them separated by comma with -i XXX. Exiting..\n";}
	if ( $parameters eq ''){ die "ERROR: the gene name  is needed. Please specify one using -p XXX. Exiting..\n";}

	my $sep = ",";
	my $gene_name = $parameters;
	
	
	my @sample_ids = split($sep,separate_input_ids($input_f,$sep));
	my $gene_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																$cfg_hash->{'db_genes_name'},"'".$gene_name."'");			
  #Print the header
  print "sample_name\tgene_cov\n";
  foreach my $sample_id (@sample_ids){
		#print "Using sample id: $sample_id\n";#DEBUGCODE
		my $gene_cov = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_genecoverage_table'},$cfg_hash->{'db_genecoverage'},
																$cfg_hash->{'db_genes_id'}.",".$cfg_hash->{'db_sample_id'},$gene_id.",".$sample_id);
		my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
																$cfg_hash->{'db_sample_id'},$sample_id);		
		print $sample_name."\t".$gene_cov."\n";
	}
}




#Given a list of sample names, returns a table in which, for each gene, there is the number of variants.
#If you give in input one of the following strings with -p, the variants are filtered:
# - DEL: deleterious mutation
# - SYN: synonymous mutations
# - NSYN: non-synonymous mutations
if ( $func eq 'VARS_ON_SAMPLES'){
	#samples ids of the sample
	die "ERROR: input $input_f does not exist\n" unless (-e $input_f);#If input file does not exist

	#if ( $parameters ne ''){ print " You want to filter the variants per $parameters..\n";}
	
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	
	#Get the sample names for which we want the table from the input file
	#print "Putting sample ids from $input_f into an array...\n";
	my @sample_ids = list_to_array($input_f,"NO_NEW_LINE");


	#Get all the variant ids and corresponding gene into an hash
	my $gene_vars;
	#print "Get all the variant ids\n";
	my $query = "SELECT ".$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'gene_field'};							
	
	#If user wants other fields can be picked
	if ($parameters ne ''){
			if ( $parameters eq 'DEL'){
					$query .= ",".$cfg_hash->{'pred_to_transform'};
			}
			if ( $parameters eq 'SYN' or $parameters eq 'NSYN' ){
				my $func_field = (split(",",$cfg_hash->{'func_cols'}))[0]; 
					$query .= ",".$func_field;
			}							
	}
	$query .= " FROM ".$cfg_hash->{'db_var_annotations_table'}.";";
	#print "Executing: $query\n";
	
	#fetch_all_rows
	my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
								$cfg_hash->{'db_pass'},$query);
	foreach my $var_field (@$res) {
		my @arr = @$var_field;
		#Get the database id of the variant
		my $varid = $arr[0];
		my $gene = $arr[1];
		my $ins_flag = 1;
		
		#If the variants must be only deleterious
		if ( $parameters eq 'DEL'){
			#Get an array with all the scores
			my @scores = ();
			my @scores_flds = split(",",$cfg_hash->{'pred_to_transform'});
			for (my $i = 2;$i<scalar(@scores_flds);$i++){
					push(@scores,$arr[$i]);
					#print $arr[$i]." ";#DEBUGCODE
			}
			#Compute the frequencies of the elements
			my %freq;
			$freq{$_}++ for @scores;
			my @hfreq = sort { $freq{$b} <=> $freq{$a} } keys %freq;
			
			#If the highest frequency of scores is not P or LP, disactivate the flag
			
			if ( $hfreq[0] !~ /P/ or $hfreq[0] eq '.'){				
				$ins_flag = 0;
				#print "$varid NO\n";#DEBUGCODE
			}
			#else{
				#print "HF: ".$hfreq[0]."\n";#DEBUGCODE
			#}
		}	

		#If the variants must be only synonymous
		if ( $parameters eq 'SYN'){
			my $func = $arr[2];
			if ( $func !~ /;syn/){
				$ins_flag = 0;
			}
		}

		#If the variants must be only nonsynonymous
		if ( $parameters eq 'NSYN'){
			my $func = $arr[2];
			if ( $func !~ /;nonsyn/){
				$ins_flag = 0;
			}
		}
						
		#The ins_flag is used to decide if the variant should be added to the count
		#it is used only if the user gives a filter
		if ( $ins_flag == 1){
			if (!(defined $gene_vars->{$gene}) ){
				$gene_vars->{$gene} = $varid;
			}else{
				$gene_vars->{$gene} .= ",$varid";
				#print "$gene:".$gene_vars->{$gene}."\n";#DEBUGCODE
			}			
		}	

	}				

	 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
		
		
	##Print the header 
	my $header = "samplename\tsampleid\t";
	foreach my $gene (keys %{$gene_vars}){
		$header .= "$gene\t";
		#print "Gene: $gene - ".scalar(split(",",$gene_vars->{$gene}))."\n";#DEBUGCODE
	}
	chop($header);
	print "$header\n";
	
	#For each sample
	foreach my $sample_id (@sample_ids){

		my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
																$cfg_hash->{'db_sample_id'},$sample_id);
																
		my $line = "$sample_name\t$sample_id\t";
		#Get all the variants that the sample has
		my @all_vars = ();
		$query = "SELECT ".$cfg_hash->{'db_var_ids'}." FROM  ".$cfg_hash->{'db_genotype_sample_table'}." WHERE ".
								" ".$cfg_hash->{'db_sample_id'}." = $sample_id AND ".lc($cfg_hash->{'db_genotype_gt'})." NOT LIKE '0%0' ".
								"AND ".lc($cfg_hash->{'db_genotype_gt'})." NOT LIKE '.%.';";	
		#print "Executing: $query\n";
		my $variant_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
												$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_var_ids'});
		
		if (scalar (keys %{$variant_ids}) == 0){die "No variant associated with: $sample_id\n";}
		#put the variants into an array, splitting by comma						
		foreach my $varids (keys %{$variant_ids}){
			push(@all_vars,split(",",$varids))
		}		
		#Get a structure for the intersection using map
		my %original =  ();
		map { $original{$_} = 1 } @all_vars;
		
		#print "Get the intersection of ".scalar(@all_vars)." vars in $sample_id...";#DEBUGCODE
		#For each gene
		foreach my $gene (keys %{$gene_vars}){		
			my @gene_vars = split(",",$gene_vars->{$gene});
		  #print " ".scalar(@gene_vars)." in $gene |";#DEBUGCODE
			my $intersected = 0;
			#If there are variants on the gene
			if ( scalar(@gene_vars) > 0 ){
				my @isect = grep { $original{$_} } @gene_vars;
				#get the number of variants in the intersection
				$intersected =  scalar(@isect);
			}
			$line .= $intersected."\t";
		}
		chop($line);
		print $line."\n";
	}
	
	#Disconnect db
  $dbh->disconnect()	

}


#Given a list of sample names, extracts a graduatory of the most involved genes with the phenotypes
if ( $func eq 'GENE_PHEN'){
	#HPO ids of the sample
	die "ERROR: input $input_f does not exist\n" unless (-e $input_f);#If input file does not exist

	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	
	#Get the sample names for which we want the table from the input file
	print "Putting sample names from $input_f into an array...\n";
	my @sample_names = list_to_array($input_f,"NO_NEW_LINE");
	
	#Get all the genes from the database
	my $genes;
	my $query = "";
	if ( length($suppl_file) < 1){
		$query = "SELECT ".$cfg_hash->{'db_genes_name'}." FROM  ".$cfg_hash->{'db_genes_table'}.";";	
		print "Executing: $query\n";
		$genes = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_genes_name'});
				
	}else{
		print "Using file $suppl_file for gene list\n";
		list_to_hash($suppl_file,\$genes);
	}

	#Print the header 
	print "Gene\t".join("\t",@sample_names)."\n";
	
	#For each gene
	foreach my $gene_name (keys %{$genes}){
		print "Using gene: $gene_name\n";
		my $line = "$gene_name\t";
		
		#Get the list of all HPO ids associated with the gene
		my $gene_hpo_ids = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_hpo_ids'},
																$cfg_hash->{'db_genes_name'},"'".$gene_name."'");			
    
    #If hpo ids for the gene are found
    if (defined $gene_hpo_ids){
			#split the hpo ids																	
			my @gene_hpo_ids = split(/\Q$mult_ann_sep\E/,$gene_hpo_ids);
			
			#For each sample
			foreach my $sample_name (@sample_names){
				#The counter saying how many phenotypes associated with the current gene are
				#observed in the sample
				my $gene_for_sample = 0;			

				#Get the sample id
				my $personid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
																$cfg_hash->{'db_sample_name'},"'".$sample_name."'");		
																
				#Get all the phenotypes of the sample
				my $query = "SELECT ".$cfg_hash->{'db_phenotype_id'}." FROM  ".$cfg_hash->{'db_samples_hpo_table'}." WHERE ".$cfg_hash->{'db_samples_personid'}." = $personid AND ".$cfg_hash->{'db_samples_hpo_observed'}." = 't';";	
				#print "Executing: $query\n";
				my $sample_phen_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_phenotype_id'});
				
				#For each phenotype
				foreach my $sample_phen_id (keys %{$sample_phen_ids}){
					#If the sample phenotype is among those associated with gene, increment
					if ( grep {/\b$sample_phen_id\b/} @gene_hpo_ids){
							$gene_for_sample++;
					}
				}
				#concatenate the value
				$line .= "$gene_for_sample\t";
			}
			#Remove the last tab and print
			chop($line);
			print "$line\n";		
		}			
	}
}


#Given a list of sample names and a list of phenotypes, filters the samples which have at least one among those phenotypes
if ( $func eq 'SAMPLES_AFFECTED'){
	#samples
	die "ERROR: input $input_f does not exist\n" unless (-e $input_f);#If input file does not exist

	#samples
	die "ERROR: hpo ids file $suppl_file does not exist and is needed\n" unless (-e $suppl_file);#If file does not exist
	
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	
	#Get the sample names for which we want the table from the input file
	print "Putting sample names from $input_f into an array...\n";
	my @sample_names = list_to_array($input_f,"NO_NEW_LINE");

	#Get the hpo ids for which wthe samples should be affected
	print "Putting HPO ids from $suppl_file into an array...\n";
	my @phenotypes = list_to_array($suppl_file,"NO_NEW_LINE");
		

	#Print the header 
	print "analysisid\tsampleid\tsamplename\thpo_ids\n";

	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
			
	#For each sample
	foreach my $sample_name (@sample_names){
		#print "Checking sample $sample_name..\n";#DEBUGCODE
		
		#The counter saying how many phenotypes associated with the current gene are
		#observed in the sample
		my $phen_for_sample = 0;			
		
		#Get the sample id obtaining those which are used for frequencies. Hence it is uniq
		my $sample_id = -1;
		my $sample_id_query =  "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}.
		" WHERE ".$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}." FROM  ".$cfg_hash->{'db_analyses_table'}." WHERE ".
		" ".$cfg_hash->{'db_analysis_infreq'}." = 1) AND ".$cfg_hash->{'db_sample_name'}." = '$sample_name' ";	
		my $sample_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$sample_id_query,$cfg_hash->{'db_sample_id'});
		if (scalar (keys %{$sample_ids}) > 1 ){
			print "ERROR I have found ".scalar (keys %{$sample_ids})." results for query $sample_id_query\n";
		}else{
			foreach my $item (keys %{$sample_ids}){
				$sample_id	= $item;
			}
		}
		
	
		my $analysis_id = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_analysis_id'},
						 $cfg_hash->{'db_sample_id'},$sample_id);
	
		my $personid = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
						 $cfg_hash->{'db_sample_name'},"'".$sample_name."'");
						 														
		#Get all the phenotypes of the sample
		my $query = "SELECT ".$cfg_hash->{'db_phenotype_id'}." FROM  ".$cfg_hash->{'db_samples_hpo_table'}." WHERE ".$cfg_hash->{'db_samples_personid'}." = $personid AND ".$cfg_hash->{'db_samples_hpo_observed'}." = 't';";	
		#print "Executing: $query\n";
		my $sample_phen_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_phenotype_id'});
		
		my $phenids_count = scalar(keys %{$sample_phen_ids});
		#print "Obtained $phenids_count hpo ids. These are present in the list:\n" unless $phenids_count==0;#DEBUGCODE
		my $line = "";
		#For each phenotype
		foreach my $sample_phen_id (keys %{$sample_phen_ids}){
			
			my $hpoid = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_hpo_id'},
						 $cfg_hash->{'db_phenotype_id'},$sample_phen_id);
												
			#If the sample phenotype is among those given in input, increment
			if ( grep {/\b$hpoid\b/} @phenotypes){
				
				if ( $phen_for_sample == 0){
					$line .= "$analysis_id\t$sample_id\t$sample_name\t$hpoid,";
					#print "$sample_phen_id  ";#DEBUGCODE
				}elsif ( $phen_for_sample > 0 ){
					$line .= $hpoid.",";
					#print "$sample_phen_id ";#DEBUGCODE
				}
				$phen_for_sample++;
			}

		}
		chop($line);
		if ( $phen_for_sample > 0 ){
			print $line."\n";
		}
	}
	#Disconnect db
	$dbh->disconnect();
}

#permits to obtain the frequencies of the variants into the database adjusting for relatedness or by subgroup of samples
						
						#To compute the allelic frequency of a given variant we first compute its:
						#	- number of times the variant is found heterozygous (V_het)
						#	- number of times the variant is found homozygous (V_hom)
						#	- number of samples for which for the variant has been computed a genotype (Tot_gen)
							
						#Then the formula to obain the allelic frequency of the variant V is
						
						#V_all_f = V_het + (V_hom * 2) / (Tot_gen * 2)
						
						#Where we multiply V_hom by two because the homozygous variant has been found in both alleles
						#and Tot_gen because the samples have two alleles. We use only those variant for which GATK
						#is able to get a genotype.
#To use this function you need to give in input:
#   -o: the output file
#   -log: a log file
#   -parameters: the filter that you want to use
if ( $func eq 'FREQUENCIES'){

	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;

	my $step = $cfg_hash->{'update_freqs_step'};

	if ( $log_file eq '' ){ die "ERROR: log file is needed. Please specify one using -log XXX. Exiting..\n";}
	if ( $out_file eq '' ){ die "ERROR: out file is needed. Please specify one using -o XXX. Exiting..\n";}
	if ( $parameters eq ''){ print "WARNING: You are not using -p to request filtering the variants for subgroup of samples.\n";}
	
	
	#Get all sequencing type and for each type do the update
	print_and_log( "Get all sequencing types\n",$log_file);
	my $query = "SELECT DISTINCT ".$cfg_hash->{'db_sequencingtype'}." FROM  ".$cfg_hash->{'db_analyses_table'}.";";	
	print_and_log( "Executing: $query\n",$log_file);
	my $seqtypes = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sequencingtype'});
				
	#For each sequencing type execute an update process of the frequencies
	if (scalar(keys %{$seqtypes}) > 1){
		foreach my $seqtype (keys %{$seqtypes}){
			print_and_log("Starting the update process of variants frequencies for $seqtype\n",$log_file);
			
			open(OUT,">$out_file.$seqtype") or die "Cannot open file $out_file.$seqtype\n";
			#################DB
			#Connect to database
			my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
			# PERL DBI CONNECT
			my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
			
				
			#Get all the variant ids
			print_and_log( "Get all the variant ids\n",$log_file);
			my $query = "SELECT ".$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_var_compid'}." FROM ".$cfg_hash->{'db_variants_table'}.";";							
			print_and_log("Executing: $query\n",$log_file);
			#fetch_all_rows
			my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
										$cfg_hash->{'db_pass'},$query);


			#Table of all variants in heterozygosity (must be in use for frequency calculation and for this type  of sequencing)
			my $het_vars_query = "SELECT ".$cfg_hash->{'db_var_ids'}." FROM ".$cfg_hash->{'db_genotype_sample_table'}.
											 " WHERE ".$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}.
											 " FROM ".$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_analysis_infreq'}."=1 AND ".$cfg_hash->{'db_sequencingtype'}."='$seqtype' ".
											 " $parameters) AND  ".$cfg_hash->{'db_genotype_gt'}." LIKE ANY ( values('1|0'),('1/0'),('0|1'),('0/1'))";

			#Table of all variants in homozygosity	(must be in use for frequency calculation and for this type  of sequencing)							 
			my $hom_vars_query = "SELECT ".$cfg_hash->{'db_var_ids'}." FROM ".$cfg_hash->{'db_genotype_sample_table'}.
											 " WHERE ".$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}.
											 " FROM ".$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_analysis_infreq'}."=1 AND ".$cfg_hash->{'db_sequencingtype'}."='$seqtype' ".
											 " $parameters) AND  ".$cfg_hash->{'db_genotype_gt'}." LIKE ANY ( values ('1|1'),('1/1'))";
																					
			#Get the count of distinct samples where the genotype has been predicted (must be in use for frequency calculation and for this type of sequencing)
			my $distinct_gen_samples_infreq =  "SELECT DISTINCT ".$cfg_hash->{'db_sample_id'}." FROM ".$cfg_hash->{'db_genotype_sample_table'}.
											 " WHERE ".$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}.
											 " FROM ".$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_analysis_infreq'}."=1 AND ".$cfg_hash->{'db_sequencingtype'}."='$seqtype'".
											 " $parameters) AND  ".$cfg_hash->{'db_genotype_gt'}." LIKE ANY ( values ('1|1'),('1/1'),('0|0'),('0/0'),('1|0'),('1/0'),('0|1'),('0/1') )";	 
			my $tot_infreq_gen_sam = get_count_from_selected($dbh,$distinct_gen_samples_infreq," ",$log_file);
			
			#If distinct genotyped samples are present...
			if ( defined $tot_infreq_gen_sam){
				#Get two lists corresponding to the variants in heterozygosity and those in homozygosity
				my $het_vars_lists = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
											$cfg_hash->{'db_pass'},$het_vars_query);
				my $hom_vars_lists = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
											$cfg_hash->{'db_pass'},$hom_vars_query);
				
				#The variants will result as lists separated by comma. I want here a full list with all the variants
				#in heterozygosity and homozygsity so that later I use the grep on it
				my $het_vars_hash;
				foreach my $het_vars_list (@$het_vars_lists){
					my @fields = @$het_vars_list;
					foreach (split(",",$fields[0])){
						$het_vars_hash->{$_}++;
					}
					#print_and_log("Splitting and saving ".$fields[0]."\n",$log_file);#DEBUGCODE
				}
				my $hom_vars_hash;
				foreach my $hom_vars_list (@$hom_vars_lists){
					my @fields = @$hom_vars_list;
					foreach (split(",",$fields[0])){
						$hom_vars_hash->{$_}++;
					}
				}
				print_and_log("Obtained ".scalar(keys(%$het_vars_hash))." het_vars and ".scalar(keys(%$hom_vars_hash))." hom vars\n",$log_file);	
				
				#Print some statistics									
				print_and_log("het_vars_query: $het_vars_query\n",$log_file); 
				print_and_log("hom_vars_query: $hom_vars_query\n",$log_file); 
				print_and_log("distinct_gen_samples_infreq: $distinct_gen_samples_infreq\n",$log_file); 
				print_and_log("Total Genotyped $seqtype samples: $tot_infreq_gen_sam\n",$log_file);
				print_and_log("Out of ".scalar(@$res)." $seqtype variants, frequencies updated:\n",$log_file);
				
				#Using some variables to print at any 10000 variants the status of the update
				my $countThresh = 10000;									
				my $var_count = 0;
				
				#For each variant get the count of occurrences from all the samples in the database
				#avoiding those whose field infreq=0. And update the variant frequency table
				foreach my $var_field (@$res) {
					$var_count++;
					#Prints the number of variants evaluated
					if (($var_count >= 0) and ($var_count % $countThresh == 0) ){
						print_and_log(" $var_count - ",$log_file);
					}
					my @arr = @$var_field;
					#Get the database id of the variant
					my $varid = $arr[0];
					my $compid = $arr[1];
					
					#print_and_log("\nComputing the frequency for $varid\n",$log_file);
					
					#Compute the number of times the variants is found hetorozygous
					my $het_var_count = "";
					if (defined $het_vars_hash->{$varid}){
						$het_var_count = $het_vars_hash->{$varid};
					}else{$het_var_count = 0;}
					#Compute the number of times the variants is found homozygous
					my $hom_var_count = "";
					if (defined  $hom_vars_hash->{$varid}){
						$hom_var_count = $hom_vars_hash->{$varid};
					}else{$hom_var_count = 0;}

																				
					#Divide heterozygosity and homozygosity by the total number of samples
					#Check if the denominator is greater than zero
					my $sc_var_freq = -1;
					
					if ( $tot_infreq_gen_sam > 0 ){
						$sc_var_freq = ($het_var_count + ($hom_var_count*2))/($tot_infreq_gen_sam*2);
					}
					##Update table
					#Truncate too much long values coming out from the division				
					my $var_freq = 	sprintf("%.10f",$sc_var_freq);			
					#Build the freq factors field that contains all the factors used in the division
					my $freq_factors = "$het_var_count".$cfg_hash->{'mult_ann_sep'}."$hom_var_count".$cfg_hash->{'mult_ann_sep'}."$tot_infreq_gen_sam";
					
					#Build fields and values to put into the db
					my $fields = $cfg_hash->{"db_$seqtype\_allele_freq"}.";".$cfg_hash->{"db_$seqtype\_freq_factors"};
					my $values = $var_freq.";'$freq_factors'";
					#Query to update
					#update_table_woconn($dbh,$cfg_hash->{'db_variants_table'},$fields,$values,$cfg_hash->{'db_var_id'},$varid);
					
					print OUT "$compid\t$var_freq\t$freq_factors\n";
					#print_and_log("Frequency of $varid: ($het_var_count + ($hom_var_count*2))/($tot_infreq_gen_sam*2)\n",$log_file);#DEBUGCODE
				}
				#Disconnect db
				$dbh->disconnect(); 		
				
				##Update the variants frequency update history table with the date
				#my $fields = $cfg_hash->{'db_sequencingtype'}.",".$cfg_hash->{'db_var_freq_update'};
				#my ($d,$m,$y) = (localtime)[3,4,5];
				#my $year = $y+1900; # In PERL year
				#$m = $m + 1; #Month in perl starts from 0 (0:january,11:december)
				#my $values = "'$seqtype','$year\_$m\_$d'";
				#insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													#$cfg_hash->{'db_pass'},$cfg_hash->{'db_var_freq_update_hist_table'},$fields,$values);
									
				$partDuration = time - $partTime;
				print_and_log("\nFinished [GLOBAL]\trecomputing variants frequencies\t$seqtype\t$partDuration\tseconds\n",$log_file);
				$partTime = time;					
			}
			close(OUT);
		}
	}	
}

=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.
  my $VERSION = 0;# Shows version number and exit.

	my $howToUse = "Use with: \nperl $program_name \n\n".
	"-c|--user_config: The configuration file with all the parameters needed for the execution. It should have been copied ".
									"into the working folder during the installation.\n".
	"\t-f|--func: functions: \n".
					"\t\tVARIANTS (given a list of variants returns the samples in which the variant is present with their genotype\n".
					"\t\t\t To use this function you need to give in input:\n".
					"\t\t\t -i: a file with a list of compid\n".
          					
					"\t\tGENE_PHEN (Given a list of sample names, returns, for each gene, the number of phenotypes which are involved with the gene.\n".
					"\t\t\t To use this function you need to give in input:\n".
					"\t\t\t -i: a file with a list of samples for which you want to investigate the phenotypes\n".
									
					"\t\tFREQUENCIES (permits to obtain the frequencies of the variants into the database adjusting for relatedness or by subgroup of samples use also -i.\n".
					"\t\t\t To use this function you need to give in input:\n".
					"\t\t\t -o: the output file\n".
					"\t\t\t -log: a log file\n".
          "\t\t\t -parameters: the filter that you want to use\n".
          
          "\t\tSAMPLES_GENE_COVERAGE: Given a list of sample identifiers and a gene name, returns for each sample the coverage of that gene".
          "\t\t\t To use this function you need to give in input:\n".
					"\t\t\t -i: a list of comma separated sample identifiers (you can use vargenius.pl -c XXX -sh_s ANALYSIS_ID\n".
          "\t\t\t -p: the name of the gene for which you want to get the coverage for the samples in -i\n".
          
          "\t\tVARS_ON_GENE: Given a list of gene names, extracts a file with all the variants present on that gene".
          "\t\t\t To use this function you need to give in input:\n".
          "\t\t\t -p: the name of the gene for which you want to get the variants\n".
          
	"\t-s|--suppl_file: supplementary file to use in the functions.\n".
	"\t-o|--out_file: output file to use in the functions .\n".
	"\t-p|--parameters: a generic argument to use paramters into the functions .\n".
	"\t-log|--log_file: log file to be used .\n".
	"\t-i|--input_file: input file to use\n";


	  
  #  Parse options
  GetOptions(
           "help" => \$HELP,
           "version" => \$VERSION,
           "c|user_config=s" => \$user_config,   #It's mandatory
           "f|func=s" => \$func,   #It's mandator
           "i|input_file=s" => \$input_f,
           "s|suppl_file=s" => \$suppl_file,
           "p|parameters=s" => \$parameters,
           "log|log_file=s" => \$log_file,
           "o|out_file=s" => \$out_file
            );

  #Print a little help
  if ( $HELP ){
    #pod2usage(1);
    print $howToUse;
    exit;
  }

  #Print version
  if ( $VERSION ){
    print "version: $program_version \n";
    exit;
  }
  
  #Func and user config are mandatory
  if ($user_config eq ''){
		print "configuration file in input is needed. Please use --user_config. Exiting... \n";
    exit;
	}

  #Func and user config are mandatory
  if ($func eq ''){
		print "function to use is needed. Please use --func. Exiting... \n";
    exit;
	}else{
		
		#Depending by the function, check:
		if ( $func eq "GENE_PHEN" ) {
			if ($input_f eq '' ){ die "Please use --input_file parameter when doing GENE_PHEN..\n";}
		}
		
	}
}
