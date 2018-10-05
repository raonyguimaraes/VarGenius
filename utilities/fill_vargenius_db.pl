#!/usr/bin/perl

#Given a list of compid returns the information about what are the samples continaing it and the sequencing type

#vargenius.pl path
use lib '/pico/work/TELET_UDP/VarGeniusBeta/';

#Using a library for database management
use LIB::db_management qw(get_id_if_exists_from_db do_query_select_all insert_into_table_locked
													insert_only_if_not_exists_locked update_table get_id_if_exists_from_db_like
													update_table_2);
#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash try_exec_command initialize_folders);

#Using a library to manage files
use LIB::files_management qw( extract_columns_from_file delete_file list_to_array get_col_index);
				
use Cwd;
use Getopt::Long;
				
####Static Variables
my $foldersFile = "folders.txt";
my $hpo_update = "";
my $hpo_family_update = "";
my $samples_info = "";
my $cand_variants = "";

####Get paarameters from the configuration files
my ($working_folder,$program_folder) = initialize_folders($foldersFile);
my $program_config = $program_folder."/CONFIGURATION/program_config.txt";
my $user_config = $program_folder."/CONFIGURATION/user_config.txt";
my $cfg_hash;
configFile2Hash($program_config,\$cfg_hash);
configFile2Hash($user_config,\$cfg_hash);

#Selection of the operation to do
parse_command_line_args();

#Parse the table with HPOs of the samples using the family id
#Needs that the samples are already present in samples_info table
if (-e $hpo_family_update  ){
	print "Parsing table $hpo_update with HPOs for samples..\n";
	parse_hpo_2_db_family($cfg_hash,$hpo_update );
}

#Parses the table with candidate variants
#Needs that the variants are present in the "variants" table
#and the analyses are present
elsif (-e $cand_variants){
	print "Parsing table $cand_variants with candidate variants for analyses..\n";
	parse_cand_variants_2_db($cfg_hash,$cand_variants );
}

#Parse the table with samples information
elsif (-e $samples_info){
	print "Parsing table $samples_info with samples information...\n";
	parse_samples_info_2_db($cfg_hash,$samples_info );
}

#Parse the table with HPOs of the samples using sample and analysisname
#Needs that the samples are already present in samples_info table
if (-e $hpo_update  ){
	print "Parsing table $hpo_update with HPOs for samples..\n";
	parse_hpo_2_db($cfg_hash,$hpo_update );
}



##########################################¢ANDIDATE VARIANTS TO DB

=head2 parse_cand_variants_2_db

 Title   : parse_cand_variants_2_db
 Usage   : parse_cand_variants_2_db( -database => 'name of the database,
                               );

 Function: Parses a table like the following
					analysis_name	model	vmodel	gene	key	fmax
					UD_NA009_TRIO	new_ar	HTm,known,rare	RARS2	chr6_88299677_T_C	0.005


 Returns : nothing

=cut 
sub parse_cand_variants_2_db {
	my $cfg_hash = shift;
	my $file = shift;

	#Used to remember analyses already inserted in db
	my $analyses_inserted;
	
	#Temporary hash for the header
	my $head_hash;
	my $sep = "\t";
	print "Opening candidate variants file $file\n";
	open (FILE, "<$file") or die "Cannot open $file\n";
		
		#Go line by line and get the field values for each sample
		my $row_num = 0;
		while (my $line = <FILE>){
			#print "Analyzing line: ".$line;#DEBUGCODE
				chomp($line);
				
				#This variable is useful for the update.
				#If it is true means that the same record has been already inserted
				my $do_not_update = -1;
				
				#######READ HEADER
				#Get the field names from the header
				if ($row_num == 0){
					#Put fields in a way to be used
					my @head = split($sep,$line);
					my $n_field = 0;
					#each field will be associated to the index
					foreach my $field (@head){
							$head_hash->{$n_field} = $field;
							$n_field++;
					}
				}else{#The row is composed of data
					my @fields = split($sep,$line);
					my $n_field = 0;
					my $ss_hash;	
					#Use an hash (ss_hash) to keep all the values for this iteration
					#put the value where the field-index says at the location of the sample id given
					foreach my $field (@fields){
							#print $head_hash->{$n_field}.":$field ";
							$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							$n_field++;
					}

					###############################
					#INSERT ANALYSIS INFO					#
					###############################
					my $sep = "]---[";
					my $analysisname = $ss_hash->{$cfg_hash->{'db_analysis_name'}};
					#print "Considering $analysisname  and variant ".$ss_hash->{$cfg_hash->{'db_var_compid'}}."\n";#DEBUGCODE
						#Get the analysisid
							my $analysisid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},
																					$cfg_hash->{'db_analysis_name'},"'".$analysisname."'");

						#Get the analysisid
							my $varid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					$cfg_hash->{'db_pass'},$cfg_hash->{'db_variants_table'},$cfg_hash->{'db_var_id'},
																					$cfg_hash->{'db_var_compid'},"'".$ss_hash->{$cfg_hash->{'db_var_compid'}}."'");
							
							#print "analysisid $analysisid  and varid $varid\n";#DEBUGCODE
							#Now try to insert the value
							if ( $analysisid > 0 and $varid > 0){

								my $fields = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_var_id'};
								my $values = "$analysisid,$varid";
								
							 #print "Starting parsing...\n";##DEBUGCODE
								########Fields that will be always updated or inserted		
								#db_analysis_name
								if ( $ss_hash->{$cfg_hash->{'db_analysis_name'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_analysis_name'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_analysis_name'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_analysis_name'}}."'";
								}	
								#db_cand_vars_avsnp
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_avsnp'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_avsnp'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_avsnp'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_avsnp'}}."'";
								}
								#qual
								if ( $ss_hash->{$cfg_hash->{'db_qual'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_qual'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_qual'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_qual'}}."";
								}	
								#uncertain samples
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_uncert_sam'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_uncert_sam'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_uncert_sam'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_uncert_sam'}}."";
								}
								#freqmax
								if ( $ss_hash->{$cfg_hash->{'db_freqmax'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_freqmax'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_freqmax'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_freqmax'}};
								}											
								#ivf
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_ivf'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_ivf'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_ivf'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_ivf'}}."";
								}
								#freq factors
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_frequency_factors_fld'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_frequency_factors_fld'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_frequency_factors_fld'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_frequency_factors_fld'}}."'";
								}
								#exac_ac_het
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_het'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_het'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_exac_ac_het'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_het'}}."";
								}
								#exac_ac_hom
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_hom'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_hom'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_exac_ac_hom'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_hom'}}."";
								}
								#exac_ac_hemi
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_hemi'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_hemi'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_exac_ac_hemi'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_exac_ac_hemi'}}."";
								}
								#gene name
								if ( $ss_hash->{$cfg_hash->{'db_genes_name'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_genes_name'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_genes_name'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_genes_name'}}."'";
								}									
								
								#refgene_aa_change
								if ( $ss_hash->{$cfg_hash->{'db_refgene_aa_change'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_refgene_aa_change'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_refgene_aa_change'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_refgene_aa_change'}}."'";
								}
								#known_variant
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_known_variant'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_known_variant'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_known_variant'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_known_variant'}}."'";
								}
								#genefunctionclass
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_genefunctionclass'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_genefunctionclass'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_genefunctionclass'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_genefunctionclass'}}."'";
								}
								#exac_pli
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_exac_pli'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_exac_pli'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_exac_pli'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_exac_pli'}}."";
								}
								#exac_prec
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_exac_prec'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_exac_prec'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_exac_prec'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_exac_prec'}}."";
								}
								#inpanel
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_inpanel'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_inpanel'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_inpanel'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_inpanel'}}."'";
								}
								#ai_selection
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_ai_selection'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_ai_selection'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_ai_selection'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_ai_selection'}}."";
								}								
								###Since the variant is automatically inserted the value is always TRUE
								#automatic_selection
								#if ( $ss_hash->{$cfg_hash->{'db_cand_vars_automatic_selection'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_automatic_selection'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_automatic_selection'};
										$values .= ",true";
								#}					

								###########################TEMPROANEO
								#print "Starting temp part...\n";##DEBUGCODE
								#Use the table solved to fill the next fields
								#0 artifact, 1 benign, 2 likely benign, 3 candidate, 4 likely pathogenic, 5 pathogenic
								my $familyid = (split("_",$ss_hash->{$cfg_hash->{'db_analysis_name'}}))[1];
								#print "Searching $varid for $familyid  \n";#DEBUGCODE
								
								my $interpretation = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																		$cfg_hash->{'db_pass'},$cfg_hash->{'db_solved_table'},$cfg_hash->{'db_cand_vars_interpretation'},
																		$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_familyid'},"$varid,'$familyid'");
								#print " interpretation $interpretation \n";#DEBUGCODE
								if (length($interpretation) > 2){
									$interpretation = lc($interpretation);
									print "The variant $varid for family: $familyid was already solved...\n";
									my $interpr_val = -1;
									if ($interpretation eq 'artifact'){
											$interpr_val = 0;
									}elsif ($interpretation eq 'benign'){
											$interpr_val = 1;
									}elsif ($interpretation eq 'likely benign'){
											$interpr_val = 2;
									}elsif ($interpretation eq 'candidate'){
											$interpr_val = 3;
									}elsif ($interpretation eq 'likely pathogenic'){
											$interpr_val = 4;
									}elsif ($interpretation eq 'pathogenic'){
											$interpr_val = 5;
									}
									print "\nCandidate ".$ss_hash->{$cfg_hash->{'db_var_compid'}}." found in the solved...\n";
									#Vmodel
											$fields .= ",".$cfg_hash->{'db_cand_vars_causative'};
											$values .= ",true";											

											#$fields .= ",".$cfg_hash->{'db_cand_vars_automatic_selection'};
											#$values .= ",true";	
											
											$fields .= ",".$cfg_hash->{'db_cand_vars_manual_selection'};
											$values .= ",true";	

											$fields .= ",".$cfg_hash->{'db_cand_vars_manual_revision'};
											$values .= ",true";	
											
											$fields .= ",".$cfg_hash->{'db_cand_vars_validated'};
											$values .= ",1";	
											
											$fields .= ",".$cfg_hash->{'db_cand_vars_interpretation'};
											$values .= ",$interpr_val";												
						
								}
								####TEMPORANEO LE ULTIME POI SE LE METTONO LORO MANUALMENTE									
																		
								#If the candidate variant for this analysis  is already present, update all the fields but
								#not interpretation and the selection criteria and the model must be appended
								my $cand_vars_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																						$cfg_hash->{'db_pass'},$cfg_hash->{'db_cand_vars_table'},$cfg_hash->{'db_cand_vars_id'},
																						$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_analysis_id'},"$varid,$analysisid");
								if ( $cand_vars_id > 0){
									print "The variant $varid for analysis: $analysisid was already existent...\n";
									#he selection criteria and the model must be appended. Then I fetch them
									#selection_criteria
									if ( $ss_hash->{$cfg_hash->{'db_cand_vars_selection_criteria'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_selection_criteria'}} ne '' ){
										my $old_sel_crit = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																						$cfg_hash->{'db_pass'},$cfg_hash->{'db_cand_vars_table'},$cfg_hash->{'db_cand_vars_selection_criteria'},
																						$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_analysis_id'},"$varid,$analysisid");											
																						
											#Here I make a check: if this candidate for this analysis has been already inserted then
											#the record must not be inserted again. I am using the selection criteria for this
											if ( $old_sel_crit eq $ss_hash->{$cfg_hash->{'db_cand_vars_selection_criteria'}})	{
													$do_not_update = 1;
											}						
											$fields .= ",".$cfg_hash->{'db_cand_vars_selection_criteria'};
											$values .= ",'$old_sel_crit$sep".$ss_hash->{$cfg_hash->{'db_cand_vars_selection_criteria'}}."'";
									}
									#Model 
									if ( $ss_hash->{$cfg_hash->{'db_cand_vars_model'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_model'}} ne '' ){
											my $old_model = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																						$cfg_hash->{'db_pass'},$cfg_hash->{'db_cand_vars_table'},$cfg_hash->{'db_cand_vars_model'},
																						$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_analysis_id'},"$varid,$analysisid");
											$fields .= ",".$cfg_hash->{'db_cand_vars_model'};
											$values .= ",'$old_model$sep".$ss_hash->{$cfg_hash->{'db_cand_vars_model'}}."'";
									}										
								}#Otherwise selection criteria, the model and interpretation are new
								else{
									#selection_criteria
									if ( $ss_hash->{$cfg_hash->{'db_cand_vars_selection_criteria'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_selection_criteria'}} ne '' ){
											$fields .= ",".$cfg_hash->{'db_cand_vars_selection_criteria'};
											$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_selection_criteria'}}."'";
									}
									#Model
									if ( $ss_hash->{$cfg_hash->{'db_cand_vars_model'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_model'}} ne '' ){
											$fields .= ",".$cfg_hash->{'db_cand_vars_model'};
											$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_model'}}."'";
									}
									if (length($interpretation) <= 2){###############TEMPORANEO
										#print "Printint interpretation $interpretation\n";
										#interpretation
										if ( $ss_hash->{$cfg_hash->{'db_cand_vars_interpretation'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_interpretation'}} ne '' ){
												$fields .= ",".$cfg_hash->{'db_cand_vars_interpretation'};
												$values .= ",".$ss_hash->{$cfg_hash->{'db_cand_vars_interpretation'}}."";
										}		
								  } ###############TEMPORANEO
								}
													


								#notes
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_notes'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_notes'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_notes'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_notes'}}."'";
								}									
					
								#The second values are used to evaluate if the db_analysis_id is present
								my $fields2 = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_var_id'};
								my $values2 = "$analysisid,$varid";		
								###If the variants is present, insert
								if ( $cand_vars_id <= 0){						
								print "Inserting $fields = $values if not $fields2 = $values2  in table  ".$cfg_hash->{'db_cand_vars_table'}."\n";#DEBUGCODE
								my ($new_group_id,$exist_group_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_cand_vars_table'},$cfg_hash->{'db_analysis_id'},
														$fields,$values,$fields2,$values2);						
								}else{						
									#We need an update, but let's check if the same variant has been already inserted before
									if ($do_not_update <= 0){
										#Update
										print "Since ".$cfg_hash->{'db_analysis_id'}."=,$analysisid".$cfg_hash->{'db_var_id'}."=$varid exist.Updating: ";
										    #" $fields = $values in table  ".$cfg_hash->{'db_cand_vars_table'}."\n";#DEBUGCODE
										
										update_table_2($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},	$cfg_hash->{'db_pass'},
										$cfg_hash->{'db_cand_vars_table'},$fields,$values,$fields2,$values2);
								  }else{
										print "WARNING: The candidate variant $varid for the analysis $analysisid with ".$cfg_hash->{'db_cand_vars_selection_criteria'}.
												"=".$ss_hash->{$cfg_hash->{'db_cand_vars_selection_criteria'}}." is already present. I won't insert it!\n";
									}
								}			
							}else{
								print "Cannot find analysisid: $analysisname or varid: ".$ss_hash->{$cfg_hash->{'db_var_compid'}}.". Line will not be inserted\n";#DEBUGCODE
							}

				}#ELSE
		$row_num++;
	}
	close(FILE);

}


###################################HPO TO DB USING FAMILY ID

=head2 parse_hpo_2_db_family

 Title   : parse_hpo_2_db_family
 Usage   : parse_hpo_2_db_family( -database => 'name of the database,
                               );

 Function: Parses a table like the following
					familyid	hpoid	observed
					AAA				HP:3939	yes
					
					gets the sampleid using the sample and analysis names 
					and executes and insert_only_if_not_exists where the hpoid and sampleid 
					must not exist yet

 Returns : nothing

=cut 
sub parse_hpo_2_db_family {
	my $cfg_hash = shift;
	my $seq_info_f = shift;

	#Used to remember groups already inserted in db
	my $samples_inserted;
	
	#Temporary hash for the header
	my $head_hash;
	my $sep = "\t";
	print "Opening info file $seq_info_f\n";
	open (FILE, "<$seq_info_f") or die "Cannot open $seq_info_f\n";
		
		#Go line by line and get the field values for each sample
		my $row_num = 0;
		while (my $line = <FILE>){
		#print "Analyzing line: ".$line;
				chomp($line);
				
				#######READ HEADER
				#Get the field names from the header
				if ($row_num == 0){
					#Put fields in a way to be used
					my @head = split($sep,$line);
					my $n_field = 0;
					#each field will be associated to the index
					foreach my $field (@head){
							$head_hash->{$n_field} = $field;
							$n_field++;
					}
				}else{#The row is composed of data
					my @fields = split($sep,$line);
					my $n_field = 0;
					my $ss_hash;	
					#Use an hash (ss_hash) to keep all the values for this iteration
					#put the value where the field-index says at the location of the sample id given
					foreach my $field (@fields){
							#print $head_hash->{$n_field}.":$field ";
							$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							$n_field++;
					}

					###############################
					#INSERT ANALYSIS INFO					#
					###############################
					#Putting the analysisname and analysis id in the analysis table
					my $familyname = $ss_hash->{'familyid'};
					print "Considering $familyname \n";
					if ( ! defined $samples_inserted->{$familyname} ){


						my $personid = get_id_if_exists_from_db_like($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
																					$cfg_hash->{'db_sample_name'},"'%".$familyname."_P'");									
						if ($personid <= 0 ){
							$personid = get_id_if_exists_from_db_like($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
																					$cfg_hash->{'db_sample_name'},"'%".$familyname."_P1'");
						}
							my $phenid = -1;
							
							if ($personid > 0 ){
								my $fields = $cfg_hash->{'db_samples_personid'};
								my $values = "$personid";									
								#Check 
								my $hpo_id = $ss_hash->{$cfg_hash->{'db_hpo_id'}};
								if ( $hpo_id ne '-' and  $hpo_id ne ''){
									#Get the hpoid
									$phenid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																						$cfg_hash->{'db_pass'},$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_phenotype_id'},
																						$cfg_hash->{'db_hpo_id'},"'".$hpo_id."'");
									if ( $phenid > 0){			
										$fields .= ",".$cfg_hash->{'db_phenotype_id'};
										$values .= ",$phenid";
									}else{
										print "WARNING: $hpo_id does not exist into the database ".$cfg_hash->{'db_name'}."\n";#DEBUGCODE
									}
								}	
								#Check 
								my $hpo_observed = $ss_hash->{$cfg_hash->{'db_samples_hpo_observed'}};
								if ( $hpo_observed ne '-' and  $hpo_observed ne '' ){
										$hpo_observed = lc($hpo_observed);
										if ( $hpo_observed eq 'yes' or $hpo_observed eq 'true' or $hpo_observed eq 'y' or $hpo_observed eq 't'){
											$hpo_observed = 'TRUE';
										}else{
											$hpo_observed = 'FALSE';	
										}
										$fields .= ",".$cfg_hash->{'db_samples_hpo_observed'};
										$values .= ",$hpo_observed";
								}	

						
								#Now try to insert the value
								if ( $phenid > 0){	
									#The second values are used to evaluate if the db_analysis_id is present
									my $fields2 = $cfg_hash->{'db_samples_personid'}.",".$cfg_hash->{'db_phenotype_id'};
									my $values2 = "$personid,$phenid";																
									print "Inserting $fields = $values in table ".$cfg_hash->{'db_samples_hpo_table'}."\n";#DEBUGCODE
									my ($new_group_id,$exist_group_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_hpo_table'},$cfg_hash->{'db_phenotype_id'},
															$fields,$values,$fields2,$values2);														
								}else{
									print "Sample $samplename has not been analyzed and will not be considered\n";#DEBUGCODE
								}
								#Save in the hash
								$samples_inserted->{$samplename} = $sample_id;									
							}
													
					}
				}#ELSE
		$row_num++;
	}
	close(FILE);
}	



=head2 parse_samples_info_2_db

 Title   : parse_samples_info_2_db
 Usage   : parse_samples_info_2_db( -database => 'name of the database,
                               );

 Function: Parses a table with multiple information for the samples
					gets the sampleid using the sample and analysis names 
					and executes an insert_only_if_not_exists where the analysisid and sampleid 
					must not exist yet

 Returns : nothing

=cut 
sub parse_samples_info_2_db {
	my $cfg_hash = shift;
	my $seq_info_f = shift;

	#Used to remember groups already inserted in db
	my $samples_inserted;
	
	#Temporary hash for the header
	my $head_hash;
	my $sep = "\t";
	print "Opening info file $seq_info_f\n";
	open (FILE, "<$seq_info_f") or die "Cannot open $seq_info_f\n";
		
			#Go line by line and get the field values for each sample
			my $row_num = 0;
			while (my $line = <FILE>){
				#print "Analyzing line: ".$line;
				chomp($line);
				
				#######READ HEADER
				#Get the field names from the header
				if ($row_num == 0){
					#Put fields in a way to be used
					my @head = split($sep,$line);
					my $n_field = 0;
					#each field will be associated to the index
					foreach my $field (@head){
							$head_hash->{$n_field} = $field;
							$n_field++;
					}
				}else{#The row is composed of data
					my @fields = split($sep,$line);
					my $n_field = 0;
					my $ss_hash;	
					#Use an hash (ss_hash) to keep all the values for this iteration
					#put the value where the field-index says at the location of the sample id given
					foreach my $field (@fields){
							#print $head_hash->{$n_field}.":$field ";
							$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							$n_field++;
					}
					#print "\n";
					###############################
					#INSERT ANALYSIS INFO					#
					###############################
					#Putting the analysisname and analysis id in the analysis table
					my $samplename = $ss_hash->{$cfg_hash->{'db_sample_name'}};
					my $familyname = $ss_hash->{$cfg_hash->{'db_familyid'}};
					print "Considering $familyname and $samplename\n";
					if ( ! defined $samples_inserted->{$samplename} ){										
						
							my $fields = $cfg_hash->{'db_sample_name'}.",".$ss_hash->{$cfg_hash->{'db_familyid'}};
							my $values = "'".$samplename."','".$familyname."'";
							#The second values are used to evaluate if the sample is present
							my $fields2 = $cfg_hash->{'db_sample_name'};
							my $values2 = "'".$samplename."'";
			
							#Get the kinship using the sample name
							my $kinship = "";
							if ($samplename =~ /_P\d*$/ ) {
								$kinship = "P";
							}elsif ($samplename =~ /_M$/) {
								$kinship = "M";
							}elsif ($samplename =~ /_F$/) {
								$kinship = "F";
							}else{
								$kinship = "P";
							}
						
							#Check 
							if ( $kinship ne ''){
									$fields .= ",".$cfg_hash->{'db_sample_kinship'};
									$values .= ",'$kinship'";
							}		
							if ( $samplename =~ /_P1$/ or $samplename =~ /_P$/  ){
									$fields .= ",".$cfg_hash->{'db_sample_affected'};
									$values .= ",1";
							}	else{
										$fields .= ",".$cfg_hash->{'db_sample_affected'};
									$values .= ",0";							
							}								
														
							#Check 
							if ( $ss_hash->{$cfg_hash->{'db_sample_dob'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_dob'}} ne ''){
									$fields .= ",".$cfg_hash->{'db_sample_dob'};
									$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_dob'}}."'";
							}
							#Check 
							if ( $ss_hash->{$cfg_hash->{'db_sample_pob'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_pob'}} ne ''){
								 my $pob = $ss_hash->{$cfg_hash->{'db_sample_pob'}};
								 #Remove strange chars
								 $pob =~ s/[\'\-\.\+]//g;;
									$fields .= ",".$cfg_hash->{'db_sample_pob'};
									$values .= ",'$pob'";
							}	
							#Check 
							if ( $ss_hash->{$cfg_hash->{'db_sample_gender'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_gender'}} ne ''){
									$fields .= ",".$cfg_hash->{'db_sample_gender'};
									$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_gender'}}."'";
							}					


							#Check 
							if ( $ss_hash->{$cfg_hash->{'db_sample_notes'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_notes'}} ne ''){
									$fields .= ",".$cfg_hash->{'db_sample_notes'};
									$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_notes'}}."'";																	
							}					
						
								print "Inserting $fields = $values in table ".$cfg_hash->{'db_samples_info_table'}."\n";#DEBUGCODE
								#my ($personid,$exist_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											#		$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
											#	$fields,$values,$fields2,$values2);														
													
							
							#Save in the hash
							$samples_inserted->{$samplename} = $samplename;																	
				}	
		}
		$row_num++;
	}
			close(FILE);	
}


###################################HPO TO DB USING ANALYSISNAME AND SAMPLENAME
=head2 parse_hpo_2_db

 Title   : parse_hpo_2_db
 Usage   : parse_hpo_2_db( -database => 'name of the database,
                               );

 Function: Parses a table like the following
					analysisname	samplename	hpoid	observed
					AAA						XXX					HP:3939	yes
					
					gets the sampleid using the sample and analysis names 
					and executes and insert_only_if_not_exists where the hpoid and sampleid 
					must not exist yet

 Returns : nothing

=cut 
sub parse_hpo_2_db {
	my $cfg_hash = shift;
	my $seq_info_f = shift;

	#Used to remember analyses already inserted in db
	my $samples_inserted;
	
	#Temporary hash for the header
	my $head_hash;
	my $sep = "\t";
	print "Opening info file $seq_info_f\n";
	open (FILE, "<$seq_info_f") or die "Cannot open $seq_info_f\n";
		
		#Go line by line and get the field values for each sample
		my $row_num = 0;
		while (my $line = <FILE>){
		#print "Analyzing line: ".$line;
				chomp($line);
				
				#######READ HEADER
				#Get the field names from the header
				if ($row_num == 0){
					#Put fields in a way to be used
					my @head = split($sep,$line);
					my $n_field = 0;
					#each field will be associated to the index
					foreach my $field (@head){
							$head_hash->{$n_field} = $field;
							$n_field++;
					}
				}else{#The row is composed of data
					my @fields = split($sep,$line);
					my $n_field = 0;
					my $ss_hash;	
					#Use an hash (ss_hash) to keep all the values for this iteration
					#put the value where the field-index says at the location of the sample id given
					foreach my $field (@fields){
							#print $head_hash->{$n_field}.":$field ";
							$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							$n_field++;
					}

					###############################
					#INSERT ANALYSIS INFO					#
					###############################
					#Putting the analysisname and analysis id in the analysis table
					my $samplename = $ss_hash->{$cfg_hash->{'db_sample_name'}};
					my $analysisname = $ss_hash->{$cfg_hash->{'db_analysis_name'}};
					print "Considering $analysisname and $samplename\n";
					if ( ! defined $samples_inserted->{$samplename} ){
						
						#my $analysisid = get_id_if_exists_from_db_like($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},
																					#$cfg_hash->{'db_analysis_name'},"'%".$analysisname."%'");
												
						#Get the analysisid
						my $analysisid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},
																					$cfg_hash->{'db_analysis_name'},"'".$analysisname."'");
						#Get the sampleid
						my $sampleid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																					$cfg_hash->{'db_sample_name'}.",".$cfg_hash->{'db_analysis_id'},"'".$samplename."',$analysisid");

						my $personid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
																					$cfg_hash->{'db_sample_name'},"'".$samplename."'");									
							my $phenid = -1;
							
							if ($personid > 0 ){
								my $fields = $cfg_hash->{'db_samples_personid'};
								my $values = "$personid";									
								#Check 
								my $hpo_id = $ss_hash->{$cfg_hash->{'db_hpo_id'}};
								if ( $hpo_id ne '-' and  $hpo_id ne ''){
									#Get the hpoid
									$phenid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																						$cfg_hash->{'db_pass'},$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_phenotype_id'},
																						$cfg_hash->{'db_hpo_id'},"'".$hpo_id."'");
									if ( $phenid > 0){			
										$fields .= ",".$cfg_hash->{'db_phenotype_id'};
										$values .= ",$phenid";
									}else{
										print "WARNING: $hpo_id does not exist into the database ".$cfg_hash->{'db_name'}."\n";#DEBUGCODE
									}
								}	
								#Check 
								my $hpo_observed = $ss_hash->{$cfg_hash->{'db_samples_hpo_observed'}};
								if ( $hpo_observed ne '-' and  $hpo_observed ne '' ){
										$hpo_observed = lc($hpo_observed);
										if ( $hpo_observed eq 'yes' or $hpo_observed eq 'true' or $hpo_observed eq 'y' or $hpo_observed eq 't'){
											$hpo_observed = 'TRUE';
										}else{
											$hpo_observed = 'FALSE';	
										}
										$fields .= ",".$cfg_hash->{'db_samples_hpo_observed'};
										$values .= ",$hpo_observed";
								}	

						
								#Now try to insert the value
								#if ( $sampleid > 0 and $analysisid > 0 and $phenid > 0){
								if ( $phenid > 0){
									#The second values are used to evaluate if the db_analysis_id is present
									#my $fields2 = $cfg_hash->{'db_sample_id'}.",".$cfg_hash->{'db_phenotype_id'};
									#my $values2 = "$sampleid,$phenid";				
									#The second values are used to evaluate if the db_analysis_id is present
									my $fields2 = $cfg_hash->{'db_samples_personid'}.",".$cfg_hash->{'db_phenotype_id'};
									my $values2 = "$personid,$phenid";																
									print "Inserting $fields = $values in table ".$cfg_hash->{'db_samples_hpo_table'}."\n";#DEBUGCODE
									my ($new_group_id,$exist_group_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_hpo_table'},$cfg_hash->{'db_phenotype_id'},
															$fields,$values,$fields2,$values2);														
								}else{
									print "Sample $samplename has not been analyzed and will not be considered\n";#DEBUGCODE
								}
								#Save in the hash
								$samples_inserted->{$samplename} = $sample_id;									
							}
													
					}
				}#ELSE
		$row_num++;
	}
	close(FILE);
}


=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.

	my $howToUse = "Use with: \nperl complete_install.pl \n\n".
	"\t-hpo|--hpo_update: Use this command if you just want to update the hpos for each sample.".
								" Needs that you have a file with HPOs information n";


	  
  #  Parse options
  GetOptions(
           "help" => \$HELP,
           "hpo|hpo_update=s" => \$hpo_update,
           "hpo_f|hpo_family_update=s" => \$hpo_family_update,
           "sam|samples_info=s" => \$samples_info,
           "cand|cand_variants=s" => \$cand_variants
            );


  #Print a little help
  if ( $HELP ){
    #pod2usage(1);
    print $howToUse;
    exit;
  }

}



sub get_sample_target_coverage{
	my $sam_cov_summ_f = shift;
	my $sample_name = shift;
	
	my $retval=-1;
	
	if ( -e $sam_cov_summ_f){
		#Define the name for the table to be printed
		my $sam_cov_summ_temp = $sam_cov_summ_f.".temp";

		#From the sample summary file remove the last row starting with "Total"
		my $command = 'grep -v '."'".'^Total'."' ".$sam_cov_summ_f.' > '.$sam_cov_summ_temp;
		#print_and_log( "Executing: $command..\n",$log_file);
		#print "Executing: $command..\n";
		try_exec_command($command) or die "Unable to execute command: $command\n";	
		
		my $sample_id_ind = get_col_index($sam_cov_summ_f,"sample_id");
		my $bases_above_ind = get_col_index($sam_cov_summ_f,"%_bases_above_".$cfg_hash->{'noncov_threshold'});
		#print_and_log( "Extracting cols $sample_id_ind,$bases_above_ind from $sam_cov_summ_temp\n",$log_file);
		#print "Extracting cols $sample_id_ind,$bases_above_ind from $sam_cov_summ_temp\n";
		extract_columns_from_file($sam_cov_summ_temp,"$sample_id_ind,$bases_above_ind",$sam_cov_summ_temp."2");
		
		my @cov_values = list_to_array($sam_cov_summ_temp."2","NO_NEW_LINE");
		
		foreach my $cov_val (@cov_values){
			my @fields = split(/\t/,$cov_val);
			if ($sample_name eq $fields[0] ){
				$retval = $fields[1];
				#print_and_log( "$sample_name coverage = $retval\n",$log_file);
				#print "$sample_name coverage = $retval\n";
			}
		}
		delete_file($sam_cov_summ_temp);
		delete_file($sam_cov_summ_temp."2");		
	}

	return $retval;
}









########################################################OLDCODE######################################################à

#=head2 parse_hpo_2_db

 #Title   : parse_hpo_2_db
 #Usage   : parse_hpo_2_db( -database => 'name of the database,
                               #);

 #Function: Parses a table like the following
					#analysisname	samplename	hpoid	observed
					#AAA						XXX					HP:3939	yes
					
					#gets the sampleid using the sample and analysis names 
					#and executes and insert_only_if_not_exists where the hpoid and sampleid 
					#must not exist yet

 #Returns : nothing

#=cut 
#sub parse_hpo_2_db {
	#my $cfg_hash = shift;
	#my $seq_info_f = shift;

	##Used to remember groups already inserted in db
	#my $samples_inserted;
	
	##Temporary hash for the header
	#my $head_hash;
	#my $sep = "\t";
	#print "Opening info file $seq_info_f\n";
	#open (FILE, "<$seq_info_f") or die "Cannot open $seq_info_f\n";
		
		##Go line by line and get the field values for each sample
		#my $row_num = 0;
		#while (my $line = <FILE>){
		##print "Analyzing line: ".$line;
				#chomp($line);
				
				########READ HEADER
				##Get the field names from the header
				#if ($row_num == 0){
					##Put fields in a way to be used
					#my @head = split($sep,$line);
					#my $n_field = 0;
					##each field will be associated to the index
					#foreach my $field (@head){
							#$head_hash->{$n_field} = $field;
							#$n_field++;
					#}
				#}else{#The row is composed of data
					#my @fields = split($sep,$line);
					#my $n_field = 0;
					#my $ss_hash;	
					##Use an hash (ss_hash) to keep all the values for this iteration
					##put the value where the field-index says at the location of the sample id given
					#foreach my $field (@fields){
							##print $head_hash->{$n_field}.":$field ";
							#$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							#$n_field++;
					#}

					################################
					##INSERT ANALYSIS INFO					#
					################################
					##Putting the analysisname and analysis id in the analysis table
					#my $familyname = $ss_hash->{'familyid'};
					#print "Considering $familyname \n";
					#if ( ! defined $samples_inserted->{$familyname} ){


						#my $personid = get_id_if_exists_from_db_like($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				#$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
																					#$cfg_hash->{'db_sample_name'},"'%".$familyname."_P'");									
						#if ($personid <= 0 ){
							#$personid = get_id_if_exists_from_db_like($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				#$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
																					#$cfg_hash->{'db_sample_name'},"'%".$familyname."_P1'");
						#}
							#my $phenid = -1;
							
							#if ($personid > 0 ){
								#my $fields = $cfg_hash->{'db_samples_personid'};
								#my $values = "$personid";									
								##Check 
								#my $hpo_id = $ss_hash->{$cfg_hash->{'db_hpo_id'}};
								#if ( $hpo_id ne '-' and  $hpo_id ne ''){
									##Get the hpoid
									#$phenid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																						#$cfg_hash->{'db_pass'},$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_phenotype_id'},
																						#$cfg_hash->{'db_hpo_id'},"'".$hpo_id."'");
									#if ( $phenid > 0){			
										#$fields .= ",".$cfg_hash->{'db_phenotype_id'};
										#$values .= ",$phenid";
									#}else{
										#print "WARNING: $hpo_id does not exist into the database ".$cfg_hash->{'db_name'}."\n";#DEBUGCODE
									#}
								#}	
								##Check 
								#my $hpo_observed = $ss_hash->{$cfg_hash->{'db_samples_hpo_observed'}};
								#if ( $hpo_observed ne '-' and  $hpo_observed ne '' ){
										#$hpo_observed = lc($hpo_observed);
										#if ( $hpo_observed eq 'yes' or $hpo_observed eq 'true' or $hpo_observed eq 'y' or $hpo_observed eq 't'){
											#$hpo_observed = 'TRUE';
										#}else{
											#$hpo_observed = 'FALSE';	
										#}
										#$fields .= ",".$cfg_hash->{'db_samples_hpo_observed'};
										#$values .= ",$hpo_observed";
								#}	

						
								##Now try to insert the value
								#if ( $phenid > 0){	
									##The second values are used to evaluate if the db_analysis_id is present
									#my $fields2 = $cfg_hash->{'db_samples_personid'}.",".$cfg_hash->{'db_phenotype_id'};
									#my $values2 = "$personid,$phenid";																
									#print "Inserting $fields = $values in table ".$cfg_hash->{'db_samples_hpo_table'}."\n";#DEBUGCODE
									#my ($new_group_id,$exist_group_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															#$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_hpo_table'},$cfg_hash->{'db_phenotype_id'},
															#$fields,$values,$fields2,$values2);														
								#}else{
									#print "Sample $samplename has not been analyzed and will not be considered\n";#DEBUGCODE
								#}
								##Save in the hash
								#$samples_inserted->{$samplename} = $sample_id;									
							#}
													
					#}
				#}#ELSE
		#$row_num++;
	#}
	#close(FILE);
#}	

	




#=head2 parse_samples_info_2_db2

 #Title   : parse_samples_info_2_db2
 #Usage   : parse_samples_info_2_db2( -database => 'name of the database,
                               #);

 #Function: Parses a table with multiple information for the samples
					#gets the sampleid using the sample and analysis names 
					#and executes an insert_only_if_not_exists where the analysisid and sampleid 
					#must not exist yet

 #Returns : nothing

#=cut 
#sub parse_samples_info_2_db2 {
	#my $cfg_hash = shift;
	#my $seq_info_f = shift;

	##Used to remember groups already inserted in db
	#my $samples_inserted;
	
	##Temporary hash for the header
	#my $head_hash;
	#my $sep = "\t";
	#print "Opening info file $seq_info_f\n";
	#open (FILE, "<$seq_info_f") or die "Cannot open $seq_info_f\n";
		
			##Go line by line and get the field values for each sample
			#my $row_num = 0;
			#while (my $line = <FILE>){
				##print "Analyzing line: ".$line;
				#chomp($line);
				
				########READ HEADER
				##Get the field names from the header
				#if ($row_num == 0){
					##Put fields in a way to be used
					#my @head = split($sep,$line);
					#my $n_field = 0;
					##each field will be associated to the index
					#foreach my $field (@head){
							#$head_hash->{$n_field} = $field;
							#$n_field++;
					#}
				#}else{#The row is composed of data
					#my @fields = split($sep,$line);
					#my $n_field = 0;
					#my $ss_hash;	
					##Use an hash (ss_hash) to keep all the values for this iteration
					##put the value where the field-index says at the location of the sample id given
					#foreach my $field (@fields){
							##print $head_hash->{$n_field}.":$field ";
							#$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							#$n_field++;
					#}
					##print "\n";
					################################
					##INSERT ANALYSIS INFO					#
					################################
					##Putting the analysisname and analysis id in the analysis table
					#my $samplename = $ss_hash->{$cfg_hash->{'db_sample_name'}};
					#my $analysisname = $ss_hash->{$cfg_hash->{'db_analysis_name'}};
					#print "Considering $analysisname and $samplename\n";
					#if ( ! defined $samples_inserted->{$samplename} ){

							##Get the analysisid
							#my $analysisid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},
																					#$cfg_hash->{'db_analysis_name'},"'".$analysisname."'");
							##Get the sampleid
							#my $sampleid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																					#$cfg_hash->{'db_sample_name'}.",".$cfg_hash->{'db_analysis_id'},"'".$samplename."',$analysisid");
																					
						
							#my $fields = $cfg_hash->{'db_sample_name'};
							#my $values = "'".$samplename."'";
							##The second values are used to evaluate if the sample is present
							#my $fields2 = $cfg_hash->{'db_sample_name'};
							#my $values2 = "'".$samplename."'";
							
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_dob'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_dob'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_dob'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_dob'}}."'";
							#}
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_pob'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_pob'}} ne ''){
								 #my $pob = $ss_hash->{$cfg_hash->{'db_sample_pob'}};
								 ##Remove strange chars
								 #$pob =~ s/[\'\-\.\+]//g;;
									#$fields .= ",".$cfg_hash->{'db_sample_pob'};
									#$values .= ",'$pob'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_kinshipinfo'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_kinshipinfo'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_kinshipinfo'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_kinshipinfo'}}."'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_gender'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_gender'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_gender'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_gender'}}."'";
							#}					
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_affected'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_affected'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_affected'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_affected'}}."'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_diseaseid'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_diseaseid'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_diseaseid'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_diseaseid'}}."'";
							#}

							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_notes'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_notes'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_notes'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_notes'}}."'";																	
							#}					
						
								#print "Inserting $fields = $values in table ".$cfg_hash->{'db_samples_info_table'}."\n";#DEBUGCODE
								##my ($personid,$exist_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													##$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
												##$fields,$values,$fields2,$values2);														
							  #my $personid = insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													#$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$fields,$values);
													
							##Now update the samples table with the coverage and the personid
							#if ( $sampleid > 0 and $analysisid > 0){

								###Get the file with the summary of samples coverage
								#my $work_fold = getcwd;	
								#my $sam_cov_summ_f = $work_fold."/".$analysisname."/".$cfg_hash->{'stats_out_f'}."/$analysisname\_DOC.".$cfg_hash->{'DOC_sam_sum'};
								#my $sample_cov = get_sample_target_coverage($sam_cov_summ_f,$samplename);		
								
								#print "Updating ".$cfg_hash->{'db_samples_personid'}." ($personid) for ".$cfg_hash->{'db_sample_id'}." $sampleid ($samplename) in table ".$cfg_hash->{'db_sample_table'}."\n";#DEBUGCODE
							 #update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
										#$cfg_hash->{'db_samples_personid'},"$personid",$cfg_hash->{'db_sample_id'},$sampleid);														
								#print "Updating ".$cfg_hash->{'db_sample_target_cov'}." ($sample_cov) for ".$cfg_hash->{'db_sample_id'}." $sampleid ($samplename) in table ".$cfg_hash->{'db_sample_table'}."\n";#DEBUGCODE
								#update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
										#$cfg_hash->{'db_sample_target_cov'},"$sample_cov",$cfg_hash->{'db_sample_id'},$sampleid);		
										
							#}else{
								#print "Sample $samplename has not been analyzed and will not be considered\n";#DEBUGCODE
							#}
							
							##Save in the hash
							#$samples_inserted->{$samplename} = $sample_id;																	
				#}	
		#}
		#$row_num++;
	#}
			#close(FILE);	
#}

			
#=head2 parse_samples_info_2_db

 #Title   : parse_samples_info_2_db
 #Usage   : parse_samples_info_2_db( -database => 'name of the database,
                               #);

 #Function: Parses a table with multiple information for the samples
					#gets the sampleid using the sample and analysis names 
					#and executes an insert_only_if_not_exists where the analysisid and sampleid 
					#must not exist yet

 #Returns : nothing

#=cut 
#sub parse_samples_info_2_db {
	#my $cfg_hash = shift;
	#my $seq_info_f = shift;

	##Used to remember groups already inserted in db
	#my $samples_inserted;
	
	##Temporary hash for the header
	#my $head_hash;
	#my $sep = "\t";
	#print "Opening info file $seq_info_f\n";
	#open (FILE, "<$seq_info_f") or die "Cannot open $seq_info_f\n";
		
			##Go line by line and get the field values for each sample
			#my $row_num = 0;
			#while (my $line = <FILE>){
				##print "Analyzing line: ".$line;
				#chomp($line);
				
				########READ HEADER
				##Get the field names from the header
				#if ($row_num == 0){
					##Put fields in a way to be used
					#my @head = split($sep,$line);
					#my $n_field = 0;
					##each field will be associated to the index
					#foreach my $field (@head){
							#$head_hash->{$n_field} = $field;
							#$n_field++;
					#}
				#}else{#The row is composed of data
					#my @fields = split($sep,$line);
					#my $n_field = 0;
					#my $ss_hash;	
					##Use an hash (ss_hash) to keep all the values for this iteration
					##put the value where the field-index says at the location of the sample id given
					#foreach my $field (@fields){
							##print $head_hash->{$n_field}.":$field ";
							#$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							#$n_field++;
					#}
					##print "\n";
					################################
					##INSERT ANALYSIS INFO					#
					################################
					##Putting the analysisname and analysis id in the analysis table
					#my $samplename = $ss_hash->{$cfg_hash->{'db_sample_name'}};
					#my $analysisname = $ss_hash->{$cfg_hash->{'db_analysis_name'}};
					#print "Considering $analysisname and $samplename\n";
					#if ( ! defined $samples_inserted->{$samplename} ){

							##Get the analysisid
							#my $analysisid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},
																					#$cfg_hash->{'db_analysis_name'},"'".$analysisname."'");
							##Get the sampleid
							#my $sampleid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																					#$cfg_hash->{'db_sample_name'}.",".$cfg_hash->{'db_analysis_id'},"'".$samplename."',$analysisid");

							##Get the file with the summary of samples coverage
							#my $work_fold = getcwd;	
							#my $sam_cov_summ_f = $work_fold."/".$analysisname."/".$cfg_hash->{'stats_out_f'}."/$analysisname\_DOC.".$cfg_hash->{'DOC_sam_sum'};
							#my $sample_cov = get_sample_target_coverage($sam_cov_summ_f,$samplename);
							
							#my $fields = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_id'};
							#my $values = "$analysisid,$sampleid";
							##The second values are used to evaluate if the analysisname is present
							#my $fields2 = $cfg_hash->{'db_analysis_id'};
							#my $values2 = "$analysisid";
						
							##Check 
							#if ( $sample_cov > 0 ){
									#$fields .= ",".$cfg_hash->{'db_sample_target_cov'};
									#$values .= ",$sample_cov";
							#}															
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_realname'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_realname'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_realname'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_realname'}}."'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_surname'}} ne '-' and   $ss_hash->{$cfg_hash->{'db_sample_surname'}} ne '' ){
									#$fields .= ",".$cfg_hash->{'db_sample_surname'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_surname'}}."'";
							#}									
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_dob'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_dob'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_dob'};
									#$values .= ",".$ss_hash->{$cfg_hash->{'db_sample_dob'}};
							#}
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_pob'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_pob'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_pob'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_pob'}}."'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_kinshipinfo'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_kinshipinfo'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_kinshipinfo'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_kinshipinfo'}}."'";
							#}					
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_tissue'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_tissue'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_tissue'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_tissue'}}."'";
							#}
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_affected'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_affected'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_affected'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_affected'}}."'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_paternity'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_paternity'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_paternity'};
									#$values .= ",".$ss_hash->{$cfg_hash->{'db_sample_paternity'}};
							#}					
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_heritability'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_heritability'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_heritability'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_heritability'}}."'";									
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_diseaseid'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_diseaseid'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_diseaseid'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_diseaseid'}}."'";
							#}
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_target_cov'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_target_cov'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_target_cov'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_target_cov'}}."'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_enrichment'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_enrichment'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_enrichment'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_enrichment'}}."'";
							#}					
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_readindex'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_readindex'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_readindex'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_readindex'}}."'";									
							#}		
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_code'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_code'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_code'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_code'}}."'";
							#}
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_cgharray'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_cgharray'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_cgharray'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_cgharray'}}."'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_phenomeid'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_phenomeid'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_phenomeid'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_phenomeid'}}."'";
							#}					
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_sequencer'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_sequencer'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_sequencer'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_sequencer'}}."'";																	
							#}
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_arrivaldate'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_arrivaldate'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_arrivaldate'};
									##$values .= ",to_date('".$ss_hash->{$cfg_hash->{'db_sample_arrivaldate'}}."', 'DD-MM-YYYYD')";
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_arrivaldate'}}."'";
							#}
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_library_prep_date'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_library_prep_date'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_library_prep_date'};
									##$values .= ",to_date('".$ss_hash->{$cfg_hash->{'db_sample_library_prep_date'}}."', 'DD-MM-YYYYD')";
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_library_prep_date'}}."'";
									
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_further_lib_prep_date'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_further_lib_prep_date'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_further_lib_prep_date'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_further_lib_prep_date'}}."'";
							#}					
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_sequencing_date'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_sequencing_date'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_sequencing_date'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_sequencing_date'}}."'";		
							#}	
								##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_further_seq_dates'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_further_seq_dates'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_further_seq_dates'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_further_seq_dates'}}."'";
							#}
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_libraryprep'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_sample_libraryprep'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_libraryprep'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_libraryprep'}}."'";
							#}	
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_othersampleid'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_othersampleid'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_othersampleid'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_othersampleid'}}."'";
							#}					
							##Check 
							#if ( $ss_hash->{$cfg_hash->{'db_sample_notes'}} ne '-'  and  $ss_hash->{$cfg_hash->{'db_sample_notes'}} ne ''){
									#$fields .= ",".$cfg_hash->{'db_sample_notes'};
									#$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_notes'}}."'";																	
							#}					


							
							#if ( $sampleid > 0 and $analysisid > 0){
								##The second values are used to evaluate if the analysisname is present
								#my $fields2 = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_id'};
								#my $values2 = "$analysisid,$sampleid";								
								##print "Inserting $fields = $values in table ".$cfg_hash->{'db_samples_info_table'}."\n";#DEBUGCODE
								##my $sample_info_id = insert_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														##$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$fields,$values);
								#my ($new_id,$exist_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													#$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_analysis_id'},
													#$fields,$values,$fields2,$values2);														
							#}else{
								#print "Sample $samplename has not been analyzed and will not be considered\n";#DEBUGCODE
							#}
							
							##Save in the hash
							#$samples_inserted->{$samplename} = $sample_id;
																						
				#}	
		#}
		#$row_num++;
	#}
			#close(FILE);	
#}
