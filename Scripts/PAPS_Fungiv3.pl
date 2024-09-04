#!/usr/bin/perl -w
#Jordi Paps, Oxford 2015; with help of Patrick Gemmell (Oxford) for some subroutines
#Beginning Perl for Bioinformatics
use DB_File;
use strict;
use AnyDBM_File;
use Fcntl;
use Term::ANSIColor;

#Script to extract groups of homology from a parsed MCL output, based on the taxonomic distribution

#Introduce file names
my $MCL_out = "Input/Orthogroups.txt"; 
my $MCL_columns_parsed = "Input/02_hogroups_gene_numbers_parsed.out";

#Check if files exists, open them
check_file ($MCL_out);
check_file ($MCL_columns_parsed);

#Create DBM files from the MCL output file for fast lookups.
my %MCL_out;
tie (%MCL_out,"AnyDBM_File", "$MCL_out.db", O_CREAT|O_RDWR, 0777) or die + "Can't open database: $!\n";
%MCL_out = MCL_output_to_hash ($MCL_out);
print "\nMCL out (showing some random lines out of ", scalar (keys %MCL_out), "):\n";
my @MCL_out_keys = keys %MCL_out;
foreach (@MCL_out_keys[0..4]){
	my $first_elements = substr ($MCL_out{$_}, 0, 80);
	print "$_ => $first_elements...\n";
} print "...\n";

#Create DBM files from the MCL parsed file for fast lookups.
my %MCL_columns_parsed;
tie (%MCL_out, "DB_File", "$MCL_out.db", O_CREAT|O_RDWR, 0777) or die "Can't open database: $!\n";
%MCL_columns_parsed = MCL_columns_parsed_to_hash ($MCL_columns_parsed);
print "\nMCL parsed (each column is one species, showing some random lines out of ", scalar (keys %MCL_columns_parsed), "):\n\n";
print "Species: Rirr Depi Odia Gcer Gros Gmar Rcla Gpyr Alep Fcal Amor Cpel Ager Dery Fmos Debu Ecan Pbras Fgeo Esp Pocc Razy Awol Aamb Lelo Pver\n\n";
my @MCL_columns_parsed_keys = keys %MCL_columns_parsed;
foreach (@MCL_columns_parsed_keys[0...4]){
	print "$_ => $MCL_columns_parsed{$_}\n";
} print "...\n";

#Create hash from the clade definition in subrout CLADES
my %spp;
my @value_list = ();
%spp = hash_spp();

#Calculate how many columns (species/terminal tip) are in the hash
walk_hash (\%spp, \@value_list);
my @total_columns = scalar @value_list;
print "\nNumber of species: ", (scalar @value_list), "\n";

#Ask user for clade/spp to check genes taxonomic distribution, perform search
my $user_input = '';	#To store the searching criteria fom user
my @arguments = ();	#To store the split of the searching criteria
my @search= ();  	#To store the taxa and options of the searching criteria
my $taxa  = '';		#To store the taxa from @search
my $option = '';	#To store options from @search (present, absent, minus, atleast,...)
my @columns = ();	#To store columns to search
my $final_search;	#To store columns to search plus the options (present, absent, minus, atleast,...)
my %final_searches = ();#To store ALL searches columns and options (present, absent, minus, atleast,...)
my @true_flags = (); 	#To store the flags that indicate if a homology group/line passes all the queries checks
my @outgroup = ();	#To store the columns left over at the end after extracting all the ingroups columns
my @good_homology_groups_spp_names = ();	#To store the groups of homology fullfilling the search criteria
my @good_homology_groups_spp_names_copy = ();	#Backup
my @good_homology_groups_columns_parsed = ();	#To store the columns of the groups of homology fullfilling the search criteria

OUTER: do {
	print "\nPlease, enter name of the clade/species to search (\"example\" for some samples, \"tree\" to print the evolutionary tree, Enter to exit):\n";
	$user_input = <STDIN>;
	chomp $user_input;
	unless ($user_input =~ /^\s*$|^exit|^quit/i) {			#Unless user wants to exit...
		if ($user_input =~ /example|help/i) {
			print_examples();			#Print examples of commands if user requests it
		} elsif ($user_input =~ /tree/i) {
			print_hash_colors(%spp);		#Print examples of commands if user requests it
		} else {
			#Here the real search starts, emptying variables for next loop and parsing the user input
			%final_searches = ();			#Empty the hash containing all the search conditions
			@good_homology_groups_spp_names = ();	#Empty the hash containing the results from previous search
			@good_homology_groups_columns_parsed = (); #Empty the hash containing the columns of the groups of homology fullfilling the search criteria
			@arguments = split (" ", $user_input);	#Decompose the user input in different arguments entered by user, each taxa query delimited by a space (" ")
			#print "\@arguments: @arguments\n";			
			foreach (@arguments) {
				#print "\$_ in \@arguments: $_\n";
				@search = split ("-", $_);	#Decompose each argument into taxa (first item) and options (second item)
				#print "\@search: @search\n";
				$taxa = $search[0];
				$option = $search[1];
				unless (defined $option && $option =~ /pre|all|abs|none|min|but|atl|onl|jus/gi){
						print "Taxa $taxa is missing valid options (present, absent, minus, atleast, only,...)\n";
						goto OUTER;
				}
				if ($taxa =~ /^out|^rest|^other/i) {		#We store outgroup conditions in the hash, as "outgroup" does not exists in the species hash of hashes
					$final_searches {$taxa."_".$option} = "outgroup_".$option;       #Store ALL the searches in a hash, keys are taxa and values the columns and options
				} else {
					@columns = obtain_taxa_columns (\$taxa, \%spp, \%MCL_columns_parsed); #Obtain the columns belonging to each taxa
					#print "Columns that will be inspected: @columns\n";
					if (scalar @columns == 0 ){
						print "Taxa $taxa not found in taxa list\n";
						goto OUTER;
					}
					unless (defined $search[1] && $search[1] =~ /pre|all|abs|none|min|but|atl|onl|jus/gi){
						print "Taxa $taxa is missing valid options (present, absent, minus, atleast, only)\n";
						goto OUTER;
					}
					$final_search = join ("_", @columns, $option); 		#Join the columns and the options to store them later in %final_search
					$final_searches {$taxa."_".$option} = $final_search;	#Store ALL the searches in a hash, keys are taxa and values the columns and options
				}	
			}
			
			#Read the hash containing the file with the MCL output parsed to columns, check the search criteria
			my @keys_MCL_columns_parsed = keys %MCL_columns_parsed;			
			#print "Colums to be inspected: @keys_MCL_columns_parsed\n";
			my $HG_counter = 0;
			foreach my $homology_group (@keys_MCL_columns_parsed){ 				#For each line of the MCL parsed file...
				#print "\n", '$homology_group: ', $homology_group, "\n";
				@true_flags = ();							#Empty the flags that indicate if a homology group/line passes all the queries checks
				my @MCL_columns = split (" ", $MCL_columns_parsed{$homology_group});	#Explode the string of the homology group to columns				
				@outgroup = @MCL_columns;						#Storing all columns in @outgroups, later the ingroups columns will be spliced out o this array
				my @keys = keys %final_searches;
				#print 'keys in %final_searches: ', "@keys" , "\n";
				#print "MCL columns_parsed:\n@MCL_columns\n";
				foreach my $query (keys %final_searches) {				#Now go query from query stored in %final_searches
					#print '$query ', $query,"\n";
					unless ($query =~ /outg/i) {					#Leave outgroup check for the end
						my @query = split ("_", $final_searches{$query});	#Explode the string of the columns to check
						#print "\@query: @query ";
						my $condition = pop @query;				#Save the options (present, absent, minus, atleast,...)
						#print "\$condition: $condition \n";
						@query = @query;				
						#print "\@query columns:  @query\n";
						my $total_taxa_in_clade = scalar @query;		#Stores total number of taxa in queried clade, used to later check options
						my $sum = 0;
						#print "Columns contens: ";
						foreach my $query_column (@query) {			#For each column corresponding to the taxa in the query, we extract the value of the column in MCL file
							#print "$MCL_columns[$query_column] ";							
							if ($MCL_columns[$query_column] != 0) {
								$sum++;
							}
							splice (@outgroup, $query_column, 1, "X"); 	#Remove column from outgroup array, to check later the remaining columns
						}
						#print "\nTotal taxa in clade:$total_taxa_in_clade\n";
						#print "Taxa of clade in this MCL: $sum\n";
						#Check if the group of homology matches the conditions requested by user for this clade
						my $check = check_conditions($condition, $total_taxa_in_clade, $sum);
						if ($check eq "wrong") {
							goto OUTER;
						} else {
							push (@true_flags, $check);
						}
					}
				}
				#OUTGROUP CHECK
				foreach (keys %final_searches) {				#Now we check if user asked for outgroup conditions and is stored in %final_searches
					if ($_ =~ /outg/i) {
						#print "OUTGROUP\n";
						#print "\@outgroup:\n@outgroup\n";
						#print "\@MCL_columns:\n@MCL_columns\n";
						my @query = split ("_", $final_searches{$_});	#Explode the string of the columns to check
						my $condition = pop @query;				#Save the options (present, absent, minus, atleast,...)
						#print '$condition in outgroup: ', $condition,"\n";
						@outgroup = grep { $_ ne "X" } @outgroup;
						#print 'MCL contens in outgroup: ', "\n@outgroup","\n";
						my $total_taxa_in_clade = scalar @outgroup;		#Stores total number of taxa in queried clade, used to later check options
						#print 'Total taxa in outgroup: ', $total_taxa_in_clade, "\n";
						my $sum = 0;
						#print "\$MCL_columns[\$query_column] in outgroup:\n";
						foreach (@outgroup) {			#For each column corresponding to the taxa in the query, we extract the value of the column in MCL file
							#print "$_ ";
							if ($_ != 0) {
								$sum++;
							}
						}
						#print 'Outgroup in this MCL: ', $sum, "\n";
						my $check = check_conditions($condition, $total_taxa_in_clade, $sum);
							if ($check eq "wrong") {
								goto OUTER;
							} else {
								push (@true_flags, $check);
							}
					}
				}
				#print "\@true_flags: @true_flags\n";
				#Check if all flags are true, then store the group of homology in %results
				my $flags = 0;
				
				foreach (@true_flags) {
					if ($_ !~ 'true' ) {
						$flags++;
					}					
				}
				if ($flags == 0) {
					push (@good_homology_groups_spp_names, "$homology_group\t$MCL_out{$homology_group}\n");
					push (@good_homology_groups_columns_parsed, "$homology_group\t$MCL_columns_parsed{$homology_group}\n");
					$HG_counter++;
				}				
			}
			@good_homology_groups_spp_names_copy = @good_homology_groups_spp_names;
			print "\nNumber of groups of homology found: $HG_counter \n";
			if ($HG_counter == 0) { goto OUTER; };
			
			#Save 4 different output files
			print "\nDo you want to see results (gene names, MCL groups, etc) and save them in files? Yes/No\n";
			my $save_files = <STDIN>;
			chomp $save_files;
			if ($save_files =~ /y/gi) {
				my $output_filename = join ("_", @arguments);
				$output_filename = "Output/".$output_filename."_$HG_counter\_HGs";
				
				# 1) Save the groups of homology that match the query, with spp names
				#print "\nShowing first sequence names for first groups of homology:\n";
				#foreach (@good_homology_groups_spp_names) {
				#	print substr($_, 0, 160), "\n\n";
				#}
				unless (open (OUTPUT1, ">$output_filename"."_MCL_genes_IDs.out")) {
					print "Can't open file to save";
				}
				print OUTPUT1 @good_homology_groups_spp_names;
				close OUTPUT1;
				
				# 2) Save the columns from MCL parsed file for groups of homology that match the query
				#print "\nShowing columns for the few first groups of homology:\n";
				#print "Species: Rirr Depi Odia Gcer Gros Gmar Rcla Gpyr Alep Fcal Amor Cpel Ager Dery Fmos Debu Ecan Pbras Fgeo Esp Pocc Razy Awol Aamb Lelo Pver\n\n";
				#foreach (@good_homology_groups_columns_parsed) {
				#	print $_;
				#}
				unless (open (OUTPUT2, ">$output_filename"."_MCL_columns_parsed.out")) {
					print "Can't open file to save";
				}
				print OUTPUT2 "HG\tRirr\tDepi\tOdia\tGcer\tGros\tGmar\tRcla\tGpyr\tAlep\tFcal\tAmor\tCpel\tAger\tDery\tFmos\tDebu\tEcan\tPbras\tFgeo\tEsp\tPocc\tRazy\tAwol\tAamb\tLelo\tPver\tn";
				print OUTPUT2 @good_homology_groups_columns_parsed;
				close OUTPUT2;
				
				# 3) Now save in return format, one taxa per line
				#@good_homology_groups_spp_names = @good_homology_groups_spp_names_copy;
				my @gene_names_return = ();
				print "\nShowing first sequence names for groups of homology:\n";
				foreach (@good_homology_groups_spp_names) {
					chomp $_;
					my @homology_group = split (" ", $_);
					my $gene_names = parse_gene_names_return_format(@homology_group);
					push (@gene_names_return, $gene_names);
					print "\n",substr($gene_names, 0, 480),"...\n";
				}
				unless (open (OUTPUT3, ">$output_filename"."_MCL_annotated_genes.out")) {
					print "Can't open file to save";
				}
				print OUTPUT3 @gene_names_return;
				close OUTPUT3;
				
				# 4) Save the names of the taxa present in the groups of homology 
				my @taxa_names = ();
				foreach (@good_homology_groups_spp_names) {
					chomp $_;
					my @homology_group = split (" ", $_);
					push (@taxa_names, parse_taxa_names(@homology_group));
				}
				#print "\nShowing the labels of all the taxa for each homology group:\n";
				#foreach (@taxa_names) {
				#	print "$_";
				#	}
				unless (open (OUTPUT4, ">$output_filename"."_taxa_names.out")) {
					print "Can't open file to save";
				}
				print OUTPUT4 @taxa_names;
				close OUTPUT4;
				print "\nResults saved to files $output_filename\n";
			}	
		}
	}
} until ($user_input =~ /^\s*$|^exit|^quit/i);

#Untie hash files
untie %MCL_out;
untie %MCL_columns_parsed;

#End of program
print color ("red"), "\nend\n\n", color ("reset");
exit;


################################################################################
#                                Subroutines                                   #
################################################################################

#CHECK_FILE: checks file properties, checks if file exists, is flat, is empty, or can be opened
#-e exists, -f flat file, -s empty
sub check_file {
    my ($filename) = @_;
    unless (-e $filename) {print "File $filename does NOT exists\n" and exit;}
    unless (-f $filename) {print "File $filename is NOT a flat file\n" and exit;}
    unless (-s $filename) {print "File $filename is empty\n" and exit;}
    unless (open (FH, $filename)) {print "File $filename can not be opened\n" and exit;}
    close FH;
    return;
}

#MCL_COLUMNS_PARSED_TO_HASH: subrout to parse a MCL output file and store it in a hash, in which the
#keys are the group homology and the values the genes found in that group
sub MCL_columns_parsed_to_hash {
	my ($filename) = @_;
	my %hash;
	my @line = '';
	my $line_counter = 0;
	
	open(FH, $filename);
	my @file = <FH>;
	foreach my $line (@file) {   	#Parse file one line at a time
		$line_counter++;
		chomp $line;
		$line =~ s/^\d*\s//;	#To remove first digits, which indicate the group of homology number/ID
		#$line =~ s/\s/\t/g;
		#$line =~ s/\t$//g;
		my $key = "HG_".$line_counter;
		$hash{$key} = "$line";	
	}
	close FH;
	return %hash;
}

#MCL_TO_HASH: subrout to parse a MCL output file and store it in a hash, in which the
#keys are the group homology and the values the genes found in that group
sub MCL_output_to_hash {
	my ($filename) = @_;
	my %hash;
	my $line_counter = 0;
	
	open(FH, $filename);
	my @file = <FH>;
	foreach my $line (@file) {   			#Parse file one line at a time
		$line_counter++;
		chomp $line;
		#$line =~ s/\s/\t/g;
		#$line =~ s/\t$//g;
		my $key = "HG_".$line_counter;
		$hash{$key} = "$line";	
	}
	close FH;
	return %hash;
}

#WALK_HASH: subroutine to traverse the hash of hashes, modified after in http://www.perlmonks.org/?node_id=116162
sub walk_hash { 
my ($hash, $value_list, $key_list) = @_;
while (my ($key, $value) = each %$hash) {
	push @$key_list, $key;
	if (ref($value) eq 'HASH') {
		walk_hash ($value, $value_list, $key_list);
	} else {
		push @$value_list, $value;
	}
	pop @$key_list;
	}
}

#QUERY: intermediate subroutine that sends the taxa search to the to recoursive subroutines PG_WALK and PRINT_EVERYTHING
sub obtain_taxa_columns { 
	my ($taxa, $spp, $MCL_parsed) = @_; 
	my @results = ();
	
	#Send query and hash      

	pg_walk (\%$spp, [], $$taxa, \@results);
	@results = sort {$a <=> $b} @results;
	return @results; 
}

#PG_WALK: subroutine to traverse the hash of hashes till finding the queried label
#Hat tip to Patrick Gemmell (Oxford) for the help, he modified the subroutine
#found in http://www.perlmonks.org/?node_id=116162
sub pg_walk {
	my ($hash, $key_list, $query_text, $localresults) = @_;
	while (my ($key, $value ) = each %$hash) {
		push @$key_list, $key;
		if ($key =~ /^$query_text/gi) {
			print "Taxa that will be searched: $key\n";
			print_everything($value , $localresults);
		} else {
			if (ref($value ) eq 'HASH') {
				pg_walk($value , $key_list, $query_text, $localresults);
			}
		}
		pop @$key_list;
    }
}

#PRINT_EVERYTHING: subroutine to traverse the hash of hashes to print the values corresponding to the query in sub PG_WALK (see above).
#Hat tip to Patrick Gemmell (Oxford) for the help, he modified the subroutine
#found in http://www.perlmonks.org/?node_id=116162
sub print_everything {
	my ($hash, $localresults, $key_list) = @_;
	while (my ($key, $value) = each %$hash) {
		push @$key_list, $key;
		if (ref($value) eq 'HASH') {
			print_everything($value, $localresults, $key_list);
		} else {
			push @$localresults, $value;
		}
		pop @$key_list;
	}
}

#CHECK_CONDITIONS: subroutine to check the conditions of presence/absence specified by the user, returns TRUE or FALSE
sub check_conditions {
	my ($condition, $total_taxa_in_clade, $sum) = @_;
	#Perform the check according to the different options (present, absent, minus, atleast, only...
	if ($condition =~ /pre|all/gi) {
		if ($sum == $total_taxa_in_clade) {
			#print "True all\n";
			return "true";
		} else {
			#print "False all\n";
			return "false";
		}
	} elsif ($condition =~ /abs|none/gi) {
		if ($sum == 0) {
			#print "True none\n";
			return "true";
		} else {
			#print "False none\n";
			return "false";
		}
	} elsif ($condition =~ /atl\D*/gi ) {
		$condition =~ s/$&//g;
		#print "\$condition: $condition\n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"at least\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum >= $condition) {
			#print "True atleast $condition\n";
			return "true";
		} else {
			#print "Atleast $condition false\n";
			return "false";
		}
	} elsif ($condition =~ /min\D*|but\D*/gi ) {
		$condition =~ s/$&//g;							
		#print "\$condition: $condition \n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"minus\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum == ($total_taxa_in_clade - $condition)) {
			#print "True minus $condition\n";
			return "true";
		} else {
			#print "Minus $condition false\n";
			return "false";
		}
	} elsif ($condition =~ /onl\D*|jus\D*/gi ) {
		$condition =~ s/$&//g;								
		#print "\$condition: $condition \n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"only\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum == $condition) {
			#print "True only $condition\n";
			return "true";
		} else {
			#print "Only $condition false\n";
			return "false";
		}
	}
}

#PARSE_GENE_NAMES_TEXT_FORMAT: subrout to extract the gene names from the annotated genomes for the list of groups of homology
# producing an output in text format
sub parse_gene_names_text_format {
	my (@homology_group) = @_;
	my @list_annotated_genomes = qw /Rirr Depi Odia Gcer Gros Gmar Rcla Gpyr Alep Fcal Amor Cpel Ager Dery Fmos Debu Ecan Pbras Fgeo Esp Pocc Razy Awol Aamb Lelo Pver/;
	#@list_annotated_genomes = reverse @list_annotated_genomes;
	my $results;

	foreach my $annotated_taxa (@list_annotated_genomes) {
		#print "\$gene: $gene\n";
		foreach my $gene (@homology_group) {
			#print "\$annotated_taxa: $annotated_taxa ";
			if ($gene =~ /$annotated_taxa/i) {
				#print color ("red"), "Match!\n", color ("reset");
				$results .= "$gene ";
			}
		}
	}
	#print "\@results in sub: $results\n";
	if ($results =~ /^HG_\d*\s*$/) {
		$results .= "This group of homology does not contain any annotated genome";
	}
	#print "Results @results\n";
	#@results = sort @results;
	$results .= "\n";
	return $results;
}

#PARSE_GENE_NAMES_RETURN_FORMAT: subrout to extract the gene names from the annotated genomes for the list of groups of homology
# producing an output in a format with returns
sub parse_gene_names_return_format {
	my (@homology_group) = @_;
	my @list_annotated_genomes = qw /Rirr Depi Odia Gcer Gros Gmar Rcla Gpyr Alep Fcal Amor Cpel Ager Dery Fmos Debu Ecan Pbras Fgeo Esp Pocc Razy Awol Aamb Lelo Pver/;
	#@list_annotated_genomes = reverse @list_annotated_genomes;
	my $results;
	my $flag = 0;

	foreach my $annotated_taxa (@list_annotated_genomes) {
		$flag = 0;
		#print "\$gene: $gene\n";
		foreach my $gene (@homology_group) {
			#print "\$annotated_taxa: $annotated_taxa ";
			if ($gene =~ /$annotated_taxa/i) {
				#print color ("red"), "Match!\n", color ("reset");
				$results .= "$gene\t";
				$flag = 1;
			}
		}
		if ($flag == 1) { $results .= "\n"; }
	}
	
	#print "\@results in sub: $results\n";
	if ($results =~ /^HG_\d*\s*$/) {
		$results .= "This group of homology does not contain any annotated genome\n";
	}
	#print "Results @results\n";
	#@results = sort @results;
	#$results .= "\n";
	return $results;
}

#PARSE_TAXA_NAMES: subrout to extract the taxa names from the list of groups of homology
sub parse_taxa_names {
	my (@homology_group) = @_;
	my %taxons;
	my $results;
	my $group_ID = splice (@homology_group, 0 ,1);

	foreach my $taxon (@homology_group) {
		$taxon =~ /^.{4}/;
		$taxons{$&} = '';
	}
	my @keys = sort keys %taxons;
	$results = join ("\t", @keys);
	$results = "$group_ID "."\t$results"."\n";
	return $results;
}

#PRINT_EXAMPLES: subroutine to print some examples to user
sub print_examples {
	print '
Clade/species names can be truncated, but the start of the clade name should match the table printed above.
Search is case insensitive.

Some search examples (first 4 digits in examples stand for rest of taxa, the other 4 for ingroup):	
  "Vertebrata-present" => genes found in ALL vertebrate species, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1111
  "Vertebrata-present Rest-absent" => genes found in ALL vertebrate species, absent in other clades/rest of taxa
	Rest of taxa 0000 Ingroup 1111
  "Vertebrata-present Rest-present" => genes found in ALL vertebrate species, present in other clades/rest of taxa
	Rest of taxa 1111 Ingroup 1111
  "Vertebrata-absent Rest-present" => genes found in rest of taxa species, absent in Vertebrata
	Rest of taxa 1111 Ingroup 0000
  "Homo-present Mus-present Rest-absent" => genes only found in humans and mice. Species can be specified one by one.

The number of species presenting/missing for a gene can be fine-tuned with minus#, atleast#, only# for both ingroup and rest of taxa:
  "Vertebrata-minus1" => found in ALL vertebrate species but one, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1110 / 1101 / 1011 / 0111
  "Vertebrata-minus2 Rest-minus1" => genes found in ALL vertebrate species but one, absent in other clades/rest of taxa
	Rest of taxa 1110 / 1101 / 1011 / 0111 Ingroup 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Vertebrata-atleast1 Rest-atleast1" => genes found in at least 1 vertebrate species and 1 rest of taxa species
	Rest of taxa 1000 / 1100 / 1110 / 1111 / 1010 / 1011 / 1001 / 1101 / 0110 / 0111 etc.
	Ingroup  1000 / 1100 / 1110 / 1111 / 1010 / 1011 / 1001 / 1101 / 0110 / 0111 etc.
  "Vertebrata-only3" => return genes found in just 3 vertebrate species, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1110 / 1101 / 1011 / 0111
  
Different criteria can be combined in a single search:
  "Vertebrata-minus1 Echinodermata-atleast2" => genes found in ALL vertebrate species but one, AND present in at least two echinoderms, absent/present in other clades/rest of taxa
	Rest of taxa ???? Vertebrata 1110 / 1101 / 1011 / 0111 Echinodermata 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Vertebrata-atleast2 Urochordata-atleast2" => genes found in 2 or more vertebrate species OR 2 or more urochordates, independently if they are present/absent in other clades/rest of taxa
	Rest of taxa ???? Vertebrata 1100 / 1010 / 1001 / 0110 / 0101 / 0011 Urochordata 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Nematoda-absent Platyhelminthes-absent Rest-present" => genes found in clades/rest of taxa, absent (convergently lost) in round worms and flatworms
	Rest of taxa 1111 Nematoda 0000 Platyhelminthes 0000
	
Carefull with nested taxa!!! Start with the greater group taking into account the conditions for the smaller group:
  To find genes in ALL chordates but missing only in humans => "Chordata-minus1 Hsap-absent"
  To find genes in ALL chordates but missing only in vertebrates => "Chordata-minus5 Vertebrata-absent"
  To find genes in at least one clade of chordates, but missing only in vertebrates => "Cephalocordata-atleast1 Urochordata-atleast1 Vertebrata-absent"
  ';
	return;
}

#HASH_SPP: a subroutine to define define the hash of hashes containing of all the clades and spp included, put here to not clutter the main program.
#Each species is assigned a numeric value, same as the column they occupy in the parsed MCL output, so Crei is the first column and Hsap occupies the last column.
#Thus, when user asks for a group, these values can be used as index to lookup lines/arrays.
#Each element shoud have the same number of levels, or the subrout "PRINT_HASH_COLORS" won't work.
sub hash_spp {
	my (%spp) = ();

#All spp, one by one
## SUPERKINGDOM ## CLADE ## KINGDOM ## PHYLUM ## SUBPHYLUM ## CLASS ## Node ## ORDER ## SPECIES ##
#{'Eukaryota'}{'Opisthokonta'}{'Holozoa'}{'Node1'}{'Node2'}{'Metazoa'}{'Eumetazoa'}{'Planulozoa'}{'Bilateria'}{'Deuterostomia'}{'Chordata'}{'Olfactores'}{'Vertebrata'}{'Homo_sapiens_(Hsap)'}{''} =  61;

	#Acaulosporaceae #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Diversisporales'}{'DiversiasporalesG1'}{'Acaulosporaceae'}{'Acaulospora_morrowiae_(Amor)'}{''} = 0;
	
	#Ambisporaceae #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1A'}{'Archaeosporales'}{'1AA'}{'Ambisporaceae'}{'Ambispora_leptoticha_(Alep)'}{''} = 1;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1A'}{'Archaeosporales'}{'1AA'}{'Ambisporaceae'}{'Ambispora_gerdemannii_(Ager)'}{''} = 2;

	##Diversisporaceae #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Diversisporales'}{'DiversiasporalesG1'}{'Diversisporaceae'}{'Diversispora_epigaea_(Depi)'}{''} = 3;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Diversisporales'}{'DiversiasporalesG1'}{'Diversisporaceae'}{'Diversispora_eburnea_(Debu)'}{''} = 4;

	#Entrophospora #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1A'}{'Entrophosporales'}{'1AB'}{'Entrophosporaceae'}{'Entrophospora_candida_(Ecan)'}{''} = 5;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1A'}{'Entrophosporales'}{'1AB'}{'Entrophosporaceae'}{'Entrophospora_spSA101_(Esp)'}{''} = 6;

	#Geosiphon #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{''}{'Archaeosporales'}{'1AA'}{'Geosiphonaceae'}{'Geosiphon_pyriformis_(Gpyr)'}{''} = 7;
	
	#Gigasporaceae #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Diversisporales'}{'DiversiasporalesG2'}{'Gigasporaceae'}{'Gigaspora_rosea_(Gros)'}{''} = 8;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Diversisporales'}{'DiversiasporalesG2'}{'Gigasporaceae'}{'Gigaspora_margarita_(Gmar)'}{''} = 9;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Diversisporales'}{'DiversiasporalesG2'}{'Gigasporaceae'}{'Cetraspora pellucida_(Cpel)'}{''} = 10;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Diversisporales'}{'DiversiasporalesG2'}{'Gigasporaceae'}{'Dentiscutata_erythropus_(Dery)'}{''} = 11;

	#Glomeraceae #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Glomerales'}{'1BA'}{'Glomeraceae'}{'Rhizophagus_irregularis_(Rirr)'}{''} = 12;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Glomerales'}{'1BA'}{'Glomeraceae'}{'Oehlia diaphana_(Odia)'}{''} = 13;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Glomerales'}{'1BA'}{'Glomeraceae'}{'Glomus_cerebriforme_(Gcer)'}{''} = 14;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Glomerales'}{'1BA'}{'Glomeraceae'}{'Rhizophagus_clarus_(Rcla)'}{''} = 15;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Glomerales'}{'1BA'}{'Glomeraceae'}{'Funneliformis_caledonius_(Fcal)'}{''} = 16;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Glomerales'}{'1BA'}{'Glomeraceae'}{'Funneliformis_mosseae_(Fmos)'}{''} = 17;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1B'}{'Glomerales'}{'1BA'}{'Glomeraceae'}{'Funneliformis_geosporum_(Fgeo)'}{''} = 18;

	#Paraglomeraceae #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1A'}{'Paraglomerales'}{'1AA'}{'Paraglomeraceae'}{'Paraglomus_brasilianum_(Pbras)'}{''} = 19;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'AMF'}{'Glomeromycotina'}{'Glomeromycetes'}{'Split1A'}{'Paraglomerales'}{'1AA'}{'Paraglomeraceae'}{'Paraglomus_occultum_(Pocc)'}{''} = 20;

	#Mucoromycotina #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'Outgroup2'}{'Mucoromycotina'}{' Mucoromycetes'}{''}{'Mucorales'}{''}{'Rhizopodaceae'}{'Rhizopus_azygosporus_(Razy)'}{''} = 21;

	#Mortierellomycotina #9
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'Outgroup1'}{'Mortierellomycotina'}{'Mortierellomycetes'}{''}{'Mortierellales'}{''}{'Mortierellaceae'}{'Actinomortierella_wolfii_(Awol)'}{''} = 22;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'Outgroup1'}{'Mortierellomycotina'}{'Mortierellomycetes'}{''}{'Mortierellales'}{''}{'Mortierellaceae'}{'Actinomortierella_ambigua_(Aamb)'}{''} = 23;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'Outgroup1'}{'Mortierellomycotina'}{'Mortierellomycetes'}{''}{'Mortierellales'}{''}{'Mortierellaceae'}{'Linnemannia_elongata_(Lelo)'}{''} = 24;
$spp {'Eukaryota'}{'Opisthokonta'}{'Fungi'}{'Mucoromycota'}{'Outgroup1'}{'Mortierellomycotina'}{'Mortierellomycetes'}{''}{'Mortierellales'}{''}{'Mortierellaceae'}{'Podila_verticillata_(Pver)'}{''} = 25;



##     DOMAIN  ## SUBDOMAIN ## SUPERGROUP ## SUPERKINGDOM # Node1 # Node2 # KINGDOM ## SUBKINGDOM1 ## SUBKINGDOM2 # SUBKINGDOM3# SUBKINGDOM4# SUPERCLADE  ##  PHYLUM ##  SUBPHYLUM ##  SPECIES ##

print_hash_colors (%spp);
return %spp;
}

#PRINT_HASH_COLORS: subroutine to traverse the hash of hashes and print the contents, after http://www.perlmonks.org/?node_id=116162
sub print_hash_colors {
    my (%spp) = @_; 
#Define variables (taxonomic ranks)for the hash of hashes

my $superkingdom;
my $clade;
my $kingdom;
my $phylum;
my $node1;
my $subphylum;
my $class;
my $node2;
my $order;
my $node3;
my $family;
my $species;


#Print all the keys of the hash of hashes with a nested FOR loop
## SUPERKINGDOM ## CLADE ## KINGDOM ## PHYLUM ## Node1 ## SUBPHYLUM ## CLASS ## Node2 ## ORDER ## Node3 ##FAMILY ## SPECIES ##
print "\nTree:";
for $superkingdom (keys %spp) {
      print "\n$superkingdom\n";
      for $clade (keys %{ $spp{$superkingdom} }) {
        unless ($clade eq '') {print color ("blue"),"$clade\n", color ("reset");}
        for $kingdom (keys %{ $spp{$superkingdom}{$clade} }) {
            unless ($kingdom eq '') { print color ("yellow"),"  $kingdom\n", color ("reset");}
            for $phylum (keys %{ $spp{$superkingdom}{$clade}{$kingdom} }) {
                unless ($phylum eq '') { print color ("green"),"    $phylum\n", color ("reset");}
			    for $node1 (keys %{ $spp{$superkingdom}{$clade}{$kingdom}{$phylum} }) {
					unless ($node1 eq '') { print color ("magenta"),"        $node1\n", color ("reset");}
					for $subphylum (keys %{ $spp{$superkingdom}{$clade}{$kingdom}{$phylum}{$node1} }){
						unless ($subphylum eq '') { print color ("bright_blue"),"		$subphylum\n", color ("reset");}
						for $class (keys %{ $spp{$superkingdom}{$clade}{$kingdom}{$phylum}{$node1}{$subphylum} }) {
							unless ($class eq '') { print color ("red"),"            $class\n", color ("reset");}
							for $node2 (keys %{ $spp{$superkingdom}{$clade}{$kingdom}{$phylum}{$node1}{$subphylum}{$class} }) {
								unless ($node2 eq '') { print color ("blue"),"              $node2\n", color ("reset");}
								for $order (keys %{ $spp{$superkingdom}{$clade}{$kingdom}{$phylum}{$node1}{$subphylum}{$class}{$node2} }) {
									unless ($order eq '') { print color ("bright_yellow"),"                $order\n", color ("reset");}
									for $node3 (keys %{ $spp{$superkingdom}{$clade}{$kingdom}{$phylum}{$node1}{$subphylum}{$class}{$node2}{$order} }) {
										unless ($node3 eq '') { print color ("bright_green"),"                  $node3\n", color ("reset");}
										for $family (keys %{ $spp{$superkingdom}{$clade}{$kingdom}{$phylum}{$node1}{$subphylum}{$class}{$node2}{$order}{$node3} }) {
											unless ($family eq '') { print color ("bright_magenta"),"                    $family\n", color ("reset");}
											for $species (keys %{ $spp{$superkingdom}{$clade}{$kingdom}{$phylum}{$node1}{$subphylum}{$class}{$node2}{$order}{$node3}{$family} }){                              
												print color ("bright_white"), "                            $species => $spp{$superkingdom}{$clade}{$kingdom}{$phylum}{$node1}{$subphylum}{$class}{$node2}{$order}{$node3}{$family}{$species}{''}\n", color ("reset");
											}
										}
									}
								}
							}
                        }
                    }
                }
            }
        }
    }
}
return;
}
