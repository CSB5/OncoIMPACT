#!/usr/bin/perl
use warnings;

#use strict;
#the last parameter are pair of directory containing sample data and the associated module file sperated by comma
#in case of multiple pairs only the last one is outputed
my ( $filter_type, $network_type, $explained_freq_file, $freq_threshold, $path_threshold, $hub_threshold, $out_file, $script_dir )  = @ARGV;

#print STDERR " *** $script_dir\n";

require "$script_dir/Construct_network.pl";

my $NB_PARAM_SCRIPT = 8;
#print STDERR " *** out_file $out_file\n";

#read the recurrently explained file
my %rec_explained = ();
open( FILE, $explained_freq_file );
while (<FILE>) {
	chop $_;
	@line = split( /\t/, $_ );

	#print STDERR "|$line[9]|"."\n";
	if (
		   $line[5] >= $freq_threshold
		&& $line[9] ne "-"    #significant pvalue
		                      #($line[8] >= $impact_threshold)
	  )
	{
		$rec_explained{ $line[0] } = $line[5];

		#print STDERR " $line[0] ********\n";
	}
}
close(FILE);

if(keys(%rec_explained) + 0 == 0){
    print STDERR " *** Aborting! There is no phenotype inferred for this data set.\n";
    exit 2;
}

#Construct the graph
my %gene_to_index;
my @index_to_gene;
my @connections;

if ( $network_type eq "NETBOX" ) {
	construct_netbox_network(
		$sample_data_dir, \@index_to_gene, \%gene_to_index,
		\@connections,    $script_dir
	);    #, "DATA_UPDATE");
}
if ( $network_type eq "DRIVER_NET" ) {
	construct_driver_net_network( \@index_to_gene, \%gene_to_index,
		\@connections, $script_dir );
}

#Collect the fold change of each gene to compute the module impact
%all_sample_gene_dysregulated = ();
for ( my $i = $NB_PARAM_SCRIPT ; $i < @ARGV ; $i++ ) {
    @tmp = split( /\,/, $ARGV[$i] );
    $sample_data_dir = $tmp[0];

    opendir( DIR, $sample_data_dir );
    @the_DATA_DIR = readdir(DIR);
    close(DIR);
    foreach my $dir_sample (@the_DATA_DIR) {
	$mutation_file_name =
	    "$sample_data_dir/$dir_sample/Genelist_Status.txt";
	
	#print STDERR " *** $mutation_file_name\n";
	if ( -e $mutation_file_name){

	    open( FILE, "$mutation_file_name" );

	    #read the file to obtain the dysregulated and mutated genes
	    my %sample_gene_dysregulated = ();
	    while (<FILE>) {
		chop($_);
		@line = split( /\t/, $_ );
		my @parts = split( /_/, $line[0] );
		my $gene_name = $parts[0];
		if ( exists $gene_to_index{$gene_name} ) {
		    my $gene_ID = get_ID( $gene_name, \%gene_to_index );
		    my $status = $parts[1];
		    if ( ( $status eq "UP" || $status eq "DOWN" ) ) {
			$fold_change = $line[1];
			$sample_gene_dysregulated{$gene_ID} = abs($fold_change);

#if($dir_sample eq "SNU119_OVARY" && $gene_ID == 3569){
#print STDERR "|".$_."|\t".$gene_ID."\t".$status."\t".$fold_change."\n";#<STDIN>
#print STDERR " *** *** WEIRD TEST ".($sample_gene_dysregulated{3569})."\n";
#}

			#$sample_gene_dysregulated{$line[0]} = abs($fold_change);
			#$sample_gene_dysregulated{$line[0]} = abs($fold_change);
			#
		    }
		}
	    }
	    close(FILE);
	    $all_sample_gene_dysregulated{$dir_sample} =
		\%sample_gene_dysregulated;

#print STDERR " *** WEIRD TEST ".($all_sample_gene_dysregulated{"SNU119_OVARY"}->{3569})."\n";
	}
    }
}

my %alteration_freq;

#the bipartite graph that reprensent the link of:
#  - altered genes
#  - phenotype genes on each samples
%alteration_to_phenotype = ();
%phenotype_to_alteration = ();
%all_sample_gene_impact  = ();

for ( my $i = $NB_PARAM_SCRIPT ; $i < @ARGV ; $i++ ) {
	@tmp = split( /\,/, $ARGV[$i] );
	$module_file = $tmp[1];
	#print STDERR " *** $i READ $module_file\n";    #<STDIN>;
	open( FILE, $module_file );
	while (<FILE>) {
		chop $_;
		next if ( $_ eq "LAST_SAMPLE" );

		@line = split( /\t/, $_ );

		$sample_name = $line[0];
		if ( !exists $all_sample_gene_impact{$sample_name} ) {
			my %map = ();
			$all_sample_gene_impact{$sample_name} = \%map;
		}

		$altered_gene = $line[1];
		$all_sample_gene_impact{$sample_name}->{$altered_gene} = 0;

		@dys_gene = split( /\;/, $line[2] );

		$find_pheno = 0;
		foreach $g (@dys_gene) {
			if ( is_phenotype( $g, $sample_name ) ) {
				$phenotype_name = $g . "_" . $sample_name;

				#initilations
				if ( !exists $alteration_freq{$altered_gene} ) {
					$alteration_freq{$altered_gene} = 0;
					my %map = ();
					$alteration_to_phenotype{$altered_gene} = \%map;
				}

				if ( !exists $phenotype_to_alteration{$phenotype_name} ) {
					my %map = ();
					$phenotype_to_alteration{$phenotype_name} = \%map;
				}

				#update the data structure
				if ( !$find_pheno ) {
					$alteration_freq{$altered_gene}++;
					$find_pheno = 1;
				}
				$alteration_to_phenotype{$altered_gene}->{$phenotype_name} = 1;
				$phenotype_to_alteration{$phenotype_name}->{$altered_gene} = 1;
			}

			#For the impact computation
			my @parts = split( /_/, $g );
			$gene_ID = get_ID( $parts[0], \%gene_to_index );
			$all_sample_gene_impact{$sample_name}->{$altered_gene} +=
			  $all_sample_gene_dysregulated{$sample_name}->{$gene_ID};

#print STDERR $altered_gene."\t".$all_sample_gene_impact{$sample_name}->{$altered_gene}."\n";<STDIN>;
		}
	}
	close(FILE);
}

#exit(0);<STDIN>;

#print STDERR
#    "BRAF: ".($alteration_freq{"BRAF"})."\n".
#    "LRP1: ".($alteration_freq{"LRP1"})."\n".
#    "IL6ST: ".($alteration_freq{"IL6ST"})."\n";#<STDIN>;

#For the average impact
my $nb_sample           = keys(%all_sample_gene_impact);
my %all_avg_gene_impact = ();
my %driver_list         = ();
foreach $s ( keys %all_sample_gene_impact ) {
	my %map = ();
	$driver_list{$s} = ();
	foreach $g ( keys %{ $all_sample_gene_impact{$s} } ) {
		if ( !exists $all_avg_gene_impact{$g} ) {
			$all_avg_gene_impact{$g} = 0;
		}
		$all_avg_gene_impact{$g} += $all_sample_gene_impact{$s}->{$g};
	}
}
foreach $g ( keys %all_avg_gene_impact ) {
	$all_avg_gene_impact{$g} = $all_avg_gene_impact{$g} / $nb_sample;
}

#Weighted set cover
my $impact_type = "AVG_IMPACT";
if ( $impact_type eq "AVG_IMPACT" ) {

#compute_driver_list(\%alteration_to_phenotype, \%phenotype_to_alteration, \%driver_list, \%all_avg_gene_impact, \%all_sample_gene_impact, "DRIVER_SAMPLE");#"DRIVER_ALL");
	compute_driver_list( \%alteration_to_phenotype, \%phenotype_to_alteration,
		\%driver_list, \%alteration_freq, \%all_sample_gene_impact,
		$filter_type );
}
else {
	my %alteration_to_phenotype_sample;
	my %phenotype_to_alteration_sample;
	foreach my $s ( keys %all_sample_gene_impact ) {

		#print STDERR " *** *** Set cover on sample $s\n";

		%alteration_to_phenotype_sample = ();
		%phenotype_to_alteration_sample = ();
		foreach my $altered_gene ( keys %alteration_to_phenotype ) {
			foreach my $phenotype_name (
				keys %{ $alteration_to_phenotype{$altered_gene} } )
			{

				#gene is altered on the sample and is a sample phenotye
				if ( exists $all_sample_gene_impact{$s}->{$altered_gene}
					&& index( $phenotype_name, $s ) != -1 )
				{
					$alteration_to_phenotype_sample{$altered_gene}
					  ->{$phenotype_name} = 1;
					$phenotype_to_alteration_sample{$phenotype_name}
					  ->{$altered_gene} = 1;
				}
			}
		}
		compute_driver_list(
			\%alteration_to_phenotype_sample,
			\%phenotype_to_alteration_sample,
			\%driver_list,
			$all_sample_gene_impact{$s},
			\%all_sample_gene_impact,
			"DRIVER_SAMPLE"
		);
	}
}

#print STDERR " *** Nb driver gene inferred: "  . ( ( keys %driver_list ) + 0 ) . "\n";

#print STDERR "sample driver list\n";
#foreach $d (keys %{$driver_list{"TCGA-10-0938-01"}}){
#    print $d."\t".$alteration_freq{$d}."\n";
#}

#compute the hub genes
my %hub_gene = ();
for ( my $gene = 0 ; $gene < @index_to_gene ; $gene++ ) {
	$nb_connection = @{ $connections[$gene] };
	if ( $nb_connection > $hub_threshold ) {
		$name = get_name( $gene, \@index_to_gene );

		#print " --- $gene $name\n";<stdin>;
		$hub_gene{ $name . "_UP" }   = 1;
		$hub_gene{ $name . "_DOWN" } = 1;
	}
}

#open(OUT, ">/home/bertrandd/pathways/Fine_Arts/background_gene_list.dat");
#for($i = 0; $i < @index_to_gene; $i++){
#    print OUT $index_to_gene[$i]->[0]."\n";
#}
#close(OUT);

open( OUT,  ">$out_file" );
open( FILE, $module_file );
my $current_sample = "lapin";
my %exp_to_mod     = ();
my %module_list    = ();
my $mod_ID;

#for the graph exploration
my @sub_graph = ();
my @belong_to_one_path;
my $nb_node = @index_to_gene;

for ( $i = 0 ; $i < $nb_node ; $i++ ) {
	$sub_graph[$i]          = 0;
	$belong_to_one_path[$i] = 0;
}

my %all_dist = ();
my @mut_dist = ();
while (<FILE>) {

	chop $_;
	@line = split( /\t/, $_ );
	$sample = $line[0];

	#print STDERR " *** $sample\n";

	#new sample module
	#FILETRING OF THE MODULE OF THE NEW SAMPLE
	if ( $current_sample ne $sample ) {
		if ( ( keys %module_list ) != 0 ) {

#need to clean the module by deleting all the explained gene that do not have a path for a utated gene to an explained gene
			foreach $m ( keys %module_list ) {

#print STDERR "clean:\n".$current_sample.".".$m."\t".(join(";", @{$module_list{$m}->[0]}))."\t".(join(";", @{$module_list{$m}->[1]}))."\n";#<STDIN>;
				%all_dist              = ();
				%sample_freq_explained = ();

				#initalize freq_explained and sub_graph
				foreach $g ( @{ $module_list{$m}->[1] } ) {
					@name = split( /\_/, $g );
					$sub_graph[ get_ID( $name[0], \%gene_to_index ) ] = 1;
					$sample_freq_explained{ get_ID( $name[0], \%gene_to_index )
					  } = 1
					  if ( is_phenotype( $g, $current_sample ) );
				}

				%sample_driver = ();
				foreach $g ( @{ $module_list{$m}->[0] } ) {
					@name = split( /\_/, $g );
					$sample_driver{ get_ID( $name[0], \%gene_to_index ) } = 1;
				}

#compute the minimum distance of each freq_explained genes to other gene in the module
				foreach $ID ( keys %sample_freq_explained ) {
					my @dist;
					for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
						$dist[$i] = -1;
					}
					compute_dist_to_gene( $ID, \@dist, $path_threshold,
						\%sample_driver );
					$all_dist{$ID} = \@dist;
				}

 #compute the minimum distance of each mutated genes to other gene in the module
				foreach $g ( @{ $module_list{$m}->[0] } ) {
					for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
						$mut_dist[$i] = -1;
					}
					$mut_gene_ID = get_ID( $g, \%gene_to_index );
					compute_dist_to_gene( $mut_gene_ID, \@mut_dist,
						$path_threshold, \%sample_freq_explained );
					$belong_to_one_path[$mut_gene_ID] =
					  1;    #the mutated gene belong to the path
					        #update belong to one path
					for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
						if (   $mut_dist[$i] != -1
							&& !$sample_freq_explained{$i}
							&& !$belong_to_one_path[$i] )
						{
							foreach my $phenotype ( keys %all_dist ) {
								$phenotype_name =
								  get_name( $phenotype, \@index_to_gene );
								$d = $all_dist{$phenotype}->[$i];

#if($phenotype_name eq "TOP2A"){print STDERR " --- dist to $phenotype_name ".(get_name($i, \@index_to_gene))."->".$d."\n";}
								if (   $phenotype != $i
									&& $d != -1
									&& $d + $mut_dist[$i] <= $path_threshold )
								{
									$belong_to_one_path[$i] = 1;
									last;
								}
							}
						}
					}
				}

#Make the difference between the recurrently explained genes and the conserved linker genes and cmpute the modlue impact
				$module_impact_unfiltered = 0;
				$module_impact            = 0;
				$nb_rec                   = 0;
				$nb_link                  = 0;

				$rec_str  = "";
				$link_str = "";
				foreach $g ( @{ $module_list{$m}->[1] } ) {
					###########################################################
					######## NEED INVESTIGATE THE NAMING ISSUES ###############
					###########################################################
					@tmp = split( /\t/, $g );
					my @parts     = split( /_/, $tmp[0] );
					my $gene_name = $parts[0];
					my $status    = $parts[1];
					my $c_gene_ID = get_ID( $gene_name, \%gene_to_index );
					#########################################################

					if ( !exists $all_sample_gene_dysregulated{$current_sample}
						->{$c_gene_ID} )
					{

						#print STDERR " --- $current_sample $g\n";<STDIN>;
					}
					$module_impact_unfiltered +=
					  $all_sample_gene_dysregulated{$current_sample}
					  ->{$c_gene_ID};

					if ( is_phenotype( $g, $current_sample ) ) {
						$rec_str .= $g . ";";
						$nb_rec++;
						$module_impact +=
						  $all_sample_gene_dysregulated{$current_sample}
						  ->{$c_gene_ID};
					}
					else {
						@name = split( /\_/, $g );
						if (
							$belong_to_one_path[ get_ID( $name[0],
								  \%gene_to_index ) ] )
						{

					 #print STDERR "------------------------ it is a link!!!\n";
							$link_str .= $g . ";";
							$nb_link++;
							$module_impact +=
							  $all_sample_gene_dysregulated{$current_sample}
							  ->{$c_gene_ID};
						}
					}
				}

				$link_str = "-;" if ( $link_str eq "" );

#print STDERR $current_sample.".".$m."\t".(join(";", @{$module_list{$m}->[0]})).";\t".$rec_str."\t".$link_str."\n";
				print OUT $current_sample . "." 
				  . $m . "\t"
				  . ( join( ";", @{ $module_list{$m}->[0] } ) ) . ";\t"
				  . $rec_str . "\t"
				  . $link_str . "\t"
				  . ( @{ $module_list{$m}->[0] } ) . "_"
				  . $nb_rec . "_"
				  . $nb_link . "\t"
				  . $module_impact_unfiltered . "\t"
				  . $module_impact . "\n";

#print STDERR $current_sample.".".$m."\t".(join(";", @{$module_list{$m}->[0]})).";\t".$rec_str."\t".$link_str."\n";

				#clean freq_explained and sub_graph

				foreach $g ( @{ $module_list{$m}->[1] } ) {
					@name = split( /\_/, $g );
					$ID = get_ID( $name[0], \%gene_to_index );
					$sub_graph[$ID]          = 0;
					$belong_to_one_path[$ID] = 0;
				}
			}

#print the module
#foreach $m (keys %module_list){
#print OUT $current_sample.".".$m."\t".(join(";", @{$module_list{$m}->[0]}))."\t".(join(";", @{$module_list{$m}->[1]}))."\n";
#}

		}
		$current_sample = $sample;
		%exp_to_mod     = ();
		%module_list    = ();
		$mod_ID         = 1;
	}
	next if ( $_ eq "LAST_SAMPLE" );
	###########################
	#next if($sample ne "JHOC5_OVARY");
	##########################

	#save the new module
	#The mutated gene should be a driver for this sample
	if ( exists $driver_list{$current_sample}->{ $line[1] } ) {
		my @exp_gene_list = split( /\;/, $line[2] );
		my %exp_map = ();

		#check if it least 1 explained gene is recurrently explained
		$flag_rec = 0;
		foreach $exp (@exp_gene_list) {
			if ( is_phenotype( $exp, $current_sample ) ) {
				$flag_rec = 1;

				#last;
			}
			$exp_map{$exp} = 1;
		}

		if ( $flag_rec == 0 ) {

	  #print STDERR "************ MODULE with no recurrently explained gene:\n";
	  #print STDERR "$_\n";
		}
		else {
			my @mut_gene_list = ( $line[1] );
			my @tab = ( \@mut_gene_list, \@exp_gene_list );
			$module_list{$mod_ID} = \@tab;

			#associte all the explained gene to the module
			foreach $exp (@exp_gene_list) {
				if ( !exists $exp_to_mod{$exp} ) {
					my %map = ( $mod_ID, 1 );
					$exp_to_mod{$exp} = \%map;
				}

				else {
					############# Could be better ????
					@all_mod      = ( keys %{ $exp_to_mod{$exp} } );
					$other_mod_ID = $all_mod[0];
					##################

					#2 different modules have an explained gene in common
					if ( $other_mod_ID != $mod_ID )
					{ #possible to assign a gene with the current ID during a previous merge
						  #merge if a gene is rec_explained (PHENOTYPE) or a HUB
						  #print STDERR " --- $exp\n";<STDIN>;
						 #if(exists $rec_explained{$exp} || exists $hub_gene{$exp} # merge using hub
						if (
							is_phenotype( $exp, $current_sample )

						#&& !exists $hub_gene{$exp} # do not allows hub to merge
						  )
						{

							#merge the 2 arrays
							my @tab_m = (
								@{ $module_list{$mod_ID}->[0] },
								@{ $module_list{$other_mod_ID}->[0] }
							);

							#to remove the genes that belong to the 2 modules
							my @tab_e = @{ $module_list{$mod_ID}->[1] };
							foreach $g ( @{ $module_list{$other_mod_ID}->[1] } )
							{
								if ( !exists $exp_map{$g} ) {
									push( @tab_e, $g );
								}
							}
							my @tab = ( \@tab_m, \@tab_e );
							$module_list{$mod_ID} = \@tab;
							delete $module_list{$other_mod_ID};

							#change the mod_ID of all the explained gene
							foreach $g ( keys %exp_to_mod ) {
								if ( exists $exp_to_mod{$g}->{$other_mod_ID} ) {
									delete $exp_to_mod{$g}->{$other_mod_ID};
									$exp_to_mod{$g}->{$mod_ID} = 1;
								}
							}
						}

#share explainend gene that is not recurrent to other sample
#the gene will belong to more than 1 module as a link => need to allow HUB genes to merge modules
#after the cleaning process the gene should belong to only 1 module as:
#  - it belong to multiple path of a phenotype => the phenotype gene belong to the same module
#  - it is removed from the module from which it do belong to a path to a phenotype
						else {
							$exp_to_mod{$exp}->{$mod_ID} = 1;
						}
					}
				}
			}

			#new module ID
			$mod_ID++;
		}
	}
}
close(FILE);
close(OUT);

sub is_phenotype {
	my ( $gene, $sample_name ) = @_;

#return (exists $rec_explained{$gene} && exists $rec_explained{$gene}->{$sample_name});
	return ( exists $rec_explained{$gene} );
}

sub compute_driver_list {
	my ( $alteration_to_phenotype, $phenotype_to_alteration, $driver_list,
		$impact_data, $impact_data_sample, $inference_type )
	  = @_;
	my @max_degree_altered_gene_list;
	my $nb_pheno_to_explain = ( keys %{$phenotype_to_alteration} );

	while ( $nb_pheno_to_explain != 0 ) {

#print STDERR "\n *** ALGO step ".((keys %driver_list)+0)." ".$nb_pheno_to_explain."\n";

		#$driver = infer_driver_degree();
		#$driver = infer_driver_freq();
		$driver = infer_driver_impact( $alteration_to_phenotype, $impact_data );
		if ( $inference_type eq "DRIVER_ALL" ) {
			foreach my $s ( keys %{$impact_data_sample} ) {
				$driver_list->{$s}->{$driver} = 1
				  if ( exists $impact_data_sample->{$s}->{$driver} );
			}
		}

		#update the data structure
		foreach $p ( keys %{ $alteration_to_phenotype{$driver} } ) {
			$p_sample = "";
			if ( $p =~ /(.+)\_(UP|DOWN)\_(.+)/ ) {
				$p_sample = $3;
			}

			if ( $inference_type eq "DRIVER_SAMPLE" ) {

				#print STDERR " *** $p $driver $p_sample\n";<STDIN>;
				$driver_list->{$p_sample}->{$driver} = 1;
			}

			foreach $a ( keys %{ $phenotype_to_alteration{$p} } ) {

				#print STDERR "\t remove edges $a $p\n";
				delete $phenotype_to_alteration->{$p}->{$a};
				delete $alteration_to_phenotype->{$a}->{$p};
			}
			delete $phenotype_to_alteration->{$p};
		}

		$nb_pheno_to_explain = ( keys %{$phenotype_to_alteration} );

		#print STDERR " *** $nb_pheno_to_explain $driver END\n\n";#<STDIN>;
	}
}

sub infer_driver_impact {
	my ( $graph, $impact_data ) = @_;
	return infer_driver_impact_degree( $graph, $impact_data );

	#return infer_driver_impact_freq($graph, $impact_data);
}

sub infer_driver_impact_freq {
	my ( $graph, $impact_data ) = @_;

	#take the node with the highest alteration frequency
	my $max_impact           = 0;
	my @max_impact_gene_list = ();
	foreach my $g ( keys %{$impact_data} ) {
		if ( get_degree( $g, $graph ) != 0 ) {
			$impact = $impact_data->{$g};
			if ( $impact > $max_impact ) {
				@max_impact_gene_list = ();
				$max_impact           = $impact;
			}

			push( @max_impact_gene_list, $g ) if ( $max_impact == $impact );
		}

	}

   #print STDERR " *** infer_driver_impact $max_impact @max_impact_gene_list\n";

	#break tie with node with highest degree
	$driver     = "lapin";
	$max_degree = 0;
	foreach my $g (@max_impact_gene_list) {
		$degree = get_degree( $g, $graph );
		if ( $degree > $max_degree ) {
			$driver     = $g;
			$max_degree = $degree;
		}
	}

#delete $alteration_freq{$driver};
#print STDERR "*** Driver found: $driver $max_freq ( ".(@max_freq_altered_gene_list+0)." ) ".$max_degree."\n";

	return $driver;

}

sub infer_driver_impact_degree {
	my ( $graph, $impact_data ) = @_;

	#take the node with the highest alteration frequency
	my $max_degree           = 0;
	my @max_degree_gene_list = ();
	foreach my $g ( keys %{$impact_data} ) {
		$degree = get_degree( $g, $graph );
		if ( $degree != 0 ) {
			if ( $degree > $max_degree ) {
				@max_degree_gene_list = ();
				$max_degree           = $degree;
			}

			push( @max_degree_gene_list, $g ) if ( $max_degree == $degree );
		}

	}

   #print STDERR " *** infer_driver_impact $max_impact @max_impact_gene_list\n";

	#break tie with node with highest degree
	$driver     = "lapin";
	$max_impact = 0;
	foreach my $g (@max_degree_gene_list) {
		$impact = $impact_data->{$g};
		if ( $impact > $max_impact ) {
			$driver     = $g;
			$max_impact = $impact;
		}
	}

#delete $alteration_freq{$driver};
#print STDERR "*** Driver found: $driver $max_freq ( ".(@max_freq_altered_gene_list+0)." ) ".$max_degree."\n";

	return $driver;

}

sub get_max_degree_node_list {
	my ( $graph, $gene_list ) = @_;

	#print STDERR "get_max_degree_node_list\n";

	my $max_degree = 0;

	foreach $a ( keys %{$graph} ) {
		$degree = get_degree( $a, $graph );

		#print STDERR $a."\t".$degree."\n";
		if ( $degree > $max_degree ) {
			$max_degree = $degree;
			@{$gene_list} = ();
		}
		push( @{$gene_list}, $a ) if ( $degree == $max_degree );
	}

	return $max_degree;

}

sub get_degree {
	my ( $node, $graph ) = @_;
	return ( 0 + ( keys %{ $graph->{$node} } ) );
}

#perform a breath first search to search from path of dysregulated gene of length <= $depth_th that do not contain hub_gene
sub compute_dist_to_gene {
	my ( $studied_gene, $dist_to_mutated_gene, $depth_th, $filtered_node ) = @_;

	#print STDERR "\nconstruct_explained_gene_set $mutated_gene\n";
	$dist_to_mutated_gene->[$studied_gene] = 0;
	my @queue = ($studied_gene);
	my @res   = ();

	while ( @queue + 0 != 0 ) {
		my $gene = shift(@queue);
		$gene_dist = $dist_to_mutated_gene->[$gene];

		#print STDERR "**** EXPLORED gene: $gene $gene_dist\n";
		if ( $gene_dist <= $depth_th ) {
			foreach $neigh ( @{ $connections[$gene] } ) {

#print STDERR "\t".$neigh."\t".@{$sample_gene_dysregulated}."\t".$sample_gene_dysregulated->[$neigh]->[0]."\t".$sample_gene_dysregulated->[$neigh]->[1]."\t".$sample_gene_mutated->[$neigh]."\t".@{$connections[$neigh]}."\n";
				$neigh_connection = @{ $connections[$neigh] };
				if (
					$dist_to_mutated_gene->[$neigh] == -1
					&&    #not already explored
					$sub_graph[$neigh] == 1
					&& $neigh_connection <= $hub_threshold
					&&    # to avoid going threw a path that contain a hub
					!exists $filtered_node->{$neigh
					} #to avoid going threw a path after a phenotype gene if a studied_gene is a mutated gene (and reversly)
				  )
				{
					$dist_to_mutated_gene->[$neigh] = $gene_dist + 1;
					push( @queue, $neigh );
				}

				#print STDERR "\n";
			}
		}
	}
}

sub search_all_path {
	my ( $gene, $current_path, $explored, $freq_explained, $sub_graph,
		$belong_to_one_path )
	  = @_;

#print STDERR "*********** search_all_path: ".$gene." ".(get_name($gene))."\n";<STDIN>;

	if ( @{$current_path} <= $path_threshold ) {

		$explored->[$gene] = 1;
		push( @{$current_path}, $gene );

		#Found a frequently explained gene
		if ( $freq_explained->[$gene] ) {
			foreach $g ( @{$current_path} ) {
				$belong_to_one_path->{$g} = 1;
			}
		}

		foreach $neigh ( @{ $connections[$gene] } ) {

	#print STDERR "\t$neigh ".$explored->[$neigh]." ".$sub_graph->[$neigh]."\n";
	#Still explored
			if (   $explored->[$neigh] == 0
				&& $sub_graph->[$neigh] == 1 )
			{

#print STDERR "\t$neigh ".$explored->[$neigh]." ".$sub_graph->[$neigh]."\n";<STDIN>;
				search_all_path(
					$neigh,          $current_path, $explored,
					$freq_explained, $sub_graph,    $belong_to_one_path
				);
			}
		}

		$explored->[$gene] = 0;
		pop( @{$current_path} );
	}

}

