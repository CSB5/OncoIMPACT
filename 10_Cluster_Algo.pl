#!/usr/bin/perl
use warnings;

my ( $data_dir, $network_type, $depth_th, $hub_th, $nb_joker,
	$fold_change_threshold, $out_dir, $script_dir )
  = @ARGV;

#$dir_TCGA_sample -> directory where the data are organize
#$network_gene_list -> ~/Software/netbox/data/Homo_sapiens.gene_info
#$network -> ~/Software/netbox/db/netbox.script
#$depth_th -> Threshold on maximum depth for defining module expansion
#$hub_th -> No Gene that is connected to more that $hub_th other gene will be included.

require "$script_dir/Construct_network.pl";

@dys_status_corress = ( "DOWN", "UP" );

`mkdir $out_dir` unless ( -d $out_dir );

#   1. Indexing
#not used due to synonyms problems/confusions

my %gene_to_index;
my @index_to_gene;
my @connections;

if ( $network_type eq "NETBOX" ) {
	construct_netbox_network(
		$data_dir,     \@index_to_gene, \%gene_to_index,
		\@connections, $script_dir
	);
}
if ( $network_type eq "DRIVER_NET" ) {
	construct_driver_net_network( \@index_to_gene, \%gene_to_index,
		\@connections, $script_dir );
}

# 3. Use the network to compute nicely $depth_th $hub_th
# Neighbour size distribution to find the the hub threshold
my $compute_threshold = 0;
if ($compute_threshold) {
	my @neigh_size = ();
	my $j          = 0;
	my $mean       = 0;
	for ( my $i = 0 ; $i < @index_to_gene ; $i++ ) {
		if ( @{ $connections[$i] } ) {
			$neigh_size[$j] = @{ $connections[$i] };
			$mean += @{ $connections[$i] };
			$j++;
		}

		#print "".(get_name($i))."\t".$neigh_size[$i]."\n";
	}
	@neigh_size_sorted = sort { $b <=> $a } @neigh_size;
	$mean = $mean / @index_to_gene;
	print STDERR "--------------- mean connexion " 
	  . $mean . "\t"
	  . $neigh_size_sorted[ int( @neigh_size_sorted / 2 ) ] . "\t"
	  . $neigh_size_sorted[ int( @neigh_size_sorted / 4 ) ] . "\t"
	  . $neigh_size_sorted[0] . "\n";

	my @tab = ();
	$infinity = 10000000;
	my $sum                        = 0;
	my $sum_avg                    = 0;
	my $gene_with_connexion        = 0;
	my $gene_with_connexion_sample = 0;
	my @explored                   = ();
	for ( my $i = 0 ; $i < @index_to_gene ; $i++ ) {

		#print STDERR "** $i\n" if($i % 50 == 0);
		if ( @{ $connections[$i] } ) {
			$gene_with_connexion_sample = 0;
			$sum                        = 0;
			for ( my $j = 0 ; $j < @index_to_gene ; $j++ ) {
				$explored[$j] = -1;
			}
			construct_explained_gene_set( $i, \@tab, \@tab, \@explored,
				$infinity, $hub_th );

			#Compute the awg path length for the gene
			for ( my $j = 0 ; $j < @index_to_gene ; $j++ ) {
				if ( $explored[$j] != -1 && $explored[$j] != 0 ) {
					$sum += $explored[$j];
					$gene_with_connexion_sample++;
				}
			}
			if ( $gene_with_connexion_sample != 0 ) {

				#if($gene_with_connexion_sample >= 6000){
				$gene_with_connexion++;

#print STDERR $i." ---- ".$sum."\t".$gene_with_connexion_sample."\t".(sprintf("%.3f", $sum / $gene_with_connexion_sample))."\n";
				$sum_avg += $sum / $gene_with_connexion_sample;
			}
		}
	}

	print STDERR "------------- " 
	  . $sum_avg . "\t"
	  . $gene_with_connexion . "\t"
	  . $sum_avg / $gene_with_connexion . "\n";
	exit;
}

#TO DO

# 4. Construct the set of explained genes for each samples
my $nb_sample = 0;

#Expression stats
my @explained_gene_frequency = ();
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	my @tab_status = ();
	for ( $j = 0 ; $j < 2 ; $j++ ) { my %hash = (); $tab_status[$j] = \%hash; }
	$explained_gene_frequency[$i] = \@tab_status;
}
my @dysregulated_gene_frequency = ();
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	my @tab_status = ();
	for ( $j = 0 ; $j < 2 ; $j++ ) { my %hash = (); $tab_status[$j] = \%hash; }
	$dysregulated_gene_frequency[$i] = \@tab_status;
}
my @explained_gene_mutatation_impact = ();
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	my @tab_status = ();
	for ( $j = 0 ; $j < 2 ; $j++ ) {
		my @tab;
		for ( $k = 0 ; $k < 2 ; $k++ ) {
			my @p = ();
			$tab[$k] = \@p;
		}
		$tab_status[$j] = \@tab;
	}

	#for($j = 0; $j < 2; $j++){my %hash = ();$tab_status[$j] = \%hash;}
	$explained_gene_mutatation_impact[$i] = \@tab_status;

#print STDERR "****************** ".$explained_gene_mutatation_impact[$i]->[1]->[0]."\n";<STDIN>;
}

#mutation stats
my @mutated_gene_with_explained_set_frequency = ();
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	my %hash = ();
	$mutated_gene_with_explained_set_frequency[$i] = \%hash;
}
my @mutated_gene_frequency = ();
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	my %hash = ();
	$mutated_gene_frequency[$i] = \%hash;
}

my %samples_explained_gene_set;
my @explored;    #array use during the breath first search
my %all_sample_gene_dysregulated
  ; #to keep trak of the dysregulation status of each gene for each samples (use at step 3 and 4)
my %sample_gene_mutated;
my @sample_gene_mutated_list = ();

opendir( DIR, $data_dir );
@the_DATA_DIR = readdir(DIR);
close(DIR);

#Output files
open( $OUT,        ">$out_dir/SAMPLE_stats.dat" );
open( $SAMPLE_EXP, ">$out_dir/SAMPLE_EXP.dat" );
open( $MODULE,     ">$out_dir/MODULE.dat" );

#open($PAIR, ">$out_dir/PAIR.dat");

open( $DIST_MAX, ">$out_dir/DIST_EXP.dat" );

#$OUT = STDERR;
my $nb_sample_to_test = -10;

my %BUG_not_present_in_network = ();

foreach my $dir_sample (@the_DATA_DIR) {
    $mutation_file_name      = "$data_dir/$dir_sample/Genelist_Status.txt";
    if ( -e $mutation_file_name ) {
	

	@sample_gene_mutated_list = ();

	#initilazed the mutational/expression gene status
	%sample_gene_mutated = ();
	my @sample_gene_dysregulated;
	my @sample_gene_explained_mutation_impact;
	my @sample_gene_explained_mutation_number;

	for ( my $i = 0 ; $i < @index_to_gene ; $i++ ) {
	    my @status_tab = ( 0, 0 );    #DOWN, UP
	    $sample_gene_dysregulated[$i] = \@status_tab;

	    my @status_tab1 = ();
	    my @status_tab2 = ();
	    for ( my $j = 0 ; $j < 2 ; $j++ ) {
		my @tab1 = ( 0, 0 )
		    ; #to store the 2 mutation impact value: sum_neigh sum_neigh_dys_value
		$status_tab1[$j] = \@tab1;

		my @tab2 = ( 0, 0 )
		    ; #to store the 2 mutation impact value: sum_neigh sum_neigh_dys_value
		$status_tab2[$j] = \@tab2;
	    }
	    $sample_gene_explained_mutation_impact[$i] = \@status_tab1;
	    $sample_gene_explained_mutation_number[$i] = \@status_tab2;

	    $sample_gene_mutated_list[$i] = 0;
	}
	my %all_explained_gene_set = ();

	$nb_sample++;
	$nb_mutated_gene               = 0;
	$nb_mutated_with_explained_set = 0;
	$nb_explained_gene             = 0;

	open( FILE, "$mutation_file_name" );
	#print STDERR " *** read file $mutation_file_name\n";    #<STDIN>;
	#read the file to obtain the dysregulated genes
	while (<FILE>) {
	    chop($_);
	    @line = split( /\t/, $_ );
	    my @parts     = split( /_/, $line[0] );
	    my $gene_name = $parts[0];
	    my $status    = $parts[1];
	    if ( exists $gene_to_index{$gene_name} )
	    { #filter out all the gene_name that do not belong to the input network
		my $gene_ID = get_ID( $gene_name, \%gene_to_index );

#if(get_name($gene_ID) ne $gene_name){
#print "********* PROBLEM with gene_ID $gene_name != ".get_name($gene_ID)."\n";<STDIN>;
#}

		if ( $status eq "MUT" || $status eq "AMPL" || $status eq "DEL" )
		{

#if($status eq "AMPL" || $status eq "DEL"){
#print STDERR "*********** $dir_sample ".(get_name($gene_ID, \@index_to_gene))." ".$status."\n" ;
#}

#if($gene_name eq "MIA3"){
#print STDERR "***********MUT $dir_sample ".(get_ID($gene_name))." ".$gene_name."\n";
#<STDIN>;
#}
		    $sample_gene_mutated{$gene_ID}                   = 1;
		    $sample_gene_mutated_list[$gene_ID]              = 1;
		    $mutated_gene_frequency[$gene_ID]->{$dir_sample} = 1;
		}
		else {
		    $fold_change = $line[1];
		    if ( ( $status eq "UP" || $status eq "DOWN" )
			 && abs($fold_change) >= $fold_change_threshold )
		    {
			$status_ID = 0;
			$status_ID = 1 if ( $status eq "UP" );
			$sample_gene_dysregulated[$gene_ID]->[$status_ID] =
			    $fold_change;
			$dysregulated_gene_frequency[$gene_ID]->[$status_ID]
			    ->{$dir_sample} = $gene_name;

#print STDERR "|".$_."|\t".$sample_gene_dysregulated[$gene_ID]->[$status_ID]."\t".$status_ID."\n";<STDIN>;
		    }
		}
	    }
	    else {
		$BUG_not_present_in_network{$gene_name} = 1;
	    }
	}
	close(FILE);

	if ( ( keys %sample_gene_mutated ) == 0 ) {
	    print STDERR "SAMPLE WITHOUT MUTATED GENES WEIRD !!!\n";   #<STDIN>;
	}

	foreach my $gene ( keys %sample_gene_mutated ) {
	    $module_str = "";
	    $nb_mutated_gene++;

	    #print "****MUT TTT $gene ".(get_name($gene))."\n";
	    #print STDERR "- construct explainend gene set for ".(get_name($gene))."\n";
	    for ( $i = 0 ; $i < @index_to_gene ; $i++ ) { $explored[$i] = -1; }

	    #my @false = ();
	    $explained_gene_set =
		construct_explained_gene_set_dijkstra( $gene,
						       \@sample_gene_dysregulated, \@sample_gene_mutated_list,
						       \@explored, $depth_th, $hub_th );

#$explained_gene_set = construct_explained_gene_set($gene, \@sample_gene_dysregulated, \@sample_gene_mutated_list, \@explored, $depth_th, $hub_th);
	    $mut_impact =
		compute_impact( $explained_gene_set, \@sample_gene_dysregulated );

	    if ( @{$explained_gene_set} > 1 ) {
		$flag_coint_dys_gene = 0;

		#compute the frequency
		foreach $eg ( @{$explained_gene_set} ) {
		    for ( my $dys_status = 0 ; $dys_status < 2 ; $dys_status++ )
		    {
			if ( $sample_gene_dysregulated[$eg]->[$dys_status] )
			{    #to remove the mutated genes and the jocker!

			    if ( !$flag_coint_dys_gene ) {
				$module_str .=
				    $dir_sample . "\t"
				    . ( get_name( $gene, \@index_to_gene ) )
				    . "\t";
				$nb_mutated_with_explained_set++;

				#print "*******MUT EXP $gene ".(get_name($gene))."\n";
				$mutated_gene_with_explained_set_frequency[$gene
				    ]->{$dir_sample} = 1;
				$flag_coint_dys_gene = 1;
			    }

			    #the module
			    $module_str .=
				( get_name( $eg, \@index_to_gene ) ) . "_"
				. $dys_status_corress[$dys_status] . ";";

#print the pair mutated gene/dys gene
#print $PAIR $dir_sample."\t".(get_name($gene, \@index_to_gene))."\t".(get_name($eg, \@index_to_gene))."_".$dys_status_corress[$dys_status]."\n";

			    if ( !exists $all_explained_gene_set{$eg} ) {
				$nb_explained_gene++;
				$all_explained_gene_set{$eg} = 1;
			    }

			    #the gene is explain at least by gene
			    #Need to be changed by an array !!!!
			    $explained_gene_frequency[$eg]->[$dys_status]
				->{$dir_sample} = $gene;

			    #to update the explained gene mutation impact
			    for ( my $k = 0 ; $k < 2 ; $k++ ) {
				if ( $sample_gene_explained_mutation_impact[$eg]
				     ->[$dys_status]->[$k] == 0 )
				{

#print STDERR $eg."\t".$dys_status."\t".$mut_impact->[$k]."\t".$sample_gene_explained_mutation_impact[$eg]->[$dys_status]->[$k]."\n";<STDIN>;
				}
				$sample_gene_explained_mutation_impact[$eg]
				    ->[$dys_status]->[$k] += $mut_impact->[$k];
				$sample_gene_explained_mutation_number[$eg]
				    ->[$dys_status]->[$k]++;
			    }
			}
		    }
		}
	    }

	    #Output the sample module
	    print $MODULE $module_str . "\n" if ( $module_str ne "" );

	}
	$samples_explained_gene_set{$dir_sample} = \%all_explained_gene_set;

	#To compute the avg of the mutation impact
	for ( my $i = 0 ; $i < @index_to_gene ; $i++ ) {
	    for ( my $j = 0 ; $j < 2 ; $j++ ) {
		for ( my $k = 0 ; $k < 2 ; $k++ ) {

#print STDERR "|".$explained_gene_mutatation_impact[$i]->[$i]."|\t|".$sample_gene_explained_mutation_impact[$eg]->[$i]."|\n";<STDIN>;
		    if ( $sample_gene_explained_mutation_impact[$i]->[$j]->[$k]
			 != 0 )
		    {

#$explained_gene_mutatation_impact[$i]->[$j]->[$k] += $sample_gene_explained_mutation_impact[$i]->[$j]->[$k];
			push(
			    @{
				$explained_gene_mutatation_impact[$i]->[$j]
				    ->[$k]
			    },
			    (
			     $sample_gene_explained_mutation_impact[$i]->[$j]
			     ->[$k] /
			     $sample_gene_explained_mutation_number[$i]
			     ->[$j]->[$k]
			    )
			    );
		    }
		}
	    }
	}

	#Output the sample avg module statistics
	print $OUT $dir_sample . "\t"
	    . $nb_mutated_gene . "\t"
	    . $nb_mutated_with_explained_set . "\t";

	#
	if ( $nb_mutated_gene != 0 ) {
	    print $OUT (
		sprintf( "%.3f",
			 $nb_mutated_with_explained_set / $nb_mutated_gene )
		) . "\t";
	}
	else {
	    print $OUT "-1" . "\t";
	}

	#
	if ( $nb_mutated_with_explained_set != 0 ) {
	    print $OUT $nb_explained_gene . "\t"
		. (
		sprintf( "%.3f",
			 $nb_explained_gene / $nb_mutated_with_explained_set )
		) . "\n";
	}
	else {
	    print $OUT $nb_explained_gene . "\t" . "-1" . "\n";
	}

	#
	last if ( $nb_sample == $nb_sample_to_test );
    }
}

print $MODULE "LAST_SAMPLE\n";

close($OUT);

#close($PAIR);
close($MODULE);
close($SAMPLE_EXP);
close($DIST_MAX);

$BUG_nb_not_present_in_network = ( keys %BUG_nb_not_present_in_network );
print STDERR
"********* $BUG_nb_not_present_in_network genes not taken into account as they do not belong to the network\n";
print STDERR "\n";
print STDERR `cat SAMPLE_stats.dat`;
print STDERR "\n";
# Step 4 extract mutated gene that contain frequently explained genes
#my @is_frequent = ();
#my $freq_th = 0.05;
#foreach $sample (keys %samples_explained_gene_set){
#    print $sample."\n";
#    foreach $set (@{$samples_explained_gene_set{$sample}}){
#	my $mut_gene = $set->[0];
#	@is_frequent = ();
#	for($i = 1; $i < @{$set}; $i++){
#print "\t".$set->[$i]."\t".(get_name($set->[$i]));<STDIN>;
#	    for(my $dys_status = 0; $dys_status < 2; $dys_status++){
#		my $freq = keys %{$explained_gene_frequency[$set->[$i]]->[$dys_status]};
#		if(1 && $freq /$nb_sample >= $freq_th){
#		    push(@is_frequent, (get_name($set->[$i]))."_".$dys_status_corress[$dys_status]);
#		    push(@is_frequent, ($freq));
#		}
#	    }
#	}
#	print "".(get_name($mut_gene))." ".(keys %{$mutated_gene_frequency[$mut_gene]})." --------- ".(join("\t", @is_frequent))."\n" if(@is_frequent);
#   }
#}

#Show mutatation frequency
open( $OUT, ">$out_dir/mut_gene_freq.dat" );
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	my $res_mut_exp = keys %{ $mutated_gene_with_explained_set_frequency[$i] };
	my $res_mut     = keys %{ $mutated_gene_frequency[$i] };
	if ( $res_mut != 0 ) {
		print $OUT ""
		  . ( get_name( $i, \@index_to_gene ) ) . "\t"
		  . ( @{ $connections[$i] } ) . "\t"
		  . $res_mut . "\t"
		  . ( sprintf( "%.3f", $res_mut / $nb_sample ) ) . "\t"
		  . $res_mut_exp . "\t"
		  . ( sprintf( "%.3f", $res_mut_exp / $nb_sample ) ) . "\t"
		  . ( sprintf( "%.3f", $res_mut_exp / $res_mut ) ) . "\n";
	}
}
close($OUT);

#show the explained gene frequency
#add the expression frequency and look at the spearman correlation !!!!
#Output the mutated gene set that explained each gene
open( $OUT,       ">$out_dir/exp_gene_freq.dat" );
open( $EXPLAINED, ">$out_dir/EXPLAINED.dat" );
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	for ( my $dys_status = 0 ; $dys_status < 2 ; $dys_status++ ) {
		my $res_exp = keys %{ $explained_gene_frequency[$i]->[$dys_status] };
		my $res_dys = keys %{ $dysregulated_gene_frequency[$i]->[$dys_status] };

		#print STDERR "".(sprintf("%.3f",$res_dys/$nb_sample))."\n";<STDIN>;

		print $OUT ""
		  . ( get_name( $i, \@index_to_gene ) ) . "_"
		  . $dys_status_corress[$dys_status] . "\t"
		  . ( @{ $connections[$i] } ) . "\t"
		  . $res_dys . "\t"
		  . ( sprintf( "%.3f", $res_dys / $nb_sample ) ) . "\t"
		  . $res_exp . "\t"
		  . ( sprintf( "%.3f", $res_exp / $nb_sample ) ) . "\t";
		if ( $res_dys != 0 ) {
			print $OUT "" . ( sprintf( "%.3f", $res_exp / $res_dys ) );
			if ( $res_exp != 0 ) {
				for ( my $k = 0 ; $k < 2 ; $k++ ) {
					@impact_sorted =
					  sort { $a <=> $b }
					  @{ $explained_gene_mutatation_impact[$i]->[$dys_status]
						  ->[$k] };
					$size =
					  @{ $explained_gene_mutatation_impact[$i]->[$dys_status]
						  ->[$k] } + 0;
					$median = $impact_sorted[0];
					if ( $size > 1 ) {
						if ( $size % 2 == 0 ) {
							$median =
							  ( $impact_sorted[ ( $size / 2 ) - 1 ] +
								  $impact_sorted[ $size / 2 ] ) / 2;
						}
						else {
							$median = $impact_sorted[ ( $size + 1 ) / 2 ];
						}
					}

					#print STDERR "**************** ".$size."\t".$median."\n";
					print $OUT "\t" . ( sprintf( "%.4f", $median ) );

#print $OUT "\t".(sprintf("%.4f", $explained_gene_mutatation_impact[$i]->[$dys_status]->[$k]/$res_exp));
				}
			}
			else {
				print $OUT "\t0\t0";
			}
		}
		else { print $OUT "-\t0\t0"; }
		print $OUT "\n";

		if ( $res_dys != 0 ) {

			#to obtain the list of gene that explained a dysregulated gene
			my %mut_gene_set = ();
			foreach
			  my $k ( keys %{ $explained_gene_frequency[$i]->[$dys_status] } )
			{
				$g = $explained_gene_frequency[$i]->[$dys_status]->{$k};
				if ( !exists $mut_gene_set{$g} ) {
					$mut_gene_set{$g} = 0;
				}
				$mut_gene_set{$g}++;
			}

			print $EXPLAINED ""
			  . ( get_name( $i, \@index_to_gene ) ) . "_"
			  . $dys_status_corress[$dys_status] . "\t"
			  . $res_exp . "\t"
			  . ( sprintf( "%.3f", $res_exp / $nb_sample ) );
			foreach $g ( keys %mut_gene_set ) {
				print $EXPLAINED "\t"
				  . ( get_name( $g, \@index_to_gene ) ) . ":"
				  . $mut_gene_set{$g};
			}
			print $EXPLAINED "\n";
		}
	}
}
close($EXPLAINED);
close($OUT);

my ( $f, $t );
for ( $i = 0 ; $i < 2 ; $i++ ) {
	if ( $i == 0 ) {
		print STDERR
		  "\n\n********** OUTPUT EXPRESSION RESULTS !!!!\n";    #<STDIN>;
		$f = "$out_dir/exp_gene_freq.dat";
		$t = "dysregulated";
	}
	else {
		print STDERR "\n********** OUTPUT MUTATION RESULTS !!!!\n";    #<STDIN>;
		$f = "$out_dir/mut_gene_freq.dat";
		$t = "mutated";
	}
	print STDERR "\n***** top 10 $t genes\n";
	print STDERR `sort -k3,3 -n -r $f | head`;

	print STDERR "\n***** top 10 $t explained genes\n";
	print STDERR `sort -k5,5 -n -r $f | head`;

	#print STDERR "\n***** top 10 ratio explained/dysregulated genes\n";
	#print STDERR `sort -k4,4 -k5,5 -n -r $f | head `;
}

#print `head FREQ.dat`;

sub compute_impact {
	my ( $explained_gene_set, $sample_gene_dysregulated ) = @_;

	#my @mut_impact = (0, 0, 0, 0, 0, 0, 0, 0);
	my @mut_impact = ( 0, 0 );
	for ( my $i = 0 ; $i < @{$explained_gene_set} ; $i++ ) {
		$gene_explained = $explained_gene_set->[$i];

		$gene_explained_dys_value =
		  abs( $sample_gene_dysregulated->[$gene_explained]->[0] +
			  $sample_gene_dysregulated->[$gene_explained]->[1] );

		$mut_impact[0]++;    #nb gene explained

		$mut_impact[1] +=
		  $gene_explained_dys_value;    #some of dysregulation of gene explained

	}

	#<STDIN>;
	return \@mut_impact;
}

#perform a breath first search to search from path of dysregulated gene of length <= $depth_th that do not contain hub_gene
sub construct_explained_gene_set_dijkstra {
	my ( $mutated_gene, $sample_gene_dysregulated, $sample_gene_mutated,
		$dist_to_mutated_gene, $depth_th, $hub_th )
	  = @_;

	#print STDERR "\nconstruct_explained_gene_set_dijkstra $mutated_gene\n";
	$dist_to_mutated_gene->[$mutated_gene] = 0;
	my @queue = ($mutated_gene);
	my @res   = ();
	$max_gene_dist = 0;

	#my $nb_joker = 1;

	while ( @queue + 0 != 0 ) {

		#search for the smallest element
		#print STDERR "***************** |@queue|"."\n";
		$min_index = 0;
		for ( my $i = 0 ; $i < @queue ; $i++ ) {
			if ( $dist_to_mutated_gene->[ $queue[$i] ] <
				$dist_to_mutated_gene->[ $queue[$min_index] ] )
			{
				$min_index = $i;
			}
		}
		$gene      = $queue[$min_index];
		$gene_dist = $dist_to_mutated_gene->[$gene];

		#remove the gene for the queue
		splice( @queue, $min_index, 1 );

		if (   $gene_dist > 100
			&& $gene_dist < 100 * ( $nb_joker + 1 )
			&& $gene_dist % 100 >= 2 )
		{
			print STDERR "remove $gene "
			  . get_name( $gene, \@index_to_gene )
			  . " $gene_dist\n";    # =========== |@queue|\n";<STDIN>;
		}

		#print STDERR "**** EXPLORED gene: $gene $gene_dist\n";
		if (   int( $gene_dist / 100 ) <= $nb_joker
			&& int( $gene_dist / 100 ) + $gene_dist % 100 <= $depth_th )
		{

#($gene_dist - ($nb_joker * 100) <= 0 || $gene_dist - ($nb_joker * 100) <= $depth_th){
#this gene is explained by the mutation
			push( @res, $gene );
			$max_gene_dist = $gene_dist if ( $gene_dist > $max_gene_dist );

			#hub could be phenotypes genes !
			if ( $mutated_gene == $gene || @{ $connections[$gene] } <= $hub_th )
			{    #not a hub
				foreach $neigh ( @{ $connections[$gene] } ) {

#print STDERR "\t".$neigh."\t".@{$sample_gene_dysregulated}."\t".$sample_gene_dysregulated->[$neigh]->[0]."\t".$sample_gene_dysregulated->[$neigh]->[1]."\t".$sample_gene_mutated->[$neigh]."\t".@{$connections[$neigh]}."\n";

					#to have an additional control on the neighborhood
					if
					  ( #$sample_gene_mutated->[$neigh] == 0 &&#the dysregulation is more likely to be explained by a mutation
						    #@{$connections[$neigh]} <= $hub_th && #not a hub
						(
							   $nb_joker != 0
							|| $sample_gene_dysregulated->[$neigh]->[0] != 0
							|| $sample_gene_dysregulated->[$neigh]->[1] != 0
						)    #at least 1 jocker OR the neighbour is dysregulated
					  )
					{

						#print STDERR "\t".$neigh." --------------> ADD\n";

					  #update the queue if it the first time we explore the gene
						push( @queue, $neigh )
						  if ( $dist_to_mutated_gene->[$neigh] == -1 )
						  ;    #not already explored;

						#update the distance
						$edge_weight = 100000;
						if (   @{$sample_gene_dysregulated} == 0
							|| $sample_gene_dysregulated->[$neigh]->[0] != 0
							|| $sample_gene_dysregulated->[$neigh]->[1] != 0 )
						{      #dysregulated
							$edge_weight = 1;
						}
						else {
							$edge_weight = 100;
						}
						$dist_to_mutated_gene->[$neigh] =
						  $gene_dist + $edge_weight
						  if ( $dist_to_mutated_gene->[$neigh] == -1
							|| $dist_to_mutated_gene->[$neigh] >
							$gene_dist + $edge_weight );

					}

					#print STDERR "\n";
				}
			}
		}
	}
	print $DIST_MAX ""
	  . ( get_name( $mutated_gene, \@index_to_gene ) ) . "\t"
	  . $max_gene_dist . "\n"
	  if ( $max_gene_dist != 0 );

	#print STDERR "---- THE RES: ".join("\t", @res)."\n";<STDIN>;
	return \@res;
}

#perform a breath first search to search from path of dysregulated gene of length <= $depth_th that do not contain hub_gene
sub construct_explained_gene_set {
	my ( $mutated_gene, $sample_gene_dysregulated, $sample_gene_mutated,
		$dist_to_mutated_gene, $depth_th, $hub_th )
	  = @_;
	print STDERR "\nconstruct_explained_gene_set $mutated_gene\n";
	$dist_to_mutated_gene->[$mutated_gene] = 0;
	my @queue = ($mutated_gene);
	my @res   = ();

	while ( @queue + 0 != 0 ) {
		my $gene = shift(@queue);
		$gene_dist = $dist_to_mutated_gene->[$gene];

		#print STDERR "**** EXPLORED gene: $gene $gene_dist\n";
		if ( $gene_dist <= $depth_th ) {

			#this gene is explained by the mutation
			push( @res, $gene );
			foreach $neigh ( @{ $connections[$gene] } ) {

#print STDERR "\t".$neigh."\t".@{$sample_gene_dysregulated}."\t".$sample_gene_dysregulated->[$neigh]->[0]."\t".$sample_gene_dysregulated->[$neigh]->[1]."\t".$sample_gene_mutated->[$neigh]."\t".@{$connections[$neigh]}."\n";

				if (
					$dist_to_mutated_gene->[$neigh] == -1
					&&    #not already explored
					 #(@{$sample_gene_dysregulated} == 0 || $sample_gene_dysregulated->[$neigh]->[0] == 1 || $sample_gene_dysregulated->[$neigh]->[1] == 1) && #dysregulated
					(
						@{$sample_gene_mutated} == 0
						|| $sample_gene_mutated->[$neigh] == 0
					)
					&& #the dysregulation is more likely to be explained by a mutation
					@{ $connections[$neigh] } <= $hub_th    #not a hub
				  )
				{

					#print STDERR "\t".$neigh." --------------> ADD\n";
					if (   @{$sample_gene_dysregulated} == 0
						|| $sample_gene_dysregulated->[$neigh]->[0] != 0
						|| $sample_gene_dysregulated->[$neigh]->[1] != 0 )
					{                                       #dysregulated
						$dist_to_mutated_gene->[$neigh] = $gene_dist + 1;
					}
					else {

					}
					push( @queue, $neigh );
				}

				#print STDERR "\n";
			}
		}
	}

	#print STDERR "---- THE RES: ".join("\t", @res)."\n";<STDIN>;
	return \@res;
}

