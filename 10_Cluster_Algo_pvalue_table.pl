#!/usr/bin/perl
use warnings;

#If FAST_CALL_FLAG is true: the method stop if at least one random sample have a higher frequency than the real frequency: $THE_GENE_FREQ_EXPLAINED
my (
	$data_dir,         $network_type,
	$depth_th,         $hub_th,
	$nb_joker,         $log2_fold_change_threshold,
	$nb_random_sample, $THE_GENE_ID,
	$THE_GENE_STATUS,  $THE_GENE_FREQ_EXPLAINED,
	$out_file, $out_file_table,   $FAST_CALL_FLAG,
	$script_dir, $seed
) = @ARGV;

#$dir_TCGA_sample -> directory where the data are organize
#$network_gene_list -> ~/Software/netbox/data/Homo_sapiens.gene_info
#$network -> ~/Software/netbox/db/netbox.script
#$depth_th -> Threshold on maximum depth for defining module expansion
#$hub_th -> No Gene that is connected to more that $hub_th other gene will be included.

require "$script_dir/Construct_network.pl";

if($seed != -1){
    srand($seed);
}

$FLAG_PRINT_TABLE = 0;
$FLAG_PRINT_TABLE = 1 if ( $FAST_CALL_FLAG == 0 );

@dys_status_corress = ( "DOWN", "UP" );

#   1. Indexing
#not used due to synonyms problems/confusions

my %gene_to_index;
my @index_to_gene;
my @connections;

print STDERR " *** CONSTRUCT NETWORK\n";
if ( $network_type eq "NETBOX" ) {

	#update not done to save computational time
	construct_netbox_network(
		$data_dir,     \@index_to_gene, \%gene_to_index,
		\@connections, $script_dir, "NO_DATA_UPDATE"
	);
}
if ( $network_type eq "DRIVER_NET" ) {
	construct_driver_net_network( \@index_to_gene, \%gene_to_index,
		\@connections, $script_dir );
}

# 4. Construct the set of explained genes for each samples
my $nb_sample = 0;

#Expression stats
#mutation stats
my @mutated_gene_frequency = ();
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	my %hash = ();
	$mutated_gene_frequency[$i] = \%hash;
}

my @explored;    #array use during the breath first search

my @sample_order                 = ();
my %all_sample_gene_dysregulated = (); #to keep track of the dysregulation status of each gene for each samples (use at step 3 and 4)
my %all_sample_gene_mutated      = ();
my %all_sample_gene_mutated_list = ();

opendir( DIR, $data_dir );
@the_DATA_DIR = readdir(DIR);
close(DIR);

my $nb_sample_to_test = -10;

my %BUG_not_present_in_network = ();

print STDERR " *** READ DIR\n";
foreach my $dir_sample (@the_DATA_DIR) {
    $mutation_file_name      = "$data_dir/$dir_sample/Genelist_Status.txt";
    if ( -e $mutation_file_name ) {

	#initilazed the mutational/expression gene status
	my @sample_gene_mutated_list = ();
	my %sample_gene_mutated      = ();
	my @sample_gene_dysregulated;

	for ( my $i = 0 ; $i < @index_to_gene ; $i++ ) {
	    my @status_tab = ( 0, 0 );    #DOWN, UP
	    $sample_gene_dysregulated[$i] = \@status_tab;

	    my @status_tab1 = ();
	    my @status_tab2 = ();
	    for ( my $j = 0 ; $j < 2 ; $j++ ) {
		my @tab1 = ( 0, 0 ); #to store the 2 mutation impact value: sum_neigh sum_neigh_dys_value
		$status_tab1[$j] = \@tab1;

		my @tab2 = ( 0, 0 ); #to store the 2 mutation impact value: sum_neigh sum_neigh_dys_value
		$status_tab2[$j] = \@tab2;
	    }

	    $sample_gene_mutated_list[$i] = 0;
	}
	my %all_explained_gene_set = ();

	$nb_sample++;

	push( @sample_order, $dir_sample );

	open( FILE, "$mutation_file_name" );

	#print STDERR "********** read file $mutation_file_name\n";#<STDIN>;
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
			 && abs($fold_change) >= $log2_fold_change_threshold )
		    {
			$status_ID = 0;
			$status_ID = 1 if ( $status eq "UP" );
			$sample_gene_dysregulated[$gene_ID]->[$status_ID] =
			    $fold_change;

#$dysregulated_gene_frequency[$gene_ID]->[$status_ID]->{$dir_sample} = $gene_name;
#print STDERR "|".$_."|\t".$sample_gene_dysregulated[$gene_ID]->[$status_ID]."\n";#<STDIN>;
		    }
		}
	    }
	    else {
		$BUG_not_present_in_network{$gene_name} = 1;
	    }
	}
	close(FILE);

	#print STDERR " ------------------------------------------ $dir_sample\n";
	$all_sample_gene_mutated_list{$dir_sample} = \@sample_gene_mutated_list;
	$all_sample_gene_mutated{$dir_sample}      = \%sample_gene_mutated;
	$all_sample_gene_dysregulated{$dir_sample} = \@sample_gene_dysregulated;

	close(FILE);
	#########################################################################################

	if ( ( keys %sample_gene_mutated ) == 0 ) {
	    #print STDERR "SAMPLE WITHOUT MUTATED GENES WEIRD !!!\n";
	    #<STDIN>;
	}

    }
}

my @all_random_freq = ();
my $gene_ID         = $THE_GENE_ID;
my $dys_status      = $THE_GENE_STATUS;

#Construct the matrix
if ($FLAG_PRINT_TABLE) {
	open( OUT_TABLE, "| gzip -c > $out_file_table" );
	print OUT_TABLE "" . ( join( "\t", @sample_order ) ) . "\n";
}
my @explained_line = ();
##################################################
#To compute the pvalue of the explained genes
my @gene_pvalue = ();

#for(my $gene_ID = 0; $gene_ID < @index_to_gene; $gene_ID++){
my @tab = ( 0, 0 );
$gene_pvalue[$gene_ID] = \@tab;

#for(my $dys_status = 0; $dys_status < 2; $dys_status++){

#print "$gene_ID $dys_status".(get_name($gene_ID, \@index_to_gene))."_".$dys_status_corress[$dys_status]."\t".$res_exp."\n";

#my $continue = 1;
my $nb_better_allowed = 2
  ; #change that by computing it using the number of replicate and the known pvalue cutoff
$nb_better_allowed = $nb_random_sample if ( $FAST_CALL_FLAG == 0 );
my $nb_better = 0;
my $nb_random_sample_tested;
for (
	$nb_random_sample_tested = 0 ;
	$nb_random_sample_tested < $nb_random_sample
	&& $nb_better <= $nb_better_allowed ;
	$nb_random_sample_tested++
  )
{

	#print STDERR "\t **** $k\n";
	#1) compute the random IDs

	print STDERR " *** $nb_random_sample_tested ";    #<STDIN>;

	#$random_ID_mut = compute_random_ID(\@index_to_gene, -12);
	#$random_ID_dys = compute_random_ID(\@index_to_gene, $gene_ID);

	print STDERR " ---\n";                            #<STDIN>;

	$random_freq    = 0;
	@explained_line = ();    #0 gene not explaine for the sample, 1 explained
	for ( my $sample_ID = 0 ; $sample_ID < @sample_order ; $sample_ID++ ) {
		$sample = $sample_order[$sample_ID];
		$explained_line[$sample_ID] = 0;

		#maintain the dysregulated genes
		#$random_ID_dys = compute_random_ID(\@index_to_gene, $gene_ID);

		#we perform the search only if the gene is dysregulated for the sample
		if (
			$all_sample_gene_dysregulated{$sample}->[$gene_ID]->[$dys_status] !=
			0 )
		{

			$random_ID_mut = compute_random_ID( \@index_to_gene, -12 );

			@random_sample_gene_dysregulated = ();
			@random_sample_gene_mutated_list = ();
			%random_sample_gene_mutated      = ();

			#print STDERR " **** $sample\n";#<STDIN>;
			#compute the random ID for that sample
			for ( my $ID = 0 ; $ID < @{$random_ID_mut} ; $ID++ ) {

				#print $ID."\n";<STDIN>;

#$random_sample_gene_dysregulated[$ID] = $all_sample_gene_dysregulated{$sample}->[$random_ID_dys->[$ID]];
				$random_sample_gene_dysregulated[$ID] =
				  $all_sample_gene_dysregulated{$sample}->[$ID]
				  ;    #maintain the dysregulated genes

				$random_sample_gene_mutated_list[$ID] =
				  $all_sample_gene_mutated_list{$sample}
				  ->[ $random_ID_mut->[$ID] ];

#print STDERR $random_sample_gene_mutated_list[$ID]."\t".$all_sample_gene_mutated_list{$sample}->[$random_ID->[$ID]]."\t".$all_sample_gene_mutated_list{$sample}->[$ID]."\n";<STDIN>;
#if(exists $all_sample_gene_mutated{$sample}->{$random_ID_mut->[$ID]}){
#if(!$random_sample_gene_mutated_list[$ID]){
#print STDERR " ************ WEIRD problem in random process $ID !\n";
#}
				if ( $random_sample_gene_mutated_list[$ID] == 1 ) {
					$random_sample_gene_mutated{$ID} = 1;
				}
			}

#print STDERR " *** dys status ".$random_sample_gene_dysregulated[$gene_ID]->[$dys_status]."\n";#<STDIN>;

	#print STDERR "- construct explainend gene set for ".(get_name($gene))."\n";
			for ( $i = 0 ; $i < @index_to_gene ; $i++ ) { $explored[$i] = -1; }

#2) run the method to search if the gene_ID could be explained by at least a mutation given the parameters
			$explained_gene_set = construct_explained_gene_set_dijkstra(
				$gene_ID,
				\@random_sample_gene_dysregulated,
				\@random_sample_gene_mutated_list,
				\@explored, $depth_th, $hub_th, 1
			);

			#3) update the frequency
			foreach $eg ( @{$explained_gene_set} ) {
				if ( exists $random_sample_gene_mutated{$eg} ) {
					$random_freq++;

#$mut_gene_ID = $random_ID_mut->[$eg];
#print STDERR
#"\t ---- $eg"."\t".
#(get_name($eg, \@index_to_gene))." ".(keys %{$mutated_gene_frequency[$eg]})." -> ".(get_name($mut_gene_ID, \@index_to_gene))." ".(keys %{$mutated_gene_frequency[$mut_gene_ID]})." **** $random_freq\n";
					$explained_line[$sample_ID] = 1;
					last;
				}
			}

			#We have already study enough sample to find a better frequency
			if ( $FAST_CALL_FLAG && $random_freq == $THE_GENE_FREQ_EXPLAINED ) {
				last;
			}
		}
	}

	#4) store the final frequency
	#print STDERR " -> $random_freq\n";
	push( @all_random_freq, $random_freq );
	$nb_better++ if ( $random_freq >= $THE_GENE_FREQ_EXPLAINED );
	print OUT_TABLE "" . ( join( "\t", @explained_line ) ) . "\n"
	  if ($FLAG_PRINT_TABLE);
}
close(OUT_TABLE) if ($FLAG_PRINT_TABLE);

#5) compute the pvalue
$pvalue = 0;
for ( my $p = 0 ; $p < @all_random_freq ; $p++ ) {
	$pvalue++ if ( $all_random_freq[$p] >= $THE_GENE_FREQ_EXPLAINED );
}
$pvalue =   $pvalue / $nb_random_sample_tested; #should be equal to the number of random test in the case of a correct solution

open(OUT, ">$out_file");
print OUT ""
    . ( get_name( $gene_ID, \@index_to_gene ) ) . "_"
    . $dys_status_corress[$dys_status] . "\t"
    . $THE_GENE_FREQ_EXPLAINED . "\t"
    . $pvalue . "\t"
    . ( $nb_random_sample - $nb_random_sample_tested ) . "\n";
close(OUT);

sub compute_random_ID {
	my ( $true_ID, $fixed_ID ) = @_;

	my @random_ID = ();
	for ( $j = 0 ; $j < @{$true_ID} ; $j++ ) {
		$random_ID[$j] = $j;
	}

	$max = @random_ID - 1;
	$max-- if ( $max == $fixed_ID );

	while ( $max != 0 ) {
		$r = int( rand($max) );

		if ( $r != $fixed_ID ) {
			$rand_gene_ID    = $random_ID[$r];
			$random_ID[$r]   = $random_ID[$max];
			$random_ID[$max] = $rand_gene_ID;

			$max--;
			$max-- if ( $max == $fixed_ID );
		}
	}

	return \@random_ID;
}

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
		$dist_to_mutated_gene, $depth_th, $hub_th, $take_mutation )
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

	  #case were only search if a mutated gene could explain a dysregulated gene
	  #we have a mutated gene
			if ( $take_mutation && $sample_gene_mutated->[$gene] ) {
				last;
			}

			$max_gene_dist = $gene_dist if ( $gene_dist > $max_gene_dist );
			if ( $gene == $mutated_gene || @{ $connections[$gene] } <= $hub_th )
			{
				foreach $neigh ( @{ $connections[$gene] } ) {

#print STDERR "\t".$neigh."\t".@{$sample_gene_dysregulated}."\t".$sample_gene_dysregulated->[$neigh]->[0]."\t".$sample_gene_dysregulated->[$neigh]->[1]."\t".$sample_gene_mutated->[$neigh]."\t".@{$connections[$neigh]}."\n";

					#to have an additional control on the neighborhood
					if (
						( $nb_joker != 0 ) ||    #at least 1 jocker
						(
							   $sample_gene_dysregulated->[$neigh]->[0] != 0
							|| $sample_gene_dysregulated->[$neigh]->[1] != 0
						)
						||    #OR the neighbour is dysregulated
						(
							   $take_mutation
							&& $sample_gene_mutated->[$neigh] != 0
						)
					  ) #OR the neighbour is a mutation and we allows to keep them
					{

						#print STDERR "\t".$neigh." --------------> ADD\n";

					  #update the queue if it the first time we explore the gene
						push( @queue, $neigh )
						  if ( $dist_to_mutated_gene->[$neigh] == -1 )
						  ;    #not already explored;

						#update the distance
						$edge_weight = 100000;
						if (
							   @{$sample_gene_dysregulated} == 0
							|| $sample_gene_dysregulated->[$neigh]->[0] != 0
							|| $sample_gene_dysregulated->[$neigh]->[1] != 0
							||    #dysregulated
							(
								$take_mutation && $sample_gene_mutated->[$neigh]
							)     #is a mutation
						  )
						{
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

#print $DIST_MAX "".(get_name($mutated_gene, \@index_to_gene))."\t".$max_gene_dist."\n" if($max_gene_dist != 0 );
#print STDERR "---- THE RES: ".join("\t", @res)."\n";<STDIN>;
	return \@res;
}

#perform a breath first search to search from path of dysregulated gene of length <= $depth_th that do not contain hub_gene
sub construct_explained_gene_set {
	my ( $mutated_gene, $sample_gene_dysregulated, $sample_gene_mutated,
		$dist_to_mutated_gene, $depth_th, $hub_th )
	  = @_;

	#print STDERR "\nconstruct_explained_gene_set $mutated_gene\n";
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

