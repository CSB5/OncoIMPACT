#!/usr/bin/perl
use warnings;

my (
	$data_dir, $network_type, $depth_th,
	$hub_th,   $nb_joker,     $fold_change_threshold,
	$out_dir,  $flag_real,    $script_dir
) = @ARGV;

#$flag_real = 0;
#$dir_TCGA_sample -> directory where the data are organize
#$network_gene_list -> ~/Software/netbox/data/Homo_sapiens.gene_info
#$network -> ~/Software/netbox/db/netbox.script
#$depth_th -> Threshold on maximum depth for defining module expansion
#$hub_th -> No Gene that is connected to more that $hub_th other gene will be included.
#flag_real in this case out_dir is the root of all the experiments, the run must be performed on the full data set and all the real directory will be upadeted at the same time at the end

require "$script_dir/Construct_network.pl";

my @length_th = ();
for ( my $i = 2 ; $i <= 20 ; $i += 2 ) {
	push( @length_th, $i );
}
push( @length_th, 1000000 );

my @dys_status_corress = ( "DOWN", "UP" );

#   1. Indexing
#not used due to synonyms problems/confusions

my %gene_to_index;
my @index_to_gene;
my @connections;

#print STDERR " *** --- $network_type\n";
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

# 4. Construct the set of explained genes for each samples
my $nb_sample = 0;

#Expression stats
my %explained_gene_frequency_all_depth = ();
foreach $depth (@length_th) {
	my @explained_gene_frequency = ();
	for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
		my @tab_status = ();
		for ( $j = 0 ; $j < 2 ; $j++ ) {
			my %hash = ();
			$tab_status[$j] = \%hash;
		}
		$explained_gene_frequency[$i] = \@tab_status;
	}
	$explained_gene_frequency_all_depth{$depth} = \@explained_gene_frequency;
}

my @dysregulated_gene_frequency = ();
for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	my @tab_status = ();
	for ( $j = 0 ; $j < 2 ; $j++ ) { my %hash = (); $tab_status[$j] = \%hash; }
	$dysregulated_gene_frequency[$i] = \@tab_status;
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

my $nb_sample_to_test = -10;

my %BUG_not_present_in_network = ();

foreach my $dir_sample (@the_DATA_DIR) {
    $mutation_file_name      = "$data_dir/$dir_sample/Genelist_Status.txt";
    if ( -e $mutation_file_name) {
	@sample_gene_mutated_list = ();
	
	#initilazed the mutational/expression gene status
	%sample_gene_mutated = ();
	my @sample_gene_dysregulated;
	
	for ( my $i = 0 ; $i < @index_to_gene ; $i++ ) {
	    my @status_tab = ( 0, 0 );    #DOWN, UP
	    $sample_gene_dysregulated[$i] = \@status_tab;
	    $sample_gene_mutated_list[$i] = 0;
	}
	my %all_explained_gene_set = ();
	
	$nb_sample++;
	$nb_mutated_gene = 0;
	
	open( FILE, "$mutation_file_name" );
	print STDERR " *** read file $mutation_file_name\n";
	
	#read the file to obtain the dysregulated and mutated genes
	while (<FILE>) {
	    chop($_);
	    @line = split( /\t/, $_ );
	    my @parts     = split( /_/, $line[0] );
	    my $gene_name = $parts[0];
	    my $status    = $parts[1];
	    if ( exists $gene_to_index{$gene_name} )
	    { #filter out all the gene_name that do not belong to the input network
		my $gene_ID = get_ID( $gene_name, \%gene_to_index );
		
		if ( $status eq "MUT" || $status eq "AMPL" || $status eq "DEL" )
		{
		    
		    $sample_gene_mutated{$gene_ID} = 1;
		    $sample_gene_mutated_list[$gene_ID] = 1;
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
			
#print STDERR "|".$_."|\t".$sample_gene_dysregulated[$gene_ID]->[$status_ID]."\t".$status_ID."\n";
		    }
		}
	    }
	    else {
		$BUG_not_present_in_network{$gene_name} = 1;
	    }
	}
	close(FILE);
	
	if ( ( keys %sample_gene_mutated ) == 0 ) {
	    print STDERR "SAMPLE WITHOUT MUTATED GENES WEIRD !!!\n";
	    <STDIN>;
	}
	
	foreach my $gene ( keys %sample_gene_mutated ) {
	    $nb_mutated_gene++;
	    
	    #print "****MUT TTT $gene ".(get_name($gene))."\n";
	    #print STDERR "- construct explainend gene set for ".(get_name($gene))."\n";
	    for ( $i = 0 ; $i < @index_to_gene ; $i++ ) { $explored[$i] = -1; }
	    
	    #my @false = ();
	    $explained_gene_set =
		construct_explained_gene_set_dijkstra( $gene,
						       \@sample_gene_dysregulated, \@sample_gene_mutated_list,
						       \@explored, $depth_th, $hub_th );
	    
	    if ( @{$explained_gene_set} > 1 ) {
		
		#compute the frequency
		foreach $eg ( @{$explained_gene_set} ) {
		    for ( my $dys_status = 0 ; $dys_status < 2 ; $dys_status++ )
		    {
			if ( $sample_gene_dysregulated[$eg]->[$dys_status] )
			{    #to remove the mutated genes and the jocker!
			    
			    #the gene is explain at least by gene at that depth
			    foreach $depth (@length_th) {
				if ( $explored[$eg] <= $depth ) {
				    $explained_gene_frequency_all_depth{$depth}
				    ->[$eg]->[$dys_status]->{$dir_sample} =
					$gene;
				}
			    }
			}
		    }
		}
	    }
	}
	last if ( $nb_sample == $nb_sample_to_test );
    }
}

$BUG_nb_not_present_in_network = ( keys %BUG_nb_not_present_in_network );


#show the explained gene frequency
#add the expression frequency and look at the spearman correlation !!!!
#Output the mutated gene set that explained each gene
my $res_dys;
my $res_exp;

if ( !$flag_real ) {
    foreach $depth_th (@length_th) {
	$f = "$out_dir/exp_gene_freq_$depth_th\_$hub_th.dat.gz";
	open( OUT, " | gzip -c >$f" );

	#$OUT = $all_out{$depth_th};
	for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
	    for ( my $dys_status = 0 ; $dys_status < 2 ; $dys_status++ ) {
		
		$res_dys = keys %{ $dysregulated_gene_frequency[$i]->[$dys_status] };
		next if($res_dys == 0);

		$res_exp = keys %{ $explained_gene_frequency_all_depth{$depth_th}->[$i]->[$dys_status] };
		
		#print STDERR "".(sprintf("%.3f",$res_dys/$nb_sample))."\n";<STDIN>;
		$g = get_name( $i, \@index_to_gene );

		#print STDERR " *** gene_name: $g $i\n";<STDIN>;
		print OUT "" 
		    #. $g . "_" . $dys_status_corress[$dys_status] . "\t"
		    #. ( @{ $connections[$i] } ) . "\t"
		    . $res_dys . "\t"
		    #. ( sprintf( "%.3f", $res_dys / $nb_sample ) ) . "\t"
		    . $res_exp 
		    #. "\t"
		    #. ( sprintf( "%.3f", $res_exp / $nb_sample ) ) 
		    . "\n";
	    }
	}
	close(OUT);
    }
}
else {
    opendir( DIR, $out_dir );
    @the_PARAM_DIR = readdir(DIR);
    close(DIR);
    print STDERR " *** SECOND STEP @the_PARAM_DIR\n";
    foreach my $real_sub_sample_dir (@the_PARAM_DIR) {
	print STDERR " *** file name $real_sub_sample_dir\n";
	if ( index( $real_sub_sample_dir, "REAL" ) != -1 ) {

	    #Get the sample of the real subsample data
	    opendir( DIR, "$out_dir/$real_sub_sample_dir" );
	    @the_real_sample_set = readdir(DIR);
	    close(DIR);
	    %the_real_sample_map = ();
	    foreach my $s (@the_real_sample_set) {
		$the_real_sample_map{$s} = 1;
	    }
	    $dir_res =
		"$out_dir/$real_sub_sample_dir/EXPLAINED_FREQ_DIFF_$fold_change_threshold/";
	    `mkdir $dir_res` if ( !-d $dir_res );

	    foreach $depth_th (@length_th) {
		$f = "$dir_res/exp_gene_freq_$depth_th\_$hub_th.dat.gz";
		print STDERR " *** res_file $f\n";
		open( OUT, " | gzip -c >$f" );
		for ( $i = 0 ; $i < @index_to_gene ; $i++ ) {
		    for ( my $dys_status = 0 ; $dys_status < 2 ; $dys_status++ )
		    {
			$res_exp = 0;
			$res_dys = 0;
			foreach my $s ( keys %the_real_sample_map ) {

			    #print STDERR " *** sss $s\n";<STDIN>;
			    $res_exp++
				if (
				    $explained_gene_frequency_all_depth{$depth_th}
				    ->[$i]->[$dys_status]->{$s} );
			    $res_dys++
				if (
				    $dysregulated_gene_frequency[$i]->[$dys_status]
				    ->{$s} );
			}

			next if($res_dys == 0);

			#print STDERR "".(sprintf("%.3f",$res_dys/$nb_sample))."\n";<STDIN>;
			$g = get_name( $i, \@index_to_gene );
			
			#print STDERR " *** gene_name: $g $i\n";<STDIN>;
			print OUT "" 
			    #. $g . "_" . $dys_status_corress[$dys_status] . "\t"
			    #. ( @{ $connections[$i] } ) . "\t"
			    . $res_dys . "\t"
			    #. ( sprintf( "%.3f", $res_dys / $nb_sample ) ) . "\t"
			    . $res_exp 
			    #. "\t"
			    #. ( sprintf( "%.3f", $res_exp / $nb_sample ) ) 
			    . "\n";
		    }
		}
		close(OUT);
	    }
	}
    }
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
	return \@res;
}
