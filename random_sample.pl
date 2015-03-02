#!/usr/bin/perl
use warnings;

my ( $data_dir, $network_type, $nb_real_sample_used, $nb_random_sample,
	$flag_mutation_fixed, $out_dir, $script_dir, $seed )
  = @ARGV;

require "$script_dir/Construct_network.pl";
my %gene_to_index;
my @index_to_gene;
my @connections;
print STDERR " *** construct network\n";
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

if($seed != -1){
    srand($seed);
}
`rm -r $out_dir`;
`mkdir $out_dir`;

opendir( DIR, $data_dir );
@the_DATA_DIR = readdir(DIR);
close(DIR);

#Compute the set of network genes that have an expression value
#gene without expression will not be swapped during the random process => the random sample will have the same graph as the true sample
my @index_to_gene_dys;
my %gene_dys_to_index;
my $ID_cmp = 0;
print STDERR " *** read real samples\n";
foreach my $dir_s (@the_DATA_DIR) {
    $mutation_file_name      = "$data_dir/$dir_s/Genelist_Status.txt";
    if ( -e $mutation_file_name ) {

	open( FILE, "$mutation_file_name" );

	#print STDERR "********** read file $mutation_file_name\n";#<STDIN>;

	while (<FILE>) {
	    chop($_);
	    @line  = split( /\t/, $_ );
	    @parts = split( /_/,  $line[0] );
	    $gene_name = $parts[0];
	    $status    = $parts[1];

	    #count the dysregulated genes that belong to the network
	    if ( ( $status eq "UP" || $status eq "DOWN" )
		 && exists $gene_to_index{$gene_name} )
	    {
		if ( !exists $gene_dys_to_index{$gene_name} ) {
		    $gene_dys_to_index{$gene_name} = $ID_cmp;
		    $index_to_gene_dys[$ID_cmp] = $gene_name;
		    $ID_cmp++;
		}
	    }
	}
	close(FILE);
    }
}

#random data set that maintain the frequencies
my @all_random_mutation_ID   = ();
my @all_random_expression_ID = ();
print STDERR " *** construct random ID\n";
for ( $i = 0 ; $i < $nb_random_sample ; $i++ ) {
    $all_random_mutation_ID[$i]   = compute_random_ID( \@index_to_gene );
    $all_random_expression_ID[$i] = compute_random_ID( \@index_to_gene_dys );
}

my @all_real_dir = ();
foreach my $dir_s (@the_DATA_DIR) {
    $mutation_file_name      = "$data_dir/$dir_s/Genelist_Status.txt";
    if ( -e $mutation_file_name ) {
	push( @all_real_dir, $dir_s );
    }
}

#Compute the number of sample used
#$nb_real_sample_used = @all_real_dir;
#if ( $nb_real_sample_usedfraction_real_sample_used ne "ALL" ) {
#    $nb_real_sample_used =
#	sprintf( "%.0f", $nb_real_sample_used * $fraction_real_sample_used );
#}

my @sample_mutated_gene_name = ();
my @sample_dys_gene_name     = ();
my @real_data_set_used       = ();
print STDERR " *** write random samples\n";

#Construct the REAL_ALL data set
`mkdir $out_dir/REAL_ALL`;
foreach my $dir_s (@all_real_dir) {
    $mutation_file = "Genelist_Status.txt";
    `mkdir $out_dir/REAL_ALL/$dir_s/`;
    `ln -s $data_dir/$dir_s/$mutation_file $out_dir/REAL_ALL/$dir_s/Genelist_Status.txt`;
}

for ( $i = 0 ; $i < $nb_random_sample ; $i++ ) {
    
    $dir_name_rand = "$out_dir/RANDOM_$i";
    `mkdir $dir_name_rand` unless ( -d $dir_name_rand );
    
#The real data dirctory are only constructed once if all the sample are analyzed
#fake directory linking to REAL_0 will be created during the optimal parameter inference
    $dir_name_real = "$out_dir/REAL_$i";
    if ( !-d $dir_name_real
	 && ( $nb_real_sample_used != @all_real_dir ) )
    {
	`mkdir $dir_name_real`;
    }
    
    #randomly choose X element of the real data
    @real_data_set_used = ();
    if ( $nb_real_sample_used == @all_real_dir ) {
	@real_data_set_used = @all_real_dir;
    }
    else {
	random_data_sampling( \@all_real_dir, \@real_data_set_used,
			      $nb_real_sample_used );
    }
    
    foreach my $dir_s (@real_data_set_used) {
	
	$mutation_file           = "Genelist_Status.txt";
	$mutation_file_name      = "$data_dir/$dir_s/Genelist_Status.txt";
		
	#directory to store the random data
	`mkdir $dir_name_rand/$dir_s`;
	
	#for the corresponding "truncated" real data set
	if ( $nb_real_sample_used != @all_real_dir ) {
	    `mkdir $dir_name_real/$dir_s`;
	    `ln -s $data_dir/$dir_s/$mutation_file $dir_name_real/$dir_s/Genelist_Status.txt`;
	}
	
	open( FILE, "$mutation_file_name" );
	
	#print STDERR "********** read file $mutation_file_name\n";#<STDIN>;
	@sample_mutated_gene_name = ();
	@sample_dys_gene_name     = ();
	
	while (<FILE>) {
	    chop($_);
	    @line  = split( /\t/, $_ );
	    @parts = split( /_/,  $line[0] );
	    $gene_name = $parts[0];
	    $status    = $parts[1];
	    
	    #count the mutated gene that belong to the network
	    
	    if ( exists $gene_to_index{$gene_name} ) {    #
		if ( $status eq "MUT" || $status eq "AMPL" || $status eq "DEL" )
		{
		    push( @sample_mutated_gene_name, $gene_name );
		    
		    #$rand_gene_ID = get_name;
		    #print OUT "".get_name($rand_gene_ID, \@index_to_gene)."_MUT\n";
		    
		}
		if ( exists $gene_dys_to_index{$gene_name}
		     && ( $status eq "UP" || $status eq "DOWN" ) )
		{
		    push( @sample_dys_gene_name, $_ );
		    
		}
	    }
	}
	close(FILE);
	
	#$nb_dys_gene_sample = @sample_dys_gene_name + 0;
	
#construct random mutation and expression gene without changing their frenquencies
	
	#Mutation
	$out_file = "$dir_name_rand/$dir_s/Genelist_Status.txt";
	open( OUT, ">$out_file" );
	
	for ( $j = 0 ; $j < @sample_mutated_gene_name ; $j++ ) {
	    $name = $sample_mutated_gene_name[$j];
	    $ID = get_ID( $name, \%gene_to_index );    #NO RANDOM ID
	    
	    $new_ID = $all_random_mutation_ID[$i]->[$ID];
	    $new_name = get_name( $new_ID, \@index_to_gene );
	    
	    if ( $flag_mutation_fixed eq "MUT_FIXED" ) {
		print OUT "" . $name . "_MUT\n";
	    }
	    else {
		print OUT "" . $new_name . "_MUT\n";
	    }
	    
	    #print STDERR "".$name."\t".$new_name."_MUT\n";<STDIN>;
	}
	
	#Dys
	for ( $j = 0 ; $j < @sample_dys_gene_name ; $j++ ) {
	    
	    #print STDERR "".$sample_dys_gene_name[$j]."\n";<STDIN>;
	    @line  = split( /\t/, $sample_dys_gene_name[$j] );
	    @parts = split( /_/,  $line[0] );
	    $name  = $parts[0];
	    $status = $parts[1];
	    
#to change to obtain 2 different random sampling
#$ID = get_ID($gene_name, \%gene_to_index);
#Only gene with an expression value in at least one real sample could have an expression value in the random samples
	    if ( exists $gene_dys_to_index{$name} ) {
		
		$new_ID = $all_random_expression_ID[$i]->[ $gene_dys_to_index{$name} ];
		$new_name = $index_to_gene_dys[$new_ID];
		print OUT ""
		    . $new_name . "_"
		    . $status . "\t"
		    . $line[1] . "\n";
		
		#print OUT "".$name."_".$status."\t".$line[1]."\n";
		
		if ( !exists $gene_to_index{$name} ) {
		    print STDERR
			" *** WEIRD random dysregulated gene that do not belong to the network "
			. $name . "\t"
			. $dir_s . "\n";
		}
		if ( !exists $gene_dys_to_index{$name} ) {
		    print STDERR
			" *** WEIRD random dysregulated gene that is never dysregulated in the real samples "
			. $name . "\t"
			. $dir_s . "\n";
		}
	    }
	    
	}
	close(OUT);
    }
}

sub random_data_sampling {
	my ( $all_real_dir, $real_data_set_used, $nb_real_sample_used ) = @_;
	my %take        = ();
	my $nb_real_dir = @{$all_real_dir};
	while ( @{$real_data_set_used} < $nb_real_sample_used ) {
		$r = int( rand($nb_real_dir) );
		if ( !exists $take{$r} ) {
			$take{$r} = 1;
			push( @{$real_data_set_used}, $all_real_dir->[$r] );
		}
	}
}

sub compute_random_ID {
	my ($true_ID) = @_;

	my @random_ID = ();
	for ( $j = 0 ; $j < @{$true_ID} ; $j++ ) {
		$random_ID[$j] = $j;
	}

	$max = @random_ID - 1;
	while ( $max != 0 ) {
		$r               = int( rand($max) );
		$rand_gene_ID    = $random_ID[$r];
		$random_ID[$r]   = $random_ID[$max];
		$random_ID[$max] = $rand_gene_ID;
		$max--;
	}

	return \@random_ID;
}

#./PIPELINE/10_Cluster_Algo_fast.pl TEST_RAND/RANDOM_0/ Homo_sapiens.gene_info netbox.script 1000000 25 0 3 TEST_RAND/
