#!/usr/bin/perl
use warnings;

my ($test_param_dir, $nb_random_sample, $min_log2_fold_change_threshold, $max_log2_fold_change_threshold, $min_hub_threshold, $max_hub_threshold) = @ARGV;

my $OUT;
my @nb_significant_true;

#my @all_random_distribution = ();
$col_dysregulated_count = 0;
$col_explained_count = 1;

#Fold change threshold possible value
my @fold_change_th = ();
for(my $i = $min_log2_fold_change_threshold; $i <= $max_log2_fold_change_threshold; $i += 0.5){
    push(@fold_change_th, $i);
}


my @hub_th;
for(my $i = $min_hub_threshold; $i <= $max_hub_threshold; $i += 5){
    push(@hub_th, $i);
}
push(@hub_th, 1000000) if($max_hub_threshold >= 100);

#path length possible value
my @length_th = ();
for(my $i = 2; $i <= 20; $i += 2){
    push(@length_th, $i);
}
push(@length_th, 1000000);

my @true_distribution = ();
my @random_distribution = ();

for(my $u = 0; $u < @fold_change_th; $u++){
    for(my $j = 0; $j < @hub_th; $j++){
	for($k = 0; $k < @length_th; $k++){
	    next if(($length_th[$k] == 1000000 && $hub_th[$j] != 1000000) ||
		    ($length_th[$k] != 1000000 && $hub_th[$j] == 1000000));
	    @true_distribution = ();
	    @random_distribution = ();
	    
	    $res_file = "EXPLAINED_FREQ_DIFF_$fold_change_th[$u]/exp_gene_freq_$length_th[$k]_$hub_th[$j].dat.gz";
	    
	    #print STDERR " *** $hub_th[$j] $length_th[$k]\n";
	    for(my $i = 0; $i < $nb_random_sample; $i++){
		construct_distribution(\@random_distribution, "$test_param_dir/RANDOM_$i/$res_file");
		construct_distribution(\@true_distribution, "$test_param_dir/REAL_$i/$res_file");
	    }
	    #print STDERR "\n\n***** nb significant genes $fold_change_th[$u] $hub_th[$j] $length_th[$k]: ".(@random_distribution+0)." ".(@true_distribution)."\n ";
	    
	    #compute the jenson-shanon divergence
	    my $P_true;
	    my $P_random;
	    
	    $P_true = contruct_discrete_probability_density(\@true_distribution);
	    $P_random = contruct_discrete_probability_density(\@random_distribution);
	    	    
	    $js = jensen_shannon( $P_true, $P_random);
	    #print STDERR " ***** jensen-shanon $hub_th[$j] $length_th[$k]\t".(sprintf("%.5f", $js))."\t".(join("\t", @{$fdr}))."\n";#<STDIN>;
	    
	    $res = "$fold_change_th[$u]\t$hub_th[$j]\t$length_th[$k]\t".(sprintf("%.5f", $js))."\n";
	    print $res;
	    #print STDERR $res;
	}
    }
}

sub construct_distribution{
    my ($distribution, $file) = @_;
    
    if(-e $file && -s $file){
	#print STDERR " *** Use sample $file\n";
	#Read $in_file\n";
	open(FILE, "gunzip -c $file | ");
	while(<FILE>){
	    #print $_;
	    @line = split(/\t/, $_);
	    chop $_;
	    $dysregulated_count = $line[$col_dysregulated_count];
	    if($dysregulated_count != 0){
		$explained_count = $line[$col_explained_count];
		push (@{$distribution}, $explained_count);
	    }
	}
	close(FILE);
    }
    else{
	if(! -e $file){
	    print STDERR " *** Aborting! The file $file does not exists\n";exit 2; 
	}
	if(! -s $file){
	    print STDERR " *** Aborting! the file $file is empty\n";exit 2; 
	}
    }
}


sub  jensen_shannon {
    my ($P, $Q) = @_;           # $P=pointeur dist P, $Q=pointeur dist Q 
    my $log_2   = log(2) ;      # Logarithme base 2
    my %deja_vu = () ;          # Identifier les symboles de la distribution P

    foreach $k (keys %{$P}){$deja_vu{$k} = 1;}
    foreach $k (keys %{$Q}){$deja_vu{$k} = 1;}
    #map{ $deja_vu{$_}++ } (keys %$P) ;
    my @symbols = sort  {$a <=> $b} keys %deja_vu ;  # Symbols de la distribution P inter Q
    my $js      = 0 ;                   # init Jensen-Shannon
    foreach my $mot (@symbols) {        # Calcul de la divergence Jensen-Shannon
	$P_v = 0;
	$P_v = $P->{$mot} if (exists $P->{$mot});

	$Q_v = 0;
	$Q_v = $Q->{$mot} if (exists $Q->{$mot});
	
	#print STDERR $mot."\t".$P_v."\t".$Q_v."\n";<STDIN>;

        $js += $P_v * log( 2 * $P_v / ( $P_v + $Q_v ) ) / $log_2 if($P_v != 0);   # P Log ( P /(P+Q)/2 )
	$js += $Q_v * log( 2 * $Q_v / ( $P_v + $Q_v ) ) / $log_2 if($Q_v != 0);   # Q Log ( Q /(P+Q)/2 )
    }
    return $js/2 ;                      # retourner divergence Jensen-Shannon
}

sub contruct_discrete_probability_density{
    my ($T) = @_;
    my %P = ();
    my $total = 0;
    foreach $x (@{$T}){
	$total++;
	$P{$x} = 0 if(! exists $P{$x} );
	$P{$x}++;
    }

    foreach $x (keys %P){
	$P{$x} = $P{$x} / $total;
    }
    #print STDERR " *************** total: $total\n";

    $total = 0;
    foreach $x (keys %P){
	$total += $P{$x};
    }

    #print STDERR " *************** total: $total\n";

    return \%P;
}

