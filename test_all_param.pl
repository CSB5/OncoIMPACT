#!/usr/bin/perl
use warnings;

#use strict;

my (
    $test_param_dir,                 $network_type,
    $min_log2_fold_change_threshold, $max_log2_fold_change_threshold,
    $min_hub_threshold, $max_hub_threshold,
    $flag_all_sample_used,            $nb_random_sample,
    $nb_process,                     $script_dir, $test_case
) = @ARGV;

my @hub_th;
for ( my $i = $min_hub_threshold ; $i <= $max_hub_threshold ; $i += 5 ) {
	push( @hub_th, $i );
}
push( @hub_th, 1000000 ) if($max_hub_threshold >= 100);

my @fold_change_th = ();
for (
	my $i = $min_log2_fold_change_threshold ;
	$i <= $max_log2_fold_change_threshold ;
	$i += 0.5
  )
{
	push( @fold_change_th, $i );
}

my $nb_max_line = 1000000;
my $nb_lines    = 1;
my $num_file    = 0;
my @all_cmds    = ();

$cmd_file_prefix = "param_cmd";
for ( $i = 0 ; $i < $nb_random_sample ; $i++ ) {
    if ( $flag_all_sample_used eq "ALL" ) {
	my $real_dir = "$test_param_dir/REAL_$i";
	if ( !-d $real_dir ) {
	    $exe = "ln -s REAL_ALL/ $test_param_dir/REAL_$i";
	    #print STDERR $exe . "\n";
	    `$exe`;
	}
    }
    for ( $k = 0 ; $k < @fold_change_th ; $k++ ) {
	$log2_fold_change_threshold = $fold_change_th[$k];
	for ( $j = 0 ; $j < @hub_th ; $j++ ) {
	    
	    #random data
	    $data_dir_name = "$test_param_dir/RANDOM_$i";
	    push_call($data_dir_name);
	    
	    #real data
	    if ( $i == 0 ) {
		$data_dir_name = "$test_param_dir/REAL_ALL/";
		push_call($data_dir_name);
	    }
	}
    }
}

close(OUT);

foreach $cmd_file (@all_cmds) {
    #$cmd = "cat $cmd_file | xargs -I cmd -P $nb_process bash -c cmd > /dev/null \n";
    $cmd = "cat $cmd_file | xargs -L 1 -P $nb_process perl &> /dev/null \n";
    print STDERR $cmd;
    print `$cmd`;
}

sub push_call {
	my ($data_dir) = @_;

	$flag_real = 0;
	$dir_res   = "$data_dir/EXPLAINED_FREQ_DIFF_$log2_fold_change_threshold";

	if ( index( $data_dir, "REAL" ) != -1 && $flag_all_sample_used ne "ALL" ) {
		$flag_real = 1;
		$dir_res   = "$data_dir/..";
	}

	`mkdir -p $dir_res` unless ( -d $dir_res );

	$file_res_pre = "$dir_res/exp_gene_freq_10_$hub_th[$j]";
	$file_res     = "$file_res_pre.dat.gz";

	if ( !-e $file_res ) {

		#print STDERR "mkdir $dir_res\n";
		#print STDERR " ********** $file_res\n";<STDIN>;

		$exe =
"$script_dir/10_Cluster_Algo_fast.pl $data_dir $network_type 1000000000 $hub_th[$j] 0 $log2_fold_change_threshold $dir_res $flag_real $script_dir 2> /dev/null";

		if ( $num_file == 0 || $nb_lines > $nb_max_line ) {
			close(OUT) if ( $num_file != 0 );
			$nb_lines = 1;
			$num_file++;
			$cmd_file = "$test_param_dir/$cmd_file_prefix\_$num_file.txt";
			push( @all_cmds, $cmd_file );
			open( OUT, ">$cmd_file" );
		}
		$nb_lines++;
		print OUT $exe . "\n";
	}
}
