#!/usr/bin/perl
use warnings;

my ( $step, $data_dir, $data_type, $nb_thread, $network_type, $script_dir, $test_case) = @ARGV;

print STDERR " **** pathway_ana step:$step data_dir:$data_dir nb_thread:$nb_thread network_type:$network_type script_dir:$script_dir test_case:$test_case\n";


$MAX_FRAC_DISREGULATED_GENE   = 0.5;
$MIN_MEDIAN_DISREGULATED_GENE = 300;

#$NB_SIMULATED_DATA_SET_PARAMETER = 100;
#$NB_SIMULATED_DATA_SET_PVALUE    = 500;

#
$SAMPLE_USED_DURING_PARAM_EXTIMATION = 50;
$NB_SIMULATED_DATA_SET_PARAMETER = 100;
#
$NB_SIMULATED_DATA_SET_PVALUE    = 500;
#
$EXPLAINED_FREQ_THRESHOLD = 0.05;
#
$SEED = -1;

#To be changed for RNA-SEQ data
my $MIN_LOG2_FOLD_CHANGE_THRESHOLD = 1;
my $MAX_LOG2_FOLD_CHANGE_THRESHOLD = 3;

if($data_type eq "RNA_SEQ"){
    $MIN_LOG2_FOLD_CHANGE_THRESHOLD = 1;
    $MAX_LOG2_FOLD_CHANGE_THRESHOLD = 1;
}

$MIN_HUB_THRESHOLD = 10;
$MAX_HUB_THRESHOLD = 100;

if(defined $test_case && $test_case eq "TEST"){
    $SAMPLE_USED_DURING_PARAM_EXTIMATION = 10;
    $NB_SIMULATED_DATA_SET_PARAMETER = 2;
    #
    $NB_SIMULATED_DATA_SET_PVALUE    = 5;
    #
    $EXPLAINED_FREQ_THRESHOLD = 0.6;
    #
    $MIN_LOG2_FOLD_CHANGE_THRESHOLD = 1.5;
    $MAX_LOG2_FOLD_CHANGE_THRESHOLD = 2;
    #
    $MIN_HUB_THRESHOLD = 40;
    $MAX_HUB_THRESHOLD = 50;
    #
    $SEED = 1;
}


$main_result_dir = "$data_dir/../ANALYSIS";
run_exe("mkdir $main_result_dir") unless ( -d $main_result_dir );
my $RUN_STATS            = 1;
my $RUN_TEST_PARAM       = 1;
my $RUN_DYS_SIGNIFICANCE = 1;

#my $RUN_DYS_IMPACT_SIGNIFICANCE = 1;
my $RUN_DRIVER_INFEREANCE = 1;


if ( $step eq "TEST_PARAM" ) {
	$RUN_DYS_SIGNIFICANCE  = 0;
	$RUN_DRIVER_INFEREANCE = 0;
}

if ( $step eq "DYS_SIGN" ) {
	$RUN_STATS      = 0;
	$RUN_TEST_PARAM = 0;
}

if ( $step eq "DRIVER_INFERENCE" ) {
	$RUN_TEST_PARAM       = 0;
	$RUN_DYS_SIGNIFICANCE = 0;
}

#script path
my $basic_stats_path      = "$script_dir/plot_basic_stats.pl";
my $sim_path              = "$script_dir/random_sample.pl";
my $test_param_path       = "$script_dir/test_all_param.pl";
my $compute_js_path       = "$script_dir/compute_all_js.pl";
my $cluster_algo_path     = "$script_dir/10_Cluster_Algo.pl";
my $phenotype_pvalue_path = "$script_dir/run_all_test_pvalue_table.pl";
my $module_inference_path = "$script_dir/11_Filter_And_Merge_Module_Impact.pl";
my $result_path           = "$script_dir/export_gene_list.pl";

#network file
# $network_gene_list = "Homo_sapiens.gene_info";
# $network = "netbox.script";

######################################################
# 0 DATA MANAGEMENT AND EXTRACTION (in construction) #
######################################################

#Give a simple overview of samples alteratation/dysregulation

$sample_stats_dir  = "$main_result_dir/SAMPLE_STATS";
$sample_stats_file = "$sample_stats_dir/basic_stats.dat";
if ($RUN_STATS) {
    
    run_exe("mkdir $sample_stats_dir") unless ( -d $sample_stats_dir );
    $exe = "$basic_stats_path $data_dir $network_type $script_dir >  $sample_stats_file 2> /dev/null";
    run_exe($exe);
}

#die
#chdir("../../..") or die;
#chdir("../..") or die;
#print STDERR " --- END";<STDIN>;
#<STDIN>;

#########################
# I PARAMETER INFERANCE #
#########################
$test_param_dir = "$main_result_dir/TEST_PARAM";
$js_file = "$test_param_dir\_js.dat";

if ($RUN_TEST_PARAM) {

	$NB_GENE_IN_NETWORK = 0;
	$NB_GENE_IN_NETWORK = 9448 if ( $network_type eq "DRIVER_NET" );
	$NB_GENE_IN_NETWORK = 9261 if ( $network_type eq "NETBOX" );

	
#to compute the min and max log fold change that will be in the interval [1 - 3] with op of 0.5
	my %median_sample_diff_gene = ();
	compute_median_sample_diff( $sample_stats_file, \%median_sample_diff_gene );
		

	#Some check for array data to avaoid getting too much dysregulated genes. 
	#Need to update for RNA_SEQ data based on normalized read depth
	if($data_type eq "ARRAY"){
	    #min should not have a sample median diff gene larger that 50% of the gene with an expression
	    for (
		$i = $MIN_LOG2_FOLD_CHANGE_THRESHOLD ;
		$i <= $MAX_LOG2_FOLD_CHANGE_THRESHOLD ;
		$i += 0.5
		)
	    {
		if ( $median_sample_diff_gene{$i} / $NB_GENE_IN_NETWORK <
		     $MAX_FRAC_DISREGULATED_GENE )
		{
		    $MIN_LOG2_FOLD_CHANGE_THRESHOLD = $i;
		    last;
		}
	    }
	    
	    #max should not have a sample median diff gene smaller than 0.01 of the gene with an expression
	    
	    for (
		$i = $MAX_LOG2_FOLD_CHANGE_THRESHOLD ;
		$i >= $MIN_LOG2_FOLD_CHANGE_THRESHOLD ;
		$i -= 0.5
		)
	    {
		if ( $median_sample_diff_gene{$i} >= $MIN_MEDIAN_DISREGULATED_GENE ) {
		    $MAX_LOG2_FOLD_CHANGE_THRESHOLD = $i;
		    last;
		}
	    }
	}
	
	$nb_sample = `wc -l $sample_stats_file | cut -d " " -f 1`;
	chop $nb_sample;
	$nb_sample--;

	#Number of sample during the parameter estimation stage
	#In case of large sample size 50 samples will be sub-sample from the whole data set
	$nb_sample_used = $SAMPLE_USED_DURING_PARAM_EXTIMATION;
	$flag_all_sample_used = "FRACTION";
	if($nb_sample < $nb_sample_used){
	    $nb_sample_used = $nb_sample;
	    $flag_all_sample_used = "ALL";
	}


	#1 run the simulation
	#the simulation are perfomrmed only if the $test_param_dir is empty
	if ( !-d $test_param_dir ) {
	    `mkdir $test_param_dir`;
	    print STDERR "Parameter Estimation\n";
	    $exe = "$sim_path $data_dir $network_type $nb_sample_used $NB_SIMULATED_DATA_SET_PARAMETER MUT_UNFIXED $test_param_dir $script_dir $SEED 2> /dev/null";
	    run_exe($exe);
	}
	
#2 Run the inference method with different parameters
#this method do not re-run any the test of random/real sample that have been previously analysed (good in case of crash)
	if ( !-e $js_file ) {
	    $exe = "$test_param_path $test_param_dir $network_type $MIN_LOG2_FOLD_CHANGE_THRESHOLD $MAX_LOG2_FOLD_CHANGE_THRESHOLD $MIN_HUB_THRESHOLD $MAX_HUB_THRESHOLD $flag_all_sample_used $NB_SIMULATED_DATA_SET_PARAMETER $nb_thread $script_dir ";
	    run_exe($exe);
	    
	    #exit;
	    
#3 compute the JS distance of the explnained frequency of the genes of random data with the real data

	    $exe = "$compute_js_path $test_param_dir $NB_SIMULATED_DATA_SET_PARAMETER $MIN_LOG2_FOLD_CHANGE_THRESHOLD $MAX_LOG2_FOLD_CHANGE_THRESHOLD $MIN_HUB_THRESHOLD $MAX_HUB_THRESHOLD > $js_file";
	    run_exe($exe);
	    
	    #compress the directory
	    $exe = "tar -zcf $test_param_dir.tgz $test_param_dir 2> /dev/null";
	    run_exe($exe);    #<STDIN>;
	}
}

#4 analyse the js file to to obtain the optimal parameters
open( FILE, $js_file );
my $best_js          = 0;
my $best_hub_value   = 0;
my $best_depth_value = 0;
my $best_log2_fold_change = 0;
while (<FILE>) {
	chop $_;
	@line = split( /\t/, $_ );

	#if($line[0] == 1 && $best_js < $line[3]){
	if ( $best_js < $line[3] ) {
		$best_js               = $line[3];
		$best_log2_fold_change = $line[0];
		$best_hub_value        = $line[1];
		$best_depth_value      = $line[2];
	}
}

#Check validity of parameter values
if($best_log2_fold_change == 0 || $best_hub_value == 0 || $best_depth_value == 0){
    print STDERR " *** Aborting! The paramater values are not correctly initialized best_log2_fold_change:$best_log2_fold_change best_hub_value:$best_hub_value best_depth_value:$best_depth_value\n";
    exit 2;
}



#print STDERR " *** $best_js $best_log2_fold_change $best_hub_value $best_depth_value END\n" ;    #<STDIN>;
$res_dir =
"$main_result_dir/RES_$best_depth_value\_$best_hub_value\_$best_log2_fold_change";

#
# $res_dir = "Ovarian_sample/STABILITY_OVER_SUBSAMPLE_ANALYSIS/SAMPLE_316/REPLICATE_1/DRIVER_NET/RES_6_100_2/";
#

if ($RUN_DYS_SIGNIFICANCE) {
	###############################
	# II.2 PHENOTYPE SIGNIFICANCE #
	###############################
	#Compute the expression/mutation frequency for the best parameters

	print STDERR "Phenotype Significance\n";

	$exe = "$cluster_algo_path $data_dir $network_type $best_depth_value $best_hub_value 0 $best_log2_fold_change $res_dir $script_dir 2> /dev/null";
	run_exe($exe);

	#Compute the significance
	$FAST_CALL_FLAG = 1;
	$exe            = "$phenotype_pvalue_path $data_dir $main_result_dir $network_type $best_depth_value $best_hub_value $best_log2_fold_change $EXPLAINED_FREQ_THRESHOLD $NB_SIMULATED_DATA_SET_PVALUE $nb_thread $FAST_CALL_FLAG $script_dir $SEED";
	run_exe($exe);
}

#if($RUN_DYS_IMPACT_SIGNIFICANCE){
###############################
# II.2 PHENOTYPE SIGNIFICANCE #
###############################
#    $test_impact_dir = "$main_result_dir/TEST_IMPACT";
#the simulation are perfomrmed only if the $test_param_dir is empty
#    if(! -d $test_impact_dir){
#	`mkdir $test_impact_dir`;
#	$exe = "$sim_path $data_dir $network_type ALL $NB_SIMULATED_DATA_SET_PVALUE MUT_FIXED $test_impact_dir";
#	run_exe($exe);
#    }



if ($RUN_DRIVER_INFEREANCE) {
	##########################################
	# III FILTER PASSENGER AND MERGE MODULES #
	##########################################

    $freq_gene = "$res_dir/exp_gene_freq_pvalue.dat";
    
    if(! -e $freq_gene){
	print STDERR " Aborting! The file $freq_gene does not exists\n";
	exit 2;
    }

    print STDERR "Construct Module [ Stringent Mode ]\n";
    $exe = "$module_inference_path DRIVER_SAMPLE $network_type $freq_gene 0.05 $best_depth_value $best_hub_value $res_dir/FINAL_MODULE_SAMPLE.dat $script_dir $data_dir,$res_dir/MODULE.dat";
    run_exe($exe);
    #
    print STDERR "Output Result -> $main_result_dir/GENE_LIST_SAMPLE\n";
    $gene_list_dir = "$main_result_dir/GENE_LIST_SAMPLE";
    run_exe("mkdir $gene_list_dir") unless ( -d $gene_list_dir );
    run_exe("cp $res_dir/FINAL_MODULE_SAMPLE.dat $gene_list_dir/FINAL_MODULE.dat");
    $exe = "$result_path $data_dir $gene_list_dir/FINAL_MODULE.dat $network_type $best_log2_fold_change NONE $gene_list_dir $script_dir 2> /dev/null";
    run_exe($exe);
    
    print STDERR "Construct Module [ Relaxed Mode ]\n";
    $exe       = "$module_inference_path DRIVER_ALL  $network_type $freq_gene 0.05 $best_depth_value $best_hub_value $res_dir/FINAL_MODULE.dat $script_dir $data_dir,$res_dir/MODULE.dat";
    run_exe($exe);
    #
    print STDERR "Output Result -> $main_result_dir/GENE_LIST \n";
    $gene_list_dir = "$main_result_dir/GENE_LIST";
    run_exe("mkdir $gene_list_dir") unless ( -d $gene_list_dir );
    run_exe("cp $res_dir/FINAL_MODULE.dat $gene_list_dir/FINAL_MODULE.dat");
    $exe = "$result_path $data_dir $gene_list_dir/FINAL_MODULE.dat $network_type $best_log2_fold_change NONE $gene_list_dir $script_dir 2> /dev/null";
    run_exe($exe);
    
}

sub compute_median_sample_diff {
	my ( $file, $median_sample_diff_gene ) = @_;
	open( FILE, $file );
	@all_value    = ();
	@diff_cut_off = ();
	while (<FILE>) {
		@line = split( /\t/, $_ );
		if ( $line[0] eq "NAME" ) {
			$cmp = 1;
			for ( $i = 1 ; $i - $cmp <= 3 ; $i++ ) {
				my @tab = ();
				push( @all_value,    \@tab );
				push( @diff_cut_off, $i - $cmp );
				$cmp += 0.5;
			}
		}
		else {
			for ( $i = 0 ; $i < @all_value ; $i++ ) {
				push( @{ $all_value[$i] }, $line[ $i + 1 ] );
			}
		}
	}
	close(FILE);

	#compute the median
	for ( $i = 0 ; $i < @all_value ; $i++ ) {
		$median_sample_diff_gene->{ $diff_cut_off[$i] } =
		  compute_median( $all_value[$i] );
	}

}

sub compute_median {
	my ($data) = @_;
	@expr_values_sorted = sort { $a <=> $b } @{$data};
	if ( $#expr_values_sorted % 2 == 0 ) {
		return ( $expr_values_sorted[ ( ($#expr_values_sorted) / 2 ) - 1 ] +
			  $expr_values_sorted[ ( ($#expr_values_sorted) / 2 ) ] ) / 2;
	}
	else {
		return $expr_values_sorted[ ( $#expr_values_sorted + 1 ) / 2 ];
	}
}

sub run_exe {
	my ($exe) = @_;
	$run = 1;

	print STDERR $exe . "\n";
	print STDERR `$exe` if ($run);
}

#sort -k7,7 -nr Ovarian_sample/GENE_LIST/ALTERATION.dat | more
#grep ZHX1 Ovarian_sample/GENE_LIST/SAMPLE_SPECIFIC_DATA/mut_*.dat | cut -f 2 | sort | uniq -c

#./PIPELINE/pathway_ana.pl ALL GBM_other_sample/ 40 20 DRIVER_NET
#./PIPELINE/pathway_ana.pl ALL GBM_FULL/ 40 20 DRIVER_NET

