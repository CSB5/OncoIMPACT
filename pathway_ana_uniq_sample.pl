#!/usr/bin/perl
use warnings;

my ($data_dir, $large_data_dir, $js_file, $large_data_explained_pvalue_file, $large_data_module_file, $large_data_set_driver_stats, $network_type, $script_dir) = @ARGV;

$main_result_dir = "$data_dir/../ANALYSIS/";
run_exe("mkdir -p $main_result_dir") unless(-d $main_result_dir);

#script path
my $cluster_algo_path ="$script_dir/10_Cluster_Algo.pl";
my $module_inference_path = "$script_dir/11_Filter_And_Merge_Module_Impact.pl";
my $result_path = "$script_dir/export_gene_list.pl";


#4 analyse the js file to to obtain the optimal parameters
open(FILE, $js_file);
my $best_js = 0;
my $best_hub_value = 0;
my $best_depth_value = 0;
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    #if($line[0] == 1 && $best_js < $line[3]){
    if($best_js < $line[3]){
	$best_js = $line[3];
	$best_log2_fold_change = $line[0];
	$best_hub_value = $line[1];
	$best_depth_value = $line[2];
    }
}

print STDERR " *** $best_js $best_log2_fold_change $best_hub_value $best_depth_value END\n";#<STDIN>;
$res_dir = "$main_result_dir/RES_$best_depth_value\_$best_hub_value\_$best_log2_fold_change";


#Compute the expression/mutation frequency for the best parameters
$exe = "$cluster_algo_path $data_dir $network_type $best_depth_value $best_hub_value 0 $best_log2_fold_change $res_dir $script_dir 2> /dev/null";
run_exe($exe);
    
##########################################
# III FILTER PASSENGER AND MERGE MODULES #
##########################################

print STDERR "Construct Module [ Stringent Mode ]\n";
$exe = "$module_inference_path DRIVER_SAMPLE $network_type $large_data_explained_pvalue_file 0.05 $best_depth_value $best_hub_value $res_dir/FINAL_MODULE_SAMPLE.dat $script_dir $large_data_dir,$large_data_module_file $data_dir,$res_dir/MODULE.dat";
run_exe($exe);
#
print STDERR "Output Result -> $main_result_dir/GENE_LIST_SAMPLE\n";
$gene_list_dir = "$main_result_dir/GENE_LIST_SAMPLE";
run_exe("mkdir $gene_list_dir") unless(-d $gene_list_dir);
run_exe("cp $res_dir/FINAL_MODULE_SAMPLE.dat $gene_list_dir/FINAL_MODULE.dat");
$exe = "$result_path $data_dir $gene_list_dir/FINAL_MODULE.dat $network_type $best_log2_fold_change $large_data_set_driver_stats $gene_list_dir $script_dir 2> /dev/null";
run_exe($exe);
    
print STDERR "Construct Module [ Relaxed Mode ]\n";
$exe = "$module_inference_path DRIVER_ALL $network_type $large_data_explained_pvalue_file 0.05 $best_depth_value $best_hub_value $res_dir/FINAL_MODULE.dat $script_dir $large_data_dir,$large_data_module_file $data_dir,$res_dir/MODULE.dat";
#
print STDERR "Output Result -> $main_result_dir/GENE_LIST \n";
run_exe($exe);
$gene_list_dir = "$main_result_dir/GENE_LIST";
run_exe("mkdir $gene_list_dir") unless(-d $gene_list_dir);
run_exe("cp $res_dir/FINAL_MODULE.dat $gene_list_dir/FINAL_MODULE.dat");
$exe = "$result_path $data_dir $gene_list_dir/FINAL_MODULE.dat $network_type $best_log2_fold_change $large_data_set_driver_stats $gene_list_dir $script_dir 2> /dev/null";
run_exe($exe);

sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";;
    print STDERR `$exe` if($run);
}



#sort -k7,7 -nr Ovarian_sample/GENE_LIST/ALTERATION.dat | more
#grep ZHX1 Ovarian_sample/GENE_LIST/SAMPLE_SPECIFIC_DATA/mut_*.dat | cut -f 2 | sort | uniq -c






