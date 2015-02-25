#!/usr/bin/perl
use warnings;

#If FAST_CALL_FLAG is true: the method stop if at least one random sample have a higher frequency than the real frequency: $THE_GENE_FREQ_EXPLAINED
my ($sample_dir, $main_result_dir, $network_type, $depth_th, $hub_th, $log2_fold_change_threshold, $explained_freq_threshold, $nb_random_test, $nb_process, $FAST_CALL_FLAG, $script_dir, $seed) = @ARGV;

require "$script_dir/Construct_network.pl";

$explained_file_body = "$main_result_dir/RES_$depth_th\_$hub_th\_$log2_fold_change_threshold/exp_gene_freq";
$sample_explained_file = "$explained_file_body.dat";

@dys_status_corress = ("DOWN", "UP");
my %gene_to_index;
my @index_to_gene;
my @connections;

#Construct the network
if($network_type eq "NETBOX"){
    construct_netbox_network($sample_dir, \@index_to_gene, \%gene_to_index, \@connections, $script_dir);
}
if($network_type eq "DRIVER_NET"){
    construct_driver_net_network(\@index_to_gene, \%gene_to_index, \@connections, $script_dir);
}

#Get the sample analysed
%sample_analysed = ();
opendir(DIR, $sample_dir);
@the_DATA_DIR = readdir(DIR);
close(DIR);
#print STDERR " *** READ DIR\n";
foreach my $sample (@the_DATA_DIR){
    $mutation_file_name = "$sample_dir/$sample/Genelist_Status.txt";
    if(-e $mutation_file_name){
	$sample_analysed{$sample} = 1;
    }
}
#<STDIN>;

#Get the frequency of the explained gene
%all_explained_freq = ();
open(FILE, $sample_explained_file);
while(<FILE>){
    chop $_;
    my @line = split(/\t/, $_);
    $all_explained_freq{$line[0]} = \@line;
}
close(FILE);
#<STDIN>;

$compute_pvalue = 1;
$compute_corrected_pvalue = 1;

$pvalue_dir = "$main_result_dir/PVALUE_TABLE/";

if($compute_pvalue){
    my $OUT;
    my $nb_max_line = 50000;
    my $nb_lines = 1;
    my $num_file = 0;
    my @all_cmds = ();
    
    `mkdir  $pvalue_dir`;
    
    for ($g = 0; $g < @index_to_gene; $g++){
	if($num_file == 0 || $nb_lines > $nb_max_line){
	    close($OUT) if($num_file != 0);
	    $nb_lines = 1;
	    $num_file++;
	    $cmd_file = "$pvalue_dir/cmd_$num_file.txt"; 
	    push(@all_cmds, $cmd_file);
	    open($OUT, ">$cmd_file");
	}
	
	for($d = 0; $d < 2; $d++){ 
	    $name = get_name($g, \@index_to_gene)."_".$dys_status_corress[$d];
	    $count_explained = $all_explained_freq{$name}->[4];
	    $freq_explained = $all_explained_freq{$name}->[5];
	    if($freq_explained >= $explained_freq_threshold){
		$pvalue_file = "$pvalue_dir/t_$g\_$d.dat";
		next if(-e $pvalue_file);

		$out_file_table = "$pvalue_dir/table_$g\_$d.dat.gz";
		$nb_lines++;
		if(-e $out_file_table){
		    #Read the table file if it exists
		    
		    #print STDERR " *** READ TABLE $out_file_table\n";#<STDIN>;
		    open(FILE, "gunzip -c $out_file_table |");
		    @sample_analysed_ID = ();
		    @all_random_freq = ();
		    $first = 1;
		    while(<FILE>){
			chop $_;
			@line = split(/\t/, $_);
			if($first == 1){
			    for($i = 0; $i < @line; $i++){
				$sample = $line[$i];
				if(exists $sample_analysed{$sample}){
				    push(@sample_analysed_ID, $i);
				}
				else{
				    #print STDERR " *** SAMPLE NOT ANALYSED $sample\n";#<STDIN>
				}
			    }
			}
			else{
			    $random_explained_freq = 0;
			    foreach $s (@sample_analysed_ID){
				#if(! defined $line[$s]){
				#print STDERR " *** $first value $s |$line[$s]|\n";#<STDIN>;
				#}
				$random_explained_freq += $line[$s];
			    }
			    push(@all_random_freq, $random_explained_freq);
			}
			$first++;
		    }
		    close(FILE);
		    #5) compute the pvalue
		    $pvalue = 0;
		    for(my $p = 0; $p < @all_random_freq; $p++){
			$pvalue++ if($all_random_freq[$p] >= $count_explained);
		    }
		    $pvalue = $pvalue / $nb_random_test;
		    open(OUT, ">$pvalue_file");
		    print OUT $name."\t".$count_explained."\t".$pvalue."\t".($nb_random_test - (@all_random_freq+0))."\n";
		    close(OUT);
		}
		else{
		    #Run the construction of the table
		    #$FAST_CALL_FLAG = 0;
		    $exe = "$script_dir/10_Cluster_Algo_pvalue_table.pl $sample_dir $network_type $depth_th $hub_th 0 $log2_fold_change_threshold $nb_random_test $g $d $count_explained $pvalue_file $out_file_table  $FAST_CALL_FLAG $script_dir $seed 2> /dev/null";
		    print $OUT $exe."\n";
		}
	    }
	}
	close(OUT);
    }
    foreach $cmd_file(@all_cmds){
	#$cmd = "cat $cmd_file | xargs -I cmd -P $nb_process bash -c cmd > /dev/null \n";
	$cmd = "cat $cmd_file | xargs -L 1 -P $nb_process perl &> /dev/null \n";
	print STDERR $cmd;#<STDIN>;
	print `$cmd`;
    }
}

if($compute_corrected_pvalue){
    #take the pvalue
    $all_pvalue_file = "$pvalue_dir/RES.dat";
    %all_pvalue = ();
    #if(! -e $all_pvalue_file){
    `cat $pvalue_dir/t* > $all_pvalue_file`;
    #}
    
    open(FILE, $all_pvalue_file);
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	#$gene_s = $line[0];

	my @parts = split(/_/,$line[0]);
	my $gene_name = $parts[0];
	my $status = $parts[1];
	if(exists $gene_to_index{$gene_name}){
	    $gene_s = get_name($gene_to_index{$gene_name}, \@index_to_gene)."_".$status;
	    #if($line[0] ne $gene_s){
	    #print STDERR " --- $line[0] -> $gene_s\n";#<STDIN>;
	    #}
	    
	    if($all_explained_freq{$gene_s}->[5] >= $explained_freq_threshold){
		$all_pvalue{$gene_s} = $line[2];
	    }
	}
    }
    close(FILE);

    #sort the pvalue and infer the corrected pvalue
    my $K = (keys %all_pvalue) + 0;
    %all_pvalue_corrected = ();
    foreach $gene (sort { $all_pvalue{$a} <=> $all_pvalue{$b} } keys %all_pvalue){
	$p_c = $all_pvalue{$gene};
	#if($p_c != 0){
	#print $p_c."\t".$K."\t".(0.05 / $K)."\n";<STDIN>;
	#}
	if ( ($FAST_CALL_FLAG && $p_c == 0) ||
	#if ( ($FAST_CALL_FLAG && $p_c <= 0.005) ||
	     (!$FAST_CALL_FLAG && $p_c <= (0.05 / $K))){
	    $K--;
	    $all_pvalue_corrected{$gene} = $p_c;
	}
	else{
	    last;
	}
    }
    
    #output the final result file
    $sample_explained_pvalue_file = "$explained_file_body\_pvalue.dat";
    #print STDERR $sample_explained_file."\n".$sample_explained_pvalue_file."\n";#<STDIN>;

    open(OUT, ">$sample_explained_pvalue_file");
    foreach $gene (keys %all_explained_freq){
	print OUT "".(join("\t", @{$all_explained_freq{$gene}}))."\t";
	if(exists $all_pvalue_corrected{$gene}){
	    print OUT $all_pvalue_corrected{$gene}."\n";
	}
	else{
	    print OUT "-"."\n";
	}
    }
    
}





