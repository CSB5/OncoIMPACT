#!/usr/bin/perl
use warnings;

my ($test_param_dir, $network_type, $depth_th, $hub_th, $log2_fold_change_threshold, $flag_all_sample_use, $nb_random_sample, $nb_process, $script_dir) = @ARGV;


my $nb_max_line = 1000000;
my $nb_lines = 1;
my $num_file = 0;
my @all_cmds = ();


$cmd_file_prefix = "param_cmd";
for($i = 0; $i < $nb_random_sample; $i++){
    #random data
    $data_dir_name = "$test_param_dir/RANDOM_$i";
    push_call($data_dir_name);
    
    #real data
    if($flag_all_sample_use eq "ALL" && $i == 0){
	$data_dir_name = "$test_param_dir/REAL_$i/";
	push_call($data_dir_name);
    }
    }
close(OUT);

foreach $cmd_file(@all_cmds){
    $cmd = "cat $cmd_file | xargs -I cmd --max-procs=$nb_process bash -c cmd > /dev/null \n";
    print STDERR $cmd;
    print `$cmd`;
    #`rm $cmd`;
}

#compute the pvalues
my $nb_random_gene = 0;
my $nb_dysregulated_gene= 1;
my $nb_gene_in_network = 1;
my %nb_dysregulated_gene_sample = ();
my %nb_phenotype_gene_sample = ();
my %nb_phenotype_random_sample = ();
for($i = -1; $i < $nb_random_sample; $i++){
    $file = "$test_param_dir/RANDOM_$i/IMPACT/impact.dat.gz";
    if($i == -1){
	$file = "$test_param_dir/REAL_0/IMPACT/impact.dat.gz";
    }
    open(FILE, "gunzip -c $file |");
    print STDERR " *** $file\n";
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	$sample = $line[0];
	$gene = $line[1];
	if($sample eq "nb-dysregulated-gene"){
	    $nb_dysregulated_gene = $line[1];
	}
	if($sample eq "nb-gene-dysregulated-sample"){
	    $nb_dysregulated_gene_sample{$line[1]} = $line[2];
	}
	if($sample eq "nb-gene"){
	    $nb_gene_in_network = $line[1];
	    
	}
	if(index($gene, "_UP") != -1 || index($gene, "_DOWN") != -1){
	    $impact = $line[2];
	    $avg_impact = $line[3];
	    
	    #print STDERR " ooooooooooooo $gene $impact\n";
	    #initialisation for the real sample
	    if($i == -1){
		if(! exists $nb_phenotype_gene_sample{$sample}){
		    $nb_phenotype_gene_sample{$sample} = 0;
		    $nb_phenotype_random_sample{$sample} = 0;
		}
		$nb_phenotype_gene_sample{$sample}++;
		#
		if(! exists $sample_gene_impact{$sample}){
		    my %map = ();
		    $sample_gene_impact{$sample} = \%map;
		    #$nb_sample++;
		}
		my @tab = ($impact, $avg_impact, 0, 0, 0, 0);
		$sample_gene_impact{$sample}->{$gene} = \@tab;
		#
	    }
	    else{
		
		#update the pvalue 
		$nb_phenotype_random_sample{$sample}++;
		foreach $g (keys %{$sample_gene_impact{$sample}}){
		    $sample_gene_impact{$sample}->{$g}->[5]++;
		    $g_impact = $sample_gene_impact{$sample}->{$g}->[0];
		    $g_avg_impact = $sample_gene_impact{$sample}->{$g}->[1];
		    if( $g_impact < $impact){
			$sample_gene_impact{$sample}->{$g}->[2]++;
		    }
		    if($g_avg_impact < $avg_impact){
			$sample_gene_impact{$sample}->{$g}->[3]++;
		    }
		    if($g_impact < $impact && $g_avg_impact < $avg_impact){
			$sample_gene_impact{$sample}->{$g}->[4]++;
		    }
		}
	    }
	}
    }
    close(FILE);
}

open (OUT, ">$test_param_dir/pvalue_impact.dat");
foreach $sample (keys %sample_gene_impact){
    foreach $gene (keys %{$sample_gene_impact{$sample}}){
	#$nb_test = $nb_dysregulated_gene_sample{$sample} * $nb_random_sample;
	$nb_test = $nb_phenotype_random_sample{$sample};
	$nb_pheno_real = $nb_phenotype_gene_sample{$sample};
	#print " *** ".$sample."\t".$gene."\t".$nb_test."\t".$nb_pheno_real."\n";
	print OUT $sample."\t".$gene."\t".$sample_gene_impact{$sample}->{$gene}->[0]."\t".$sample_gene_impact{$sample}->{$gene}->[1];
	for(my $i = 2; $i <= 4; $i++){
	    print OUT "\t".($sample_gene_impact{$sample}->{$gene}->[$i]/($nb_test));
	}
	for(my $i = 2; $i <= 4; $i++){
	    print OUT "\t".($sample_gene_impact{$sample}->{$gene}->[$i]/($nb_test)*$nb_pheno_real);
	}
	print OUT "\n";
    }
}

close(OUT);

sub push_call{
    my ($data_dir) = @_;
    
    $dir_res = "$data_dir/IMPACT";
    `mkdir -p $dir_res` unless (-d $dir_res);
    
    $file_res = "$dir_res/impact.dat.gz";
    
    if(! -e $file_res){
	#print STDERR "mkdir $dir_res\n";
	#print STDERR " ********** $file_res\n";<STDIN>;
	$exe = "$script_dir/10_Cluster_Algo_impact.pl $data_dir $network_type $depth_th $hub_th 0 $log2_fold_change_threshold $dir_res $script_dir 2> /dev/null";
	
	if($num_file == 0 || $nb_lines > $nb_max_line){
	    close(OUT) if($num_file != 0);
	    $nb_lines = 1;
	    $num_file++;
	    $cmd_file = "$test_param_dir/$cmd_file_prefix\_$num_file.txt"; 
	    push(@all_cmds, $cmd_file);
	    open(OUT, ">$cmd_file");
	}
	$nb_lines++;
	print OUT $exe."\n";
    }
}
