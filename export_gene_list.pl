#!/usr/bin/perl
use warnings;

my ($data_dir, $module_file, $network_type, $fold_change_threshold, $large_data_set_driver_stats_file, $out_dir, $script_dir) = @ARGV;

require "$script_dir/Construct_network.pl";

#cancer gene census
my $cancer_gene_census_file = "$script_dir/cancer_gene_census.csv";
my %cancer_census = ();
open(FILE, $cancer_gene_census_file);
while(<FILE>){
    @line = split(/\s+/, $_);    
    #chomp($line[0]);
    $cancer_census{$line[0]} = 1;
}
close(FILE);

#pan cancer data
my $pan_cancer_file = "$script_dir/pancancer_driver_list.csv";
my %pan_cancer = ();
open(FILE, $pan_cancer_file);
while(<FILE>){
    chop($_);
    @line = split(/\t/, $_);
    #chomp($line[0]);
    $pan_cancer{$line[0]} = 1 if($line[7] eq "High_Confidence_Driver");
}
close(FILE);

#Read the large data set info in case of the use of a data base
my %large_data_set_driver_stats = ();
if($large_data_set_driver_stats_file ne "NONE"){
    open(FILE, $large_data_set_driver_stats_file);
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	my @stats = ($line[14], $line[6], $line[1]);#IMPACT DRIVER_FREQUENCY  MUTATION_FREQUENCY
	$large_data_set_driver_stats{$line[0]} = \@stats;
    }
    close(FILE);
}



my %gene_to_index;
my @index_to_gene;
my @connections;

if($network_type eq "NETBOX"){
    #update not done to save computational time
    construct_netbox_network($data_dir, \@index_to_gene, \%gene_to_index, \@connections, $script_dir);
}
if($network_type eq "DRIVER_NET"){
    construct_driver_net_network(\@index_to_gene, \%gene_to_index, \@connections, $script_dir);
}


my %sample_gene_mutated = ();
my %sample_gene_dysregulated = ();

my %sample_gene_driver = ();
my %sample_gene_phenotype = ();

my %sample_gene_expression = ();

opendir(DIR, $data_dir);
@the_DATA_DIR = readdir(DIR);
close(DIR);


my %naive_mut_list = ();
my %naive_dys_list = ();

my %infer_mut_list = ();
my %infer_dys_list = ();
my %infer_linker_list = ();

@dys_type_order = ("ALL", "DOWN", "UP");
@mut_type_order = ("ALL", "MUT", "DEL", "AMPL", "BOTH");

$nb_samples = 0;
my %gene_impact;
my %gene_impact_filtered;
my %sample_gene_impact;
my %sample_gene_impact_filtered;

#Produce the naive gene lists 
foreach my $dir_sample (@the_DATA_DIR){
    $mutation_file_name = "$data_dir/$dir_sample/Genelist_Status.txt";
    if(-e $mutation_file_name){
	
	my %map = ();
	$sample_gene_dysregulated{$dir_sample} = \%map;
	
	my %map2 = ();
	$sample_gene_expression{$dir_sample} = \%map2;
	
	#print STDERR " *** ".$mutation_file_name."\n";#<STDIN>;
	$nb_samples++;
	open(FILE, "$mutation_file_name");
	while(<FILE>){
	    chop ($_);
	    @line = split(/\t/, $_);
	    my @parts = split(/_/,$line[0]);
	    my $gene_name = $parts[0];
	    my $status = $parts[1];
	    #initaillaze gene impact sample
	    my %map1;
	    $sample_gene_impact{$dir_sample} = \%map1;
	    my %map2;
	    $sample_gene_impact_filtered{$dir_sample} = \%map2;

	    if(exists $gene_to_index{$gene_name}){#filter out all the gene_name that do not belong to the input network
		
		$gene_ID = get_ID($gene_name, \%gene_to_index );
		
		#gene impact initalization
		if(! exists $gene_impact{$gene_ID}){
		    my %map = ("MUT", 0, "UP", 0, "DOWN", 0);
		    $gene_impact{$gene_ID} = \%map;
		    #
		    my %map1 = ("MUT", 0, "UP", 0, "DOWN", 0);
		    $gene_impact_filtered{$gene_ID} = \%map1;

		    my %map2 = ("MUT", 0, "UP", 0, "DOWN", 0);
		    $sample_gene_impact_filtered{$dir_sample}->{$gene_ID} = \%map2;
		    #
		    my %map3 = ("MUT", 0, "UP", 0, "DOWN", 0);
		    $sample_gene_impact_filtered{$dir_sample}->{$gene_ID} = \%map3;
		}

		if (($status eq "UP" || $status eq "DOWN")){
		    
		    if(! exists $infer_linker_list{$gene_ID}){
			my %map = ("ALL", 0, "UP", 0,"DOWN", 0);
			$infer_linker_list{$gene_ID} = \%map;
		    }
		    
		    $g_exp = $line[1];
		    $sample_gene_expression{$dir_sample}->{$gene_ID} = $g_exp;
		    if(abs($g_exp) >= $fold_change_threshold){
			
			if(! exists $sample_gene_dysregulated{$dir_sample}->{$gene_ID}){
			    my %map = ("UP", 0,"DOWN", 0);
			    $sample_gene_dysregulated{$dir_sample}->{$gene_ID} = \%map;
			    
			    my %map2 = ("UP", 0,"DOWN", 0);
			    $sample_gene_phenotype{$dir_sample}->{$gene_ID} = \%map2;

			}
			
			if(! exists $naive_dys_list{$gene_ID}){
			    my %map = ("ALL", 0, "UP", 0,"DOWN", 0);
			    $naive_dys_list{$gene_ID} = \%map;
			    
			    my %map2 = ("ALL", 0, "UP", 0,"DOWN", 0);
			    $infer_dys_list{$gene_ID} = \%map2;
			    
			}
			
			$sample_gene_dysregulated{$dir_sample}->{$gene_ID}->{$status} = 1;
			$naive_dys_list{$gene_ID}->{$status}++;
			$naive_dys_list{$gene_ID}->{"ALL"}++;
		    }
		       
		    
		    
		}
	    }
	}
	close(FILE);
	

	#Read for alteration, allows to filter alteration that do not correlate with dyregulated genes
	open(FILE, "$mutation_file_name");
	while(<FILE>){
	    chop ($_);
	    @line = split(/\t/, $_);
	    my @parts = split(/_/,$line[0]);
	    my $gene_name = $parts[0];
	    my $status = $parts[1];
	    
	    if(exists $gene_to_index{$gene_name}){#filter out all the gene_name that do not belong to the input network
		$gene_ID = get_ID($gene_name, \%gene_to_index );
		
		if ($status eq "MUT" ||
		    ($status eq "AMPL" && (1 || $sample_gene_dysregulated{$dir_sample}->{$gene_name}->{"UP"} > 0)) ||
		    ($status eq "DEL" && (1 || $sample_gene_dysregulated{$dir_sample}->{$gene_name}->{"DOWN"} < 0))){

			if(! exists $sample_gene_mutated{$dir_sample}->{$gene_ID}){
			    my %map = ("MUT", 0, "AMPL", 0, "DEL", 0, "BOTH", 0);
			    $sample_gene_mutated{$dir_sample}->{$gene_ID} =\%map;
			    
			    my %map2 = ("MUT", 0, "AMPL", 0, "DEL", 0, "BOTH", 0);
			    $sample_gene_driver{$dir_sample}->{$gene_ID} =\%map2;
			}
			
			if(! exists $naive_mut_list{$gene_ID}){
			    my %map = ("ALL", 0, "MUT", 0, "AMPL", 0, "DEL", 0, "BOTH", 0);
			    $naive_mut_list{$gene_ID} =\%map;

			    my %map2 = ("ALL", 0, "MUT", 0, "AMPL", 0, "DEL", 0, "BOTH", 0);
			    $infer_mut_list{$gene_ID} = \%map2;
			}
		
			$naive_mut_list{$gene_ID}->{$status}++;
			$sample_gene_mutated{$dir_sample}->{$gene_ID}->{$status} = 1;
			
			#CNV and SNV
			if($sample_gene_mutated{$dir_sample}->{$gene_ID}->{"MUT"} && 
			   ($sample_gene_mutated{$dir_sample}->{$gene_ID}->{"AMPL"} || $sample_gene_mutated{$dir_sample}->{$gene_ID}->{"DEL"})){
			    $sample_gene_mutated{$dir_sample}->{$gene_ID}->{"BOTH"} = 1;
			    $naive_mut_list{$gene_ID}->{"BOTH"}++;
			}
			else{
			    $naive_mut_list{$gene_ID}->{"ALL"}++;
			}
		}
	    }
	}
	close(FILE);
    }
}

#Produce the driver/phenotype gene lists 
#Read the results file
open(FILE, $module_file);
my @dys_gene;
my @mut_gene;
my @linker_gene;
my @all_mod_gene;
my @all_mod_gene_type;
`rm -r $out_dir/MODULE_LIST`;
`mkdir $out_dir/MODULE_LIST`;

my %sample_driver = ();
my %sample_pheno = ();
foreach $s (keys %sample_gene_mutated){
    my %map = ();
    $sample_pheno{$s} = \%map;
    my %map1 = ();
    $sample_driver{$s} = \%map1;
}
#open all the file

#######################
#Read the module file to obtain the sample module list
#######################
my %sample_gene_module = ();
my %sample_gene_driver_pheno = ();
my %all_gene_in_module = ();
my @module_order = ();

while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $module_ID = $line[0];
    push(@module_order, $module_ID);
    
    @temp = split(/\./, $module_ID);
    $sample_name = $temp[0];
    
    

    @dys_gene = split(/\;/, $line[2]);
    @mut_gene = split(/\;/, $line[1]);
    @linker_gene = ();
    if($line[3] ne "-;"){
	@linker_gene = split(/\;/, $line[3]);
    }
    
    #$sample_data_dir = "$out_dir/SAMPLE_SPECIFIC_DATA/$sample_name";
    #`mkdir $sample_data_dir

    @all_mod_gene = @mut_gene;
    @all_mod_gene_type = @mut_gene;
    $module_impact_filtered = $line[6];
    $module_impact = $line[5];
    
    #print STDERR " --- size = ".(@all_mod_gene+0)."\n";#<STDIN>;
    
    foreach $g (@mut_gene){
	
	#print STDERR " *** MUT"."\t".$sample_name."\t".$g."\n";
	
	$gene_ID = get_ID($g, \%gene_to_index );
	$sample_gene_impact{$sample_name}->{$gene_ID}->{"MUT"} += $module_impact;
	$sample_gene_impact_filtered{$sample_name}->{$gene_ID}->{"MUT"} += $module_impact_filtered;
	#
	$gene_impact{$gene_ID}->{"MUT"} += $module_impact;
	$gene_impact_filtered{$gene_ID}->{"MUT"} += $module_impact_filtered;

	if(exists $sample_gene_mutated{$sample_name}->{$gene_ID}){

	    foreach $type (keys %{$sample_gene_mutated{$sample_name}->{$gene_ID}}){
		if($sample_gene_mutated{$sample_name}->{$gene_ID}->{$type} != 0){
		    $infer_mut_list{$gene_ID}->{$type}++;
		    $infer_mut_list{$gene_ID}->{"ALL"}++;

		    if(!exists $sample_driver{$sample_name}->{$gene_ID}){
			$sample_driver{$sample_name}->{$gene_ID} = $type;
		    }
		    else{
			$sample_driver{$sample_name}->{$gene_ID} .= "_".$type if($type ne "BOTH");
		    }
		}
	    }
	    $infer_mut_list{$gene_ID}->{"ALL"} -= 2 if($sample_gene_mutated{$sample_name}->{$gene_ID}->{"BOTH"} != 0);
	}
	else{
	    print STDERR " *** $g DRIVER GENE PREDICTED in $sample_name BUT NOT MUTATED -_-'\n";
	}
    }

    foreach $g (@dys_gene){
	
	my @parts = split(/_/,$g);
	my $gene_name = $parts[0];
	push(@all_mod_gene, $gene_name);
	push(@all_mod_gene_type, $g);
	
	$gene_ID = get_ID($gene_name, \%gene_to_index );
	#print STDERR " *** DYS"."\t".$sample_name."\t".$gene_name."\n";
	
	if(exists $sample_gene_dysregulated{$sample_name}->{$gene_ID}){
	    foreach $type (keys %{$sample_gene_dysregulated{$sample_name}->{$gene_ID}}){
		if($sample_gene_dysregulated{$sample_name}->{$gene_ID}->{$type} != 0){
		    $infer_dys_list{$gene_ID}->{$type}++ ;
		    $infer_dys_list{$gene_ID}->{"ALL"}++ ;
		    
		    $sample_pheno{$sample_name}->{$gene_ID} = $type;
		    
		    $sample_gene_impact{$sample_name}->{$gene_ID}->{$type} += $module_impact;
		    $sample_gene_impact_filtered{$sample_name}->{$gene_ID}->{$type} += $module_impact_filtered;
		    #
		    $gene_impact{$gene_ID}->{$type} += $module_impact;
		    $gene_impact_filtered{$gene_ID}->{$type} += $module_impact_filtered;
		}

	    }
	}
	else{
	    print STDERR " *** $gene_name PHENOTYPE GENE PREDICTED in $sample_name BUT NOT DYSREGULATED -_-'\n";
	}
    }

    #print STDERR " --- size = ".(@all_mod_gene+0)."\n";#<STDIN>;


    if(@linker_gene != 0){
	foreach $g (@linker_gene){
	    my @parts = split(/_/,$g);
	    my $gene_name = $parts[0];
	    push(@all_mod_gene, $gene_name);
	    push(@all_mod_gene_type, $g);
	    $gene_ID = get_ID($gene_name, \%gene_to_index);
	    if(exists $sample_gene_dysregulated{$sample_name}->{$gene_ID}){
		$infer_linker_list{$gene_ID}->{"ALL"}++ ;
	    }
	    else{
		print STDERR " *** $g |$gene_name| LINKER GENE PREDICTED in $sample_name BUT NOT DYSREGULATED -_-'\n"; 
	    }
	}
    }
    
    #print STDERR " --- size = ".(@all_mod_gene+0)."\n";<STDIN>;


    #out a gene list for each module to be prossess by DAVID
    open(OUT_M, ">$out_dir/MODULE_LIST/$module_ID.mod");
    print OUT_M "".(join("\n", @all_mod_gene))."\n";
    close(OUT_M);

    #Update the data structure to output the sample gene in module matrix
    #the expression type of dys regulated gene is taken into account, the type of mutation not
    if(! exists $sample_gene_module{$sample_name}){
	my %map = ();
	$sample_gene_module{$sample_name} = \%map;
	my %map2;
	$sample_gene_driver_pheno{$sample_name} = \%map2;
    }
    foreach my $gene (@all_mod_gene_type){
	$sample_gene_module{$sample_name}->{$gene} = 1;
	if(! exists $all_gene_in_module{$gene}){
	    my %map = ();
	    $all_gene_in_module{$gene} = \%map;
	}
	$all_gene_in_module{$gene}->{$module_ID} = 1;
    }
}
close(FILE);

#divide the gene impcat by the number of samples
foreach $gene_ID (%gene_impact){
    foreach $type (keys %{$gene_impact{$gene_ID}}){
	#print STDERR $gene_impact{$gene_ID}->{$type}."\t".$nb_samples."\n";<STDIN>;
	$gene_impact{$gene_ID}->{$type} = $gene_impact{$gene_ID}->{$type}/$nb_samples;
	$gene_impact_filtered{$gene_ID}->{$type} = $gene_impact_filtered{$gene_ID}->{$type}/$nb_samples;
    }
}


#to write the sample specific driver phenotye gene file
my $sample_out_dir = "$out_dir/SAMPLE_SPECIFIC_DATA";
`rm -r $sample_out_dir`;
`mkdir $sample_out_dir`;

#my $SAMPLE_DATA_HEADER = "#GENE\tTYPE\tSAMPLE_IMPACT\tDATA_SET_IMPACT\tDRIVER_FREQUENCY\tMUTATION_FREQUENCY\tCANCER_CENSUS\tPAN_CANCER";
my $SAMPLE_DATA_HEADER = "#GENE\tTYPE\tSAMPLE_IMPACT\tDATA_SET_IMPACT\tDRIVER_FREQUENCY\tMUTATION_FREQUENCY\tCANCER_CENSUS\tPAN_CANCER";
if($large_data_set_driver_stats_file ne "NONE"){
    $SAMPLE_DATA_HEADER = "#GENE\tTYPE\tSAMPLE_IMPACT\tDATA_BASE_IMPACT\tDATA_BASE_DRIVER_FREQUENCY\tDATA_BASE_MUTATION_FREQUENCY\tCANCER_CENSUS\tPAN_CANCER";
}


#read the annotation file
my %all_annotation = ();
#read_annotation_file();

foreach $s (keys %sample_driver){
    

    #for driver
    open(OUT, ">$sample_out_dir/temp.dat");
    foreach $g (keys %{$sample_driver{$s}}){
	$g_name = get_name($g, \@index_to_gene);
	$str = $g_name."\t".
	    #(write_coord($g_name))."\t".
	    $sample_driver{$s}->{$g}."\t".
	    $sample_gene_impact_filtered{$s}->{$g}->{"MUT"}."\t";
	
	if($large_data_set_driver_stats_file eq "NONE"){
	    $str .= $gene_impact_filtered{$g}->{"MUT"}."\t".(sprintf("%.3f", $infer_mut_list{$g}->{"ALL"}/$nb_samples))."\t".(sprintf("%.3f", $naive_mut_list{$g}->{"ALL"}/$nb_samples))."\t";
	}
	else{
	    if(exists $large_data_set_driver_stats{$g_name}){
		$str .= join("\t", @{$large_data_set_driver_stats{$g_name}})."\t";
	    }
	    else{
		$str .= "0\t0\t0\t";
	    }
	}
	
	$str .= (get_gene_annotation_status($g_name))."\n";
	
	print OUT $str;
	
    }
    close(OUT);
    
    #sorting based on impact
    #in case of same value (possible if in same module) order according to average gene impact in the whole samples
    
    $file = "$sample_out_dir/$s.txt";
    `echo -e \"$SAMPLE_DATA_HEADER\" > $file`;
    `sort -k3,3 -nr -k4,4 -nr $sample_out_dir/temp.dat >> $file`;
    
    #$file = "$sample_out_dir/driver_$s.dat";
    #`echo -e \"$SAMPLE_DATA_HEADER\" > $file`;
    #`sort -k4,4 -nr $sample_out_dir/temp.dat >> $file`;
    

    #for phenotype
    #open(OUT, ">$sample_out_dir/temp.dat");
    #foreach $g (keys %{$sample_pheno{$s}}){
	#$g_name = get_name($g, \@index_to_gene);
	#$type =  $sample_pheno{$s}->{$g};
	#print OUT $g_name."\t".
	    #(write_coord($g_name))."\t".
	#$type."\t".(sprintf("%.3f", $infer_dys_list{$g}->{$type}/$nb_samples))."\t".(sprintf("%.3f", $naive_dys_list{$g}->{$type}/$nb_samples))."\t".(get_gene_annotation_status($g_name))."\n";	
    #}
    #close(OUT);
    #$file = "$sample_out_dir/phenotype_$s.dat";
    #`echo -e \"$SAMPLE_DATA_HEADER\" > $file`;
    #`sort -k4,4 -nr $sample_out_dir/temp.dat >> $file`;
}
`rm $sample_out_dir/temp.dat`;


#Write the output files
my $naive_data;
my $infer_data;
my $out_file_name;


#Frequency gene list
for($i = 0; $i < 4; $i++){
    if($i <= 1){
	$naive_data = \%naive_mut_list;
	$infer_data = \%infer_mut_list;
	$type_order = \@mut_type_order;
	$out_file_name = "ALTERATION";
	$out_file_name = "NUM_".$out_file_name if($i == 1); 
    }
    else{
	$naive_data = \%naive_dys_list;
	$infer_data = \%infer_dys_list;
	$type_order = \@dys_type_order;
	$out_file_name = "DYSREGULATION";
	$out_file_name = "NUM_".$out_file_name if($i == 3); 
    }
    
    $out_file_name .= ".dat";
    open(OUT, ">$out_dir/$out_file_name");
    
    foreach $g (keys %{$naive_data}){
	$g_name = get_name($g, \@index_to_gene);
	$res = $g_name;
	
	foreach $type (@{$type_order}){
	    if($i == 0 || $i == 2){
		$res .= "\t".(sprintf("%.3f", $naive_data->{$g}->{$type}/$nb_samples));
	    }
	    else{
		$res .= "\t".$naive_data->{$g}->{$type};
	    }
	}

	foreach $type (@{$type_order}){
	    if($i == 0 || $i == 2){
		$res .= "\t".(sprintf("%.3f", ($infer_data->{$g}->{$type}/$nb_samples)));
	    }
	    else{
		$res .= "\t".$infer_data->{$g}->{$type};
	    }
	}
	
	$res .= get_gene_annotation_status($g_name);
	
	#the impact
	if($i <= 1){
	    $res .= "\t".$gene_impact{$g}->{"MUT"}."\t".$gene_impact_filtered{$g}->{"MUT"};
	}
	else{
	   $res .= "\t".($gene_impact{$g}->{"DOWN"}+$gene_impact{$g}->{"UP"})."\t".$gene_impact{$g}->{"DOWN"}."\t".$gene_impact{$g}->{"UP"}; 
	   $res .= "\t".($gene_impact_filtered{$g}->{"DOWN"}+$gene_impact_filtered{$g}->{"UP"})."\t".$gene_impact_filtered{$g}->{"DOWN"}."\t".$gene_impact_filtered{$g}->{"UP"}; 
	}
	print OUT $res."\n";
	
    }
    
    close(OUT);
}



open(OUT_ALL_EXPR, ">$out_dir/ALL_EXPR_MATRIX");
open(OUT_EXPR, ">$out_dir/EXPR_MATRIX");
open(OUT_PHENO, ">$out_dir/PHENO_MATRIX");
#all alterated genes are reported 
open(OUT_ALL_ALT, ">$out_dir/ALL_ALT_MATRIX");
#only driver altered genes are reported 
open(OUT_ALT, ">$out_dir/ALT_MATRIX");
#only patient specific drivers are reported
open(OUT_DRIVER, ">$out_dir/DRIVER_MATRIX");
#open(OUT_PHENO_ALT, ">$out_dir/PHENO_ALT_MATRIX");
open(OUT_MODULE, ">$out_dir/MODULE_MATRIX");
open(OUT_SAMPLE_MODULE, ">$out_dir/SAMPLE_MODULE_MATRIX");
#HEADER
my @expr = ();
my @sample_order = ();
foreach $sample (keys %sample_gene_expression){
    if(keys %{$sample_pheno{$sample}} != 0){
	push(@sample_order, $sample);
    }
    else{
	print STDERR " *** WARNING SAMPLE WITHOUT PHENOTYPE $sample\n";
    }
}
print OUT_ALL_EXPR (join("\t", @sample_order))."\n";
print OUT_EXPR (join("\t", @sample_order))."\n";
print OUT_PHENO (join("\t", @sample_order))."\n";
#print OUT_PHENO_ALT (join("\t", @sample_order))."\n";
print OUT_ALL_ALT (join("\t", @sample_order))."\n";
print OUT_ALT (join("\t", @sample_order))."\n";
print OUT_DRIVER (join("\t", @sample_order))."\n";
print OUT_SAMPLE_MODULE (join("\t", @sample_order))."\n";
print OUT_MODULE (join("\t", @module_order))."\n";

my $matrix_abscence_value = 1;
my $matrix_precense_value = 2;
#expression value matrix for heat map
foreach $gene (keys %naive_dys_list){
    #For the matrix of all the genes
    @expr = ();
    foreach $sample (@sample_order){
	$e = 0;
	if(exists $sample_gene_expression{$sample}->{$gene}){
	    $e = $sample_gene_expression{$sample}->{$gene};
	}
	push(@expr, $e);
    }
    print OUT_ALL_EXPR (join("\t", @expr))."\n";
    
    #FOR THE INFER GENE
    if($infer_dys_list{$gene}->{"ALL"} != 0 
       	){
	#The expression based list
	@expr = ();
	foreach $sample (@sample_order){
	    $e = 0;
	    if(exists $sample_gene_expression{$sample}->{$gene}){
		$e = $sample_gene_expression{$sample}->{$gene};
	    }
	    push(@expr, $e);
	}
	print OUT_EXPR (join("\t", @expr))."\n";
	
	#The PHENO BASED list
	for($i = 1; $i < @dys_type_order; $i++){
	    $type = $dys_type_order[$i];
	    
	    if($infer_dys_list{$gene}->{$type} != 0){
		@expr = ();#$sample;
		foreach $sample (@sample_order){
		    #print STDERR $gene."\t".$type."\n";<STDIN>;
		    $e = $matrix_abscence_value;
		    
		    if(exists $sample_pheno{$sample}->{$gene} &&
		       $sample_pheno{$sample}->{$gene} eq $type){
			$e = $matrix_precense_value;
		    }
		    push(@expr, $e);
		}
		print OUT_PHENO (join("\t", @expr))."\n";
		#print OUT_PHENO_ALT (join("\t", @expr))."\n";
	    }
	}
    }
}
close(OUT_ALL_EXPR);
close(OUT_EXPR);

foreach $gene (sort {$gene_impact_filtered{$b}->{"MUT"} <=> $gene_impact_filtered{$a}->{"MUT"}} keys %naive_mut_list){
    #print STDERR " *** $gene ".$gene_impact_filtered{$gene}->{"MUT"}."\n";<STDIN>;
    @expr_all = ();
    @expr = ();
    @expr_driver = ();
    foreach $sample (@sample_order){
	#Altered
	$e = $matrix_abscence_value;
	if(exists $sample_gene_mutated{$sample}->{$gene}){
	    $e = $matrix_precense_value;
	}
	push(@expr_all, $e);
	
	if($infer_mut_list{$gene}->{"ALL"} != 0){
	    #Driver 
	    $e = $matrix_abscence_value;
	    #print STDERR $gene."\t".$type."\n";<STDIN>;
	    if(exists $sample_driver{$sample}->{$gene}){
		$e = $matrix_precense_value;
	    }
	    push(@expr_driver, $e);
	    
	    #Altered but driver in one sample
	    $e = $matrix_abscence_value;
	    if(exists $sample_gene_mutated{$sample}->{$gene}){
		$e = $matrix_precense_value;
	    }
	    push(@expr, $e);
	}
    }
    
    if($infer_mut_list{$gene}->{"ALL"} != 0){
	print OUT_DRIVER (join("\t", @expr_driver))."\n";
	print OUT_ALT "".get_name($gene, \@index_to_gene)."\t".(join("\t", @expr))."\n";
    }
    print OUT_ALL_ALT  "".get_name($gene, \@index_to_gene)."\t".(join("\t", @expr_all))."\n";
    #print OUT_PHENO_ALT (join("\t", @expr))."\n";
}

close(OUT_ALT);
close(OUT_DRIVER);
#close(OUT_PHENO_ALT);


#module value matrix for heat map
foreach $gene (keys %all_gene_in_module){
    @expr = ();
    #all module of a sample or on the same line
    foreach $sample (@sample_order){
	$e = $matrix_abscence_value;
	if(exists $sample_gene_module{$sample}->{$gene}){
	    $e = $matrix_precense_value;
	}
	push(@expr, $e);
    }
    print OUT_SAMPLE_MODULE (join("\t", @expr))."\n";
    if(@expr+0 != @sample_order){
	print STDERR " *** MODULE_MATRIX WRONG SAMPLE NUMBER !!!!\n";<STDIN>;
    }
    
    #all module are on separate lines
    @expr = ();
    foreach $module_ID (@module_order){
	$e = $matrix_abscence_value;
	if(exists $all_gene_in_module{$gene}->{$module_ID}){
	   $e = $matrix_precense_value;
	}
	push(@expr, $e);
    }
    print OUT_MODULE (join("\t", @expr))."\n";
    
}
close(OUT_SAMPLE_MODULE);
close(OUT_MODULE);

run_exe("$out_dir/MATRIX_DIR");
run_exe("mv $out_dir/*MATRIX $out_dir/MATRIX_DIR");


sub get_gene_annotation_status{
    my ($g_name) = @_;
    my $res = "";
    if(exists $cancer_census{$g_name}){
	$res .= "\tCC";
    }
    else{$res .= "\t-";}

    if(exists $pan_cancer{$g_name}){
	$res .= "\tPC";
    }
    else{$res .= "\t-";}
	
    return $res;
}

#S100A10


sub write_coord{
    my ($gene_name) = @_;
    $str = "-";
    if( exists $all_annotation{$gene_name} ){
	$str = $all_annotation{$gene_name}->{"CHROM"}.":".$all_annotation{$gene_name}->{"TX_S"}."-".$all_annotation{$gene_name}->{"TX_E"};
    }
    #print STDERR $gene_name."\t".$str."\n";<STDIN>;
    return $str;
}



sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";
    print STDERR `$exe` if($run);
}
