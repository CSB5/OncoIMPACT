#!/usr/bin/perl
use warnings;

my ($data_dir, $network_type, $script_dir) = @ARGV;

require "$script_dir/Construct_network.pl";
my %gene_to_index;
my @index_to_gene;
my @connections;

if(@ARGV > 1){
    if($network_type eq "NETBOX"){
	construct_netbox_network($data_dir, \@index_to_gene, \%gene_to_index, \@connections, $script_dir);
    }
    if($network_type eq "DRIVER_NET"){
	construct_driver_net_network(\@index_to_gene, \%gene_to_index, \@connections, $script_dir);
    }

    print STDERR "------ nb gene in network ".(@index_to_gene)."\n";
}



opendir(DIR, $data_dir);
@the_DATA_DIR = readdir(DIR);
close(DIR);


%sample_gene_mutated = ();
%sample_gene_CNV = ();
%sample_gene_dysregulated = ();

my $nb_fold_change = 13;
my $inc_fold_change = 0.5;
my @fold_change_value = ();
for($i = 0; $i < $nb_fold_change; $i++){
    push(@fold_change_value, $i*$inc_fold_change);
}
foreach my $dir_sample (@the_DATA_DIR){
    $mutation_file_name = "$data_dir/$dir_sample/Genelist_Status.txt";
    if(-e $mutation_file_name){
	open(FILE, $mutation_file_name);
	my @fold_change_stats = ();
	
	for($i = 0; $i < $nb_fold_change; $i++){
	    push(@fold_change_stats, 0);
	}
	$sample_gene_dysregulated{$dir_sample} = 0;
	$sample_gene_mutated{$dir_sample} = 0;
	$sample_gene_CNV{$dir_sample} = 0;
	while(<FILE>){
	    #print $_;
	    chop $_;
	    @line = split(/\t/, $_);

	    @line = split(/\t/, $_);
	    my @parts = split(/_/,$line[0]);
	    my $gene_name = $parts[0];
	    my $status = $parts[1];

	    if(@ARGV == 1 || exists $gene_to_index{$gene_name}){
		
		if ($status eq "MUT" || $status eq "AMPL" || $status eq "DEL"){
		    if($status eq "MUT"){
			$sample_gene_mutated{$dir_sample}++;
		    }
		    else{
			$sample_gene_CNV{$dir_sample}++;
		    }
		}
		
		if ($status eq "UP" || $status eq "DOWN"){
		    $fold_change = abs($line[1]);
		    for($i = 0; $i < @fold_change_stats; $i++){
			if($fold_change != 0 && $fold_change >= $fold_change_value[$i]){
			    $fold_change_stats[$i]++;
			}
		    }
		}
	    }
	}
	$sample_gene_dysregulated{$dir_sample} = \@fold_change_stats;
    }
}

$header = "NAME";
for($i = 0; $i < $nb_fold_change; $i++){
    $h = "DYS_$fold_change_value[$i]";
    $h =~ s/\./_/;
    $header .= "\t$h";
}
$header .= "\tMUT\tCNV\n";
print $header;
foreach $s (keys %sample_gene_dysregulated){
    print $s."\t".join("\t", @{$sample_gene_dysregulated{$s}})."\t".$sample_gene_mutated{$s}."\t".$sample_gene_CNV{$s}."\n";
}



