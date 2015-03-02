use warnings;
#use strict;

sub construct_driver_net_network{
    my ($index_to_gene, $gene_to_index, $connections, $script_dir) = @_;
    #print STDERR " *** BEGIN DRIVER_NET NETWORK CONSTRUCTION\n";
    $network = "$script_dir/network_FIsInGene_041709.txt";
    construct_network_pair($network, $index_to_gene, $gene_to_index, $connections);
    #print STDERR " *** END DRIVER_NET NETWORK CONSTRUCTION\n";
}


sub construct_netbox_network{
    my ($data_dir, $index_to_gene, $gene_to_index, $connections, $flag_update, $script_dir) = @_;
    $network_gene_list = "$script_dir/Homo_sapiens.gene_info";
    $network = "$script_dir/netbox.script";
    construct_gene_ID($network_gene_list, $index_to_gene, $gene_to_index);
    #to obtain the gene that belong to the network
    $updated_index_to_gene = construct_network($network, $index_to_gene, $gene_to_index, $connections);

    #print STDERR "size: ".(@{$index_to_gene})."\t".(@{$updated_index_to_gene})."\n";
    
    @{$index_to_gene} = @{$updated_index_to_gene};#simply copy the array
    
    #print STDERR "size: ".(@{$index_to_gene})."\t".(@{$updated_index_to_gene})."\n";

    if(defined $flag_update && $flag_update ne "NO_DATA_UPDATE"){
	update_gene_index_using_data($data_dir, $index_to_gene, $gene_to_index);
    }
}

sub construct_network_pair{
    my ($network_file, $index_to_gene, $gene_to_index, $connections) = @_;
    
    #print STDERR " *** construct_network_pair\n";

    open(FILE, $network_file);
    my $index = @{$index_to_gene}+0;
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	for($i = 0; $i < @line; $i++){
	    $gene = $line[$i];
	    #print STDERR " *** $i $gene\n";
	    
	    #a new gene
	    if(! exists $gene_to_index->{$gene}){
		$gene_to_index->{$gene} = $index;
		
		my @tab = ($gene);
		push(@{$index_to_gene}, \@tab);
		
		my @tab2 = ();
		push(@{$connections}, \@tab2);
		
		$index++;
	    }
	}
	
	$connection_exist = 0;
	$index1 = get_ID($line[0], $gene_to_index);
	$index2 = get_ID($line[1], $gene_to_index);
	
	#a new connection
	if(!belong($connections->[$index1], $index2)){
	    push(@{$connections->[$index1]}, $index2);
	    push(@{$connections->[$index2]}, $index1);
	}
	#print STDERR " *** $index1 $index2\t".@{$connections->[$index1]}."\t".@{$connections->[$index2]}."\n";<STDIN>;
    }
    #print STDERR " *** end construction\n";
}



sub construct_gene_ID{
    my ($network_gene_list, $index_to_gene, $gene_to_index) = @_;
    print STDERR " *** construct gene name table\n";
    #Read the gene file 2 times to avoid problem with synonyms !!!
    #1) compute the index of the main gene name
    my $index = 0;
    open(GEN, $network_gene_list);
    while (<GEN>){
	my $line = $_;
	chomp($line);
	my @parts = split (/\t/, $line);
	my $gene = $parts[2];
	$gene_to_index->{$gene} = $index;
	my @tab = ();
	$index_to_gene->[$index] = \@tab;
	push(@{$index_to_gene->[$index]}, $gene);
	$index++;
    }
    close(GEN);
    
    #2) Add the index of the synonys that do not have any index to avoid confusion as in the following example:
    #9606	26	ABP1	-	ABP|AOC1|DAO|DAO1|KAO
    #9606	1610	DAO	-	DAAO|DAMOX|MGC35381|OXDA
    open(GEN, $network_gene_list);
    while (<GEN>){
	my $line = $_;
	chomp($line);
	my @parts = split (/\t/, $line);
	my $gene = $parts[2];
	$index = $gene_to_index->{$gene};
	if (index($parts[4], "|") != -1){
	    my @synonyms = split(/\|/,$parts[4]);
	    foreach my $syn (@synonyms){
		if(! exists $gene_to_index->{$syn}){
		    $gene_to_index->{$syn} = $index ;
		    push(@{$index_to_gene->[$index]}, $syn);
		}
	    }
	}
    }
    close(GEN);
}

sub construct_network{
    my ($network, $index_to_gene, $gene_to_index, $connections) = @_;
    print STDERR " *** construct network \n";#<STDIN>;
    #  2. read network info into matrix
    #my @connections = ();
    read_network($connections, $index_to_gene, $gene_to_index, $network);
    #remove all the genes that do not belong to the network
    my @updated_index_to_gene = ();
    $updated_index = 0;
    for(my $i = 0; $i < @{$index_to_gene}; $i++){
	if(@{$connections->[$i]}){
	    $updated_index_to_gene[$updated_index] = $index_to_gene->[$i];
	    $updated_index++;
	}
	else{
	    for($j = 0; $j < @{$index_to_gene->[$i]}; $j++){
		delete $gene_to_index->{$index_to_gene->[$i]->[$j]};
	    }
	}
    }
    #@index_to_gene = ();
    
    $index_to_gene = \@updated_index_to_gene;
    for(my $i = 0; $i < @{$index_to_gene}; $i++){
	foreach $g (@{$index_to_gene->[$i]}){
	    if(get_name($i, $index_to_gene) eq "MMP9"){
		#print STDERR $g."\t".(get_ID($g))."\t".$i."\n";<STDIN>;
	    }
	    
	    #print STDERR $g."\t".$i."\n";<STDIN>;
	    $gene_to_index->{$g} = $i;
	}
    }
    #read the network a 2nd time to obtain construct the edges with the updated index
    #print STDERR "********** construct network 2\n";<STDIN>;
    undef @{$connections};
    read_network($connections, $index_to_gene, $gene_to_index, $network);

    return \@updated_index_to_gene;

}

sub get_ID{
    my ($gene, $gene_to_index) = @_;
    return $gene_to_index->{$gene};
}

sub get_name{
    my ($id, $index_to_gene) = @_;
    if(!exists $index_to_gene->[$id]){
	print ("Sorry there is NO gene-name for the ID $id \n");
	<STDIN>;#exit();
    }
    return $index_to_gene->[$id]->[0];
}


sub read_network{
    my ($connections, $index_to_gene, $gene_to_index, $net) = @_;
    for (my $i=0; $i < @{$index_to_gene}; $i++){
	my @tab = ();
	push(@{$connections}, \@tab);
    }
    open(NET, $net);

    my @line_parts;
    my @genes;
    while(<NET>){
	my $line = $_;
	chomp($line);
	
	if (index($line, "INSERT INTO INTERACTION VALUES") != -1){
	    @line_parts = split(/,/,$line);
	    @genes = split(/'/,$line_parts[2]);
	    $gene1 = $genes[1];
	    @genes = split(/'/,$line_parts[3]);
	    $gene2 = $genes[1];
	    $index1 = get_ID($gene1, $gene_to_index);
   	    $index2 = get_ID($gene2, $gene_to_index); 

	    if(!belong($connections->[$index1], $index2)){
		push(@{$connections->[$index1]}, $index2);
		push(@{$connections->[$index2]}, $index1);
	    }
	    else{
		#print STDERR "relation present multiple time $gene1 $gene2 ** $_";<STDIN>;
	    }
	}
    }
    close(NET);
}

#to remove weird conflict due to gene with synonymous name is the input data sets from the ahsh table index_to_gene
#only one of the synonymous gene (look code for tie break) is used during the analysis
sub update_gene_index_using_data{
    my ($data_dir, $index_to_gene, $gene_to_index) = @_;

    opendir(DIR, $data_dir);
    @the_DATA_DIR = readdir(DIR);
    close(DIR);
    
    
    #To track problem in gene naming, espessialy synonymous gene names in
    my %gene_ID_name_update = ();

    foreach my $dir_sample (@the_DATA_DIR){
	$mutation_file_name = "$data_dir/$dir_sample/Genelist_Status.txt";
	if(-e $mutation_file_name){
	    
	    $nb_sample++;
	    $nb_mutated_gene = 0;
	
	    open(FILE, "$mutation_file_name");
	    #print STDERR " *** read file $mutation_file_name\n";#<STDIN>;
	    #read the file to obtain the dysregulated and mutated genes
	    while(<FILE>){
		chop ($_);
		@line = split(/\t/, $_);
		my @parts = split(/_/,$line[0]);
		my $gene_name = $parts[0];
		my $status = $parts[1];
		if(exists $gene_to_index->{$gene_name}){#filter out all the gene_name that do not belong to the input network
		    my $gene_ID = get_ID($gene_name, $gene_to_index );
		    
		    #to make sure that the outputed gene name is the one present in the input files
		    $expected_gene_name = get_name($gene_ID, $index_to_gene);
		    if($expected_gene_name ne $gene_name){
			
			if(exists $gene_ID_name_update{$gene_ID}){
			    $used_gene_name = $gene_ID_name_update{$gene_ID};
			    $index_to_gene->[$gene_ID]->[0] = $used_gene_name;
			    
			    print STDERR " *** WARNING GENE $gene_ID WITH SYNONYMOUS NAMES DETECTED: $gene_name and ".$expected_gene_name. " -> ONLY DATA FROM GENE ".$used_gene_name." WILL BE USED/REPORTED DURING THE ANALYSIS\n";#<STDIN>;
			    #remove the gene for which we add a conflict
			    if($used_gene_name ne $gene_name){
				delete $gene_to_index->{$gene_name};
			    }
			    else{
				delete $gene_to_index->{$expected_gene_name};
			    }
			    #if(exists $gene_to_index->{$gene_name}){
			    #print STDERR " *** WEIRD !!!!!\n";<STDIN>;
			    #}
			    
			}
			
			else{

			    #save the 1st gene name of the synomym list will be used in case of further conflict
			    $gene_ID_name_update{$gene_ID} = $expected_gene_name;
			    #To output the gene name present in the input file
			    $index_to_gene->[$gene_ID]->[0] = $gene_name;
			    
			    print STDERR " *** SWAP GENE NAME ACCORDING TO DATA: $gene_ID: $expected_gene_name => $gene_name\n";#.($index_to_gene->[$gene_ID]->[0])."\n";#<STDIN>;	

			}
			
		    }
		    
		    else{
			$gene_ID_name_update{$gene_ID} = $expected_gene_name;
		    }
		    
		}
	    }
	    close(FILE);
	    #last if($nb_sample == $nb_sample_to_test);
	}
    }   
}

sub is_connected{
    my ($g1, $g2, $connections) = @_;
    return belong($connections->[$g1], $g2);
}

sub belong{
    my ($tab, $val) = @_;
    my $res = 0;
    my $comp = 0;
    while($comp < @{$tab} && !$res){
	$res = ($tab->[$comp] == $val);
	$comp++;
    }
    return $res;
}

return 1;
