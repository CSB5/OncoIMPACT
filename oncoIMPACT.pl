#!/usr/bin/perl
use warnings;

my ( $configFile, %config, $flag_debug );

my $help_message = "
This script prepares the data for oncoIMPACT run before launching oncoIMPACT.

Usage:
	oncoIMPACT.pl <config file> <optional: 1 to enable debug mode>

Version:
	0.9.3

Author:
	Burton Chia - chiakhb\@gis.a-star.edu.sg
	Denis Bertrandd - bertrandd\@gis.a-star.edu.sg\n";

if ( @ARGV == 0 ) {
	print $help_message;
	exit 0;
}


( $configFile, $flag_debug ) = @ARGV;

$flag_debug = 0 if(! defined $flag_debug);

#print STDERR " *** flag_debug $flag_debug\n";

# Sanity check on user provided parameters
unless(-s $configFile){
	print STDERR "Aborting! Config file does not exist or is empty. Please check the config file and try again.\n";
	exit 3;
}
if($flag_debug != 0 && $flag_debug != 1){
	print STDERR "Aborting! Debug flag contains an invalid option. Please check the parameter provided and try again.\n";
	exit 1;
}


# Check that dependent system programs are present
#print STDERR "[Dependencies] Checking the presence and version of required system programs\n" if $flag_debug;
#my $programPath;
# xargs
#chomp($programPath = `command -v xargs`);
#if($programPath eq ""){
#	print STDERR "Aborting! System command 'xargs' not found! Please ensure you are running this programme in a suitable environment.\n";#
#	exit 2;
#} elsif($flag_debug){
#	print STDERR "[xargs] Path: $programPath\n";
#	print STDERR "[xargs] Version:" . `xargs --version | head -n 1`;
#}


#Initialize the default parameters
$config{'testMode'} = 0;
$config{'dataType'} = "ARRAY";
# Read config file
print "Reading config file. Please wait...";
read_config( $configFile, \%config );
print "done.\n";


# Prep output directory
system("mkdir $config{'outDir'}") unless ( -d $config{'outDir'} );


# Check validity of parameters in config file
unless(-e $config{'outDir'} && -w $config{'outDir'}){
	print STDERR "Aborting! Output directory does not exist or is unwritable. Please check the config file and try again.\n";
	exit 1;
}
unless(-s $config{'scriptDir'}){
	print STDERR "Aborting! Scripts directory does not exist or is empty. Please check the config file and try again.\n";
	exit 1;
}
unless ($config{'numThreads'} =~ /^\d+?$/) {
    print STDERR "Aborting! numThreads is not given as an integer. Please check the config file and try again.\n";
    exit 1;
}
unless(-s $config{'cnv'}){
	print STDERR "Aborting! CNV file does not exist or is empty. Please check the config file and try again.\n";
	exit 1;
}
unless(-s $config{'exp'}){
	print STDERR "Aborting! Expression file does not exist or is empty. Please check the config file and try again.\n";
	exit 1;
}
unless(-s $config{'snp'}){
	print STDERR "Aborting! SNP file does not exist or is empty. Please check the config file and try again.\n";
	exit 1;
}
#
unless($config{'testMode'} == 0 || $config{'testMode'} == 1){
	print STDERR "Aborting! testMode option is not a valid flag - valid options are '0' or '1'. Please check the config file and try again.\n";
	exit 1;
}
unless($config{'dataType'} eq "ARRAY" || $config{'dataType'} eq "RNA_SEQ"){
	print STDERR "Aborting! dataType option is not a valid - valid options are 'ARRAY' or 'RNA_SEQ'. Please check the config file and try again.\n";
	exit 1;
}

# Prep data
unless ( -s $config{'outDir'} . "/COMPLETE_SAMPLES" ) {
	print "Preparing CNV data. Please wait...";
	prep_cnv($config{'cnv'}, $config{'outDir'});
	print "done.\n";

	print "Preparing SNP data. Please wait...";
	prep_snp($config{'snp'}, $config{'outDir'});
	print "done.\n";

	print "Preparing Expression data. Please wait...";
	prep_exp($config{'exp'}, $config{'outDir'});
	print "done.\n";

	print "Merging and cleaning output directory. Please wait...";
	merge_and_clean($config{'outDir'}, "COMPLETE_SAMPLES", "INCOMPLETE_SAMPLES");
	print "done.\n";
}


# Run oncoIMPACT
if(exists $config{'dataBase'}){
    print "\nImport data base $config{'dataBase'} Please wait...";
    import_data_base();
    print "done.\n";
    #
    print "\nRunning oncoIMPACT discovery mode using the data base: $config{'dataBase'}. Please wait...";
    run_oncoIMPACT_discovery();
}
else{
    print "\nRunning oncoIMPACT. Please wait...";
    run_oncoIMPACT();
    export_data_base() if(exists $config{'databaseExport'});
}

# system("rm -r $config{'outDir'}/COMPLETE_SAMPLES $config{'outDir'}/INCOMPLETE_SAMPLES");
print "done.\n";

#Output the results
output_final_result();


### Sub-routines ###
sub import_data_base{
    my $data_base_dir = "$config{'outDir'}/TMP_DATABASE";
    if(! -d $data_base_dir){
	#Uncompress the data
	system("mkdir $data_base_dir");
	system("tar -C $data_base_dir -zxvf  $config{'dataBase'} &> /dev/null");
	#print STDERR " *** Decompression done\n";<STDIN>;
	#
	$config{'dataBaseDir'}=$data_base_dir;
	#
	#Construct the directory structure for the data set
	prep_cnv("$data_base_dir/CNV.dat", $data_base_dir);
	#print STDERR " *** CNV Done\n";<STDIN>;
	prep_snp("$data_base_dir/SNP.dat", $data_base_dir);
	prep_exp("$data_base_dir/EXPR.dat", $data_base_dir);
	merge_and_clean("$data_base_dir/", "SAMPLE_INFO", "INCOMPLETE_SAMPLES");
	#
	system("rm -r $data_base_dir/INCOMPLETE_SAMPLES")
    }
}

sub export_data_base{
    $data_base_export_dir = $config{'databaseExport'};
    $analysis_dir = "$config{'outDir'}/ANALYSIS";
    system("mkdir $data_base_export_dir");
    #Copy the oncoIMPACT prediction
    system("cp $analysis_dir/TEST_PARAM_js.dat $data_base_export_dir/JS.dat");
    system("cp $analysis_dir/RES_*/MODULE.dat $data_base_export_dir/MODULE.dat");
    system("cp $analysis_dir/RES_*/exp_gene_freq_pvalue.dat $data_base_export_dir/PHENO.dat");

    #Get the gene frequencies 
    $result_dir = "GENE_LIST";
    $result_dir = "GENE_LIST_SAMPLE" if($config{'dataType'} eq "RNA_SEQ");
    system("cp $analysis_dir/$result_dir/ALTERATION.dat $data_base_export_dir/");

    #Copy the data set
    system("cp $config{'cnv'} $data_base_export_dir/CNV.dat");
    system("cp $config{'snp'} $data_base_export_dir/SNP.dat");
    system("cp $config{'exp'} $data_base_export_dir/EXPR.dat");
    #Compress the data
    system("tar -C $data_base_export_dir -zcf $data_base_export_dir.tgz .");#To get a flat archive
    system("rm -r $data_base_export_dir");
}

sub prep_cnv {

    my ($file, $outDir ) = @_;

    #print STDERR " *** prep_cnv $file $outDir\n";

    my ( %ht, @samples, @temp, $sample, $gene, $outDir_s);
    open( MATRIX, "$file" );
    print STDERR "[prep_cnv] Opening file:$file\n" if $flag_debug;

    ## Process header
    chomp( @samples = split( /\t/, <MATRIX> ) );
    for ( my $i = 1 ; $i < @samples ; $i++ ) {
	my %cnv;
	$ht{ $samples[$i] } = \%cnv;
    }

    ## Process genes
    while (<MATRIX>) {
	chop $_;
	@temp = split( /\t/, $_ );
	for ( my $i = 1 ; $i < @temp ; $i++ ) {
	    $ht{ $samples[$i] }->{ $temp[0] } = $temp[$i];
	    print STDERR
		"[prep_cnv] Pushing:$samples[$i](sample),$temp[0](gene),$temp[$i](value)\n"
		if $flag_debug;
	}
    }

    close(MATRIX);

    ## Generate results
    foreach $sample ( sort keys %ht ) {
	$outDir_s = "$outDir/$sample";
	system("mkdir $outDir_s") unless ( -s "$outDir_s" );

	print STDERR "[prep_cnv] Writing results to file: $outDir_s/CNV_Data.txt\n" if $flag_debug;
	open( OUT, "> $outDir_s/CNV_Data.txt" );
	foreach $gene ( sort keys %{ $ht{$sample} } ) {
	    print OUT "${gene}_AMPL\t$ht{$sample}->{$gene}\n"
		if ( $ht{$sample}->{$gene} == 1 );

	    #
	    print OUT "${gene}_DEL\t$ht{$sample}->{$gene}\n"
		if ( $ht{$sample}->{$gene} == -1 );
	}
	close(OUT);
    }
}    # end prep_cnv


sub prep_snp {

    my ($file, $outDir) = @_;

    #print STDERR " *** prep_snp $file $outDir\n";

    my ( %ht, @samples, @temp, $sample, $gene, $outDir_s);
    open( MATRIX, "$file" );
    print STDERR "[prep_snp] Opening file:$file\n" if $flag_debug;

    ## Process header
    chomp( @samples = split( /\t/, <MATRIX> ) );
    for ( my $i = 1 ; $i < @samples ; $i++ ) {
	my %snp;

	#print STDERR " **** |$samples[$i]|\n";
	$ht{ $samples[$i] } = \%snp;
    }

    ## Process genes
    while (<MATRIX>) {
	chomp( @temp = split( /\t/, $_ ) );
	for ( my $i = 1 ; $i < @temp ; $i++ ) {
	    $ht{ $samples[$i] }->{ $temp[0] } = $temp[$i];
	    print STDERR
		"[prep_snp] Pushing:$samples[$i](sample),$temp[0](gene),$temp[$i](value)\n"
		if $flag_debug;
	}
    }

    close(MATRIX);

    ## Generate results
    foreach $sample ( sort keys %ht ) {
	$outDir_s = "$outDir/$sample";
	system("mkdir $outDir_s") unless ( -d "$outDir_s" );

	print STDERR "[prep_snp] Writing results to file:$outDir_s/SNP_Data.txt\n" if $flag_debug;
	open( OUT, "> $outDir_s/SNP_Data.txt" );
	foreach $gene ( sort keys %{ $ht{$sample} } ) {
	    print OUT "${gene}_MUT\n" if ( $ht{$sample}->{$gene} == 1 );
	}
	close(OUT);
    }
}    # end prep_snp


sub prep_exp {

    my ($file, $outDir) = @_;

    #print STDERR " *** prep_exp $file $outDir\n";

    my ( %ht, @samples, @temp, $sample, $gene, $outDir_s);
    open( MATRIX, "$file" );
    print STDERR "[prep_exp] Opening file: $file\n" if $flag_debug;  

    ## Process header
    chomp( @samples = split( /\t/, <MATRIX> ) );
    for ( my $i = 1 ; $i < @samples ; $i++ ) {

	#print STDERR " *** $i $samples[$i]\n";<STDIN>;
	my %exp;
	$ht{ $samples[$i] } = \%exp;
    }

    ## Process genes
    while (<MATRIX>) {
	chomp( @temp = split( /\s+/, $_ ) );

	#print STDERR " *** $temp[0] (gene)\n";#<STDIN>;
	for ( my $i = 1 ; $i < @temp ; $i++ ) {
	    $ht{ $samples[$i] }->{ $temp[0] } = $temp[$i];
	}
    }

    close(MATRIX);

    ## Generate results
    foreach $sample ( sort keys %ht ) {
	$outDir_s = "$outDir/$sample";
	system("mkdir $outDir_s") unless ( -d "$outDir_s" );

	print STDERR  "[prep_exp] Writing results to file: $outDir_s/EXPR_Data.txt\n" if $flag_debug;
	open( OUT, "> $outDir_s/EXPR_Data.txt" );
	foreach $gene ( sort keys %{ $ht{$sample} } ) {
	    print OUT "${gene}_UP\t$ht{$sample}->{$gene}\n"
		if ( $ht{$sample}->{$gene} > 0 );
	    print OUT "${gene}_DOWN\t$ht{$sample}->{$gene}\n"
		if ( $ht{$sample}->{$gene} < 0 );
	}
	close(OUT);
    }
}    # end prep_exp


sub merge_and_clean {
    my ($out_dir, $complete_sample_dir, $incomplete_sample_dir) = @_;
    #print STDERR " *** merge_and_clean $out_dir $complete_sample_dir $incomplete_sample_dir\n";
    my $sysCall;
    my $dir = $out_dir . "/$complete_sample_dir";
    system("rm -r $dir") if ( -d $dir );
    system("mkdir $dir");

    $dir = $out_dir . "/$incomplete_sample_dir";
    system("rm -r $dir") if ( -d $dir );
    system("mkdir $dir");

    opendir( DIR, "$out_dir" );
    my @samples = readdir(DIR);
    close(DIR);

    foreach my $dir (@samples) {
	$sample_dir = "$out_dir/$dir";
	
	next if (! -d $sample_dir || substr($dir,0,1) eq "." || $dir eq "$complete_sample_dir" || $dir eq "$incomplete_sample_dir");

	$cnv_file   = "$sample_dir/CNV_Data.txt";
	$snv_file   = "$sample_dir/SNP_Data.txt";
	$expr_file  = "$sample_dir/EXPR_Data.txt";
	if (   -e $cnv_file
	       && -e $snv_file
	       && -e $expr_file )
	{

	    $out_file_name = "$sample_dir/Genelist_Status.txt";
	    $sysCall = "cat $cnv_file $snv_file $expr_file > $out_file_name";
	    print STDERR "[System]$sysCall\n" if $flag_debug;
	    system($sysCall);
	    $sysCall = "mv $sample_dir $out_dir/$complete_sample_dir/";
	    print STDERR " **** [System]$sysCall\n" if $flag_debug;
	    system($sysCall);
	}
	elsif( -e $cnv_file
	       || -e $snv_file
	       || -e $expr_file )
	{
	    $sysCall = "mv $sample_dir $out_dir/$incomplete_sample_dir/";
	    print STDERR " **** [System]$sysCall\n" if $flag_debug;
	    system($sysCall);
	}
    }
    #print STDERR "** end merge_and_clean\n";<STDIN>;
}    # end merge_and_clean


sub run_oncoIMPACT {
	my ($sysCall, @temp);
	$sysCall = "$config{'scriptDir'}/pathway_ana.pl ALL $config{'outDir'}/COMPLETE_SAMPLES $config{'dataType'} $config{'numThreads'} DRIVER_NET $config{'scriptDir'}";
	$sysCall .= " TEST " if $config{'testMode'}; 
	$sysCall .= "&> $config{'outDir'}/run.log";
	print STDERR "$sysCall\n";
	system($sysCall);
}    # end run_oncoIMPACT

sub run_oncoIMPACT_discovery {
	my ($sysCall, @temp);
	my $data_base_dir =  $config{'dataBaseDir'};
	$sysCall = "$config{'scriptDir'}/pathway_ana_uniq_sample.pl $config{'outDir'}/COMPLETE_SAMPLES  $data_base_dir/SAMPLE_INFO $data_base_dir/JS.dat $data_base_dir/PHENO.dat $data_base_dir/MODULE.dat $data_base_dir/ALTERATION.dat DRIVER_NET $config{'scriptDir'} &> $config{'outDir'}/run.log";
	print STDERR "$sysCall\n";
	system($sysCall);
}

sub output_final_result{
    
    $result_dir = "GENE_LIST";
    $result_dir = "GENE_LIST_SAMPLE" if($config{'dataType'} eq "RNA_SEQ");
    
    #The data set driver gene prediction file
    $final_res_file = "$config{'outDir'}/driver_list.txt";
    print STDERR "Output final driver list in $final_res_file\n\n";
    open(OUT, "> $final_res_file");
    print OUT "GENE\tDRIVER_FREQUENCY\tDRIVER_SNV_FREQUENCY\tDRIVER_DELTION_FREQUENCY\tDRIVER_AMPLIFICATION_FREQUENCY\tCANCER_CENSUS\tPAN_CANCER\tIMPACT\tMUTATION_FREQUENCY\tSNV_FREQUENCY\tDELTION_FREQUENCY\tAMPLIFICATION_FREQUENCY\n";
    
    open(IN, "sort -k15,15 -nr $config{'outDir'}/ANALYSIS/$result_dir/ALTERATION.dat |");
    
    while(<IN>){
	chomp(@temp = split(/\t/, $_));
	print OUT "$temp[0]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[11]\t$temp[12]\t$temp[14]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\n" if($temp[14] != 0);
    }
    
    close(OUT);
    close(IN);

    #The sample specific driver gene prediction file
    $final_dir = "$config{'outDir'}/sample_driver_list";
    $intermidate_dir = "$config{'outDir'}/ANALYSIS/$result_dir/SAMPLE_SPECIFIC_DATA";
    #system("mkdir $final_dir");
    print STDERR "Output final sample specific driver lists in $final_dir\n\n";
    system("mv $intermidate_dir $final_dir");
    
}    # end run_oncoIMPACT



sub read_config {
	my ( $file, $hashRef ) = @_;
	my @temp;
	open( FILE, $file );

	while (<FILE>) {
		chop $_;

		@temp = split( /=/, $_ );
		if ( @temp == 2 ) {
			$hashRef->{ $temp[0] } = $temp[1];
		}
	}
	close(FILE);
}    # end read_config

