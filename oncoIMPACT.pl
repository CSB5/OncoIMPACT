#!/usr/bin/perl
use warnings;

my ( $configFile, %config, $flag_debug );

my $help_message = "
This script prepares the data for oncoIMPACT run before launching oncoIMPACT.

Usage:
	oncoIMPACT.pl <config file> <subsample size> <optional: 1 to enable debug mode>
	
* indicates required parameters	


Version:
	0.9.2

Author:
	Burton Chia - chiakhb\@gis.a-star.edu.sg
	Denis Bertrandd - bertrandd\@gis.a-star.edu.sg\n";

if ( @ARGV == 0 ) {
	print $help_message;
	exit 0;
}

$flag_debug = 0;
( $configFile, $subsampleSize, $flag_debug ) = @ARGV;


# Sanity check on user provided parameters
unless(-s $configFile){
	print STDERR "Aborting! Config file does not exist or is empty. Please check the config file and try again.\n";
	exit 3;
}
if($subsampleSize > 1){
	print STDERR "Aborting! Variable subsample size is greater than 1. Please ensure that this variable is a fraction and try again.\n";
	exit 1;
}
if($flag_debug != 0 || $flag_debug != 1){
	print STDERR "Aborting! Debug flag contains an invalid option. Please check the parameter provided and try again.\n";
	exit 1;
}


# Check that dependent system programs are present
print STDERR "[Dependencies] Checking the presence and version of required system programs\n" if $flag_debug;
my $programPath;
# awk
chomp($programPath = `command -v awk`);
if($programPath eq ""){
	print STDERR "Aborting! System command 'awk' not found! Please ensure you are running this programme in a suitable environment.\n";
	exit 2;
} elsif($flag_debug){
	print STDERR "[awk] Path: $programPath\n";
	print STDERR "[awk] Version:" . `awk --version | head -n 1`;
}
# xargs
chomp($programPath = `command -v xargs`);
if($programPath eq ""){
	print STDERR "Aborting! System command 'xargs' not found! Please ensure you are running this programme in a suitable environment.\n";
	exit 2;
} elsif($flag_debug){
	print STDERR "[xargs] Path: $programPath\n";
	print STDERR "[xargs] Version:" . `xargs --version | head -n 1`;
}


# Read config file
print "Reading config file. Please wait...";
read_config( $configFile, \%config );
print "done.\n";


# Prep output directory
system("mkdir $config{'outDir'}") unless ( -s $config{'outDir'} );


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
	print STDERR "Aborting! cnv file does not exist or is empty. Please check the config file and try again.\n";
	exit 1;
}
unless(-s $config{'exp'}){
	print STDERR "Aborting! exp file does not exist or is empty. Please check the config file and try again.\n";
	exit 1;
}
unless(-s $config{'snp'}){
	print STDERR "Aborting! snp file does not exist or is empty. Please check the config file and try again.\n";
	exit 1;
}
unless($config{'testMode'} == 0 || $config{'testMode'} == 1){
	print STDERR "Aborting! testMode option is not a valid flag - valid options are '0' or '1'. Please check the config file and try again.\n";
	exit 1;
}


# Prep data
unless ( -s $config{'outDir'} . "/COMPLETE_SAMPLES" ) {
	print "Preparing CNV data. Please wait...";
	prep_cnv();
	print "done.\n";

	print "Preparing SNP data. Please wait...";
	prep_snp();
	print "done.\n";

	print "Preparing Expression data. Please wait...";
	prep_exp();
	print "done.\n";

	print "Merging and cleaning output directory. Please wait...";
	merge_and_clean();
	print "done.\n";
}


# Run oncoIMPACT
print "\nRunning oncoIMPACT. Please wait...";
run_oncoIMPACT();
print "done.\n";


### Sub-routines ###
sub prep_cnv {
	my ( %ht, @samples, @temp, $sample, $gene, $outDir );
	open( MATRIX, "$config{'cnv'}" );
	print STDERR "[prep_cnv] Opening file:$config{'cnv'}\n" if $flag_debug;

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
		$outDir = "$config{'outDir'}/$sample";
		system("mkdir $outDir") unless ( -s "$outDir" );

		print STDERR "[prep_cnv] Writing results to file: $outDir/CNV_Data.txt\n" if $flag_debug;
		open( OUT, "> $outDir/CNV_Data.txt" );
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
	my ( %ht, @samples, @temp, $sample, $gene, $outDir );
	open( MATRIX, "$config{'snp'}" );
	print STDERR "[prep_snp] Opening file:$config{'snp'}\n" if $flag_debug;

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
		$outDir = "$config{'outDir'}/$sample";
		system("mkdir $outDir") unless ( -d "$outDir" );

	  	print STDERR "[prep_snp] Writing results to file:$outDir/SNP_Data.txt\n" if $flag_debug;
		open( OUT, "> $outDir/SNP_Data.txt" );
		foreach $gene ( sort keys %{ $ht{$sample} } ) {
			print OUT "${gene}_MUT\n" if ( $ht{$sample}->{$gene} == 1 );
		}
		close(OUT);
	}
}    # end prep_snp


sub prep_exp {
	my ( %ht, @samples, @temp, $sample, $gene, $outDir );
	open( MATRIX, "$config{'exp'}" );
	print STDERR "[prep_exp] Opening file: $config{'exp'}\n" if $flag_debug;  

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
		$outDir = "$config{'outDir'}/$sample";
		system("mkdir $outDir") unless ( -d "$outDir" );

		print STDERR  "[prep_exp] Writing results to file: $outDir/EXPR_Data.txt\n" if $flag_debug;
		open( OUT, "> $outDir/EXPR_Data.txt" );
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
	my $sysCall;
	my $dir = $config{'outDir'} . "/COMPLETE_SAMPLES";
	system("rm -r $dir") if ( -s $dir );
	system("mkdir $dir");

	$dir = $config{'outDir'} . "/INCOMPLETE_SAMPLES";
	system("rm -r $dir") if ( -s $dir );
	system("mkdir $dir");

	opendir( DIR, "$config{'outDir'}" );
	my @samples = readdir(DIR);
	close(DIR);

	foreach my $dir (@samples) {
		next if (substr($dir,0,1) eq ".");
		$sample_dir = "$config{'outDir'}/$dir";
		$cnv_file   = "$sample_dir/CNV_Data.txt";
		$snv_file   = "$sample_dir/SNP_Data.txt";
		$expr_file  = "$sample_dir/EXPR_Data.txt";
		if (   -s $cnv_file
			&& -s $snv_file
			&& -s $expr_file )
		{

			$out_file_name = "$sample_dir/Genelist_Status.txt";
			$sysCall = "cat $cnv_file $snv_file $expr_file > $out_file_name";
			print STDERR "[System]$sysCall\n" if $flag_debug;
			system($sysCall);
			$sysCall = "mv $sample_dir $config{'outDir'}/COMPLETE_SAMPLES/";
			print STDERR "[System]$sysCall\n" if $flag_debug;
			system($sysCall);
		}
		else {
			$sysCall = "mv $sample_dir $config{'outDir'}/INCOMPLETE_SAMPLES/";
			print STDERR "[System]$sysCall\n" if $flag_debug;
			system($sysCall);
		}
	}
}    # end merge_and_clean


sub run_oncoIMPACT {
	my ($sysCall, @temp);
	$sysCall = "$config{'scriptDir'}/pathway_ana.pl ALL $config{'outDir'}/COMPLETE_SAMPLES $subsampleSize $config{'numThreads'} DRIVER_NET $config{'scriptDir'} &> $config{'outDir'}/run.log";
	$sysCall = "$config{'scriptDir'}/pathway_ana.pl ALL $config{'outDir'}/COMPLETE_SAMPLES $subsampleSize $config{'numThreads'} DRIVER_NET $config{'scriptDir'} TEST &> $config{'outDir'}/run.log" if $config{'testMode'};
	print STDERR "[System]$sysCall\n" if $flag_debug;
	system($sysCall);

	$final_res_file = "$config{'outDir'}/driver_list.txt";
	open(OUT, "> $final_res_file");
	print OUT "GENE\tDRIVER_FREQUENCY\tDRIVER_SNV_FREQUENCY\tDRIVER_DELTION_FREQUENCY\tDRIVER_AMPLIFICATION_FREQUENCY\tCANCER_CENSUS\tPAN_CANCER\tIMPACT\tMUTATION_FREQUENCY\tSNV_FREQUENCY\tDELTION_FREQUENCY\tAMPLIFICATION_FREQUENCY\n";
	open(IN, "sort -k15,15 -nr $config{'outDir'}/oncoIMPACT_analysis/GENE_LIST/ALTERATION.dat |")
	
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		print OUT "$temp[0]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[11]\t$temp[12]\t$temp[14]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\n" if($temp[14] != 0);
	}
	
	close(OUT);
	close(IN);
	print "Final result in $final_res_file\n\n";
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
