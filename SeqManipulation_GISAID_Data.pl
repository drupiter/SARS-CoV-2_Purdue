#!/usr/bin/perl

# The scripts written below were produced to facilitate research into the spread
# of SARS-CoV-2 through a university community.  They were created to perform specific
# functions, such as combining and reformatting sample information and sequence data
# from GISAID and calculating simple distances between aligned target sequences with
# a collection of test sequences.  This code is very likely sub-otpimal in many ways,
# which we acknowledge.  Please use as you see fit and without reservation. 

$inSeqsFile = $ARGV[0];
$inDataFile = $ARGV[1];
$outFASTA = $ARGV[2];

## Set of functions to properly format sequence files
# %GISAID_seqs = getSeqs($inSeqsFile);
# %GISAID_data = getData($inDataFile);
# %GISAID_names = getNames(%GISAID_data);
#
## Print a FASTA file with the newly named sequences
# open(OUTFASTA,">".$outFASTA);
# foreach $key (keys %GISAID_seqs) {
# print OUTFASTA ">".$GISAID_names{$key}."\n";
# print OUTFASTA $GISAID_seqs{$key}."\n";}
# close(OUTFASTA);


## Set of functions to read in seqs and create hash of distances between target and test
## sequences
# print "READ IN TARGET SEQS\n\n";
# %target = getSeqs($inSeqsFile);
# print "READ IN TEST SEQS\n\n";
# %test = getSeqs($inDataFile);
# foreach $key (keys %target) {$target{$key} = $test{$key};}
#
# %distances = getMinDistances(\%target,\%test);
# @dist_dist = sumDistances(\%distances);
#
# $count = 0;
# print "\nDISTRIBUTION OF DISTANCES BETWEEN TARGET SEQS AND TEST SEQS\n";
# foreach $dis (@dist_dist) {print $count++,"\t".$dis."\n";}
#
# open(OUTFILE,">".$outFASTA);
# foreach $key (keys %distances) {print OUTFILE $key."\t".$distances{$key}."\n";}
# close(OUTFILE);


%my_distances = readDistHash("temp.distances.txt");
%my_close_distances = returnHashOfShorterDist(2,\%my_distances);
%all_seqs = getSeqs($inSeqsFile);
open(OUTFILE,">",$outFASTA);
foreach $key (keys %my_close_distances) {print OUTFILE ">".$key."\n".$all_seqs{$key}."\n";}
close(OUTFILE);

sub returnHashOfShorterDist {
	my ($min_dist,$inHash) = @_;
	
	my %out_hash = ();
	
	my %temp_hash = %{$inHash};
	
	foreach $key (keys %temp_hash) {
		if ($temp_hash{$key} <= $min_dist) {$out_hash{$key} = $temp_hash{$key};}
		else {}
	}
	
	return %out_hash;}

sub readDistHash {
	my ($inHashFile) = @_;
	
	my %out_hash = ();
	
	open(INFILE,$inHashFile);
	foreach $line (<INFILE>) {
		chomp($line);
		($key,$value) = split(/\t/,$line);
		$out_hash{$key} = $value;}
	close(INFILE);
	
	return %out_hash;
	}


sub combineDistHashes {
	my (@hashes) = @_;
	
	my %out_hash = ();
	
	foreach $hash (@hashes) {
		%temp_hash = %{$hash};
		foreach $key (keys %temp_hash) {
			if (exists $out_hash{$key}) {}
			else {$out_hash{$key} = $temp_hash{$key};}
		}
	}
	return %out_hash;
	}

sub getOnDistances {
	my ($max_dist,$distance_hash) = @_; 
	
	my %dist_seqs = %{$distance_hash};
	my %targ_seqs = ();
	
	foreach $key (keys %dist_seqs) {
		if ($dist_seqs{$key} <= $max_dist) {$targ_seqs{$key} = $dist_seqs{$key};}
		else {}
		}
	
	return %targ_seqs;}

sub sumDistances {
	my ($target_hash) = @_;
	
	#        0...             ...10...             ...20...             ...30,>30
	my @a = (0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0, 0);
	
	my %dist_hash = %{$target_hash};

	print "\nSUMMARIZING DISTANCES\n";	
	foreach $key (keys %dist_hash) {
		$dist = $dist_hash{$key};
		if ($dist == 0) {$a[0]++;}
		elsif ($dist == 1) {$a[1]++;}
		elsif ($dist == 2) {$a[2]++;}
		elsif ($dist == 3) {$a[3]++;}
		elsif ($dist == 4) {$a[4]++;}
		elsif ($dist == 5) {$a[5]++;}
		elsif ($dist == 6) {$a[6]++;}
		elsif ($dist == 7) {$a[7]++;}
		elsif ($dist == 8) {$a[8]++;}
		elsif ($dist == 9) {$a[9]++;}
		elsif ($dist == 10) {$a[10]++;}
		elsif ($dist == 11) {$a[11]++;}
		elsif ($dist == 12) {$a[12]++;}
		elsif ($dist == 13) {$a[13]++;}
		elsif ($dist == 14) {$a[14]++;}
		elsif ($dist == 15) {$a[15]++;}
		elsif ($dist == 16) {$a[16]++;}
		elsif ($dist == 17) {$a[17]++;}
		elsif ($dist == 18) {$a[18]++;}
		elsif ($dist == 19) {$a[19]++;}
		elsif ($dist == 20) {$a[20]++;}
		elsif ($dist == 21) {$a[21]++;}
		elsif ($dist == 22) {$a[22]++;}
		elsif ($dist == 23) {$a[23]++;}
		elsif ($dist == 24) {$a[24]++;}
		elsif ($dist == 25) {$a[25]++;}
		elsif ($dist == 26) {$a[26]++;}
		elsif ($dist == 27) {$a[27]++;}
		elsif ($dist == 28) {$a[28]++;}
		elsif ($dist == 29) {$a[29]++;}
		elsif ($dist == 30) {$a[30]++;}
		else {$a[31]++;}
	}
	
	return @a;}

# Function that returns MINIMUM distances between sequences in the target set and
# sequences in the test set.
# This function takes hash inputs for the target and test sets of sequences.
sub getMinDistances {
	my ($target_hash,$test_hash) = @_;
	my %a = ();
	my $dist = 0;
	
	my %target_seqs = %{$target_hash};
	my %test_seqs = %{$test_hash};
	
	foreach my $key (keys %test_seqs) {$a{$key} = 99999;} 
	
	print "CALCULATING DISTANCES\n";
	foreach $target_key (keys %target_seqs) {
		foreach $test_key (keys %test_seqs) {
			$dist = getDist($target_seqs{$target_key},$test_seqs{$test_key});
			if ($dist > $a{$test_key}) {}
			else {$a{$test_key} = $dist;}
			}
		}
	
	return %a;}


# Function that reads in two sequences as strings and returns the number of differences
# between them as an integer.
# This function ignores positions with 'N' or '-' values in either sequence. 
sub getDist {
	my ($seq1,$seq2) = @_;
	my $seq_dist = 0;
	my $seq_length = -99;
	my $count = 0;
	
	@seq1 = split(//,$seq1);
	@seq2 = split(//,$seq2);
	$seq_length = @seq1 - 1;
	
	foreach $nuc1 (@seq1) {
		if ($nuc1 eq $seq2[$count])	{}
		elsif ($nuc1 eq 'N') {}
		elsif ($nuc1 eq 'n') {}
		elsif ($nuc1 eq '-') {}
		elsif ($seq2[$count] eq 'N') {}
		elsif ($seq2[$count] eq 'n') {}
		elsif ($seq2[$count] eq '-') {}
		else {
			$seq_dist = $seq_dist + 1;}
		$count++;}

	return $seq_dist;}

# Function that reads in sequences in FASTA format, interleaved or no.
# Function returns a hash with sequence names as keys and sequences as values.
sub getSeqs {
	my ($inFile) = @_;
	my %a = ();
	my $hash_key = 'default';
	
	open(INFILE,$inFile);
	foreach my $line (<INFILE>) {
		if ($line =~ m/>/) {
			chomp($line);
			($junk,$line) = split(/>/,$line);
			$hash_key = $line;
			$line = <INFILE>;
			chomp($line);
			$a{$hash_key} = $line;}
		else {
			chomp($line);
			$a{$hash_key} = $a{$hash_key}.$line;}
	}
	close(INFILE);
	
	return %a;}

# Function that reads in metadata in tab spaced values (TSV) form.
# Function returns a hash with sequence names as keys and metadata as a long string.
sub getData {
	my ($inFile) = @_;
	my %a = ();

	open(INFILE,$inFile);
	foreach $line (<INFILE>) {
		if ($line =~ m/hCoV-19/) {
			chomp($line);
			($key,@line) = split(/\t/,$line);
			$a{$key} = join('___',@line);}
		else {}
		}
	close(INFILE);
	
	return %a;}

# Function that reads in a hash of metadata and creates taxa names with embedded data.
# Function returns a hash with original sequence IDs as keys and data rich names as a long
# string.
# The fields that are embedded in the taxa names are hard-coded in this function; this
# should be changed - it would not be hard to do.
sub getNames {
	my (%seqData_temp) = @_;
	my %a = ();
	my @b = ();
	
	# GISAID order of information in metadata TSV file
	# 0 = virus
	# 1 = gisaid_epi_isl
	# 2 = genbank_accession
	# 3 = date
	# 4 = region
	# 5 = country
	# 6 = division
	# 7 = location
	# 8 = region_exposure
	# 9 = country_exposure
	# 10 = division_exposure
	# 11 = segment
	# 12 = length
	# 13 = host
	# 14 = age
	# 15 = sex
	# 16 = Nextstrain_clade
	# 17 = pangolin_lineage
	# 18 = GISAID_clade
	# 19 = originating_lab
	# 20 = submitting_lab
	# 21 = authors
	# 22 = url
	# 23 = title
	# 24 = paper_url
	# 25 = date_submitted
	# 26 =purpose_of_sequencing

	foreach $key (keys %seqData_temp) {
		my $b = $seqData_temp{$key};
		@b_array = split(/___/,$b);
		# hard-coded in the fields to include in the taxon name, alas
		$a{$key} = $key."|".$b_array[1]."|".$b_array[5]."|".$b_array[3];
		}
	
	return %a;}
