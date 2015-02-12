#!/usr/bin/perl -w

###################################################
# This code is written to map barcodes in FW & RC #
# fastq files and keep the ones which meets the   #
# criteria.					  #	
#						  #
#	Written By: Manish Boolchandani		  #
#	Part of: PaRFumS			  #
#						  #
###################################################

#Quick Use Cmd: perl ./preprocessing_v1.pl --bcfile barcode_files.txt --fwd_file ../lane1_NoIndex_FW.fastq --rev_file ../lane1_NoIndex_RC.fastq --bol --exact --dir ProcessedFiles

use strict;	
use IO::Handle;
use Getopt::Long;

# Global flags and arguments,
# Set by command line argumens
my $barcode_file ;
my $barcodes_at_eol = 0 ;
my $barcodes_at_bol = 0 ;
my $exact_match = 0 ;
my $unmatch = 0;
my $allowed_mismatches = 1;
my $fwdseq_file ;
my $revseq_file  ;
my $outputDir;

# Global variables
# Populated by 'create_output_files'
my %filenames;
my %trackFiles;
my %counts;
my $barcodes_length;
my %barcodes;
my $input_file_io;
my $flagFile = 0;

#Global file variables
my ($fwdID, $revID, $fwdSeq, $revSeq, $fwdQual, $revQual);

#START OF PROGRAM

#Subroutine to get the arguments from command line
&parse_command_line();
#Setting up variables for including unmatch entries
if ($unmatch) {
	$counts{'FW_unmatch'} = 0;
	$counts{'RC_unmatch'} = 0;
}
#Subroutine to read and load barcode file in hash with following structure:
# $barcodes{BARCODE} = IDENTIFIER
&load_barcode_file($barcode_file);

#Open the forward and reverse fastq file.
open (my $fwdIN, '< :raw', $fwdseq_file) or die "Error: failed to open FW file\n"; #:raw allows to be read the files in ASCII.
open (my $revIN, '< :raw', $revseq_file) or die "Error: failed to open RC file\n";

#Subroutine to read both fastq files and match barcodes.
&match_sequences();

#Subroutine to display summary on STDOUT.
&writeOutput();

exit 0;
#END OF PROGRAM


sub match_sequences {
	#Subroutine to initialize all output identifier files.
	&create_output_files();
	#Reading fastQ files for analysis
	while(defined($fwdID = <$fwdIN>) && defined($revID = <$revIN>)) {
		#Initializing variables with reads information
		if (($fwdID =~ /^\@/) && ($revID =~ /^\@/)) {
			my @fwID = split(/\s+/, $fwdID); my @reID = split(/\s+/, $revID);
			chomp($fwdSeq = <$fwdIN>); chomp($revSeq = <$revIN>);
			<$fwdIN>; <$revIN>;
			chomp($fwdQual = <$fwdIN>); chomp($revQual = <$revIN>);
			#Getting barcodes from reads
			my ($fwdBar, $revBar);
			if ($barcodes_at_eol) {
				$fwdBar = substr($fwdSeq, -$barcodes_length);
				$revBar = substr($revSeq, -$barcodes_length);
			}
			else {
				$fwdBar = substr($fwdSeq, 0, $barcodes_length);
				$revBar = substr($revSeq, 0, $barcodes_length);
			}
			#Matching barcodes from FW and RC reads and to the ones on list
			if ((exists $barcodes{$fwdBar}) && ($fwdBar eq $revBar)) {
				if ($barcodes_at_eol) {
					$fwdSeq = substr($fwdSeq, 0, -$barcodes_length);
					$fwdQual = substr($fwdQual, 0, -$barcodes_length);
					$revQual = substr($revQual, 0, -$barcodes_length);
					$revSeq = substr($revSeq, 0, -$barcodes_length);
				}
				else {
					$fwdSeq = substr($fwdSeq, $barcodes_length);
					$fwdQual = substr($fwdQual, $barcodes_length);
					$revQual = substr($revQual, $barcodes_length);
					$revSeq = substr($revSeq, $barcodes_length);
				}
				my $fileTrack = $trackFiles{$barcodes{$fwdBar}};
				$counts{$barcodes{$fwdBar}}++;				
				&write_record($fileTrack, 1);
			}
			else {
				if ($unmatch) {
					my $fileT= $trackFiles{"FW_unmatch"};
					$counts{"FW_unmatch"}++;
					&write_record($fileT, 2);
					$fileT= $trackFiles{"RC_unmatch"};
					$counts{"RC_unmatch"}++;
					&write_record($fileT, 3);
				}
			}
		}
	}		
}	

sub create_output_files {
	foreach my $barcode (keys %barcodes) {
		$counts{$barcodes{$barcode}} = 0;
	}
	foreach my $fileName (keys %counts) {
		open (my $fileOut, "> ./$outputDir/$fileName.txt") or die "$!";
		$trackFiles{$fileName} = $fileOut;
	}
}

sub write_record {
	my $file = shift;
	my $flag = shift;
	if ($flag == 1) {
		print $file "$fwdID\n$fwdSeq\n+\n$fwdQual\n$revID\n$revSeq\n+\n$revQual\n";
	}
	elsif ($flag == 2) {
		print $file "$fwdID\n$fwdSeq\n+\n$fwdQual\n";
	}
	elsif ($flag == 3) {
		print $file "$revID\n$revSeq\n+\n$revQual\n";
	}
}	

sub writeOutput {
	foreach my $k (keys %counts) {
		print "$k\t$counts{$k}\n";
	}
}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "bcfile=s" => \$barcode_file,
				  "eol"  => \$barcodes_at_eol,
				  "bol"  => \$barcodes_at_bol,
				  "exact" => \$exact_match,
				  "fwd_file=s" => \$fwdseq_file,
				  "rev_file=s" => \$revseq_file,
				  "mismatches=i" => \$allowed_mismatches,
				  "dir=s" => \$outputDir,
				  "unmatch" => \$unmatch,
				  "help" => \$help
				  ) ;

	usage() if ($help);

	die "Error: barcode file not specified (use '--bcfile [FILENAME]')\n" unless defined $barcode_file;

	if ($barcodes_at_bol == $barcodes_at_eol) {
		die "Error: can't specify both --eol & --bol\n" if $barcodes_at_eol;
		die "Error: must specify either --eol or --bol\n" ;
	}

	$allowed_mismatches = 0 if $exact_match;

	die "Error: invalid value for mismatches (valid values are 0 or more)\n" if ($allowed_mismatches<0);

	if (not(-d $outputDir)) {
		system("mkdir", "$outputDir");
		print "$outputDir made\n";
	}
	else {
		print "$outputDir exists: Overwriting files\n";
	}
	exit unless $result;
}

sub load_barcode_file ($) {
	my $filename = shift;

	open BCFILE,"<$filename" or die "Error: failed to open barcode file ($filename)\n";
	while (<BCFILE>) {
		next if m/^#/;
		chomp;
		my ($barcode, $ident) = split ;

		$barcode = uc($barcode);

		# Sanity checks on the barcodes
		die "Error: bad data at barcode file ($filename) line $.\n" unless defined $barcode;
		die "Error: bad barcode value ($barcode) at barcode file ($filename) line $.\n"
			unless $barcode =~ m/^[AGCT]+$/;

		die "Error: bad identifier value ($ident) at barcode file ($filename) line $. (must be alphanumeric)\n"
			unless $ident =~ m/^\w+$/;

		die "Error: badcode($ident, $barcode) is shorter or equal to maximum number of " .
		    "mismatches ($allowed_mismatches). This makes no sense. Specify fewer  mismatches.\n"
		    	if length($barcode)<=$allowed_mismatches;

		$barcodes_length = length($barcode) unless defined $barcodes_length;
		die "Error: found barcodes in different lengths. this feature is not supported yet.\n"
			unless $barcodes_length == length($barcode);

		$barcodes{$barcode} = $ident;
	}
	close BCFILE;
}

sub usage {
print<<EOF;
PaRFuMs preprocessing file, by Manish Boolchandani (manish\@wustl.edu), 11feb2015

This program reads forward and reverse FASTQ files and splits it into several 
smaller files, based on barcode matching. 
Output files will be writen to disk
Summary will be printed to STDOUT.

usage: $0 --bcfile FILE --fwd_file FORWARD_FASTQ_FILE --rev_file REVERSE_FASTQ_FILE 
	--dir OUTPUT_DIR [--bol|--eol] [--unmatch] [--help]

Arguments:

--bcfile FILE			- Barcodes file name. (Barcodes and Identifiers separated
				by tab)
--fwd_file FW_FASTQ_FILE	- Forward FastQ file.
--rev_file RC_FASTQ_FILE 	- Reverse FASTQ file.
--dir OUTPUT_DIR		- Directory in which all output files will be saved.
--bol				- Try to match barcodes at the BEGINNING of sequences.
		  		(What biologists would call the 5' end, and programmers
		  		would call index 0.)
--eol				- Try to match barcodes at the END of sequences.
		  		(What biologists would call the 3' end, and programmers
		  		would call the end of the string.)
		  		NOTE: only one of --bol, --eol must be specified; default
				is --bol.
--unmatch			- Make the files of unmatched entries too. (FW_unmatch & RC_unmatch)
		  		(Default is not to store unmatch entries)
--help				- This helpful help screen.

EOF
exit 1;
}
