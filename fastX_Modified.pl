#!/usr/bin/perl -w

#Use Options: perl ./fastX_Modified.pl --bcfile barcode_files.txt --fwd_file ../lane1_NoIndex_FW.fastq --rev_file ../lane1_NoIndex_RC.fastq --bol --exact --dir ProcessedFiles
use strict;	
use IO::Handle;
use Getopt::Long;

# Global flags and arguments,
# Set by command line argumens
my $barcode_file ;
my $barcodes_at_eol = 0 ;
my $barcodes_at_bol = 0 ;
my $exact_match = 0 ;
my $allowed_mismatches = 1;
my $fwdseq_file ;
my $revseq_file  ;
my $outputDir;

# Global variables
# Populated by 'create_output_files'
my %filenames;
my %files;
my %counts = ( 'unmatched' => 0 );
my $barcodes_length;
my %barcodes;
my $input_file_io;


#Start of the Program

&parse_command_line();
&load_barcode_file($barcode_file);
&match_sequences($fwdseq_file, $revseq_file);
&writeOutput();

exit 0;

sub match_sequences {
	my $fwd_filename = shift;
	my $rev_filename = shift;
	print "$fwd_filename\n$rev_filename\n";
	open (my $fwdIN, "<$fwd_filename") or die "Error: failed to open FW file\n";
	open (my $revIN, "<$rev_filename") or die "Error: failed to open RC file\n";
	#my ($fline, $rline);
	while(defined(my $fline = <$fwdIN>) && defined(my $rline = <$revIN>)) {
		if (($fline =~ /^\@/) && ($rline =~ /^\@/)) {
			my @fwID = split(/\s+/, $fline); my @revID = split(/\s+/, $rline);
			#print "$fwID[0]\t$revID[0]\n";
			if ($fwID[0] eq $revID[0]) {
			chomp(my $fSeq = <$fwdIN>); chomp(my $rSeq = <$revIN>);
				<$fwdIN>; <$revIN>;
				chomp(my $fQual = <$fwdIN>); chomp(my $rQual = <$revIN>);
				my ($fwdBar, $revBar);
				if ($barcodes_at_eol) {
					$fwdBar = substr($fSeq, -$barcodes_length);
					$revBar = substr($rSeq, -$barcodes_length);
				}
				else {
					$fwdBar = substr($fSeq, 0, $barcodes_length);
					$revBar = substr($rSeq, 0, $barcodes_length);
				}
				if (exists $barcodes{$fwdBar}) {
					if ($fwdBar eq $revBar) {
						if ($barcodes_at_eol) {
							$fSeq = substr($fSeq, 0, -$barcodes_length);
							$fQual = substr($fQual, 0, -$barcodes_length);
							$rQual = substr($rQual, 0, -$barcodes_length);
							$rSeq = substr($rSeq, 0, -$barcodes_length);
						}
						else {
							$fSeq = substr($fSeq, $barcodes_length);
							$fQual = substr($fQual, $barcodes_length);
							$rQual = substr($rQual, $barcodes_length);
							$rSeq = substr($rSeq, $barcodes_length);
						}				
						$filenames{$barcodes{$fwdBar}}{$fwID[0]} = [ $fSeq, $fQual, $rSeq, $rQual ];
					}
					else {
						$filenames{"RC_unmatch"}{$fwID[0]} = [ $fSeq, $fQual ];
					}
				}
				elsif (exists $barcodes{$revBar}) {
					$filenames{"FW_unmatch"}{$revID[0]} = [ $rSeq, $rQual ];
				}
				else {
					$filenames{"unmatch"}{$revID[0]} = [ $fSeq, $fQual, $rSeq, $rQual ];
				}
			}
		}
	}		
}	

sub writeOutput {
	my $sum = 0;
	print "Identifier\tCount\n";
	foreach my $k (sort keys %filenames) {
		open (my $out, "> $outputDir/$k.txt") or die "$!";
		my $idCount = scalar(keys %{$filenames{$k}});
		print "$k\t$idCount\n";
		$sum = $sum + $idCount;
		foreach my $kx (keys %{$filenames{$k}}) {
			my @outArr = @{$filenames{$k}{$kx}};
			if (scalar(@outArr) > 2) {
				print $out "$kx FW\n$outArr[0]\n+\n$outArr[1]\n";
				print $out "$kx RC\n$outArr[2]\n+\n$outArr[3]\n";
			}
			else {
				print $out "$kx\n$outArr[0]\n+\n$outArr[1]\n";
			}
		}
	}
	print "Total\t$sum\n";
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


