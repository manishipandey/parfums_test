#!/usr/bin/perl -w

# Written by Alejandro Reyes and Kevin Forsberg

##Version 3 requires that the index be a perfect match (though calculates the stats for imperfect matches) and that the Forward read be

use strict;

my $usage = "Get_BC_Gautam.pl Table.txt <MappingFile> <IlluminaOutputForReads>  <IlluminaOutputRevReads> >statsfile.txt";

if (@ARGV != 3) {
	die "\n$usage \n\n";
}

my $infile = shift; #Mapping File
my $wordCount = `wc -l $infile`; #number of lines in a mapping file
my @result = split /\s+/, $wordCount;#++
my $numSamples = $result[0]; #stores number of lines in mapping file
my $j = 0;
my %pm;
my $mm = 0;
my %mmCnt;
my @readNames;
my $mismatches = 0;
my $readCount = 0;
my %qualRemoved = ();
my $forName = "";
my $revName = "";
my $forSeq = "";
my $revSeq = "";
my $forQual = "";
my $revQual = "";
my $prntcnt = 0;

open (IN, "<$infile") or die ("Couldn't open file: $infile\n");

my $file1=shift; # forward sequence read file
my $file2=shift; # reverse sequence read file

#This code identify all the barcodes in use in Table and send the appropriate Grep command.

my %BC1=();
my %forHandle;
my %revHandle;
# read mapping file
while (my $line=<IN>){
	chomp $line;

	my @entry=split /\s+/, $line;#++
	my $bc1 = $entry[0];
	$BC1{$bc1}=$entry[2]; # stores barcode and sampleid in a hash from mapping file
}

close IN;

open (OUT3, ">DistributionStats_$infile") or die ("Couldn't open file: Stats\n");
open (BClog,">GetBC_log_$infile") or die ("Coudn't open file: Stats\n");
print OUT3 "Sample\tPerfectMatches\t1_Mismatch\t2_Mismatches\n";

print BClog "File Variables\tValues\n";
print BClog "\$file1\t$file1\n";
print BClog "\$file2\t$file2\n";


open(IN, "<$file1") || die ("Couldn't open file: $file1\n"); # forward sequence read file
open (IN2, "<$file2") or die ("Couldn't open $file2\n"); # reverse sequence read file

while ((my $line=<IN>) && (my $line2 =<IN2>)){
  my $current;
  my $check = 0;
  my $cnt = 0;
  chomp $line;
  chomp $line2;
  if (($line =~ m/^\@HW/) && ($line2 =~ m/^\@HW/)){
  	my @nameTable = split / /, $line;#/
  	my @nameTable2 = split / /, $line2;#/
    my $tableLen = scalar(@nameTable);
    my $tlen2 = scalar(@nameTable2);

  	print BClog "\@nameTable2\t@nameTable2\n";
    print BClog "\@nameTable\t@nameTable\n";
  	print BClog "\$tableLen\t$tableLen\n";
    print BClog "\$tlen2\t$tlen2\n";

    if (($tableLen == 2) && ($tlen2 == 2)){
  		$readCount++;
  		$nameTable[0] =~ s/\@//;
  		$nameTable2[0] =~ s/\@//;
  		$forName = $nameTable[0]."\#0\/1";
  		$revName = $nameTable2[0]."\#0\/3";
  		$forSeq = <IN>;
  		$revSeq = <IN2>;
		chomp $forSeq;
		chomp $revSeq;

        print BClog "\$readCount\t$readCount\n";
        print BClog "\$forName\t$forName\n";
  	    print BClog "\$revName\t$revName\n";
        print BClog "\$forSeq\t$forSeq\n";
        print BClog "\$revSeq\t$revSeq\n";

    }
	elsif (($tableLen == 1) && ($tlen2 == 1)){
  		$readCount++;
  		if ($prntcnt == 0){
  			print "The Standard FastQ Header-Format was not used. A workaround was invoked, if the format is as occurrs following Nicole's phiX removal script all is fine. If not, errors may ensue. Double Check to make sure there were no errors.\n";
  			$prntcnt++;
  		}
  		$nameTable[0] =~ s/\@//;
  		$nameTable2[0] =~ s/\@//;
  		$nameTable[0] =~ s/\/1//;
  		$nameTable2[0] =~ s/\/2//;
  		$forName = $nameTable[0]."\#0\/1";
  		$revName = $nameTable2[0]."\#0\/3";
  		$forSeq = <IN>;
  		$revSeq = <IN2>;
		chomp $forSeq;
		chomp $revSeq;
  	}
  }
  if (($line =~ m/^\+$/) && ($line2 =~ m/^\+$/)){
  	$forQual = <IN>;
  	$revQual = <IN2>;
  	chomp $forQual;
  	chomp $revQual;
  	my $seqLength = length($forSeq);
  	if ($seqLength == 108){
  	    my $tempSeq = substr($forSeq, 0, 101);
  	    $forSeq = $tempSeq;
    }

    print BClog "\$forQual\t$forQual\n";
    print BClog "\$revQual\t$revQual\n";

    print BClog "\$seqLength\t$seqLength\n";
    print BClog "\$forSeq\t$forSeq\n";

    my $revSeqLength = length($revSeq);
  	foreach my $key (keys(%BC1)){
		if (($revSeq =~ m/^$key/i) && ($forSeq =~ m/^$key/i)){
			$current = $key;
			$pm{$key}++;
			$check = 1;
			last;
		}elsif ($mismatches > 0){
			my @seqArray = split //, $forSeq;
			my @revSeqArray = split //, $revSeq;
                my @bcArray = split //, $key;
                my $seqLength = length($forSeq);
                if (($revSeqLength == 101) && ($seqLength == 101)){
                	for (my$i=0;$i<8;$i++){
                       	if ($revSeqArray[$i] =~ m/$bcArray[$i]/){
                                $cnt++;
                        }
                        if ($seqArray[$i] =~ m/$bcArray[$i]/){
                        	$cnt++;
                        }
                	}
                }
         }
	}

    print BClog "\$check\t$check\n";
    print BClog "\$cnt\t$cnt\n";
    print BClog "\$mismatches\t$mismatches\n";

	if ($check == 0){
		next;
	}
	if (($cnt >= (16-$mismatches)) || ($check == 1)){
		my $forRead = substr($forSeq, 8, 93);
		my $revRead = substr($revSeq, 8, 93);
		$forQual = substr ($forQual, 8, 93);
		$revQual = substr ($revQual, 8, 93);
			if (exists($forHandle{$current})){
				my $fh = $forHandle{$current};
				my $fh2 = $revHandle{$current};
				print $fh "$forName\:\n$forRead\n+\n$forQual\n";
				print $fh2 "$revName\:\n$revRead\n+\n$revQual\n";
			} else {
				my $handle;
				my $handle2;
				open ($handle, ">Forward_$BC1{$current}.txt") or die ("Couldn't open file: >Forward_$BC1{$current}.txt\n");
				open ($handle2, ">Reverse_$BC1{$current}.txt") or die ("Couldn't open file: >Reverse_$BC1{$current}.txt\n");
				print $handle "$forName\:$forRead\:$forQual\n";
				print $handle2 "$revName\:$revRead\:$revQual\n";
				$forHandle{$current} = $handle;
				$revHandle{$current} = $handle2;
			}

	}
	if ($cnt >= 14){
      $mm = 16 - $cnt;
      ${$mmCnt{$current}}[$mm]++;
    }
  }
}

foreach my $key (keys(%BC1)){
	print OUT3  "$BC1{$key}\t$pm{$key}\t${$mmCnt{$key}}[1]\t${$mmCnt{$key}}[2]\n";
	close($forHandle{$key});
	close($revHandle{$key});
}
print OUT3 "There were $readCount reads from this lane, in total.\n";
close(IN);
close(IN2);
close(OUT3);
close(BClog);
