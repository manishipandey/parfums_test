#!/usr/bin/perl -w

# Written by Alejandro Reyes and Kevin Forsberg
## For a given prefix will run the whole assembly pipeline

use strict;
use constant WAIT_TIME => '300' ; # 300 seconds = 5 minutes
use constant NO_CHANGE_IN_FILE_FOR_X_SECONDS => '60';
my $MAX_READS = 4000000;

if ((@ARGV < 2) || (@ARGV > 3)) {
	die "\n\nUsage: PARFuMS.pl <Prefix> <MappingFile> <OptionalFlag> >& FileWithStdErrStdOut\n\tIt needs to have a file in the working directory that is called <prefix>.fasta\n\n";
}
my $type = "a";
if ($ARGV[2] =~ m/^-a$/i){
	print "\n\nYou have chosen your function of interest as: ANTIBIOTIC RESISTANCE.  Classification Keyword List(s) will be chosen accordingly. \n\n";
	$type = "a";
}
elsif ($ARGV[2] =~ m/^-b$/i){
	print "\n\nYou have chosen your function of interest as: BIOFUEL-RELATED.  Classification Keyword List(s) for this type do not yet exist, PARFuMS will give only its standard output. \n\n";
	$type = "b";
}
elsif ($ARGV[2] =~ m/^-c$/i){
	print "\n\nYou have chosen your function of interest as: CARBOHYDRATE-RELATED.  Classification Keyword List(s) for this type do not yet exist, PARFuMS will give only its standard output. \n\n";
	$type = "c";
}
else{
	print "\n\nPARFuMS does not understand your flag indicating your function of interest, it will assume the default function: ANTIBIOTIC RESISTANCE.  Classification Keyword List(s) will be chosen accordingly. \n\n";
	$type = "a";
}

my $home="/home/gdlab/shared/old-shared";
my $return_file="";

my $dir=`pwd`;
chomp $dir;
my $prefix = shift;
my $nd= "$dir/$prefix";

my $numReads = `grep _1 $prefix.fasta | wc -l`;
my $tempFasta = "$prefix"."_short.fasta";
my $loop = 0;
my $x = 0;
my $srcc_command = "";
my $srcc_in = "";
my $srcc_suf = "";
my $returnCommand = "";

if (-e "$nd/$prefix.fasta.gz" && -e "$nd/$prefix\_clean.fna.gz"){
  system("cd $nd/; gzip -d $prefix.fasta.gz; gzip -d $prefix\_clean.fna.gz; mv $prefix.fasta ../.; mv $prefix\_clean.fna ../.; rm *; mv ../$prefix.fasta .; mv ../$prefix\_clean.fna .");
}
elsif (-e "$prefix.fasta" or -e "$nd/$prefix.fasta.gz" ){
  if (-e "$nd/$prefix.fasta.gz"){
    system("cd $nd/; gzip -d $prefix.fasta.gz; mv $prefix.fasta ../.; cd ..; rm -r $nd");
  }
  my $file = shift;

  my $adapterA="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGGTTG";
  my $adapterB="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTG";
  my $endVector1="GACCTCGAGGGGGGGCCCGGTACCTTTCTCCTCTTTAATGAATTCGGXXXXAGATCGGAAGAGCG";
  my $endVector1R="CGCTCTTCCGATCTYYYYCCGAATTCATTAAAGAGGAGAAAGGTACCGGGCCCCCCCTCGAGGTC";
  my $endVector2="GACGGTATCGATAAGCTTGATATCGAXXXXAGATCGGAAGAGCG";
  my $endVector2R="CGCTCTTCCGATCTYYYYTCGATATCAAGCTTATCGATACCGTC";

  my $bc_line=`grep $prefix $file`;
  my @bc_arr = split /\s+/, $bc_line;#++
  my $bc=reverse($bc_arr[0]);
  $bc =~ tr/ATCG/TAGC/;

  system ("mkdir $nd");

  chomp $numReads;
  if ($numReads > $MAX_READS){
	my $percent = int(($MAX_READS / $numReads) * 100);
	print"\nThere Exist $numReads reads in the given FASTA file, whereas this pipeline is only verified for up to $MAX_READS reads. Thus, approximately $percent percent of the original reads were maintained in a shortened fasta file, for use in the remainder of the assembly pipeline\n";
	$percent = (1 - ($percent / 100));
	system ("mv $prefix.fasta $nd/$prefix.fasta");
	system ("cd $nd/; perl $home/scripts/PARFuMS_scripts/GrabRandomReads.pl $prefix.fasta $percent");
	$return_file=&wait_for_file("$nd/$tempFasta");
	system ("cd $nd/; mv $prefix.fasta $prefix.fasta.original");
	system ("cd $nd/; mv $tempFasta $prefix.fasta");
	system ("cd $nd/; gzip $prefix.fasta.original");
	$loop = 1;
   }

  # First will make its own copy of the adapters that will include the BC.
  $endVector1 =~ s/XXXX/$bc/;
  $endVector1R =~ s/YYYY/$bc_arr[0]/;
  $endVector2 =~ s/XXXX/$bc/;
  $endVector2R =~ s/YYYY/$bc_arr[0]/;

  $adapterA=$bc.$adapterA;
  $adapterB=$bc.$adapterB;
  open (OUT, ">$nd/Adapters_$prefix.fna") or die ("Couldn't open file $nd/Adapters_$prefix.fna\n");
  print OUT ">Vector1\n$endVector1\n>Vector1R\n$endVector1R\n>Vector2\n$endVector2\n>Vector2R\n$endVector2R\n>Adapter1\n$adapterA\n>Adapter2\n$adapterB\n";
  close OUT;

  if ($loop < 1){
  system("cd $nd/; mv ../$prefix.fasta .");
  }

  ## will split cross_match in jobs of 200,000 seqs and parse them
  system ("cd $nd/; $home/scripts/PARFuMS_scripts/split_run_check_combine_mod.pl -i $prefix.fasta -s crossMatch -n 200000 -p 'cross_match INCLUDE_INFILE Adapters_$prefix.fna -gap1_only -minmatch 6 -minscore 10 -gap_init -3 | $home/scripts/PARFuMS_scripts/CleanAdapter-Illumina_PE_mod.pl INCLUDE_INFILE'");
  $srcc_command = "cross_match INCLUDE_INFILE Adapters_$prefix.fna -gap1_only -minmatch 6 -minscore 10 -gap_init -3 | $home/scripts/PARFuMS_scripts/CleanAdapter-Illumina_PE_mod.pl INCLUDE_INFILE";
  $srcc_in = "$prefix.fasta";
  $srcc_suf = "crossMatch";
  $returnCommand=&check_srcc($srcc_command, $srcc_in, $srcc_suf);
  $return_file=&wait_for_file_crossMatch("$nd/$prefix.fasta.crossMatch");
  ##### HEREHAVE TO WAIT UNTIL $prefix.crossMatch is generated

  #Joins the cleaned files and remove the trash
  system("cd $nd/; cat $prefix.*.clean > $prefix\_clean.fna; rm $prefix.*.clean; rm $prefix.*.log; rm $prefix.*crossMatch; rm -r SGE");
}else{
  die ("Does not exist file $prefix.fasta\n");
}

my $counts=`grep ">" $nd/$prefix.fasta -c`;
chomp($counts);
system("cd $nd/; echo 'Total\t$counts' >> $prefix\_reads.txt; gzip --best $prefix.fasta;");
$counts=`grep ">" $nd/$prefix\_clean.fna -c`;
chomp($counts);
system("cd $nd; echo 'NoAdapter\t$counts' >> $prefix\_reads.txt;");

# Now is going to map against the PhiX vector
system("cd $nd; fr-hit -a $prefix\_clean.fna -d $home/datafiles/PARFuMS_datafiles/PhiX.fna -o $prefix\_Map_PhiX.txt -c 70 -m 40 -r 0" );
system("cd $nd; $home/scripts/PARFuMS_scripts/FR-Hit_clean_Gautam.pl $prefix\_Map_PhiX.txt NonPhiX_$prefix.txt");
system("cd $nd; $home/scripts/PARFuMS_scripts/GetSequences_inverted.pl NonPhiX_$prefix.txt $prefix\_clean.fna > $prefix\_NoPhiX.fna");


# Now is going to map to the vector and remove it. Is going to use crossmatch so is faster and more sensitive. Is going to remove anything that has a hit within the initial or last bases.

system ("cd $nd/; $home/scripts/PARFuMS_scripts/split_run_check_combine_mod.pl -i $prefix\_NoPhiX.fna -s vector -n 200000 -p 'cross_match INCLUDE_INFILE $home/datafiles/PARFuMS_datafiles/RevVector.fas -gap1_only -minmatch 6 -minscore 10 -gap_init -3 | $home/scripts/PARFuMS_scripts/RemoveVector_Gautam.pl INCLUDE_INFILE'");
$srcc_command = "cross_match INCLUDE_INFILE $home/datafiles/PARFuMS_datafiles/RevVector.fas -gap1_only -minmatch 6 -minscore 10 -gap_init -3 | $home/scripts/PARFuMS_scripts/RemoveVector_Gautam.pl INCLUDE_INFILE";
$srcc_in = "$prefix\_NoPhiX.fna";
$srcc_suf = "vector";
$returnCommand=&check_srcc($srcc_command, $srcc_in, $srcc_suf);
$return_file=&wait_for_file_NoPhiX("$nd/$prefix\_NoPhiX.fna.vector");

system("cd $nd/; cat $prefix*.clean > $prefix\_NoVector.fna; rm $prefix*.clean; rm $prefix*.log; rm $prefix*vector; rm -r SGE;");
$counts=`grep ">" $nd/$prefix\_NoVector.fna -c`;
chomp($counts);
system("cd $nd/; echo 'NoVector\t$counts' >> $prefix\_reads.txt;");

system("cd $nd/; $home/scripts/PARFuMS_scripts/split_run_check_combine_mod.pl -n 10000 -s test -i $prefix\_NoVector.fna -p 'velveth assemble_INCLUDE_INFILE 31 -shortPaired INCLUDE_INFILE; velvetg assemble_INCLUDE_INFILE -cov_cutoff 10 -ins_length 100 -min_contig_lgth 100'");
$srcc_command = "velveth assemble_INCLUDE_INFILE 31 -shortPaired INCLUDE_INFILE; velvetg assemble_INCLUDE_INFILE -cov_cutoff 10 -ins_length 100 -min_contig_lgth 100";
$srcc_in = "$prefix\_NoVector.fna";
$srcc_suf = "test";
$returnCommand=&check_srcc($srcc_command, $srcc_in, $srcc_suf);


# Creates the file $prefix_reads.txt that should have some statistics of the assembly, also the folders assemble_* all to be deleted

$return_file=&wait_for_file_NoVectorTest("$nd/$prefix\_NoVector.fna.test");
##### HERE I HAVE TO WAIT UNTIL $prefix.test is generated


### Join all contigs and filter redundand with cd-hit 90% id
system("cd $nd/; cat assemble_$prefix\_NoVector.fna*/contigs.fa > $prefix.velvet_contigs.fna");
## Test to see contigs have passed initial velvet filters
my $velvetline = `ls -l $nd/$prefix.velvet_contigs.fna`;
my @velvArray = split /\s+/, $velvetline;#++
if ($velvArray[4] == 0){
	die "\n\nNo initial velvet contigs were created. This is likely because nothing passed the user-inputted filters.  PARFuMS ended after the first velvet-run for this sample.\n\n";
}
system("cd $nd/; rm -r assemble_$prefix*; rm $prefix\_NoVector.fna.test; rm -r SGE");
system("cd $nd/; cd-hit-est -i $prefix.velvet_contigs.fna -o $prefix\_cd-hit.fna -g 1 -r 1");
system("cd $nd/; rm *.clstr");

### Map the raw reads to the contigs and remove potential chimeras formed
system("cd $nd/; fr-hit -d $prefix\_cd-hit.fna -o $prefix\_Map.txt -a $prefix\_NoVector.fna -m 30");
system("cd $nd/; $home/scripts/PARFuMS_scripts/FR-Hit_cleanChimera.pl $prefix\_Map.txt $prefix\_cd-hit.fna > $prefix\_NoChimera.fna");
system("cd $nd/; cd-hit-est -i $prefix\_NoChimera.fna -o $prefix\_NoChimera_cd-hit.fna -g 1 -r 1");

### Map the raw reads to the contigs and get unmapped reads round 1 & 2.
system("cd $nd/; fr-hit -d $prefix\_NoChimera_cd-hit.fna -a $prefix\_NoVector.fna -m 30 -o $prefix\_Map2.txt");
system("cd $nd/; $home/scripts/PARFuMS_scripts/FR-Hit_get_unmapped.pl $prefix\_Map2.txt $prefix");
system("cd $nd/; cat $prefix\_NoChimera_cd-hit.fna $prefix\_missing.contigs.fna > $prefix\_ForCD-Hit.fna");
system("cd $nd/; cd-hit-est -i $prefix\_ForCD-Hit.fna -o $prefix\_cd-hit2.fna -g 1 -r 1");
system("cd $nd/; fr-hit -d $prefix\_cd-hit2.fna -a $prefix\_NoVector.fna -m 30 -o $prefix\_Map2a.txt");
system("cd $nd/; $home/scripts/PARFuMS_scripts/FR-Hit_get_unmapped2.pl $prefix\_Map2a.txt $prefix");
system("cd $nd/; cat $prefix\_cd-hit2.fna $prefix\_missing.contigs2.fna > $prefix\_ForCD-Hit2.fna");
system("cd $nd/; makeblastdb -in Adapters_$prefix.fna -dbtype nucl; /srv/cgs/local/ncbi-blast/latest/bin/blastn -db Adapters_$prefix.fna -query $prefix\_ForCD-Hit2.fna -outfmt 6 -dust no -evalue 1e-4 | $home/scripts/PARFuMS_scripts/CleanAdapterContigs.pl $prefix\_ForCD-Hit2.fna > $prefix.ForCD-Hit2.fna");
system("cd $nd/; cd-hit-est -i $prefix.ForCD-Hit2.fna -o $prefix\_cd-hit3.fna -g 1 -r 1");

### Clean potential chimeras and phrap the contigs from the 3 assemblies
system("cd $nd/; fr-hit -d $prefix\_cd-hit3.fna -o $prefix\_Map3.txt -a $prefix\_NoVector.fna -m 30");
system("cd $nd/; $home/scripts/PARFuMS_scripts/FR-Hit_cleanChimera.pl $prefix\_Map3.txt $prefix\_cd-hit3.fna > $prefix\_ForPhrap1.fna");
system("cd $nd/; phrap -minmatch 25 -maxmatch 40 -bandwidth 1 -minscore 30 -penalty -5 -gap_init -4 -gap_ext -3 $prefix\_ForPhrap1.fna");
system("cd $nd/; cat $prefix\_ForPhrap1.fna.contigs $prefix\_ForPhrap1.fna.singlets $prefix\_ForPhrap1.fna.problems > $prefix\_phrapContigs.fna");
system("cd $nd/; rm $prefix\_ForPhrap1.fna.*");

### Map reads to the contigs and clean for potential new chimeras
system("cd $nd/; cd-hit-est -i $prefix\_phrapContigs.fna -o $prefix\_phrap_cd-hit.fna -g 1 -r 1");
system("cd $nd/; rm *.clstr");

### Map reads again and link between contigs, phrap again.
system("cd $nd/; fr-hit -d $prefix\_phrap_cd-hit.fna -o $prefix\_Map4.txt -a $prefix\_NoVector.fna -m 30");
system("cd $nd/; $home/scripts/PARFuMS_scripts/FR-Hit_link_light.pl $prefix\_Map4.txt $prefix\_phrap_cd-hit.fna $prefix; gzip --best $prefix\_clean.fna;");
system("cd $nd/; phrap -minmatch 25 -maxmatch 40 -bandwidth 1 -minscore 30 $prefix\_ForPhrap.fna");
system("cd $nd/; cat $prefix\_ForPhrap.fna.contigs $prefix\_ForPhrap.fna.singlets $prefix\_ForPhrap.fna.problems > $prefix\_phrap.fna");
system("cd $nd/; cd-hit-est -i $prefix\_phrap.fna -o $prefix\_phrap_cd-hit1.fna -g 1 -r 1");
system("cd $nd/; rm *.clstr");
system("cd $nd/; rm $prefix\_ForPhrap.fna.*");

#### Do the stiching by function.
system("cd $nd/; $home/scripts/PARFuMS_scripts/split_run_check_combine_mod.pl -n 80 -s blastCOG -i $prefix\_phrap_cd-hit1.fna -M 6 -p '/srv/cgs/local/ncbi-blast/latest/bin/blastx -outfmt 6 -evalue 1E-4 -seg no -db $home/datafiles/PARFuMS_datafiles/COG_with_len.fa -query INCLUDE_INFILE'");
$srcc_command = "/srv/cgs/local/ncbi-blast/latest/bin/blastx -outfmt 6 -evalue 1E-4 -seg no -db $home/datafiles/PARFuMS_datafiles/COG_with_len.fa -query INCLUDE_INFILE";
$srcc_in = "$prefix\_phrap_cd-hit1.fna";
$srcc_suf = "blastCOG";
$returnCommand=&check_srcc($srcc_command, $srcc_in, $srcc_suf);
$return_file=&wait_for_file_blastCOG("$nd/$prefix\_phrap_cd-hit1.fna.blastCOG");
system("cd $nd/; $home/scripts/PARFuMS_scripts/Preparse_Cog_Blast1_Gautam.pl $home/datafiles/PARFuMS_datafiles/COG.mappings.v8.2.txt $prefix\_phrap_cd-hit1.fna $prefix\_phrap_cd-hit1.fna.blastCOG > $prefix.preparse_COG; $home/scripts/PARFuMS_scripts/Posparse_Cog_Blast1_Gautam.pl $prefix.preparse_COG > $prefix.annotation.txt; ");
system("cd $nd/; makeblastdb -in $prefix\_phrap_cd-hit1.fna -dbtype nucl; $home/scripts/PARFuMS_scripts/split_run_check_combine_mod.pl -n 10000 -s BlastForStich -i $prefix\_NoVector.fna -p '/srv/cgs/local/ncbi-blast/latest/bin/blastn -query INCLUDE_INFILE -db $prefix\_phrap_cd-hit1.fna -outfmt 6 -evalue 1e-7 -dust no'");
$srcc_command = "/srv/cgs/local/ncbi-blast/latest/bin/blastn -query INCLUDE_INFILE -db $prefix\_phrap_cd-hit1.fna -outfmt 6 -evalue 1e-7 -dust no";
$srcc_in = "$prefix\_NoVector.fna";
$srcc_suf = "BlastForStich";
$returnCommand=&check_srcc($srcc_command, $srcc_in, $srcc_suf);
$return_file=&wait_for_file_blastSTITCH("$nd/$prefix\_NoVector.fna.BlastForStich");
system("cd $nd/; $home/scripts/PARFuMS_scripts/MakeGraph_Gautam.pl $prefix\_phrap_cd-hit1.fna $prefix\_NoVector.fna.BlastForStich > $prefix.graph.dot");
system("cd $nd/; $home/scripts/PARFuMS_scripts/Stich_Contigs_Gautam.pl $prefix\_phrap_cd-hit1.fna $prefix.graph.dot $prefix.annotation.txt > $prefix.stich.fna");
system("cd $nd/; cd-hit-est -i $prefix.stich.fna -o $prefix.cd-hit-last.fna -g 1 -r 1");

##### Rename Final contigs adding meadian coverage and length
system("cd $nd/; fr-hit -d $prefix.cd-hit-last.fna -a $prefix\_NoVector.fna -m 30 -o $prefix\_Map5.txt");
system("cd $nd/; $home/scripts/PARFuMS_scripts/FR-Hit_RenameFinalContigs.pl $prefix\_Map5.txt $prefix.cd-hit-last.fna  > $prefix.lastContigs.fna");

##### Blast against COG to do functional annotation
system("cd $nd/; /srv/cgs/local/ncbi-blast/latest/bin/blastx -outfmt 6 -evalue 1E-4 -seg no -db $home/datafiles/PARFuMS_datafiles/COG_with_len.fa -query $prefix.lastContigs.fna -out $prefix.lastContigs.fna.blastCOG_Final");
system("cd $nd/; $home/scripts/PARFuMS_scripts/Preparse_Cog_Blast1_Gautam.pl $home/datafiles/PARFuMS_datafiles/COG.mappings.v8.2.txt $prefix.lastContigs.fna $prefix.lastContigs.fna.blastCOG_Final > $prefix.preparse_COG.final.txt; $home/scripts/PARFuMS_scripts/Posparse_Cog_Blast1_Gautam.pl $prefix.preparse_COG.final.txt > $prefix.annotation.final.txt; ");

##### Molly Gibson's Python-based GeneFinder + HMM-annotater, supplementing COG-based annotation #####
system ("cd $nd; /home/gdlab/mgibson/fmg_tools/scripts/annotate_functional_selections.py -contigs $prefix.lastContigs.fna --resfams -o $prefix.annotation -f");

#### Now do some cleaning. In particular big Blast files
system("cd $nd/; rm -r SGE *.nhr *.nin *.nsq *.clstr");
#system("cd $nd/; rm $prefix.clean2.* $prefix.FirstRoundMapped.txt $prefix.ForSplit.blastn $prefix.outLinkAssembler.txt $prefix.linkedContigs.fna.blastCOG $prefix.preparse_COG");

if ($type =~ m/a/){
	print "\n\nYou have either chosen your type as ANTIBIOTIC RESISTANCE, or it was chosen by default.\nAll annotations will now be classified into one of several resistance-related categories, according to the antibiotic-resistance keyword file in the PARFuMS datafile-directory, and outputted to an additional tab-delimited file.\n*Pro-tip: For easier reading, view the file in excel*\n\n";
#	system ("cd $nd/; python $home/scripts/PARFuMS_scripts/ProcessingPARFuMSData-Resistance.py $prefix.annotation.final.txt $prefix.lastContigs.fna $prefix");
}
elsif ($type =~ m/b/){
	print "\n\nYou would like to have annotations classified into BIOFUEL-RELATED categories, and PARFuMS would like to fulfill your desire. Unfortunately, no biofuel-related keyword file exists in the PARFuMS datafile-directory, and your wish will be ignored. No additional analysis will be performed.\n\n";
}
elsif ($type =~ m/c/){
	print "\n\nYou would like to have annotations classified into CARBOHYDRATE-RELATED categories, and PARFuMS would like to fulfill your desire. Unfortunately, no biofuel-related keyword file exists in the PARFuMS datafile-directory, and your wish will be ignored. No additional analysis will be performed.\n\n";
}
else{
	print "\n\nPARFuMS cannot figure out a type of functional selection performed, in its confusion it refuses to classify any annotations into meaningful categories. No additional analysis will be performed\n";
}

if ($loop == 1){
	system ("cd $nd/; gunzip $prefix.fasta.gz");
	system ("cd $nd/; gunzip $prefix.fasta.original.gz");
	system ("cd $nd/; mv $prefix.fasta $tempFasta");
	system ("cd $nd/; mv $prefix.fasta.original $prefix.fasta");
	system ("cd $nd/; gzip $prefix.fasta");
	system ("cd $nd/; gzip $tempFasta");
}

#Running Mira to prouduce visualization of contigs/reads for easy hand-checking of assembled contigs
#my $gunzipfile = $prefix."_clean.fna.gz";
#my $cleanFile = $prefix."_clean.fna";
#my $MiraIn = $prefix."_in.solexa.fasta";
#my $MiraBack = $prefix."_backbone_in.fasta";
#my $MiraOut = $prefix."_mapping";
#my $MiraLog = $prefix."_mappingLog.txt";
#my $f1 = $MiraOut."_assembly";
#my $f2 = $MiraOut."_d_results";
#my $f3 = $MiraOut."_out.ace";
#system ("cd $nd/; gunzip $gunzipfile");
#system ("cd $nd/; mv $cleanFile $MiraIn");
#system ("cd $nd/; cp $prefix.lastContigs.fna $MiraBack");
#system ("cd $nd/; mira -projectin=$prefix -projectout=$MiraOut --job=mapping,genome,accurate,solexa  -fasta -SB:bft=fasta -AS:urd=no -MI:sonfs=no SOLEXA_SETTINGS -LR:wqf=no:mxti=no -AS:epoq=no  > $MiraLog");
#$return_file=&wait_for_file("$nd/$f1/$f2/$f3");
#system ("cd $nd/; mv $f1/$f2/$f3 $f3");
#system ("cd $nd/; mv $MiraIn $cleanFile");
#system ("cd $nd/; rm $MiraBack");
#system ("cd $nd/; gzip $cleanFile");




#Checking to see if there were any errors when the alignment was run
my $runSum = $prefix."_runSummary.txt";
my $error1 = `grep ERROR $runSum | wc -l`;
my $error2 = `grep error $runSum | wc -l`;
my $error3 = `grep Error $runSum | wc -l`;
my $error4 = `grep Illegal $runSum | wc -l`;
my $error5 = `grep illegal $runSum | wc -l`;
my $error6 = `grep ILLEGAL $runSum | wc -l`;
my $error7 = `grep "ALIGNMENT ERROR" $runSum | wc -l`;
my $work = `grep Workaround $runSum | wc -l`;
chomp $error1;
chomp $error2;
chomp $error3;
chomp $error4;
chomp $error5;
chomp $error6;
chomp $error7;
chomp $work;
my $numErrors = $error1 + $error2 + $error3 + $error4 + $error5 + $error6;
if ($numErrors > 0){
	open (OUT, ">$prefix.HAS_AN_ERROR");
	print OUT "A total of $numErrors errors are present in $prefix\n\n";
	print OUT "$error1 of them were found due to the keyword ERROR (case-sensitive)\n$error2 of them were found due to the keyword error (case-sensitive)\n$error3 of them were found due to the keyword Error (case-sensitive)\n$error4 of them were found due to the keyword Illegal (case-sensitive)\n$error5 of them were found due to the keyword illegal (case-sensitive)\n$error6 of them were found due to the keyword ILLEGAL (case-sensitive)\n";
	print OUT "\nOf the $numErrors errors, $error7 are 'Alignment Errors' with particular reads in Phrap. As far as I can tell (Kevin's opinion, not gospel), these errors don't seem to have major effects...especially if they are relatively low in overall number\n";
	print OUT "Use the grep command to view errors in the runSummary file\n";
}
if ($work > 0){
	open (OUT2, ">$prefix.INVOKED_A_WORKAROUND");
	print OUT2  "A total of $work workarounds were invoked, use grep and the keyword 'Workaround' to view these invocations\n";
}



################### Subroutines #####################

sub wait_for_file{
  my($fileName) = @_;

  # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

  print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

  # Check for the finished file every WAIT_TIME minutes
  until (-e "$fileName") {
    print "$fileName does not exist yet at ".localtime()."\n";
    sleep WAIT_TIME;
  }

  print "File $fileName found at ".localtime()."\n";

  # make sure the file is no longer being modified (changing size)
  my$time = 0;
  my$size = -s $fileName;

  until ($time == NO_CHANGE_IN_FILE_FOR_X_SECONDS) {
    if ($size == -s $fileName) {
      $time++;
      sleep 1;
      if ($time%5 == 0) {
	print "No change in file size for $time seconds\n";
      }
    } else {
      $time = 0;
      $size = -s $fileName;
      print "The file size changed, sleeping for 1 minute\t";
      sleep 60;
      print"Waking up, try again\n";
    }
  }
  print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
  return $fileName;
}

sub wait_for_file_blastCOG{
  my($fileName) = @_;

  # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

  print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

  # Check for the finished file every WAIT_TIME minutes
  until ((-e "$fileName") || $x > 13) {
    print "$fileName does not exist yet at ".localtime()."\n";
    $x++;
    sleep WAIT_TIME;
  }
 if ($x <= 13){
  print "File $fileName found at ".localtime()."\n";

  # make sure the file is no longer being modified (changing size)
  my$time = 0;
  my$size = -s $fileName;

  until ($time == NO_CHANGE_IN_FILE_FOR_X_SECONDS) {
    if ($size == -s $fileName) {
      $time++;
      sleep 1;
      if ($time%5 == 0) {
	print "No change in file size for $time seconds\n";
      }
    } else {
      $time = 0;
      $size = -s $fileName;
      print "The file size changed, sleeping for 1 minute\t";
      sleep 60;
      print"Waking up, try again\n";
    }
    $x = 0;
  }
  print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
  return $fileName;
 }
 else{
 	my $y = (($x * WAIT_TIME)/60);
 	print "\nWorkaround: The blastCOG File $fileName was not found after $y minutes of wait-time. There was likely an problem during PARFuMS with split_run_check_combine_mod.pl. We will try running the blast command without first splitting the file...\n\n";
 	system("cd $nd/; /srv/cgs/local/ncbi-blast/latest/bin/blastx -db $home/datafiles/PARFuMS_datafiles/COG_with_len.fa -query $prefix\_phrap_cd-hit1.fna -evalue 1E-4 -out $prefix\_phrap_cd-hit1.fna\.blastOUT -outfmt 6 -seg no");
 	system("cd $nd/; mv $prefix\_phrap_cd-hit1.fna\.blastOUT $prefix\_phrap_cd-hit1.fna.blastCOG");
 	system("cd $nd/; rm $prefix\_phrap_cd-hit1.fna.c*");
 	$x = 0;
 }
}

sub wait_for_file_NoVectorTest{
  my($fileName) = @_;

  # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

  print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

  # Check for the finished file every WAIT_TIME minutes
  until ((-e "$fileName") || $x > 13) {
    print "$fileName does not exist yet at ".localtime()."\n";
    $x++;
    sleep WAIT_TIME;
  }
 if ($x <= 13){
  print "File $fileName found at ".localtime()."\n";

  # make sure the file is no longer being modified (changing size)
  my$time = 0;
  my$size = -s $fileName;

  until ($time == NO_CHANGE_IN_FILE_FOR_X_SECONDS) {
    if ($size == -s $fileName) {
      $time++;
      sleep 1;
      if ($time%5 == 0) {
	print "No change in file size for $time seconds\n";
      }
    } else {
      $time = 0;
      $size = -s $fileName;
      print "The file size changed, sleeping for 1 minute\t";
      sleep 60;
      print"Waking up, try again\n";
    }
    $x = 0;
  }
  print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
  return $fileName;
 }
 else{
 	my $y = (($x * WAIT_TIME)/60);
 	my $summaryFile = $prefix."_runSummary.txt";
 	my $totalJobs = 0;
 	my $failedJobs = 0;
 	my $c_num = "";
 	my $z = 0;
 	open (SUMMARY, $summaryFile) || print "Could not open $summaryFile. You'll have to open the summary and check the % reads lost manually...sorry brah!";
 		while (my $line = <SUMMARY>){
			chomp $line;
			if ($line =~ m/sequence $prefix\_NoVector.fna/){
				my @lineArray = split /\s+/, $line;
				$totalJobs = $lineArray[5];
				my $nline = <SUMMARY>;
				my @arr1 = split /\s+/, $nline;
				my @arr2 = split /\.sh/, $arr1[1];
				$c_num = $arr2[0];
				last;
			}
 		}
 	close SUMMARY;
 	my $c_file = $prefix."/"."$c_num".".check_status";
 	$failedJobs = `grep -c "complete not found" $c_file`;
 	chomp $failedJobs;
 	my $percentageFailed = 0;
 	if ($failedJobs =~ m/^\d/){
 		if ($failedJobs =~ m//){
 			$percentageFailed = ($failedJobs/$totalJobs)*100;
 		}
 	}
 	else{
 		my $com_filePresent = 0;
 		for (my$i=0;$i<$totalJobs;$i++){
 			my $com_file = $prefix."/"."$c_num"."\."."$i"."\."."complete";
 			my $temp = `ls $com_file`;
 			chomp $temp;
 			if ($temp =~ m/^$prefix/){
 				$com_filePresent++;
 			}
 		}
 		$failedJobs = $totalJobs - $com_filePresent;
 		if ($failedJobs =~ m//){
 			$percentageFailed = ($failedJobs/$totalJobs)*100;
 		}
 	}
 	my $NoVectorTestSkipFile = $fileName;
 	open (OUT, ">$NoVectorTestSkipFile") or die ("Couldn't open file: $NoVectorTestSkipFile\n");
 	print OUT "Temp Text To Initialize File";
 	print "\nWorkaround: The File $fileName was not found after $y minutes of wait-time. There was likely an problem during PARFuMS with split_run_check_combine_mod.pl. $failedJobs of $totalJobs jobs failed, or $percentageFailed\% of the reads, were not considered in the assembly\n\n";
 	system("cd $nd/; rm $prefix\_NoVector.fna.c*");
 	$x = 0;
 }
}

sub wait_for_file_blastSTITCH{
  my($fileName) = @_;

  # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

  print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

  # Check for the finished file every WAIT_TIME minutes
  until ((-e "$fileName") || $x > 13) {
    print "$fileName does not exist yet at ".localtime()."\n";
    $x++;
    sleep WAIT_TIME;
  }
 if ($x <= 13){
  print "File $fileName found at ".localtime()."\n";

  # make sure the file is no longer being modified (changing size)
  my$time = 0;
  my$size = -s $fileName;

  until ($time == NO_CHANGE_IN_FILE_FOR_X_SECONDS) {
    if ($size == -s $fileName) {
      $time++;
      sleep 1;
      if ($time%5 == 0) {
	print "No change in file size for $time seconds\n";
      }
    } else {
      $time = 0;
      $size = -s $fileName;
      print "The file size changed, sleeping for 1 minute\t";
      sleep 60;
      print"Waking up, try again\n";
    }
    $x = 0;
  }
  print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
  return $fileName;
 }
 else{
 	my $y = (($x * WAIT_TIME)/60);
 	print "\nWorkaround: The blastSTITCH File $fileName was not found after $y minutes of wait-time. There was likely an problem during PARFuMS with split_run_check_combine_mod.pl. We will try running the blast command without first splitting the file...\n\n";
 	system("cd $nd/; /srv/cgs/local/ncbi-blast/latest/bin/blastn -db $prefix\_phrap_cd-hit1.fna -query $prefix\_NoVector.fna -evalue 1e-7 -out $prefix\_NoVector.fna\.blastSTICHOUT -outfmt 6 -dust no");
 	system("cd $nd/; mv $prefix\_NoVector.fna\.blastSTICHOUT $prefix\_NoVector.fna.BlastForStich");
 	system("cd $nd/; rm $prefix\_NoVector.fna.c*");
 	$x = 0;
 }
}

sub wait_for_file_crossMatch{
  my($fileName) = @_;

  # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

  print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

  # Check for the finished file every WAIT_TIME minutes
  until ((-e "$fileName") || $x > 13) {
    print "$fileName does not exist yet at ".localtime()."\n";
    $x++;
    sleep WAIT_TIME;
  }
  if ($x <= 13){
  print "File $fileName found at ".localtime()."\n";

  # make sure the file is no longer being modified (changing size)
  my$time = 0;
  my$size = -s $fileName;

  until ($time == NO_CHANGE_IN_FILE_FOR_X_SECONDS) {
    if ($size == -s $fileName) {
      $time++;
      sleep 1;
      if ($time%5 == 0) {
	print "No change in file size for $time seconds\n";
      }
    } else {
      $time = 0;
      $size = -s $fileName;
      print "The file size changed, sleeping for 1 minute\t";
      sleep 60;
      print"Waking up, try again\n";
    }
    $x = 0;
  }
  print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
  return $fileName;
  }
  else{
 	my $y = (($x * WAIT_TIME)/60);
	print "\nWorkaround: The crossmatch File $fileName was not found after $y minutes of wait-time. There was likely an problem during PARFuMS with split_run_check_combine_mod.pl. We will try running the crossmatch command without first splitting the file...\n\n";
 	system("cd $nd/; cross_match $prefix.fasta Adapters_$prefix.fna -gap1_only -minmatch 6 -minscore 10 -gap_init -3 | $home/scripts/PARFuMS_scripts/CleanAdapter-Illumina_PE_mod.pl $prefix.fasta");
	my $summaryFile = $prefix."_runSummary.txt";
	my $fail = 0;
 	open (SUMMARY, $summaryFile) || print "Could not open $summaryFile. You'll have to open the summary and check the % reads lost manually...sorry broseph!";
 		while (my $line = <SUMMARY>){
			chomp $line;
			if ($line =~ m/FATAL ERROR/){
				$fail = 1;
				last;
			}
			elsif (($line =~ m/Printing matches/) && ($line =~ m/Sorting pairs/)){
				$fail = 2;
			}
 		}
 		if ($fail==1){
 			print "Re-running cross-match without splitting took up too much memory, and therefore failed. This strategy did not work, and we will now take the portion of reads that passed the initial split_run_check_run_combine.pl cross-match step, and move forward...\n\n";
 			my $initialReads = `grep -c "_1" $prefix\/$prefix.fasta`;
 			chomp $initialReads;
 			my $consideredFastas = `ls $prefix\/*fasta.*.clean`;
 			chomp $consideredFastas;
 			my @conArray = split /\s+/, $consideredFastas;
 			my @usedFastas;
 			foreach my $usedFasta (@conArray){
 				chop $usedFasta;
 				chop $usedFasta;
 				chop $usedFasta;
 				chop $usedFasta;
 				chop $usedFasta;
 				chop $usedFasta;
 				my @usedArray = split /\//, $usedFasta;
 				$usedFasta = $usedArray[1];
 				push (@usedFastas, $usedFasta);
 			}
 			my $catCommand = join(' ', @usedFastas);
 			chomp $catCommand;
 			system ("cd $nd/; cat $catCommand > $prefix.fasta.consideredReads.tmp");
 			system("cd $nd/; mv $prefix.fasta.consideredReads.tmp $prefix.fasta.consideredReads");
 			my $consideredReads = `grep -c "_1" $prefix\/$prefix.fasta.consideredReads`;
 			chomp $consideredReads;
 			my $perc = (($consideredReads / $initialReads) * 100);
 			print "Due to a split_run_check problem, and the fact you had many initial reads, only $consideredReads reads of a total of $initialReads reads successfully made it to cross_match to be considered for adapter-removal, $perc\% of the reads. \n\n";
 			system("cd $nd/; rm $prefix.fasta.clean");
 			system("cd $nd/; rm $prefix.fasta.consideredReads");
			open(FH,">$fileName") or die "Can't create $fileName\n";
 			print FH "Temp Text To Initialize File";
 			close (FH);
 		}
 		elsif ($fail==2){
 			print "Re-running cross-match without splitting worked!\n\n";
 			system("cd $nd/; rm $prefix.fasta.c*clean");
			open(FH,">$fileName") or die "Can't create $fileName\n";
 			print FH "Temp Text To Initialize File";
 			close (FH);
 		}
 		else{
 			system("cd $nd/; rm $prefix.fasta.clean");
 			die "split_run_check_combine_mod.pl failed at the initial cross_match step, the sample will need to be re-run. Oy Vey!\n\n";
 		}
 close (SUMMARY);
 	$x = 0;
 }
}

sub wait_for_file_NoPhiX{
  my($fileName) = @_;

  # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

  print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

  # Check for the finished file every WAIT_TIME minutes
  until ((-e "$fileName") || $x > 13) {
    print "$fileName does not exist yet at ".localtime()."\n";
    $x++;
    sleep WAIT_TIME;
  }
  if ($x <= 13){
  print "File $fileName found at ".localtime()."\n";

  # make sure the file is no longer being modified (changing size)
  my$time = 0;
  my$size = -s $fileName;

  until ($time == NO_CHANGE_IN_FILE_FOR_X_SECONDS) {
    if ($size == -s $fileName) {
      $time++;
      sleep 1;
      if ($time%5 == 0) {
	print "No change in file size for $time seconds\n";
      }
    } else {
      $time = 0;
      $size = -s $fileName;
      print "The file size changed, sleeping for 1 minute\t";
      sleep 60;
      print"Waking up, try again\n";
    }
    $x = 0;
  }
  print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
  return $fileName;
  }
  else{
 	my $y = (($x * WAIT_TIME)/60);
	print "\nWorkaround: The crossmatch File $fileName was not found after $y minutes of wait-time. There was likely an problem during PARFuMS with split_run_check_combine_mod.pl. We will try running the crossmatch command without first splitting the file...\n\n";
 	system("cd $nd/; cross_match $prefix\_NoPhiX.fna $home/datafiles/PARFuMS_datafiles/RevVector.fas -gap1_only -minmatch 6 -minscore 10 -gap_init -3 | $home/scripts/PARFuMS_scripts/RemoveVector_Gautam.pl $prefix\_NoPhiX.fna");
	my $summaryFile = $prefix."_runSummary.txt";
	my $fail = 0;
 	open (SUMMARY, $summaryFile) || print "Could not open $summaryFile. You'll have to open the summary and check the % reads lost manually...sorry broseph!";
 		while (my $line = <SUMMARY>){
			chomp $line;
			if ($line =~ m/FATAL ERROR/){
				$fail = 1;
				last;
			}
			elsif (($line =~ m/Printing matches/) && ($line =~ m/Sorting pairs/)){
				$fail = 2;
			}
 		}
 		if ($fail==1){
 			print "Re-running cross-match without splitting took up too much memory, and therefore failed. This strategy did not work, and we will now take the portion of reads that passed the initial split_run_check_run_combine.pl cross-match step, and move forward...\n\n";
 			my $initialReads = `grep -c "_1" $prefix\/$prefix\_NoPhiX.fna`;
 			chomp $initialReads;
 			my $consideredFastas = `ls $prefix\/$prefix\_NoPhiX.c*.clean`;
 			chomp $consideredFastas;
 			my @conArray = split /\s+/, $consideredFastas;
 			my @usedFastas;
 			foreach my $usedFasta (@conArray){
 				chop $usedFasta;
 				chop $usedFasta;
 				chop $usedFasta;
 				chop $usedFasta;
 				chop $usedFasta;
 				chop $usedFasta;
 				my @usedArray = split /\//, $usedFasta;
 				$usedFasta = $usedArray[1];
 				push (@usedFastas, $usedFasta);
 			}
 			my $catCommand = join(' ', @usedFastas);
 			chomp $catCommand;
 			system ("cd $nd/; cat $catCommand > $prefix\_NoPhix.fna.consideredReads.tmp");
 			system("cd $nd/; mv $prefix\_NoPhix.fna.consideredReads.tmp $prefix\_NoPhix.fna.consideredReads");
 			my $consideredReads = `grep -c "_1" $prefix\/$prefix\_NoPhix.fna.consideredReads`;
 			chomp $consideredReads;
 			my $perc = (($consideredReads / $initialReads) * 100);
 			print "Due to a split_run_check problem, and the fact you had many reads in $prefix\_NoPhiX.fna, only $consideredReads reads of a total of $initialReads reads successfully made it to cross_match to be considered for vector-removal, $perc\% of the reads. \n\n";
 			system("cd $nd/; rm $prefix\_NoPhiX.fna.clean");
 			system("cd $nd/; rm $prefix\_NoPhix.fna.consideredReads");
			open(FH,">$fileName") or die "Can't create $fileName\n";
 			print FH "Temp Text To Initialize File";
 			close (FH);
 		}
 		elsif ($fail==2){
 			print "Re-running cross-match without splitting worked!\n\n";
 			system("cd $nd/; rm $prefix\_NoPhiX.c*.clean");
			open(FH,">$fileName") or die "Can't create $fileName\n";
 			print FH "Temp Text To Initialize File";
 			close (FH);
 		}
 		else{
 			system("cd $nd/; rm $prefix\_NoPhiX.fna.clean");
 			die "split_run_check_combine_mod.pl failed at the Vector-Removal cross_match step, the sample will need to be re-run. Oy Vey!\n\n";
 		}
 close (SUMMARY);
 	$x = 0;
 }
}

sub check_srcc{
	my($srcc_com, $srcc_in, $srcc_suff) = @_;
	my $summaryFile = $prefix."_runSummary.txt";
	my $newCommand = "---";
	open (SUMMARY, $summaryFile) || print "Could not open $summaryFile. You need to have a properly named summary file for PARFuMS to work correctly, this run probably will fail\n\n";
	while (my $line = <SUMMARY>){
		chomp $line;
		last;
	}
	close (SUMMARY);
	my $tailText = `tail $summaryFile`;
		if ($tailText =~ m/PARFuMS is checking for split_run_check_combine submissions of only one job/){
			$newCommand = $srcc_com;
			$newCommand =~ s/INCLUDE_INFILE/$srcc_in/g;
			my $tempfile = "$prefix\/"."$srcc_in"."."."$srcc_suff";
			if ($newCommand =~ m/blast/i){
				my @temArray = split /\//, $tempfile;
				my $tempCommand = $newCommand." -out $temArray[1]";
				$newCommand = $tempCommand;
			}
			else{
			open (FH, ">$tempfile");
			print FH "temp text to initialize temp file\n";
			close (FH);
			}
			print "\n\nPARFuMS is directly re-submitting the following command:\n\n$newCommand\n\n";
			system ("cd $nd/; $newCommand");
		}
return $newCommand;
}
