#!/usr/bin/perl
# taken from DK 111031
# modified by Ray Wang for making pdb2vall.py
# modified by YoshitakaM for BLAST+ and psiblast 20181004

use FindBin qw($Bin);

my $psiblast_num_cpus = 8;

# possible blastpgp locations
my @psiblast = (
	"$Bin/../../../../../../src/blast/bin/blastpgp", # Robetta location
	"$Bin/../../blast/bin/blastpgp", # fragment_tools location
	"$Bin/../pdb_scripts/blastpgp"	# pdb2vall location
);

# possible nr database locations
my @nr = (
	"$Bin/../../../../../databases/nr/nr_pfilt", # Robetta location
	"$Bin/../../databases/nr_pfilt", # fragment_tools location
	"$Bin/../database/ncbi/nr_pfilt", # pdb2vall location
	"$Bin/../database/ncbi/nr_segfilt_x",
	"$Bin/../../../../../databases/nr/nr",
	"$Bin/../database/ncbi/nr",
	"$Bin/../../databases/nr"
);

# find blastpgp
my $psiblast;
foreach my $b (@psiblast) {
	if (-s $b) {
		$psiblast = $b;
		last;
	}
}

# find nr database
my $DB;
foreach my $d (@nr) {
	if (-s "$d.pal") {
		$DB = $d;
		last;
	}
}

# read config
open(F, "$Bin/../blastplus_pdb2vall.cfg") or die "ERROR! cannot open config file $Bin/../blastplus_pdb2vall.cfg: $!\n";
while (<F>) {
  if (/^\s*psiblast\s*=\s*(\S+)\s*/) {
    $psiblast = $1;
  } elsif (/^\s*psiblast_num_cpus\s*=\s*(\d+)\s*/) {
    $psiblast_num_cpus = $1;
  } elsif (/^\s*nr\s*=\s*(\S+)\s*/) {
    $DB = $1;
  }
}
close(F);

(-s $psiblast) or die "ERROR! BLAST+ psiblast $psiblast does not exist\n";

my $psiblast = "$psiblast -num_threads $psiblast_num_cpus";

# (-s $DB) or die "ERROR! NR database $DB does not exist\n";


my $DEBUG = 0;

my $round1       = '0.0000000001';
my $round2       = '0.00001';

if (scalar@ARGV != 1) {
  print "USAGE: $0 <fasta>\n";
  exit(1);
}

my ($fasta) = @ARGV;
chomp $fasta;

my $fastaname = $fasta;
$fastaname =~ s/^.*\/([^\/]+)\s*$/$1/gs;
my $fastaprefix = $fastaname;
$fastaprefix =~ s/\.fasta$//;

# get the sequence from the fasta file
my $sequence;
open(SEQFILE, "$fasta") or die "ERROR! cannot open $fasta: $!\n";
my $has_comment = 0;
my $has_eof     = 0;
my $eof;
while (<SEQFILE>) {
  $eof = $_;
  s/\s//g;
  if (/^>/) { $has_comment = 1; next; }
  chomp;
  $sequence .= $_;
}
close(SEQFILE);
$has_eof = 1 if ($eof =~ /\n$/);
($has_comment && $has_eof) or die "fasta file must have a comment and end with a new line!\n";
(!$DEBUG) || print "Sequence: $sequence\n";

if (-s "$fastaprefix.checkpoint") {
  print "$fastaprefix.checkpoint exists\n";
  print "Done!\n";
  exit(0);
}

my $shell;

if (!-s "$fastaprefix.1.chk") {

	#sleep int(rand ( $sleep_secs )) + 1;

# -comp_based_stats 1   # 1 : Composition-based statistics as in NAR  29:2994--3005, 2001
# -num_iterations 2     # Number of iterations to perform (0 means run until convergence)
# -inclusion_ethresh    # E-value inclusion threshold for pairwise alignments
# -evalue               # Expectation value (E) threshold for saving hits
# -num_alignments       # Number of database sequence to show alignments for (B)
  $shell = "$psiblast -comp_based_stats 1 -num_iterations 2 -inclusion_ethresh $round1 -evalue $round1 -num_alignments 0 -db $DB ";
  $shell .= " -out $fastaprefix.1.blast -out_pssm $fastaprefix.1.chk -out_ascii_pssm $fastaprefix.1.matrix ";
  $shell .= " -query $fastaname";

	#print "rundir: $rundir\n";
  print "$shell\n";
  system($shell);
}

if (!-s "$fastaprefix.2.chk") {

	#sleep int(rand ( $sleep_secs )) + 1;

  $shell = "$psiblast -comp_based_stats 1 -num_iterations 2 -inclusion_ethresh $round2 -evalue $round2 -num_alignments 0 -db $DB ";
  $shell .= " -out $fastaprefix.2.blast -out_pssm $fastaprefix.2.chk -out_ascii_pssm $fastaprefix.2.matrix ";
  $shell .= " -in_pssm $fastaprefix.1.chk";

  print "$shell\n";
  system($shell);
}


# parse the checkpoint matrix using sequence_profile_scripts/fasta.py.
system("python3 $Bin/convert_asciichk_to_checkpoint.py -p ${fastaprefix}");

print "done!\n";
