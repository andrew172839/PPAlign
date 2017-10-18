#!/usr/bin/perl
use strict;

my %RES = ('ALA' => 'A', 'VAL' => 'V', 'PHE' => 'F', 'ILE' => 'I', 'LEU' => 'L', 'ASP' => 'D', 'GLU' => 'E', 'LYS' => 'K', 'SER' => 'S', 'THR' => 'T', 'TYR' => 'Y', 'CYS' => 'C', 'ASN' => 'N', 'GLN' => 'Q', 'PRO' => 'P', 'MET' => 'M', 'ARG' => 'R', 'HIS' => 'H', 'TRP' => 'W', 'GLY' => 'G', 'MSE' => 'M');
my $seq = $ARGV[0];
my $filelist = $ARGV[1];
my $Wprof = 0.;
my $DB = '~/PPAlign';
open FILE, "<$filelist";
open LOG, ">$seq.log";
my @lists = <FILE>;
close FILE;
chomp @lists;
my $isblast = 1;
my $temp_freq = $seq.'.temp.freq';
chomp $seq;
if ($isblast == 1) {
	print LOG "start running psipred\n";
	&psipred($seq);
	&prepare_seq($seq);
}
my $list;
my @output;
my %out;
open SEQ, "<$seq.seq";
my @fs = <SEQ>;
close SEQ;
chomp @fs;
shift @fs;
my $tt = join('', @fs);
my $seqlen = length($tt);
my (%r_l, %r_e);
@output = `./threadall.bin $filelist $seq 4.5 0.5 0.5 3 $Wprof`;
chomp @output;
my ($tl, $energy, $align_len);
foreach (@output) {
	s/^\s+//;
	($tl, $energy, $align_len) = split;
	next if ($align_len <= 0);
	my $en2 = $energy / $align_len;
	$r_l {$tl} = $align_len;
	$r_e {$tl} = $energy;
	print LOG "$tl\t$energy\t$align_len\n";
}
my @sort_list = sort {$r_e {$a} <=> $r_e {$b}} keys(%r_e);
open OUT, ">$seq.result";
foreach $list (@sort_list) {
	print OUT "$list\t$r_l{$list}\t$r_e{$list}\n";
}
close OUT;
my @model_name = ();
for (0..9) {
	my $ii = $_ + 1;
	my $best = $sort_list[$_];
	push @model_name, $best;
	print LOG "best is $best\n";
	&gen_model($best, $seq);
}
my $i;
my @energy_rank = `ls $seq-*.pdb | ./wupo.bin | sort -n -k2`;
open MOD, ">$seq.model.list";
print MOD @energy_rank;
close MOD;
foreach (@energy_rank) {
	chomp;
	my ($tl_pdb, $tl_en) = split;
	push @model_name, $tl_pdb;
}
for $i (1..10) {
	my $i_pos = index($model_name[$i - 1], '.pdb');
	my $parent_id = substr($model_name[$i - 1], $i_pos - 5, 5);
	open FPDB, ">$seq"."_model_00$i.pdb";
	open IN, "<$model_name[$i-1]";
	print FPDB "PFRMAT TS\n";
	print FPDB "TARGET $seq\n";
	print FPDB "AUTHOR 1775-8745-8849\n";
	print FPDB "METHOD OPUS\n";
	print FPDB "MODEL $i\n";
	print FPDB "PARENT $parent_id\n";
	while (<IN>) {
		chomp;
		if ($_ =~ /^ATOM/) {
			substr($_, 56, 4) = '1.00';
			substr($_, 66, length($_) - 66) = '';
			print FPDB "$_\n";
		}
	}
	print FPDB "TER\n";
	print FPDB "END\n";
	close FPDB;
	close IN;
}
close LOG;

sub fastaprint() {
	chomp;
	my $seq = $_[0];
	my $name = $_[1];
	my $fn = $_[2];
	my $i = 0;
	my $line;
	open FILE, ">$fn";
	print FILE ">$name\n";
	while (($i * 80 +80) <= length($seq)) {
		$line = substr($seq, $i * 80, 80);
		print FILE "$line\n";
		$i++;
	}
	my $a = length($seq) - $i * 80 ;
	if ($a != 0) {
		$line = substr($seq, $i * 80, $a);
		print FILE "$line\n";
	}
	close FILE;
}

sub getseq() {
	my $pdb = substr($_[0], 0, 4);
	my $chain = substr($_[0], length($_[0]) - 1, 1);
	my $seq = ();
	$pdb = "$pdb.pdb" if (length($pdb) == 4);
	$chain = "\U$chain";
	$pdb = "\L$pdb";
	my $lib = $DB.'/pdb';
	open FILE, "<$lib/$pdb";
	while (<FILE>) {
		my $head = substr($_, 0, 6);
		my $chid = substr($_, 21, 1);
		if ($head eq 'ATOM' || $head eq 'HETATM') {
			if (($chain eq  '_') || ($chain eq ' ') || ( $chain eq '')) {
				$chain = $chid;
			}
			if (($chain eq $chid) && (substr($_, 13, 3) eq 'CA '))  {
				my $res = substr($_, 17, 3);
				$seq . =  $RES {$res};
			}
		}
		elsif ($head eq 'TER') {
			if (($chain eq  '_') || ($chain eq ' ') || ( $chain eq '')) {
				$chain = $chid;
			}
			if ($chain eq $chid)  {
				last;
			}
		}
	}
	close FILE;
	return $seq;
}

sub prepare_seq() {
	my $seq = $_[0];
	open IN, "<$temp_freq";
	open OUT, ">$seq.freq";
	my $line;
	while ($line = <IN>) {
		print OUT $line if (substr($line, 0, 5) =~ /\s*[0-9]+/);
	}
	close OUT;
	close IN;
	open IN, "<$seq.horiz";
	my @lines = <IN>;
	close IN;
	open OUT, ">$seq.ss";
	my ($line, $conf, $pred, $AA) = ("", "", "", "");
	my $len;
	foreach $line (@lines) {
		chomp $line;
		$len = length($line);
		$conf = $conf.substr($line, 6, $len) if (substr($line, 0, 4) eq 'Conf');
		$pred = $pred.substr($line, 6, $len) if (substr($line, 0, 4) eq 'Pred');
		$AA = $AA.substr($line, 6, $len) if (substr($line, 0, 4) eq 'AA');
	}
	my $a;
	$len = length($pred);
	print OUT "$len\n";
	for (my $i = 0; $i < $len; $i++) {
		$a = substr($pred, $i, 1);
		print OUT "$a\n";
	}
	close OUT;
}

sub psipred() {
	my $seq = $_[0];
	my ($dbname, $ncbidir, $execdir, $datadir);
	$dbname = 'nr';
	$ncbidir = '~//blast-2.2.10/bin';
	$execdir = '~/psipred24/bin';
	$datadir = '~/psipred24/data';
	system "cp -f $seq.seq $seq.psitmp.fasta";
	system "$ncbidir/blastpgp -b 0 -j 3 -h 0.001 -d $dbname -i $seq.psitmp.fasta -Q $temp_freq -C $seq.psitmp.chk >& $seq.blast";
	system "echo $seq.psitmp.chk > $seq.psitmp.pn";
	system "echo $seq.psitmp.fasta > $seq.psitmp.sn";
	system "$ncbidir/makemat -P $seq.psitmp";
	system "$execdir/psipred $seq.psitmp.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $seq.ss";
	system "$execdir/psipass2 $datadir/weights_p2.dat 1 1.0 1.3 $seq.ss2 $seq.ss > $seq.horiz";
	print LOG "finishing psipred\n";
}

sub gen_model() {
	my $tn = $_[0];
	my $seq = $_[1];
	my @output = `./thread.bin $tn $seq 4.5 0.5 0.5 3 1 $Wprof`;
	print @output;
	my $tmpalign = $seq.".tmpalign";
	my $tmpscript = $seq.".tmpscript.py";
	&gen_align(\@output, $tn, $seq, $tmpalign);
	&gen_script($tmpscript, $tn, $seq, $tmpscript);
	system "mod8v1 $tmpscript";
	system "cat $seq.B99990001.pdb >> $seq-$tn.pdb";
}

sub gen_script {
	my $pdbpath = $DB."/pdb";
	my $pdbname = substr($_[1], 0, 4);
	my $tmpalign = $seq.".tmpalign";
	open SCRIPT, ">$_[0]";
	print SCRIPT "from modeller.automodel import *\n";
	print SCRIPT "log.verbose()\n";
	print SCRIPT "env = environ()\n";
	print SCRIPT "env.io.atom_files_directory = \'$pdbpath\'\n";
	print SCRIPT "a = automodel(env, alnfile = \'$tmpalign\', knowns = \'$pdbname\', sequence = \'$_[2]\')\n";
	print SCRIPT "a.starting_model = 1\n";
	print SCRIPT "a.ending_model = 1\n";
	print SCRIPT "a.make()\n";
	close SCRIPT;
}

sub gen_align {
	my $libpath = $DB.'/seq';
	my $line_len = 60;
	my $tmpseqpdb = $seq.".tmpseqpdb";
	my $pdbseq = &getseq($_[1]);
	&fastaprint($pdbseq, $_[1], $tmpseqpdb);
	my @modseq = `./seqalign.bin $tmpseqpdb $libpath/$_[1].seq`;
	print @modseq;
	my %st;
	my ($ttt, @ttt2);
	chomp @modseq;
	foreach $ttt (@modseq) {
		@ttt2 = split /\t/, $ttt;
		$st {$ttt2[1]} = $ttt2[0];
	}
	my @align = @{$_[0]};
	print @align;
	chomp @align;
	my $pdbname = substr($_[1], 0, 4);
	my $pdbchain = substr($_[1], 4, 1);
	$pdbchain = '@'if ($pdbchain eq '_');
	open STRUCT, "<$libpath/$_[1].seq";
	my @struct = <STRUCT>;
	close STRUCT;
	chomp @struct;
	shift @struct;
	my $temp = join('', @struct);
	@struct = split('', $temp);
	open SEQ, "<$_[2].seq";
	my  @seq = <SEQ>;
	close SEQ;
	chomp @seq;
	shift @seq;
	$temp = join('', @seq);
	@seq = split('', $temp);
	print @seq;
	my ($i_seq, $i_str, @seq_align, @str_align);
	my $i;
	for ($i = 0; $i <= $#align; $i++) {
		my ($i_str, $i_seq) = split /   /, $align[$i];
		if (($i_str > 0) && ($st {$i_str} == 0)) {
			print "$i  $i_str    $i_seq    $struct[$i_str - 1] \n";
			if ($i_seq != 0) {
				push @str_align, '-';
				push @seq_align, $seq[$i_seq - 1];
			}
			next;
		}
		print "$i_str   $i_seq  fg $struct[$i_str - 1]\n";
		if ($i_str == 0) {
			push @str_align, '-';
		}
		elsif ($i_str > 0) {
			push @str_align, $struct[$i_str - 1];
		}
		if ($i_seq ==  0) {
			push @seq_align, '-';
		}
		elsif ($i_seq > 0) {
			push @seq_align, $seq[$i_seq - 1];
		}
	}
	push  @seq_align, '*';
	push  @str_align, '*';
	print "@seq_align\n";
	print "@str_align\n";
	open OUT, ">$_[3]";
	print OUT ">P1;
	$pdbname\n";
	print OUT "structure:$pdbname", ':.:', $pdbchain, ':.:', $pdbchain, ':::0.00:0.00', "\n";
	my($i, $j) = (0, 0);
	my $buffer;
	while ($i <= $#str_align) {
		$buffer .= $str_align[$i];
		$i++;
		$j++;
		if ($j == $line_len) {
			print OUT "$buffer\n";
			$j = 0;
			$buffer = '';
			next;
		}
	}
	print OUT "$buffer\n\n" if ($j != 0);
	my($i, $j) = (0, 0);
	my $buffer;
	print OUT ">P1;$_[2]\n";
	print OUT "sequence:$_[2]:::::::0.00:0.00\n";
	while ($i <= $#seq_align) {
		$buffer .= $seq_align[$i];
		$i++;
		$j++;
		if ($j == $line_len) {
			print OUT "$buffer\n";
			$j = 0;
			$buffer = '';
			next;
		}
	}
	print OUT "$buffer\n" if ($j != 0);
	close OUT;
}
