#!/usr/bin/env perl
# Script to produce Met-restricted S2S BSJ fragments from
#  human_hg19_circRNAs_putative_spliced_sequence.fa.gz (circBase; http://www.circbase.org/)
#  v1.0 230823
#  v1.1 231123

use strict;
use warnings;
use Getopt::Long;

my $trim2MET;
my $fdrGroup;
my $help;

GetOptions("trim2MET!" => \$trim2MET,
	   "addFDRgroup!" => \$fdrGroup,
	   "help!" => \$help);

die &usage if ($help);
die &usage unless ($ARGV[0]);

# Lookup hash for translate sub
my %code = (TTT=>'F', TTC=>'F', TTA=>'L', TTG=>'L',
            TCT=>'S', TCC=>'S', TCA=>'S', TCG=>'S',
            TAT=>'Y', TAC=>'Y', TAA=>'*', TAG=>'*',
            TGT=>'C', TGC=>'C', TGA=>'*', TGG=>'W',
            CTT=>'L', CTC=>'L', CTA=>'L', CTG=>'L',
            CCT=>'P', CCC=>'P', CCA=>'P', CCG=>'P',
            CAT=>'H', CAC=>'H', CAA=>'Q', CAG=>'Q',
            CGT=>'R', CGC=>'R', CGA=>'R', CGG=>'R',
            ATT=>'I', ATC=>'I', ATA=>'I', ATG=>'M',
            ACT=>'T', ACC=>'T', ACA=>'T', ACG=>'T',
            AAT=>'N', AAC=>'N', AAA=>'K', AAG=>'K',
            AGT=>'S', AGC=>'S', AGA=>'R', AGG=>'R',
            GTT=>'V', GTC=>'V', GTA=>'V', GTG=>'V',
            GCT=>'A', GCC=>'A', GCA=>'A', GCG=>'A',
            GAT=>'D', GAC=>'D', GAA=>'E', GAG=>'E',
            GGT=>'G', GGC=>'G', GGA=>'G', GGG=>'G');

if ($trim2MET) {
    &trim2MET(@ARGV);
    exit 0;
}

if ($fdrGroup) {
    &addFDRgroup(@ARGV);
    exit 0;
}

my $file = shift;

$/ = "\n>"; # fasta sequences
    
if ($file =~ /gz$/) { open IN, "zcat $file |" or die "Can't open $file: $!"; }
elsif ($file =~ /bz2$/) { open IN, "bzcat $file |" or die "Can't open $file: $!"; }
else { open IN, $file or die "Can't open $file: $!"; }

my $fxSeq;  # multiple units of circRNA (4x or 2x)
my $fxSeqL; # length
my $cid;    # circRNA ID
my @h;      # header elements
my @oh;     # overhangs
my @frm;    # frames
my @bsj;    # single (OH-zero) or three BSJs (codon coordinates)
my @bsjP;   # single (OH-zero) or three BSJ AA(s)
my @atg;    # ATG starts
my @stp;    # STOPs
my %all;    # collect all the coords with descriptive hybrid keys

my $cnt = 0;   # circRNA count
my %lDup = (); # redundant protein sequences
    
while (<IN>) {
    @oh = ();
    @frm = ();
    @bsj = ();
    @bsjP = ();
    @atg = ();
    @stp = ();
    %all = ();
    $cnt++;
    
    s/>//g;           # remove '>'s
    s/^(.+)$//m;      # remove header
    my $header = $1;  # save header
    @h = split /\|/, $header;
    $cid = $h[0];
    $cid =~ s/^>//;
    s/\n//g;          # remove all line breaks
    my $seq = $_;
    $seq = uc $seq;   # to uppercase (just in case)
    my $seqL = length($seq);
    if ($seqL < 28) {
	print STDERR "Input circRNA seq too short: $header\n\n";
	next;
    }
    
    my $seqOH = $seqL % 3; # overhang (OH) of frame 0 in seq
    my $unitL = $seqL - $seqOH;

    my $bsjNseq = join('|', substr($seq, $seqL-5), substr($seq, 0, 5)); # graphic position of BSJ in nucletoides (log only)

    $fxSeq = join('', $seq, $seq, $seq, $seq); # 4x unit circRNA seq
    $fxSeqL = length($fxSeq);
       
    if ($seqOH == 0) { # different from other cases; manually create frames 1 and 2 [OH is 0, 2, 1 - but all stay in frame]
	$fxSeq = join('', $seq, $seq); # just two units
	$fxSeqL = length($fxSeq);
	$bsj[0] = $unitL - 2; # this coord is same for all three frames since fxSeq starts on codon boundary!
	# Three separate frames 0 1 2 on 1nt cycled sequence
	for my $f (0, 1, 2) { # three frames
	    print STDERR "$cid OH-zero in frame 0 $bsjNseq $seqL $fxSeqL\n" if ($f == 0); # just once (OH is 2 for frame 1 and 1 for frame 2)
	    %all = (); # empty the hash
	    if ($f == 0) { @oh = (0, 0); }
	    elsif ($f == 1) { @oh = (2, 2); }
	    elsif ($f == 2) { @oh = (1, 1); }
	    @frm = ($f, $f);
	    if ($f != 0) {
		my $bp = substr($fxSeq, 0, 1, '');
		$fxSeq .= $bp; # add first nuc to end
	    }
	    my $zero = 1; # set only for this special case
	    &processCR($zero);
	}
    }
    elsif ($seqOH == 1) { # frame sequence is 0 2 1 0; need only 0 2 1 [OH: 1 2 0]
	@oh = (1, 2, 0);
	@frm = (0, 2, 1); # the frame *before* the BSJ
	@bsj = ($unitL + 1, (2 * $unitL) + 1, (3 * $unitL) + 1);
	
	print STDERR "$cid $bsjNseq $seqL $fxSeqL\n";
	my $zero = 0;
	&processCR($zero);
    }
    elsif ($seqOH == 2) { # frame sequence is 0 1 2 0; need only 0 1 2 [OH: 2 1 0]
	@oh = (2, 1, 0);
	@frm = (0, 1, 2); # the frame *before* the BSJ
	@bsj = ($unitL + 1, (2 * $unitL) + 4, (3 * $unitL) + 4);
	
	print STDERR "$cid $bsjNseq $seqL $fxSeqL\n";
	my $zero = 0;
	&processCR($zero);
    }
    print STDERR "\n"; # circRNA block separator in log
}
close IN;

exit 0;

sub usage {
    "Usage: circRNA_MS_ref_fasta.pl [<options>] <input file> [output to STDOUT and STDERR]
where options are
-trim2MET       trim back globalApp seq if first ATG is seq-internal
-addFDRgroup    add or change an FDR group in the fasta header (PE=1 for UniProt/contaminants; PE=4 for circRNAs)
-help           print this\n";
}

sub processCR {
    my $zero = shift; # OH-zero (1) or not (0)
    # Gather data - elements are BSJ-# BSJS-# STP-# ATG-#
    for (my $j=0;$j<@bsj;$j++) {
	my $nt = ($oh[$j] == 0) ? 6 : 3; # zero OHs have 2aa BSJ; other OHs have 1aa BSJ
	my $jaa = &translate(substr($fxSeq, $bsj[$j] - 1, $nt));
	push @bsjP, $jaa; # before possible OH-zero rejection
	my $hk;
	if ($jaa =~ /\*/) {
	    if ($zero) {
		print STDERR "REJECT: STOP at BSJ $jaa frame $frm[0] $cid\n";
		return; # skip this sequence
	    }
	    $hk = 'BSJS-' . $bsj[$j];
	    print STDERR "INFO: STOP at $hk $jaa [$j] $cid\n"; # don't reject yet!
	}
	else { $hk = 'BSJ-' . $bsj[$j]; } # 'normal' BSJ hybrid key
	$all{$hk} = $bsj[$j];
	
    }
    my $hasATG;
    @atg = &getATGs(\$fxSeq); # triplets of coordN coordN.Kozak coordP.Kozak
    for (my $i=0;$i<@atg;$i+=3) {
	my $hk = 'ATG-' . $atg[$i];
	$all{$hk} = $atg[$i];
	$hasATG = 1;  # at least one ATG
    }
    if ($zero) {
	if (! $hasATG) {
	    print STDERR "REJECT: no ATGs in frame $frm[0] $cid\n";
	    return;
	}
    }
    @stp = &getSTOPs(\$fxSeq);
    my $hasSTOP;
    for my $stp (@stp) {
	$hasSTOP = 1; # at least one STOP
	next if (exists $all{"BSJS-$stp"}); # skip STOPs created at BSJs - already in %all
	my $hk = 'STP-' . $stp;
	$all{$hk} = $stp;
    }
    # Order the hybrid keys by ascending coords
    my @all = sort { $all{$a}<=>$all{$b} } keys %all;
    
    # Create lookup for BSJ elements -> position in fxSeq (0, 1 or 2)  [not used for OH-zero]
    my %bs2p = (); # hash or
    my @bs2p = (); # array
    my $n = 0;
    for my $e (@all) {
	if ($e =~ /^BSJ/) {
	    $bs2p{$e} = $n;
	    push @bs2p, $e;
	    $n++;
	}
    }
    
    # Isolate the minimal s2s in the all list that contains each BSJ (or several BSJs)
    my @s2s = (); # list of lists
    # Get position (not coord!) of each BSJ in @all and expand 5' and 3' to STOP or first/last element
    my @i = (); # single BSJ oposition for OH-zero; three BSJ positions for others
    for (my $i=0;$i<@all;$i++) {
	push @i, $i if ($all[$i] =~ /^BSJ/); # will also catch BSJS
    }
    if ($zero) {
	my $i = $i[0]; # single BSJ
	my @sa = reverse(@all[0..$i]); # back; include BSJ
	for my $e (@sa) {
	    unshift @s2s, $e;
	    last if ($e =~ /^STP/);
	}
	@sa = (@all[$i+1..$#all]); # forward; exclude BSJ
	for my $e (@sa) {
	    push @s2s, $e;
	    last if ($e =~ /^STP/);
	}
    }
    else {
	for my $j (0, 1, 2) {
	    my $i = $i[$j]; # to simplify list of lists in @s2s
	    my @sa = reverse(@all[0..$i]); # back; include BSJ
	    for my $e (@sa) {
		unshift @{ $s2s[$j] }, $e;
		last if ($e =~ /^STP/); # no STOP on /BSJS/
	    }
	    @sa = (@all[$i+1..$#all]); # forward; exclude BSJ
	    for my $e (@sa) {
		push @{ $s2s[$j] }, $e;
		last if ($e =~ /^STP/); # no STOP on /BSJS/
	    }
	    print STDERR "S2S $j @{ $s2s[$j] }\n";
	}
    }
    
    # Now output the sequences containing one or more BSJs; identify and skip BSJs with STOPs
    my $skipTwo; # set to skip ($p = 2) after fusion
    my %seen = (); # log BSJ position (0, 1 or 2) as processed
    for my $p (0, 1, 2) {
	if ($zero) {
	    next unless ($p == $frm[0]);
	}
	if ($p == 2 and ($skipTwo)) {
	    print STDERR "SKIP [2] since already processed [2] + [0] fusion $cid\n";
	    next;
	}
	if ($zero) { print STDERR "S2S [$p] @s2s\n"; }
	else { print STDERR "S2S [$p] @{ $s2s[$p] }\n"; }
	my $cidF = $cid . '_'. $p;
	
	my $beg;
	my $end;
	my $bsjf;
	my @fusion = (); # merged [2] + [0]
	my $fusion;
	my @atgF = (); # 'extended' coords for fusion ATGs

	if ($zero) {
	    if ($s2s[0] =~ /^STP/) { $beg = $all{$s2s[0]}; }
	    else { $beg = 1; }
	    if ($s2s[-1] =~ /^STP/) { $end = $all{$s2s[-1]} + 2; }
	    else { $end = $fxSeqL; }
	    $bsjf = substr($fxSeq, $beg - 1, ($end - $beg + 1));
	}
	elsif (($p == 0 and $s2s[0][0] !~ /^STP/) or ($p == 2 and $s2s[2][-1] !~ /^STP/)) {
	    # Merge [2] + [0] into @fusion and update the ATG coords in the downstream [0] part into @atgF
	    $fusion = 1; # only set for this $p
	    $skipTwo = 1; # set until next circRNA
	    
	    for my $e (@{ $s2s[2] }) {
		if ($bsjP[2] =~ /\*/) {
		    next unless ($e =~ /^BSJS/); # take only the BSJS
		}
		push @fusion, $e; # include first BSJ
		if ($e =~ /^ATG/) {
		    push @atgF, $all{$e}; # just the coord
		}
		last if ($e eq $bs2p[2]); # include [2] BSJ|BSJS - check w/ @bs2p
	    }
	    my $cd = $bsj[2] + 2; # coord difference [2] -> [0]
	    for my $e (@{ $s2s[0] }) {
		if ($e =~ /^ATG/) {
		    push @atgF, ($all{$e} + $cd); # just the updated coord
		}
		last if ($e =~ /^BSJ/); # only ATGs upstream of BSJ
	    }
	    for my $e (@{ $s2s[0] }) {
		push @fusion, $e; # all of the [0] downstream of BSJ
		last if ($e =~ /^BSJS/); # stop at BSJS if present
	    }
	    
	    print STDERR "FUSION: [$p] merged @fusion; ATGs: @atgF\n";
	    
	    # Assemble the BSJF seq: first and last elements are STP or BSJS
	    my $fs = $fusion[-1]; # last element
	    $end = $all{$fs} + 2;
	    my $zeroSeq = substr($fxSeq, 0, $end);
	    
	    $fs = $fusion[0]; # first element
	    $beg = $all{$fs};
	    $end = $bsj[2] + 2; # bsj[2] first codon end
	    my $twoSeq = substr($fxSeq, $beg - 1, $end - $beg + 1);
	    $bsjf = $twoSeq . $zeroSeq;
	}
	else { # default - no fusion, just [0], [1] and/or [2]; may/may not have STOPs
	    # Check for STOP at BSJ
	    if ($bsjP[$p] =~ /\*/) {
		print STDERR "REJECT: STOP at BSJ $bsjP[$p] frame $p $cid\n";
		next;
	    }
	    if ($bsjP[$p-1] =~ /\*/) { # Check for (non-$p) BSJS in element list
		my @e = @{ $s2s[$p] };
		my @f = reverse(@e);
		@e = ();
		for my $e (@f) {
		    unshift @e, $e;
		    last if ($e =~/^BSJS/);
		}
		@{ $s2s[$p] } = @e; # replace
	    }
	    $beg = $all{$s2s[$p][0]};
	    $end = $all{$s2s[$p][-1]} + 2;
	    $bsjf = substr($fxSeq, $beg - 1, ($end - $beg + 1));
	}
	
	my $bsjfT = &translate($bsjf);
	my $bsjfTL = length($bsjfT);
	
	my @info = ();
	push @info, 'No_US_STOP', 'No_DS_STOP' unless ($hasSTOP);
	
	# Trim to 23/24aa surrounding BSJ and reject BSJfrags w/o upstream MET
	my @e = ();
	if ($zero) {
	    @e = @s2s;
	}
	else {
	    @e = @{ $s2s[$p] };
	    @e = @fusion if ($fusion);
	}
	my @fr = (); # 0, 1 or 2; can have multiple BSJs in s2s; generate dynamically for all
	
	# Convert @e to %s2sT with AA coords; ATG-aaCord > aaCoord
	my %s2sT = (); # s2s elements with AA coords
	my @bsjAAc; # BSJ AA coords
	my $ncd;
	$ncd = $bsj[2] + 2 if (! $zero); # nuc coord difference [2] -> [0]; add before converting to AAc

	if ($zero) {
	    if ($bsjfT =~ /^\*/) { # has upstream STP
		my $cd = $all{$s2s[0]} + 3; # start coord of first codon AFTER STP-
		for my $e (@s2s) {
		    my ($n, $c) = split /-/, $e;
		    $c = ($c - $cd + 3) / 3; # convert to AAcoord
		    my $ne = join('-', $n, $c);
		    $s2sT{$ne} = $c;
		    $bsjAAc[$p] = $c if ($e =~ /^BSJ/); # single BSJ
		}
	    }
	    else {
		for my $e (@s2s) { # no upstream STP
		    my ($n, $c) = split /-/, $e;
		    $c = ($c + 2) / 3; # convert to AAcoord
		    my $ne = join('-', $n, $c); # direct nucCoord->aaCoord
		    $s2sT{$ne} = $c;
		    $bsjAAc[$p] = $c if ($e =~ /^BSJ/); # single BSJ
		}
	    }
	    @fr = ($p); # single BSJ
	    $seen{$p}++; # log as processed
	}
	elsif ($bsjfT =~ /^\*/) { # has upstream STP or BSJS
	    my $ts = $all{$e[0]}; # [2] start (to ID [0] elements)
	    my $cd = $ts + 2; # one less than start coord of first codon AFTER STP-
	    for my $e (@e) {
		next if ($e =~ /^STP/); # skip the STOPs
		my ($n, $c) = split /-/, $e;
		$c += $ncd if ($c < $ts); # update coord in [0]
		$c = ($c - $cd + 2) / 3; # convert to AAcoord
		my $ne = join('-', $n, $c);
		$s2sT{$ne} = $c;
		if ($e =~ /^BSJ/) {
		    my $x = $bs2p{$e};
		    $bsjAAc[$x] = $c;
		    push @fr, $x; # dynamically create (fr) array
		    $seen{$x}++; # log as processed
		}
	    }
	}
	else {
	    for my $e (@e) { # no upstream STP; assume seq starts at 1
		next if ($e =~ /^STP/); # skip the STOPs
		my ($n, $c) = split /-/, $e;
		$c += $ncd if ($c < $bsj[2]); # update coord
		$c = ($c + 2) / 3; # convert to AAcoord
		my $ne = join('-', $n, $c); # direct nucCoord->aaCoord
		$s2sT{$ne} = $c;
		if ($e =~ /^BSJ/) {
		    my $x = $bs2p{$e};
		    $bsjAAc[$x] = $c;
		    push @fr, $x; # dynamically create (fr) array
		    $seen{$x}++; # log as processed
		}
	    }
	}
	
	my @s2sT = sort { $s2sT{$a}<=>$s2sT{$b} } keys %s2sT;
	print STDERR "s2sT: @s2sT\nfr: @fr\n";
	
	# Go thru protein sequence as many times as there are BSJs (fusion is 2, 0)
	for my $f (@fr) {
	    if ($seen{$f} > 1) {
		print STDERR "Already processed BSJ [$f] in $cid\n";
		next;
	    }
	    print STDERR "[$f] BSJ:$bsjAAc[$f] BSJP:$bsjP[$f] @s2sT\n";
	    
	    if ($bsjP[$f] =~ /\*/) {
		print STDERR "REJECT: STOP present in [$f] $bsjP[$f] $cid\n";
		next;
	    }
	    
	    $cidF =~ s/_\d+$//;
	    $cidF .= "_$f";
	    
	    # Check for ATG before BSJ in @s2sT
	    my $bad;
	    my $atg;
	    for my $e (@s2sT) {
		if ($e =~ /^ATG/) {
		    $atg = 1;
		}
		if ($e eq 'BSJ-'.$bsjAAc[$f]) {
		    $bad = 1 unless ($atg); # ATG must be before bsjAAc[$f]
		    last;
		}
	    }
	    if ($bad) {
		print STDERR "REJECT: no ATG between 5' STOP and BSJ-$bsjAAc[$f] in [$f] $cid\n";
		next;
	    }
	    
	    # Trim bsjfT to 23aa-BSJ(2aa)-23aa
	    $bsjfT =~ s/^\*//; # delete N-term STOP
	    $bsjfT =~ s/\*$//; # delete C-term STOP
	    $bsjfTL = length($bsjfT);
	    
	    my $lhs;
	    my $lhsC; # start coord for LHS (1-based)
	    my $rhs;
	    my $bsjL = length($bsjP[$f]);
	    
	    if ($bsjAAc[$f] < (26 - $bsjL)) { # not enough for 23aa prefix
		$lhs = substr($bsjfT, 0, $bsjAAc[$f] - 1); # all pep before BSJ
		$lhsC = 1;
	    }
	    else {
		$lhs = substr($bsjfT, $bsjAAc[$f] - 26 + $bsjL, 25 - $bsjL);
		$lhsC = $bsjAAc[$f] - 26 + $bsjL + 1;
	    }
	    if ($bsjfTL < ($bsjAAc[$f] + $bsjL + 25 - $bsjL)) { # not enough PEP for 23aa suffix
		print STDERR "Not enough PEP seq for RHS\n";
	    }
	    $rhs = substr($bsjfT, $bsjAAc[$f] - 1 + $bsjL, 25 - $bsjL);
	    
	    print STDERR "LHS: $lhs  START: $lhsC\n";
	    print STDERR "BSJ: $bsjP[$f] $bsjAAc[$f]\n";
	    print STDERR "RHS: $rhs\n";
	    
	    my $bsjfTP = $lhs.$bsjP[$f].$rhs; # distinct from $bsjfT
	    
	    if (length($bsjfTP) < 8) {
		print STDERR "REJECT short peptide: $bsjfTP\n";
		next;
	    }
	    
	    # Now get final coords of BSJ and ATGs
	    my @info2 = ();
	    my @atgR = (); # coords ATG elements in 23-2-23 region of bsjfT
	    my $newPos = $bsjAAc[$f] - $lhsC + 1; # for BSJ
	    for my $e (@s2sT) {
		next unless ($e =~ /^ATG/); # ATGs only
		my ($n, $c) = split /-/, $e;
		$c = $c - $lhsC + 1;
		$c-- if ($c <= 0); # correct for zero-less scale :-)
		push @atgR, $c; # list all ATGs upstream of BSJ
		last if ($e =~ /^BSJ/); # only ATGs upstream of BSJ; woops never reached!
	    }
	    push @info2, 'BSJ:'.$bsjP[$f].$newPos;
	    push @info2, join(':', 'ATG', join('|', @atgR));
	    
	    # Test if seen sequence before
	    if (exists $lDup{$bsjfTP}) {
		print STDERR "REJECT: $cid [$f] has identical sequence to [$lDup{$bsjfTP}]\n";
		next;
	    }
	    else {
		$lDup{$bsjfTP} = $cidF;
	    }
	    my $info;
	    if ($#info < 0) { # empty
		$info = join('|', @info2); # only BSJ/ATG info
	    }
	    else { $info = join(' ', join('|', @info), join('|', @info2)); } # separate STOP info from BSJ/ATG info
	    my $vH = join('|', $cidF, @h[1..$#h]);
	    
	    $bsjfTP =~ s/(.{60})/$1\n/g;
	    $bsjfTP =~ s/\s+$//s;
	    print ">$vH $info\n$bsjfTP\n";
	}
    }
}

sub trim2MET {
    $/ = "\n>"; # fasta sequences
    my $file = shift;
    
    if ($file =~ /gz$/) { open IN, "zcat $file |" or die "Can't open $file: $!"; }
    elsif ($file =~ /bz2$/) { open IN, "bzcat $file |" or die "Can't open $file: $!"; }
    else { open IN, $file or die "Can't open $file: $!"; }

    my %lDup = (); # duplication of protein sequences (rare)
    
    while (<IN>) {
	s/>//g;           # remove '>'s
	s/^(.+)$//m;      # remove header
	my $header = $1;  # save header
	my $pep = $_;
	$pep =~ s/\n//g;  # remove all line breaks
	
	# Is the first ATG in the seq or upstream?
	my @h  = split " ", $header; # BSJ|ATG are [-1]
	my @x = split /\|/, $h[-1];
	my ($bsjAA, $bsjAAc) = $x[0] =~ /^BSJ:(\D+)(\d+)$/;
	my @atg = @x[1..$#x];
	die "Problem with $header\n" unless ($atg[0] =~ /^ATG:/);
	$atg[0] =~ s/^ATG://;

	if (($atg[0] <= 1) or ($header =~ /No_US_STOP/)) { # no need to trim seq
	    
	    if (exists $lDup{$pep}) {
		print STDERR "REJECT: $header $pep has identical sequence to [$lDup{$pep}]\n";
		next;
	    }
	    else { $lDup{$pep} = $header; }
	    
	    $pep =~ s/(.{60})/$1\n/g;
	    $pep =~ s/\s+$//s;
	    
	    my @nc = ();
	    for my $c (@atg) {
		push @nc, $c;
	    }
	    my $new = join('|', $x[0], join(':', 'ATG', join('|', @nc)));
	    if ($#h == 1) { print ">$h[0] $new\n$pep\n"; }
	    else { print ">$h[0] $h[1] $new\n$pep\n"; }
	}
	elsif ($atg[0] > 1) { # trim seq
	    my $delta = $atg[0] - 1;
	    my $nPep = substr($pep, $delta);
	    
	    if (length($nPep) < 8) {
		print STDERR "REJECT short peptide: $nPep ($h[0])\n";
		next;
	    }
	    
	    if (exists $lDup{$nPep}) {
		print STDERR "REJECT: $header $nPep has identical sequence to [$lDup{$nPep}]\n";
		next;
	    }
	    else { $lDup{$nPep} = $header; }
	    
	    $nPep =~ s/(.{60})/$1\n/g;
	    $nPep =~ s/\s+$//s;
	    
	    $bsjAAc -= $delta;
	    my $nb = 'BSJ:'.$bsjAA.$bsjAAc;
	    
	    my @nc = ();
	    for my $c (@atg) {
		$c -= $delta;
		push @nc, $c;
	    }
	    
	    my $new = join('|', $nb, join(':', 'ATG', join('|', @nc)));
	    
	    if ($#h == 1) { print ">$h[0] $new\n$nPep\n"; }
	    else { print ">$h[0] $h[1] $new\n$nPep\n"; }
	}
	else { die "ATG of $header after BSJ\n"; }
    }
    close IN;
}

sub addFDRgroup {
    $/ = "\n>"; # fasta sequences
    my $file = shift;
    
    if ($file =~ /gz$/) { open IN, "zcat $file |" or die "Can't open $file: $!"; }
    elsif ($file =~ /bz2$/) { open IN, "bzcat $file |" or die "Can't open $file: $!"; }
    else { open IN, $file or die "Can't open $file: $!"; }

    my %lDup = (); # duplication of protein sequences (rare)
    
    while (<IN>) {
	s/>//g;           # remove '>'s
	s/^(.+)$//m;      # remove header
	my $header = $1;  # save header
	my @h = split " ", $header;
	my $pep = $_;
	$pep =~ s/\n//g;  # remove all line breaks

	my $newH;
	my $pe;
	if ($h[0] =~ /^sp/ or $h[0] =~ /^rev_sp/ ) { # UniProt, contamination or decoy
	    for my $e (@h) {
		if ($e =~ /^PE=\d+$/) {
		    $e = "PE=1";
		    $pe = 1;
		}
	    }
	    if ($pe) { $newH = join(" ", @h); }
	    else { $newH = join(" ", @h, 'PE=1'); }
	}
	elsif ($h[0] =~ /^hsa_circ/) { # circRNA; remodel UniProt style
	    if ($header =~ /^PE=\d+$/) { die "circRNA already has PE assignment: $header\n"; }
	    my @cid = split /\|/, $h[0];
	    $newH = join(' ', join('|', $cid[0], $cid[0]), "circ$cid[3]", 'OS=Homo sapiens OX=9696', "GN=$cid[3]", 'PE=4');
	}
	else { print STDERR "Unrecognised header: $header\n"; }
	
	$pep =~ s/(.{60})/$1\n/g;
	$pep =~ s/\s+$//s;
	print ">$newH\n$pep\n";
    }
}

sub translate {
    my $seq = shift;
    my $p = '';
    my $offset = 0;
    while ($offset < length $seq) {
        my $codon = substr($seq, $offset, 3);
        if (length($codon) == 2) {
            $codon .= 'N';
        }
        my $aa;
        if (exists $code{$codon}) { $aa = $code{$codon}; }
        else { $aa = 'X'; }
        $p .= $aa;
        $offset += 3;
    }
    $p =~ s/X$//; # remove terminal incomplete AA
    return $p;
}

sub getATGs { # single frame [0] search; returned coords are base 1
    my $s = shift; # ref to sequence
    # Find ATGs in single forward frame [0]
    my @pos = ();
    my $p = 0;
    while (my $c = substr($$s, $p, 3)) {
	last if (length($c) < 3);
	if ($c eq 'ATG') {
	    push @pos, ($p + 1); # 1 based coord
	    my $kSeq;
	    if ($p >= 3) { $kSeq = substr($$s, $p - 3, 7); } # position of Kozak consensus
	    else { $kSeq = substr($$s, 0, $p + 4); } # incomplete consensus
	    if ($kSeq =~ /[AG]..ATGG/) { # test core consensus
		push @pos, ($p + 1).".$kSeq.".'(K)'; # convert to base 1
		push @pos, (($p + 3)/3).".$kSeq.".'(K)'; # convert to base 1
	    }
	    else {
		push @pos, ($p + 1).".$kSeq";
		push @pos, (($p + 3)/3).".$kSeq";
	    }
	}
	$p += 3; # next codon
    }
    return @pos;
}

sub getSTOPs { # single frame [0] search; returned coords are base 1
    my $s = shift; # ref to sequence
    # Find STOPs in single forward frame [0]
    my @pos = ();
    my $p = 0;
    while (my $c = substr($$s, $p, 3)) {
	last if (length($c) < 3);
	if ($c =~ /T(?:AA|AG|GA)/) {
	    push @pos, ($p + 1); # 1 based coord
	}
	$p += 3; # next codon
    }
    return @pos;
}

