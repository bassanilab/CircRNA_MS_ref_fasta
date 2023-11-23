# CircRNA_MS_ref_fasta
Create an MS reference fasta for circRNAs

## Introduction

CircRNA_MS_ref_fasta is a Perl script to produce mass spectrometry reference fasta sequences for immunopeptidomics that span the back-splice junction (BSJ) of circular RNAs (circRNAs). The script was designed to work with sequences from the circBase database (http://www.circbase.org/) in which the BSJ is located at the 3' end of the nucleotide sequence entry. Translated BSJ sequence fragments containing one or more BSJs are generated using a stop-to-stop (s2s) strategy in the three forward frames of each circRNA. The search for BSJs flanked by stop codons is initiated on concatenated circRNA sequences: 2x copies for circRNAs that contain an integral number of codons; and 4x copies for circRNAs that have an incomplete 3' terminal codon. This ensures that all BSJ derived AAs can be identified by scanning the concatenated sequences. If the first BSJ fragment located at the 5' end of the 4x concatenated circRNA does not have an upstream stop codon, it is fused with the previous/upstream (i.e. third) BSJ fragment in the concatemer. Where possible, the BSJ translated peptide sequences are truncated to 48 or 49 amino acids with the BSJ AA(s) centrally positioned to optimize immunopeptide searching. Peptide sequences are discarded if there is no methionine (MET) upstream of the BSJ AA(s), and sequences can be N-terminally truncated to start at MET, using the -trim2MET option, if the first MET is located within the shortened BSJ fragment.


## Usage

Running the script with no options or input file gives:

```
Usage: circRNA_MS_ref_fasta.pl [<options>] <input file> [output to STDOUT and STDERR]
where options are
-trim2MET       trim back globalApp seq if first ATG is seq-internal
-addFDRgroup    add or change an FDR group in the fasta header (PE=1 for UniProt/contaminants; PE=4 for circRNAs)
-help           print this
```
A typical usage would be:

```
$ circRNA_MS_ref_fasta.pl circBase_sequences.fa.gz > circBase_s2s_bsjf_pep.fa 2> circBase_s2s_bsjf_pep.log
```

Example sequence output:
```
>hsa_circ_0000021_0|chr1:16044387-16047883+|NM_015164|PLEKHM2 BSJ:KA24|ATG:-6|6|10|59
HLVRIMILEMCFQQCRLYPAQTGKAVPGCTWPSTRTPWRATCGCSRRT
>hsa_circ_0000021_2|chr1:16044387-16047883+|NM_015164|PLEKHM2 No_US_STOP|No_DS_STOP BSJ:G25|ATG:-48|98
PSSEDYDFGDVFPAVPSVPSTDWEGRAWLYLALNENSLESYLRLFQENL
```
The header of the fasta sequences is composed of 2 or 3 elements separated by a space. The first element contains the original circBase ID, the chromosomal location of the circRNA plus its Refseq ID and host gene name. The suffix added to the circBase ID distinguishes BSJ fragments originating from distinct reading frames of the circRNA sequence. The second element is present only for BSJ fragments that do not contain an upstream (US) and/or downstream (DS) stop codon, and may signal the presence of infinite translation in the given reading frame. The last element contains the local sequence amino acid coordinates of the BSJ(s) and the start (ATG) codons relative to the beginning of the peptide sequence. The BSJ AAs are always internal to the local sequence, whereas the ATGs can be located outside of the sequence (so upstream ATGs will have negative coordinates).


Example log output:
```
hsa_circ_0000021 OH-zero in frame 0 GGAAG|GCCGT 435 870
S2S [0] STP-202 ATG-346 ATG-379 ATG-391 BSJ-433 ATG-538 STP-565
s2sT: STP-0 ATG-48 ATG-59 ATG-63 BSJ-77 ATG-112 STP-121
fr: 0
[0] BSJ:77 BSJP:KA STP-0 ATG-48 ATG-59 ATG-63 BSJ-77 ATG-112 STP-121
LHS: HLVRIMILEMCFQQCRLYPAQTG  START: 54
BSJ: KA 77
RHS: VPGCTWPSTRTPWRATCGCSRRT
REJECT: no ATGs in frame 1 hsa_circ_0000021
S2S [2] ATG-217 BSJ-433 ATG-652
s2sT: ATG-73 BSJ-145 ATG-218
fr: 2
[2] BSJ:145 BSJP:G ATG-73 BSJ-145 ATG-218
LHS: PSSEDYDFGDVFPAVPSVPSTDWE  START: 121
BSJ: G 145
RHS: RAWLYLALNENSLESYLRLFQENL
```
Each circRNA processed by the script has an entry in the log output (to STDERR), and contains information about the circRNA and the derived S2S BSJ fragments derived from the three reading frames: circBase ID; BSJ location at the nucleotide level; circRNA length; 2x or 4x concatenated length. CircRNAs that contain an integral number of codons are denoted "OH-zero". BSJ fragments lacking a MET or containing a peptide sequence that is identical to previous output fasta are logged and skipped.

Truncation to the first MET in the S2S fragment in the fasta peptide sequences can be achieved with the -trim2MET option as follows:
```
$ circRNA_MS_ref_fasta.pl -trim2MET circBase_s2s_bsjf_pep.fa > circBase_s2s_bsjf_pep_MET.fa 2> circBase_s2s_bsjf_pep_MET.log
```
Only those fasta peptides containing a MET that is the first MET in the S2S fragment will be truncated.

The -addFDRgroup option adds or changes an FDR group in the fasta headers of, for example, the circBase_s2s_bsjf_pep_MET.fa file for use in group-specific FDR mass spectrometry analysis, compatible with FragPipe v20.0 onwards (https://github.com/Nesvilab/FragPipe). Two groups are implemented: one for circRNAs, designated 'PE=4'; and one for UniProt and contaminant fasta, designated 'PE=1'. Addition of the 'PE=4' group designation to a circRNA fasta modifies the header to a UniProt-like format:
```
$ circRNA_MS_ref_fasta.pl -addFDRgroup circBase_s2s_bsjf_pep_MET.fa > circBase_s2s_bsjf_pep_MET_gsFDR.fa 2> circBase_s2s_bsjf_pep_MET_gsFDR.log

>hsa_circ_0000021_0|chr1:16044387-16047883+|NM_015164|PLEKHM2 BSJ:KA24|ATG:-6|6|10|59
HLVRIMILEMCFQQCRLYPAQTGKAVPGCTWPSTRTPWRATCGCSRRT

>hsa_circ_0000021_0|hsa_circ_0000021_0 circPLEKHM2 OS=Homo sapiens OX=9696 GN=PLEKHM2 PE=4
HLVRIMILEMCFQQCRLYPAQTGKAVPGCTWPSTRTPWRATCGCSRRT
```

## License

Copyright (C) LICR - Ludwig Institute for Cancer Research, Lausanne, Switzerland

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation.


## Developer

Brian Stevenson, CHUV/SIB, Lausanne, Switzerland (brian.stevenson@unil.ch)

