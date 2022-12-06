#!/usr/bin/perl5.30.0 -w

use strict;

#print "1 \n";
#my $pre_dir = '/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/database/Wheat/';
#my $sp='IWGSC';
#my $gene='TraesCS5A02G542800';
#my $target_ch='chr5A';
#my $Working_dir='/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/User/69_166_59_242_TraesCS5A02G542800';
#my $p_size = 30000;
#my $a_size = 0;
#my $flag1 = 2;
#my $flag2 = 2;

my $pre_dir = $ARGV[0];
my $sp=$ARGV[1];
my $gene=$ARGV[2];
my $target_ch=$ARGV[3];
my $Working_dir=$ARGV[4];
my $p_size = $ARGV[5];
my $a_size = $ARGV[6];
my $flag1 = $ARGV[7];
my $flag2 = $ARGV[8];

my $ref_gb=$sp;
my $sp_dir = $pre_dir.$sp.'/'; 
my $sp_gb_dir = $pre_dir;

mkdir $sp_dir.'Candidate_genes/' unless -e $sp_dir.'Candidate_genes/';

my $gene_dir = $sp_dir.'Candidate_genes/'.$gene.'/'; # modifify this as the directory you can write into.

mkdir $gene_dir unless -e $gene_dir;

my $gene_mrna_file = $Working_dir.'/'.$gene.'_ref';

opendir (SpDir, $sp_gb_dir) || die;

my @genomes = ($ref_gb);

foreach (readdir SpDir) {
    
	push @genomes, $_ unless $_ =~ /\./ || $_ eq $ref_gb;
     
	}

my $md = '"m D"';



#push @genomes, 'Parent1', 'Parent2', 'query';
 push @genomes, 'Parent1', 'Parent2', $gene;


if ($flag1 == 1) {
	my $gene_syn_file_raw = $gene_dir.$gene.'_Haplotype_syn';
	open (Syn, '>'.$gene_syn_file_raw) || die;
	print Syn "query\tqry_S\tqry_E\tGenome\tch\tsbj_St\tsbj_E\tSize\tSimilarity\n";

	foreach my $g (@genomes) {
  	my $g_fas = $sp_gb_dir.$g.'/'.$g.'_'.$target_ch.'.fa';
  	my $g_fas_db = $g_fas.'.nin';

		if ($g =~ /Parent1/) {
			$g_fas = $Working_dir.'/'.$g.'/'.$g.'_'.$target_ch.'.fa';
			$g_fas_db = $g_fas.'.nin';
			}
		if ($g =~ /Parent2/) {
			$g_fas = $Working_dir.'/'.$g.'/'.$g.'_'.$target_ch.'.fa';
			$g_fas_db = $g_fas.'.nin';
			}
	#	if ($g =~ /query/) {
        if ($g =~ $gene) {

			$g_fas = $Working_dir.'/'.$g.'/'.$g.'_'.$target_ch.'.fa';
			$g_fas_db = $g_fas.'.nin';
			}
    
		my $Ref_mrna_out = $gene_dir.$gene.'-'.$g.'gb_out_m8';

		if (-e $g_fas_db) {
      system("/usr/local/bin/blastn -db $g_fas -query $gene_mrna_file -out $Ref_mrna_out -outfmt 6 -evalue 20");
		
		Extract_syn($Ref_mrna_out, $g, $target_ch, \*Syn);
		unlink $Ref_mrna_out  unless  $g eq $ref_gb;
	}

    }

	close Syn;	
}

if ($flag2 == 2) {
	print "s\n";
	my $gene_syn_file = $gene_dir.$gene.'_Haplotype_syn';
	my $gene_syn_fa = $gene_dir.$gene.'_Haplotype.fa';
	my $blast_self_o = $gene_dir.$gene.'_Haplotype-Self_out_m8';
	my $blast_self_o2 = $gene_dir.$gene.'_Haplotype-Self_out';
	my $blast_mrna_o = $gene_dir.$gene.'_ref_mRNA-Haplotype_out_m8';
	my $blast_mrna_o2 = $gene_dir.$gene.'_ref_mRNA-Haplotype_out';
	my $region_anno_file = $gene_dir.$gene.'_Haplotype_anno';
	my $seq_gap = $gene_dir.$gene.'_Haplotype_N_Gaps';
	my $var_file = $gene_dir.$gene.'_Variations';

     my $grass_repeats_fa = '/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/database/Repeat/grasrep.fa'; ## 05/17/2022

     my $blast_repeats_o = $gene_dir.$gene.'_repMask2'; # Repeats


	open (Gap, '>'.$seq_gap) || die;
	print Gap "Genome\tGAP_Start\tGAP_End\n";
	open (Anno, '>'.$region_anno_file) || die;
	print Anno "Genome\tch\tType\tStart\tEnd\tStrand\trefID\n";
	my ($syn_block_hashref, $syn_strand_hashref) = Parse_Syn_file ($gene_syn_file, $ref_gb);

	open (Syn, '>'.$gene_syn_fa) || die;   

	foreach my $g (keys %$syn_block_hashref) {
		my $g_strand = $$syn_strand_hashref{$g};
		my $g_fas = $sp_gb_dir.$g.'/'.$g.'_'.$target_ch.'.fa.gz'; ## 09/19/22

			if ($g =~ /Parent1/) {
			$g_fas = $Working_dir.'/'.$g.'/'.$g.'_'.$target_ch.'.fa.gz';

			}

			if ($g =~ /Parent2/) {
			$g_fas = $Working_dir.'/'.$g.'/'.$g.'_'.$target_ch.'.fa.gz';

			}

		#	if ($g =~ /query/) {
			if ($g =~ $gene) {
				
			$g_fas = $Working_dir.'/'.$g.'/'.$g.'_'.$target_ch.'.fa.gz';

			}

		next unless -e $g_fas;
		my @chs = keys %{ $$syn_block_hashref{$g} };
		my $ch = $chs[0];
		my @bps;
		push @bps, split /\t/, $_ foreach (@{ $$syn_block_hashref{$g}{$ch} });
		@bps = sort { $a <=> $b} @bps;
		my ($s1, $e1);

		if ($g_strand > 0) {
			  $s1 = $bps[0]  - $p_size; $e1 = $bps[-1] + $a_size;

			} else {  
				$s1 = $bps[-1] + $p_size; $e1 = $bps[0]  - $a_size;

				}
	
		print "$s1 == > $e1\n" if $g eq 'mace';		

		my $region = $ch.':'.$s1.'-'.$e1;            
		my $region0 = '>'.$region;     
		print "samtools faidx $g_fas $region\n";
		my $fas;
		if ($s1 < $e1) {
			$fas = `samtools faidx $g_fas $region`;  
		} else {
			($s1, $e1) = ($e1, $s1);
			$region = $ch.':'.$s1.'-'.$e1;
			$region0 = '>'.$region;
			print "samtools faidx $g_fas $region\n";
			$fas = `samtools faidx $g_fas $region -i --mark-strand no`;  
		}
		  
		my $newID = '>'.$g;
              
    $fas =~ s/$region0/$newID/g;   

		my $gap_hashref = Ns_position_size($fas);
		foreach my $gap (sort { $a <=> $b} keys %$gap_hashref) {
			my $size =  $$gap_hashref{$gap};
			next unless $size > 100;
			my $end = $gap + $size;
			print Gap $g."\t".$gap."\t".$end."\n";
			}

		print Syn $fas;                            

		my $gb_gff_file = $sp_gb_dir.$g.'/'.$g.'_gb.gff3';
		if (-e $gb_gff_file) {
			my $awk_para = "'".'$1=='.$target_ch.'&&$4>='.$s1.'&&$5<='.$e1."'";
			my @CDSs = `awk $awk_para $gb_gff_file`;
			foreach my $info (@CDSs) {
				my @t = split /\s+/, $info;
				next unless $t[2] eq 'CDS' || $t[2] eq 'gene';
				my $arrow = $t[6] eq '-' ? 1 : 2;
				my ($gID, $x) = ('.', '.'); 
				if ($g eq $ref_gb && $t[2] eq 'gene') { ($gID, $x) = $t[-1] =~ /ID\=(\S+)\;(N)ame/}
				print Anno join "\t", ($g, $t[0], $t[2], $t[3] - $s1, $t[4] - $s1, $arrow, $gID);
				print Anno "\n";
				}
			}
		}

	system("/usr/local/bin/makeblastdb -in $gene_syn_fa -dbtype nucl");

	system("/usr/local/bin/blastn -db $gene_syn_fa -query $gene_syn_fa -out $blast_self_o -outfmt 6 -evalue 10");  ## -W 11 or -W 9 ??
	system("/usr/local/bin/blastn -db $gene_syn_fa -query $gene_mrna_file -out $blast_mrna_o -outfmt 6 -evalue 10");
	system("/usr/local/bin/blastn -db $gene_syn_fa -query $gene_mrna_file -out $blast_mrna_o2 -outfmt 0 -evalue 10");

    system("/usr/local/bin/blastn -db $grass_repeats_fa -query $gene_syn_fa -out $blast_repeats_o -outfmt 6 -evalue 10"); #Repeats

	}
##################################################

sub Extract_syn {
	my ($f, $g, $ch, $o) = @_;
	open (F, $f) || die;
	my %hash;
	while (<F>) {
		chomp;
		my @t = split /\t/;
		next unless $t[0] =~ /\_mRNA/;
		next unless $t[1] eq $ch; ## if $t[1] =~ /^\d/ && 
		next unless $t[2] > 90;
#		next unless $t[9] < 1e-10
		my $info = $t[0]."\t".$t[6]."\t".$t[7]."\t".$g."\t".$t[1]."\t".$t[8]."\t".$t[9]."\t".$t[3]."\t".$t[2]."\n";
		print $o $info; 
	}
	close F;

 }
sub Parse_Syn_file {
	my ($f, $ref_gb) = @_;
	my (%hash, %hash2);
	open (F, $f) || die;

	my $Ref_strand;

	while (<F>) { 
		chomp;
		my @t = split /\t/;
		next if $t[0] eq 'query';
		
		my ($genome, $ch, $s_s, $s_e) = @t[3..$#t];
		  $Ref_strand =  $s_e - $s_s if $genome eq $ref_gb;
		my $g_strand = $s_e - $s_s;
		if ($g_strand * $Ref_strand > 0) {
			$hash2{$genome} =  1;
		} else {
			$hash2{$genome} =  -1;
		}
		print $genome.' '.$g_strand."\n"; 
		push @ { $hash{$genome}{$ch} } , $s_s."\t".$s_e;
		}
	close F;
	return (\%hash, \%hash2);	
	}
sub Ns_position_size {
	my ($fas) = @_;
	my %hash;
	my @t = split //, $fas;
	my $n_0 = 0;
	my $n_tag = 0;
	for (my $i = 0; $i <= $#t; $i ++) {
		if ($t[$i] eq 'N') {
			$n_tag = $i;
			$hash{$n_0} ++
			}
		my $diff = $i - $n_tag;
		$n_0 = $i if $diff > 1;
		}
	return \%hash;	
	}