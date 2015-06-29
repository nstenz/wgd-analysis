#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $max_forks = 40;

my %rev_codon_table = (
	S => qr/((AG[C|T])|(TC.))/,
	F => qr/TT[T|C]/,
	L => qr/((TT[A|G])|(CT.))/,
	Y => qr/TA[T|C]/,
	C => qr/TG[T|C]/,
	W => qr/TGG/,
	P => qr/CC./,
	H => qr/CA[T|C]/,
	Q => qr/CA[A|G]/,
	R => qr/((AG[A|G])|(CG.))/,
	I => qr/AT[T|C|A]/,
	M => qr/ATG/,
	T => qr/AC./,
	N => qr/AA[T|C]/,
	K => qr/AA[A|G]/,
	V => qr/GT./,
	A => qr/GC./,
	D => qr/GA[T|C]/,
	E => qr/GA[A|G]/,
	G => qr/GG./,
	#X => qr/((TA[A|G])|TGA)/,
	X => qr/.../,
	'*' => qr/((TA[A|G])|TGA)/,
);

my %codon_table = (
	TTT => "F", TTC => "F", TTA => "L", TTG => "L",
	TCT => "S", TCC => "S", TCA => "S", TCG => "S",
	#TAT => "Y", TAC => "Y", TAA => "X", TAG => "X",
	TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
	#TGT => "C", TGC => "C", TGA => "X", TGG => "W",
	TGT => "C", TGC => "C", TGA => "*", TGG => "W",
	CTT => "L", CTC => "L", CTA => "L", CTG => "L",
	CCT => "P", CCC => "P", CCA => "P", CCG => "P",
	CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
	CGT => "R", CGC => "R", CGA => "R", CGG => "R",
	ATT => "I", ATC => "I", ATA => "I", ATG => "M",
	ACT => "T", ACC => "T", ACA => "T", ACG => "T",
	AAT => "N", AAC => "N", AAA => "K", AAG => "K",
	AGT => "S", AGC => "S", AGA => "R", AGG => "R",
	GTT => "V", GTC => "V", GTA => "V", GTG => "V",
	GCT => "A", GCC => "A", GCA => "A", GCG => "A",
	GAT => "D", GAC => "D", GAA => "E", GAG => "E",
	GGT => "G", GGC => "G", GGA => "G", GGG => "G",
	# ambiguities where we can still determine amino acid
	TTY => "F", 
	TTR => "L", 
	CTN => "L", CTM => "L", CTR => "L", CTW => "L", 
	CTS => "L", CTY => "L",	CTY => "L", CTK => "L", 
	CTV => "L", CTH => "L", CTD => "L", CTB => "L",
	TCY => "S",
	TCN => "S", TCM => "S", TCR => "S", TCW => "S", 
	TCS => "S", TCY => "S",	TCY => "S", TCK => "S", 
	TCV => "S", TCH => "S", TCD => "S", TCB => "S",
	TAY => "Y", 
	#TAR => "X",
	TAR => "*",
	TGT => "C",
	CCN => "P", CCM => "P", CCR => "P", CCW => "P", 
	CCS => "P", CCY => "P",	CCY => "P", CCK => "P", 
	CCV => "P", CCH => "P", CCD => "P", CCB => "P",
	CAY => "H",
	CAR => "Q",
	AGR => "R",
	CGN => "R", CGM => "R", CGR => "R", CGW => "R", 
	CGS => "R", CGY => "R",	CGY => "R", CGK => "R", 
	CGV => "R", CGH => "R", CGD => "R", CGB => "R",
	ATY => "I", ATM => "I", ATW => "I", ATH => "I",
	ACN => "T", ACM => "T", ACR => "T", ACW => "T", 
	ACS => "T", ACY => "T",	ACY => "T", ACK => "T", 
	ACV => "T", ACH => "T", ACD => "T", ACB => "T",
	AAY => "N",
	AAR => "K",
	GTN => "V", GTM => "V", GTR => "V", GTW => "V", 
	GTS => "V", GTY => "V",	GTY => "V", GTK => "V", 
	GTV => "V", GTH => "V", GTD => "V", GTB => "V",
	GCN => "A", GCM => "A", GCR => "A", GCW => "A", 
	GCS => "A", GCY => "A",	GCY => "A", GCK => "A", 
	GCV => "A", GCH => "A", GCD => "A", GCB => "A",
	GAY => "D",
	GAR => "E",
	GGN => "G", GGM => "G", GGR => "G", GGW => "G", 
	GGS => "G", GGY => "G",	GGY => "G", GGK => "G", 
	GGV => "G", GGH => "G", GGD => "G", GGB => "G",
);

# check that dependencies are present in user's PATH
my $blat = check_path_for_exec("blat");
#my $blast = check_path_for_exec("blastp");
my $blastp = check_path_for_exec("blastp");
my $blastn = check_path_for_exec("blastn");
my $muscle = check_path_for_exec("muscle");
my $clustalw2 = check_path_for_exec("clustalw2");
my $makeblastdb = check_path_for_exec("makeblastdb");
my $kaks_calc = check_path_for_exec("KaKs_Calculator");
my $protein_ortho = check_path_for_exec("proteinortho5.pl");

my %ks;
my $lower_ks;
my $upper_ks;
my $model = "YN";
my $threshold = 300;
my $min_contig_length = 300;
my @transcriptomes;
GetOptions(
	"ks|k:s"                 => \&get_ks_range,
	"transcriptomes|t:s{2,}" => \@transcriptomes,
	"help|h"                 => \&help,
);

#die "You must specify a Ks range to investigate.\n".&usage if (!defined($lower_ks) || !defined($upper_ks));
die "You need to specify at least two fasta files containing transcriptomes.\n".&usage if (scalar(@transcriptomes) < 2);
foreach my $transcriptome (@transcriptomes) {
	die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);
}

# translate each transcriptome into best frame
my @pids;
my @gene_pairs;
#my @proteomes;
#foreach my $transcriptome (@transcriptomes) {
foreach my $index (0 .. $#transcriptomes) {

	$lower_ks = $ks{"$index-LOWER"};
	$upper_ks = $ks{"$index-UPPER"};
	$lower_ks = (sort {$a <=> $b} keys %ks)[0] if (!defined($lower_ks)); 
	$upper_ks = (sort {$a <=> $b} keys %ks)[1] if (!defined($upper_ks)); 
	die "You must specify a Ks range to investigate.\n".&usage if (!defined($lower_ks) || !defined($upper_ks));

	my $transcriptome = $transcriptomes[$index];
	(my $base_name = $transcriptomes[$index]) =~ s/(.*)\..*/$1/;
	#(my $name = $transcriptome) =~ s/.*\/|\..*//g;
	#(my $proteome = $transcriptome) =~ s/(.*\/?)\Q$name\E(\..*)/$1$name-translated$2/;
	#my $proteome = $base_name."-translated.fasta";
	(my $proteome = $transcriptome) =~ s/\.cds$/.pep/;
	push(@gene_pairs, $base_name."_k$lower_ks-$upper_ks.fasta");

	#push(@proteomes, $proteome);

	my $pid = fork();
	if ($pid == 0) {

		my %align = parse_fasta($transcriptome);

#		open(my $out, ">", $proteome);
#		foreach my $contig (keys %align) {
#			print {$out} ">$contig\n";
#			print {$out} get_best_frame($align{$contig}),"\n";
#		}
#		close($out);

		# perform intraspecies blat for base transcriptome (FIRST in list)
		#system("$blat $proteomes[$index] $proteomes[$index] -prot -out=pslx -noHead $proteomes[$index].pslx >/dev/null");
		#system("$blat $proteome $proteome -prot -out=pslx -noHead $proteome.pslx >/dev/null");
		#system("$blat $proteome $proteome -prot -out=pslx -noHead $proteome.pslx");
		system("$blat $proteome $proteome -prot -out=pslx -noHead $proteome.pslx") if (!-e "$proteome.pslx");

		#system("$makeblastdb -in $proteome -input_type fasta -dbtype prot");

		#system("$blast -db $proteome -query $proteome -out=$proteome.out -outfmt '6 qseqid sseqid qseq sseq' -num_threads=20 -evalue 1e-05");
		#system("$blast -db $proteome -query $proteome -out=$proteome.out -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen slen' -num_threads=30 -evalue 1e-05");
		#system("$blast -db $proteome -query $proteome -out=$proteome.out -outfmt 6 -num_threads=30 -evalue 1e-05");

#		die;

		# parse blat output 
		#open(my $blat_output, "<", "$proteomes[0].pslx");
		open(my $blat_output, "<", "$proteome.pslx");
		#open(my $blat_output, "<", "$proteome.pslx");
		#open(my $blat_output, "<", "$proteome.out");
#		chomp(my @data = <$blat_output>);
#		close($blat_output);

		my $id = 0;
		my (%queries, %matches);
		#foreach my $line (@data) {
		#while (chomp(my $line = <$blat_output>)) {
		while (my $line = <$blat_output>) {
			chomp($line);

			# split the line on tabs and extract the information we want
			my @line = split("\t", $line);
			my ($query_name, $match_name, $query_align, $match_align) = 
				($line[9], $line[13], $line[21], $line[22]);
			#my ($query_name, $match_name, $query_align, $match_align, $query_len, $match_len) = 
	#		my ($query_name, $match_name, $query_start, $query_end, $match_start, $match_end) = 
	#			#($line[0], $line[1], $line[2], $line[3], $line[4], $line[5]);
	#			($line[0], $line[1], $line[6], $line[7], $line[8], $line[9]);
#			my ($query_name, $match_name, $query_align, $match_align, $query_len, $match_len) = 
#				($line[0], $line[1], $line[12], $line[13], $line[14], $line[15]);

#			$query_start--; $query_end--; # blast uses 1-based indexing and we don't want that
#			$match_start--; $match_end--; # blast uses 1-based indexing
#
#			my $query_range = "$query_start-$query_end";
#			my $match_range = "$match_start-$match_end";

			# exclude self match
			if ($query_name ne $match_name) {

				my @query_align =  split(",", $query_align);
				my @match_align =  split(",", $match_align);
				my $query_align_prot = join("", @query_align);

				my $query_nuc_align = $align{$query_name};
				my $match_nuc_align = $align{$match_name};

				my $trans_query_align;
				foreach my $align (@query_align) {
					$trans_query_align .= reverse_translate({"DNA" => $query_nuc_align, "PROT" => $align});
				}

				my $trans_match_align;
				foreach my $align (@match_align) {
					$trans_match_align .= reverse_translate({"DNA" => $match_nuc_align, "PROT" => $align});
				}

				$query_align = $trans_query_align;
				$match_align = $trans_match_align;

#				$query_align = reverse_translate({"DNA" => $query_nuc_align, "PROT" => $query_align});
#				$match_align = reverse_translate({"DNA" => $match_nuc_align, "PROT" => $match_align});

#				my $hit_frac = 0.3;

			#	if (length($query_align) >= $threshold && length($match_align) >= $threshold
			#		&& length($query_align) / $query_len >= $hit_frac && length($match_align) / $match_len >= $hit_frac) {
				if (length($query_align) >= $threshold && length($match_align) >= $threshold) {
					#&& length($query_align) / $query_len >= $hit_frac && length($match_align) / $match_len >= $hit_frac) {
					# check that the match hasn't already been extracted in reverse order
					if (!exists($queries{"$match_name-$query_name"})) {

						my $name = "q_$query_name"."_t_$match_name";
						my $pair = {'QUERY_ALIGN' => $query_align,
									'QUERY_ALIGN_PROT' => $query_align_prot,
									'MATCH_ALIGN' => $match_align,
									'LENGTH' => length($query_align)};
					#	my $pair = {'QUERY_RANGE', => $query_range,,
					#				'MATCH_RANGE' => $match_range};

						# check if there is already a match between these two sequences
						# if there is a match, the longer length one will be output
						if (exists($matches{$name})) {
							my $current_length = $matches{$name}->{'LENGTH'};
							if ($current_length <= length($query_align)) {
								$matches{$name} = $pair;
							}

		#					$matches{$name}->{'LENGTH'} += length($query_align);
		#					$matches{$name}->{'QUERY_ALIGN'} .= $query_align;
		#					$matches{$name}->{'MATCH_ALIGN'} .= $match_align;

		#					$matches{$name}->{'QUERY_RANGE'} = union($query_range, $matches{$name}->{'QUERY_RANGE'});
		#					$matches{$name}->{'MATCH_RANGE'} = union($match_range, $matches{$name}->{'MATCH_RANGE'});
						}
						else {
							$matches{$name} = $pair;
						}
						$queries{"$query_name-$match_name"}++;
					}
				}
			}
		}
		close($blat_output);
		#undef(@data);

		my @output;
		foreach my $key (sort { $a cmp $b} keys %matches) {
			my $pair = $matches{$key};

			(my $query_name = $key) =~ s/q_(.*)_t_.*/$1/;
			(my $match_name = $key) =~ s/q_.*_t_(.*)/$1/;

			push(@output, ">$key\n");
			push(@output, "$pair->{QUERY_ALIGN}\n");
			push(@output, "$pair->{MATCH_ALIGN}\n\n");
		#	my $meow = get_partial_seq($pair->{QUERY_RANGE}, $align{$query_name});
		#	my $meow2 = get_partial_seq($pair->{MATCH_RANGE}, $align{$match_name});
		#	push(@output, get_partial_seq($pair->{QUERY_RANGE}, $align{$query_name}),"\n");
		#	push(@output, get_partial_seq($pair->{MATCH_RANGE}, $align{$match_name}),"\n");
		#	print "$query_name $match_name: $pair->{QUERY_RANGE} $pair->{MATCH_RANGE} (".length($meow)." ".length($meow2).")\n";
			$id++;
		}

		die "No homologs meeting the nucleotide length threshold found.\n" if ($id == 0);
		print "$id pair(s) met the requirements.\n";

		#open(my $blat_atx, ">", "$base_transcriptome.atx");
		#open(my $blat_atx, ">", "$name.atx");
		open(my $blat_atx, ">", "$base_name.atx");
		print {$blat_atx} @output;
		close($blat_atx);


		# use KaKs_Calculator to determine Ks between homologus pairs
		#system("$kaks_calc -i $base_transcriptome.atx -o $base_transcriptome"."_kaks -m $model >/dev/null");
		#system("$kaks_calc -i $name.atx -o $name"."_kaks -m $model >/dev/null");
		#system("$kaks_calc -i $base_name.atx -o $base_name".".kaks -m $model >/dev/null");
		system("$kaks_calc -i $base_name.atx -o $base_name.kaks -m $model >/dev/null") if (!-e "$base_name.kaks");

		# for base transcriptome, subselect gene pairs with Ks between $lower_ks and $upper_ks
		#open(my $kaks_output, "<", $base_transcriptome."_kaks");
		#open(my $kaks_output, "<", $name."_kaks");
		open(my $kaks_output, "<", $base_name.".kaks");
#		chomp(@data = <$kaks_output>);
#		close($kaks_output);

		#shift(@data);

		my %genes;
		open(my $gene_subset, ">", $base_name."_k$lower_ks-$upper_ks.fasta");
		#foreach my $line (@data) {
		#while (chomp(my $line = <$kaks_output>)) {
		while (my $line = <$kaks_output>) {
			chomp($line);

			next if ($line =~ /^Sequence/);

			my @line = split(/\s+/, $line);
			my $ks = $line[3];
			$ks = 0 if ($ks eq "NA");

			if ($ks >= $lower_ks && $ks <= $upper_ks) {
				#my $pair_name = $line[0];
				(my $pair_name = $line[0]) =~ s/^>//;

				if ($pair_name =~ /q_(.*)_t_(.*)/) {
					my ($query_name, $match_name) = ($1, $2);
					print $query_name, " ", $match_name,"\n";
					#print {$gene_subset} ">$pair_name\n";
					#print {$gene_subset} get_best_frame($align{$query_name}),"\n";

					my $found_match = 0;
					#foreach my $key (sort { $a cmp $b } keys %matches) {
					foreach my $key (keys %matches) {
						#if ($key =~ /\Q$pair_name\E/) {
						if ($key eq $pair_name) {
							my $pair = $matches{$key};

							print {$gene_subset} ">$pair_name\n";
							#print {$gene_subset} "$pair->{QUERY_ALIGN}\n";
							#print {$gene_subset} get_best_frame($pair->{QUERY_ALIGN}),"\n";
							print {$gene_subset} $pair->{QUERY_ALIGN_PROT},"\n";
							$found_match++;
							last;
						}
					}
					die "Could not locate match: '$pair_name'\n" if (!$found_match);
				}
			}
		}
		close($kaks_output);
		#open(my $gene_subset, ">", $base_transcriptome."_k$lower_ks-$upper_ks.fasta");
		#open(my $gene_subset, ">", $name."_k$lower_ks-$upper_ks.fasta");
#		open(my $gene_subset, ">", $base_name."_k$lower_ks-$upper_ks.fasta");
#		foreach my $gene (keys %genes) {
#			print {$gene_subset} ">$gene\n";
#			print {$gene_subset} get_best_frame($align{$gene}),"\n";
#		}
		close($gene_subset);
#		unlink("$base_name-qt.pslx", "$base_name-q.fasta", "$base_name-t.fasta");

		exit(0);
	}
	else {
		push(@pids, $pid);
	}
}
foreach my $pid (@pids) {
	waitpid($pid, 0);
}
unlink(@pids);

#die;

# need to make it so gene names in proteinortho output have q_t_ format for easier parsing

# run proteinortho using subset of genes of base transcriptome against all other transcriptomes 
print "$protein_ortho @gene_pairs\n";
system("$protein_ortho @gene_pairs -clean");

# parse proteinortho output and create a new multisequence alignment for each hit with proteins common to all transcriptomes
my %families;
my $count = 0;
my @species_order;
open(my $ortho_output, "<", "myproject.proteinortho");
#chomp(my @ortho_output = <$ortho_output>);
#close($ortho_output);

#shift(@ortho_output);

#foreach my $line (@ortho_output) {
while (my $line = <$ortho_output>) {
	chomp($line);
	my @line = split("\t", $line);
	#my $num_species = $line[0];
	my ($num_species, $num_genes) = ($line[0], $line[1]);

	if ($. == 1) {
		@species_order = @line[3 .. $#line]; 
		next;
	}

	# only want matches which include all species
	#if ($num_species == scalar(@transcriptomes)) {
	if ($num_species == scalar(@transcriptomes) && $num_genes == scalar(@transcriptomes)) {
		#my @species = @line[3, $#line];	
#		my @species = @line[3 .. $#line];	
#		foreach my $index (0 .. $#species) {
#			my $species = $species[$index];
#			my @pairs = split(",", $species);
#			foreach my $pair (@pairs) {
#				if ($pair =~ /q_(.*)_t_(.*)/) {
#					my ($contig1, $contig2) = ($1, $2);
#
#					$families{$count}[$index]->{$contig1}++;
#					$families{$count}[$index]->{$contig2}++;
#				}
#			}
#		}
#		$count++;
		my @contigs = @line[3 .. $#line];	

		my $num_genes = 0;
		foreach my $index (0 .. $#contigs) {
			#my $contig = $contigs[$index];

			my $contigs = $contigs[$index];
			if ($contigs =~ /q_(.*)_t_(.*)/) {
				my ($contig1, $contig2) = ($1, $2);

				$contigs[$index] = [$contig1, $contig2];
			}

#			my @split_contig = split(",", $contig);
			#$num_genes += scalar(@split_contig);

			#$contigs[$index] = \@split_contig;
		}
		$families{$count} = \@contigs;
		$count++;

		print "$count: ";
		foreach my $contig (@contigs) {
			print "@{$contig} ";
		}
		print "\n";
	}
	else {
		#print "Ignoring match which was not found in all given transcriptomes.\n";
	}
}
close($ortho_output);

#reorient transcriptome array to match protein ortho output
#foreach my $species (@species_order) {
foreach my $index (0 .. $#species_order) {
	#(my $file = $species_order[$index]) =~ s/_k$lower_ks-$upper_ks\.fasta$/.cds/;
	(my $file = $species_order[$index]) =~ s/_k.*?-.*?\.fasta$/.cds/;
	$transcriptomes[$index] = $file;
}

mkdir("families") if (!-e "families");

my @unlink;
$SIG{CHLD} = 'IGNORE';
foreach my $family (sort { $a <=> $b } keys %families) {
#foreach my $family (1 .. 1) {
	my $members = $families{$family};
	#my @quartets = reduce_family_to_quartets($members);

#	foreach my $index (0 .. $#quartets) {
#		my $quartet = $quartets[$index];
		until(okay_to_run()) {};

		my $pid = fork();
		if ($pid == 0) {
			setpgrp();
			analyze_family({'ID' => $family, 'MEMBERS' => $members});
			#analyze_family({'ID' => $family."_$index", 'MEMBERS' => $quartet});
			exit(0);
		}
		else {
			push(@pids, $pid);
		}

#	}
	#die;
}

#foreach my $family (sort { $a <=> $b } keys %families) {
#	my @family = @{$families{$family}};
#	print "family: \n";
#	#foreach my $species (@family) {
#	foreach my $index (0 .. $#family) {
#		my $species = $family[$index];
#		print " species: \n";
#		foreach my $key (keys %{$species}) {
#			print "  ",$key,"\n";
#		}
#	}
#	print "\n";
#}

# align with clustalw2
# reverse translate alignments if originally nucleotide
# determine supported tree topology of the alignment 

sub combine {
	my ($list, $n) = @_;
	die "Insufficient list members" if ($n > @$list);

	return map [$_], @$list if ($n <= 1);

	my @comb;

	for (my $i = 0; $i+$n <= @$list; ++$i) {
		my $val  = $list->[$i];
		my @rest = @$list[$i + 1 .. $#$list];
		push(@comb, [$val, @$_]) for combine(\@rest, $n - 1);
	}

	return @comb;
}

sub analyze_family {
	my $args = shift;

	my $id = $args->{ID};
	my @members = @{$args->{MEMBERS}};
	#my $members = @{$args->{MEMBERS}};

	my %sequences;
#	foreach my $index (0 .. $#members) {
#		#(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\..*/$2/;
#		#print "base_name: $base_name\n";
#		my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $members[$index]});
#		#$sequences{$base_name} = $sequence;
#		#print ">$base_name\n$sequence\n\n";
#	}

	my $total_members = 0;
	foreach my $index (0 .. $#members) {

		my $member = $members[$index];
		(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\.fasta\.transdecoder\.cds/$2/;
		#foreach my $contig (@{$member}) {
		foreach my $index2 (0 .. scalar(@{$member}) - 1) {
			my $contig = @{$member}[$index2];
			my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $contig});
			$sequences{$base_name."_$index2"} = $sequence;
			#$names{$base_name."_$index2"} = $contig;
			$total_members++;
		}
	}

#	my $family_raw = "$id-raw.fasta";
#	my $family_aligned = "$id-aligned.nex";
#	my $family_reduced = "$id-reduced.fasta";
	my $family_raw = "families/$id-raw.fasta";
	my $family_aligned = "families/$id-aligned.nex";
	my $family_reduced = "families/$id-reduced.fasta";

	push(@unlink, $family_raw, $family_reduced, $family_raw.".nhr", $family_raw.".nin", $family_raw.".nsq", $family_raw.".out");
#	foreach my $run (1 .. $nruns) {
#		push(@unlink, "$family_aligned.run$run.t");
#		push(@unlink, "$family_aligned.run$run.p");
#	}
	$SIG{INT} = sub { unlink(@unlink); exit(0) };
	$SIG{KILL} = sub { unlink(@unlink); exit(0) };

	open(my $raw, ">", $family_raw);
	foreach my $taxon (sort { $a cmp $b } keys %sequences) {
		print $taxon,"\n";
		print {$raw} ">$taxon\n";
		print {$raw} "$sequences{$taxon}\n";
	}
	close($raw);

	#print "[$id] Reducing sequence to homologous sites.\n";

	reorient_contigs($id, \%sequences);
	%sequences = parse_fasta($family_raw);

	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

	my $has_minus;

	# remove matches in reverse direction
	foreach my $index (reverse(0 .. $#lines)) {
		my @line = split("\t", $lines[$index]);
		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		if ($s_strand eq "minus") {
			splice(@lines, $index, 1);	
			$has_minus++;
		}
	}
	#exit(0) if (!defined($has_minus));

	my %counts;
	#my %strands;
	my %matches;
	my %quartet;
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];

		my @line = split("\t", $line);
		#my ($query, $match, $q_start, $q_end) = ($line[0], $line[1], $line[2], $line[3]);
		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		next if ($query eq $match);

#		next if ($s_strand eq "minus"); ###################

		$counts{$query}++;

	#	if ($s_strand eq "minus" && !exists($strands{$query})) {
	#		$strands{$match}++;
	#	}

		my $next_line = " \t \t";
		if ($index + 1 <= $#lines) {
			$next_line = $lines[$index + 1];
		}
		my @next_line = split("\t", $next_line);
		my ($next_query, $next_match) = ($next_line[0], $next_line[1]);

		$q_start--; $q_end--; # blast uses 1-based indexing and we don't want that

		my $q_range = "$q_start-$q_end";

		if (exists($matches{"$query-$match"})) {
			$q_range = union($q_range, $matches{"$query-$match"});
		}

	#	print "$query-$match: $q_range ($s_strand)\n";
	#	print "  next : $next_query $next_match\n";
		if ("$query-$match" ne "$next_query-$next_match") {
			if (exists($quartet{$query})) {
				$quartet{$query} = intersection($q_range, $quartet{$query});	
			}
			else {
				$quartet{$query} = $q_range;
			}
		}
		#print "\n" if $next_query ne $query;
		$matches{"$query-$match"} = $q_range;
	}
	if (keys %quartet < 4) {
		unlink(@unlink);
		die "[$id] Could not identify homologous sites for a species in this family.\n";
	}
	foreach my $taxon (keys %counts) {
		my $count = $counts{$taxon};
		#if ($count != scalar(@transcriptomes) - 1) {
		#if ($count < scalar(@transcriptomes) - 1) {
		#if ($count < scalar(@members) - 1) {
		if ($count < $total_members - 1) {
			die "[$id] Some taxa did not share homology.\n";
			unlink(@unlink);
			exit(0);
		}
	}

#	foreach my $contig (keys %strands) {
#		print "$contig: $strands{$contig}\n";
#	}

	open(my $out, ">", $family_reduced);
	#foreach my $contig (keys %quartet) {
	foreach my $contig (sort { $a cmp $b } keys %quartet) {
		#print "frame: $frames{$contig}\n";

#		my $seq;
#		if ($strands{$contig}) {
#			$seq = rev_comp($sequences{$contig});
#		}
#		else {
#			$seq = $sequences{$contig};
#		}
		my $seq = $sequences{$contig};
		my $range = $quartet{$contig};

		print {$out} ">$contig\n";
#		print "[$id] $contig : $range\n";

		my $reduced_length = 0;
		foreach my $segment (split(",", $range)) {
			#print {$out} get_partial_seq($segment, $seq),"\n";
			my $sequence = get_partial_seq($segment, $seq);
			print {$out} "$sequence\n";

			$reduced_length += length($sequence);
		}
		if ($reduced_length < $min_contig_length) {
			unlink(@unlink);
			#die "[$id] A sequence did not meet the minimum length requirements.\n";
			die "[$id] A sequence did not meet the minimum length requirements ($reduced_length).\n";
		}
		print {$out} "\n";
	}
	close($out);

	print "[$id] Aligning reduced sites with muscle.\n";
	system("$muscle -in $family_reduced -out $family_aligned >/dev/null 2>&1") || die;
	#system("$clustalw2 $family_reduced -output=nexus -outfile=$family_aligned >/dev/null 2>&1") || die;

	my %family_aligned = parse_fasta($family_aligned);
	#write_nexus({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned});
	#write_phylip({'OUT' => "$family-aligned.phy", 'ALIGN' => \%aligned_quartet});
	write_phylip({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned});

	chdir("families");
	system("raxmlHPC -f a -m GTRGAMMA -p 123123 -x 123123 -#100 -s ../$family_aligned -n $id");

	unlink(@unlink);
}

sub reorient_contigs { 
	#my ($family_raw, $sequences) = (@_);
	my ($id, $sequences) = (@_);

	#my $family_raw = "$id-raw.fasta";
	my $family_raw = "families/$id-raw.fasta";
	#my $family_aligned = "$id-aligned.nex";
	#my $family_reduced = "$id-reduced.fasta";

	my %sequences = %{$sequences};
	#(my $id = $family_raw) =~ s/^(\d+)-.*/$1/;

	system("$makeblastdb -in $family_raw -input_type fasta -dbtype nucl >/dev/null") || die;
	#system("$blastn -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstrand' -evalue 1e-05 -perc_identity=50 >/dev/null") || die;
	system("$blastn -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstrand' -evalue 1e-05 >/dev/null") || die;

	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

	#my %matches;
	#my %counts;
	my %hits;
	my %quartet;
	my %strands;
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];

		my @line = split("\t", $line);
		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		next if ($query eq $match);

		$quartet{$query}++;

		if ($s_strand eq "minus" && !exists($strands{$query}) && !exists($hits{"$query-$match"})) {
		#if ($s_strand eq "minus" && !exists($strands{$query})) {
			$strands{$match}++;
		}
		$hits{"$query-$match"}++;
	}
#	if (keys %quartet < 4) {
#		#unlink(@unlink);
#		die "[$id] Could not identify homologous sites for a species in this family.\n";
#	}
	foreach my $taxon (keys %quartet) {
		my $count = $quartet{$taxon};
		#if ($count != scalar(@transcriptomes) - 1) {
		if ($count < scalar(@transcriptomes) - 1) {
			unlink(@unlink);
			#die "[$id] Some taxa did not share homology.\n";
			exit(0);
		}
	}

	if (keys %strands > 0) {
#		print "[$id] Reorienting contig(s).\n";
		open(my $raw, ">", $family_raw);
		foreach my $contig (keys %quartet) {

			my $seq;
			if ($strands{$contig}) {
				$seq = rev_comp($sequences{$contig});
			}
			else {
				$seq = $sequences{$contig};
			}

			print {$raw} ">$contig\n";
			print {$raw} "$seq\n";
		}
		close($raw);
		reorient_contigs($id, \%sequences);
	}
}

sub write_phylip {
	my $settings = shift;

	my $out_name = $settings->{'OUT'};
	my %align = %{$settings->{'ALIGN'}};

	my $ntaxa = scalar(values %align);
	my $nchar = length((values %align)[0]);

	open(my $out, ">", $out_name);
	print {$out} "$ntaxa $nchar\n";
	foreach my $taxon (sort {$a cmp $b} keys %align) {
		print {$out} "$taxon $align{$taxon}\n";
	}
	close($out);
}


sub get_best_frame {

	my $for_frame1 = shift;
	my $rev_frame1 = reverse($for_frame1);

	my $for_frame2 = substr($for_frame1, 1, length($for_frame1) - 1);
	my $for_frame3 = substr($for_frame2, 1, length($for_frame2) - 1);

	my $rev_frame2 = substr($rev_frame1, 1, length($rev_frame1) - 1);
	my $rev_frame3 = substr($rev_frame2, 1, length($rev_frame2) - 1);

	chop($for_frame1) foreach (1 .. length($for_frame1) % 3);
	chop($for_frame2) foreach (1 .. length($for_frame2) % 3); chop($for_frame3) foreach (1 .. length($for_frame3) % 3);

	chop($rev_frame1) foreach (1 .. length($rev_frame1) % 3);
	chop($rev_frame2) foreach (1 .. length($rev_frame2) % 3);
	chop($rev_frame3) foreach (1 .. length($rev_frame3) % 3);

	$for_frame1 = translate($for_frame1);
	$for_frame2 = translate($for_frame2);
	$for_frame3 = translate($for_frame3);

	$rev_frame1 = translate($rev_frame1);
	$rev_frame2 = translate($rev_frame2);
	$rev_frame3 = translate($rev_frame3);

#	print $for_frame1,"\n\n",$for_frame2,"\n\n",$for_frame3,"\n\n";
#	print $rev_frame1,"\n\n",$rev_frame2,"\n\n",$rev_frame3,"\n\n";

	my $for_frame1_longest_orf = get_longest_orf($for_frame1);
	my $for_frame2_longest_orf = get_longest_orf($for_frame2);
	my $for_frame3_longest_orf = get_longest_orf($for_frame3);

	my $rev_frame1_longest_orf = get_longest_orf($rev_frame1);
	my $rev_frame2_longest_orf = get_longest_orf($rev_frame2);
	my $rev_frame3_longest_orf = get_longest_orf($rev_frame3);

#	print $for_frame1_longest_orf,"\n",$for_frame2_longest_orf,"\n",$for_frame3_longest_orf,"\n";
#	print $rev_frame1_longest_orf,"\n",$rev_frame2_longest_orf,"\n",$rev_frame3_longest_orf,"\n";

	my %orfs;
	determine_longest({"ORFS" => \%orfs, "FRAME" => $for_frame1, "FRAME_LONGEST_ORF" => $for_frame1_longest_orf});
	determine_longest({"ORFS" => \%orfs, "FRAME" => $for_frame2, "FRAME_LONGEST_ORF" => $for_frame2_longest_orf});
	determine_longest({"ORFS" => \%orfs, "FRAME" => $for_frame3, "FRAME_LONGEST_ORF" => $for_frame3_longest_orf});

	determine_longest({"ORFS" => \%orfs, "FRAME" => $rev_frame1, "FRAME_LONGEST_ORF" => $rev_frame1_longest_orf});
	determine_longest({"ORFS" => \%orfs, "FRAME" => $rev_frame2, "FRAME_LONGEST_ORF" => $rev_frame2_longest_orf});
	determine_longest({"ORFS" => \%orfs, "FRAME" => $rev_frame3, "FRAME_LONGEST_ORF" => $rev_frame3_longest_orf});

	my $longest = (sort {$b <=> $a} keys %orfs)[0];
	#print "$longest ",$orfs{$longest},"\n\n";

	#substitute out * for X to allow for reverse translation from blat
	(my $best_frame = $orfs{$longest}) =~ s/\*/X/g;

	#return $orfs{$longest};
	return $best_frame;

#	print translate($for_frame1),"\n\n",translate($for_frame2),"\n\n",translate($for_frame3),"\n\n";
#	print translate($rev_frame1),"\n\n",translate($rev_frame2),"\n\n",translate($rev_frame3),"\n\n";
}

sub determine_longest {
	my $settings = shift;

	my %orfs = %{$settings->{"ORFS"}};
	my $frame = $settings->{"FRAME"};
	my $frame_longest_orf = $settings->{"FRAME_LONGEST_ORF"};

	if (!exists($orfs{$frame_longest_orf})) {
		#$orfs{$frame_longest_orf} = $frame;
		$settings->{"ORFS"}{$frame_longest_orf} = $frame;
	}
	else {
		#print "entry exists for $frame_longest_orf, breaking tie\n";
		my $new_num_stops = 0;
		$new_num_stops++ while ($frame =~ /\*/g);
		#$new_num_stops++ while ($frame =~ /X/g);

		my $old_num_stops = 0;
		$old_num_stops++ while ($orfs{$frame_longest_orf} =~ /\*/g);
		#$old_num_stops++ while ($orfs{$frame_longest_orf} =~ /X/g);

		#print "new frame has $new_num_stops, old frame has $old_num_stops\n";

		$settings->{"ORFS"}{$frame_longest_orf} = $frame if ($new_num_stops < $old_num_stops);
		#$orfs{$frame_longest_orf} = $frame if ($new_num_stops < $old_num_stops);
	}
}

sub get_longest_orf {
	my ($frame) = @_;

	$frame =~ s/\*+/*/g;
	#$frame =~ s/X+/X/g;

	my $longest = 0;
	#while ($frame =~ /[^.|\*](.*?)\*/g) {
	while ($frame =~ s/^(.*?\*)//g) {
	#while ($frame =~ s/^(.*?X)//g) {
		my $length = length($1) - 1;
		$longest = $length if ($length > $longest);
	}
	my $length = length($frame) - 1;
	$longest = $length if ($length > $longest);
	#print $longest,"\n";

	return $longest;
}

sub translate {
	my $frame = shift;

	my $translated;
	for (my $i = 0; $i < length($frame); $i += 3) {
		my $codon = substr($frame, $i, 3);
		#print $codon,"\n" if (!exists($codon_table{uc($codon)}));

		if (exists($codon_table{uc($codon)})) {
			$translated .= $codon_table{uc($codon)};
		}
		else {
			$translated .= "X";
		}
		#print $codon,"\n" if (!exists($codon_table{uc($codon)}));
		#$translated .= $codon_table{uc($codon)};
	}

	return $translated;
}

sub get_ks_range {
	my ($arg, $ks_range) = @_;

#	my @ks_range = split(",", $ks_range);
#	@ks_range = sort { $a <=> $b } @ks_range;
#
#	$lower_ks = $ks_range[0];
#	$upper_ks = $ks_range[1];

	my @ks_ranges = split(",", $ks_range);
	#foreach my $ks_range (@ks_ranges) {
	foreach my $index (0 .. $#ks_ranges) {
		my $ks_range = $ks_ranges[$index];
		my @ks_range = split("-", $ks_range);
		@ks_range = sort { $a <=> $b } @ks_range;

		$ks{"$index-LOWER"} = $ks_range[0];
		$ks{"$index-UPPER"} = $ks_range[1];

#		$lower_ks = $ks_range[0];
#		$upper_ks = $ks_range[1];
	}

	return;
}

sub parse_fasta {
	my $filename = shift;

	my %align;
	open(my $alignment_file, '<', $filename) 
		or die "Could not open '$filename': $!\n";
	chomp(my @data = <$alignment_file>);
	close($alignment_file);
	
	my $taxon;
	foreach my $line (@data) {
		if ($line =~ /^>(\S+)/) {
			$taxon = $1;
		}
		else {
			$taxon =~ s/-/_/g;
			#$align->{$taxon} .= $line;
			$align{$taxon} .= $line;
		}
	}
	return %align;
}

sub check_path_for_exec {
	my $exec = shift;
	
	my $path = $ENV{PATH}.":."; #include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec);
	}
	die "Could not find the following executable: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path));
	return $exec_path;
}

sub is_nucl_align {
	my $align = shift;
	my %align = %{$align};

	my $total = 0;
	my $count = 0;
	foreach my $align (values %align) {
		$count++ while ($align =~ /A|T|C|G|-|\?|N/g);
		$total += length($align);
	}

	print "",($count / $total),"\n";

	return ($count / $total > 0.9);
}

#sub reverse_translate {
#	my $settings = shift;
#
#	my $dna = $settings->{'DNA'};
#	my $prot = $settings->{'PROT'};
#
#	my $regex;
#	foreach my $index (0 .. length($prot) - 1) {
#		my $char = substr($prot, $index, 1);
#		$regex .= $rev_codon_table{$char};
#	}
#
#	my $translation;
#	if ($dna =~ /($regex)/) {
#		$translation = $1;
#	}
#	elsif (reverse($dna) =~ /($regex)/) {
#		$translation = $1;
#	}
#	else {
#		die "Protein sequence could not be reverse translated.\n";
#	}
#
#	return $translation;
#}

sub reverse_translate {
	my $settings = shift;

	my $dna = $settings->{'DNA'};
	my $prot = $settings->{'PROT'};

	my @gaps;

	my $regex;
	foreach my $index (0 .. length($prot) - 1) {
		my $char = substr($prot, $index, 1);
		if ($char ne "-") {
			$regex .= $rev_codon_table{$char};
		}
		else {
			push(@gaps, $index * 3);
		}
	}

	my $translation;
	if ($dna =~ /($regex)/) {
		$translation = $1;
	}
	elsif (reverse($dna) =~ /($regex)/) {
		$translation = $1;
	}
	else {
		die "Protein sequence could not be reverse translated.\n";
	}

	my $count = 0;
	foreach my $gap (@gaps) {
		#$translation = substr($translation, 0, $gap)."---".substr($translation, $gap, length($translation) - 1);
		$translation = substr($translation, 0, $gap)."---".substr($translation, $gap, length($translation) - $gap);
		#print $translation,"\n";
		$count++;
		#die if $count > 5;
	}

#	print "\nPROT:  $prot\n\n";
#	print "DNA:   $dna\n\n";
#	print "TRANS: $translation\n\n";

	return $translation;
}

sub union {
	my ($add, $old) = (@_);
	
	my $union;

	my @add = split(",", $add);
	my @old = split(",", $old);

	my %old;
	foreach my $segment (@old) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}
	foreach my $segment (@add) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}

	my @sites = sort { $a <=> $b } keys %old;

	my @ends;
	my @starts = ($sites[0]);
	my $previous = $sites[0] - 1;
	#foreach my $site (@sites) {
	foreach my $index (0 .. $#sites) {
		my $site = $sites[$index];
		if ($previous != $site - 1) {
			push(@ends, $sites[$index - 1]);
			push(@starts, $site);
		}
		$previous = $site;
	}
	push(@ends, $sites[$#sites]);

	#print "starts: @starts\n";
	#print "ends:   @ends\n";
	foreach my $index (0 .. $#starts) {
		$union .= $starts[$index]."-".$ends[$index].",";
	}
	chop($union);

	#print "$add + $old = $union\n";

	return $union;
}

sub intersection {
	my ($new, $old) = (@_);

	my @new = split(",", $new);
	my @old = split(",", $old);

	my $intersection;

	my %old;
	foreach my $segment (@old) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}
	foreach my $segment (@new) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}

	foreach my $site (keys %old) {
		if ($old{$site} != 2) {
			delete($old{$site});
		}
	}

	my @sites = sort { $a <=> $b } keys %old;

	return "0-0" if (scalar(@sites) == 0);

	my @ends;
	my @starts = ($sites[0]);
	my $previous = $sites[0] - 1;
	#foreach my $site (@sites) {
	foreach my $index (0 .. $#sites) {
		my $site = $sites[$index];
		if ($previous != $site - 1) {
			push(@ends, $sites[$index - 1]);
			push(@starts, $site);
		}
		$previous = $site;
	}
	push(@ends, $sites[$#sites]);

	foreach my $index (0 .. $#starts) {
		$intersection .= $starts[$index]."-".$ends[$index].",";
	}
	chop($intersection);

	#print "$new intersect $old = $intersection\n";

	return $intersection;
}

sub get_partial_seq {
	my ($range, $seq) = (@_);

	my $partial_seq;

	my @range = split(",", $range);
	foreach my $segment (@range) {
		my ($start, $end) = split("-", $segment);
		$partial_seq .= substr($seq, $start, $end - $start + 1);	
	}
	
	return $partial_seq;
}

sub get_seq_from_fasta {
	my $args = shift;

	my $file = $args->{'FILE'};
	my $taxon = $args->{'TAXON'};
	#$taxon =~ s/-/_/g;

	my $seq;
	my $is_desired_seq;

	open(my $align, "<", $file);
	while (my $line = <$align>) {
		chomp($line);
		if ($line =~ /^>\Q$taxon\E/) {
			$is_desired_seq++;
		}
		elsif ($line =~ /^>/ && $is_desired_seq) {
			last;
		}
		elsif ($is_desired_seq) {
			$seq .= $line;
		}
	}
	close($align);

	return $seq;
}

sub rev_comp {
	my $seq = shift;
	
	my %comp = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C');

	my $rev_comp;
	my $rev = reverse($seq);
	foreach my $index (0 .. length($rev) - 1) {
		my $char = substr($rev, $index, 1);
		$rev_comp .= $comp{$char};
	}

	return $rev_comp;
}

sub okay_to_run {
	my $running_forks = 0;
	foreach  my $index (reverse(0 .. $#pids)) {
		my $pid = $pids[$index];
		if (kill 0, $pid) {
			$running_forks++;
		}
		else {
			splice(@pids, $index, 1);
		}
	}
	return ($running_forks < $max_forks);
}

sub help {
	return "";
}

sub usage {
	return "";
}

