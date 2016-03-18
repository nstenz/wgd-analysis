#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Fcntl qw(:flock SEEK_END);

# For reverse translating peptides
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

# For multithreading
my $free_cpus = get_free_cpus();

# Cutoff values
my $cutoff = 70; # bootstrap support for RAxML
#my $pid_cutoff = 80; # amino acid level
my $pid_cutoff = 60; # amino acid level
my $max_forks = 40;

my $output_dir = "families.wgd-plot";

# Check that dependencies are present in user's PATH
my $blat = check_path_for_exec("blat");
my $blastp = check_path_for_exec("blastp");
my $muscle = check_path_for_exec("muscle");
my $makeblastdb = check_path_for_exec("makeblastdb");
my $kaks_calc = check_path_for_exec("KaKs_Calculator");
my $protein_ortho = check_path_for_exec("proteinortho5.pl");

# Check that a version of Transdecoder is present in user's PATH
my $transdecoder = check_path_for_exec("TransDecoder", 1);
my $transdecoder_orfs = check_path_for_exec("TransDecoder.LongOrfs", 1);
my $transdecoder_predict = check_path_for_exec("TransDecoder.Predict", 1);
if (!defined($transdecoder) || (!defined($transdecoder_orfs) && !defined($transdecoder_predict))) {
	die "Could not locate required TransDecoder executables (TransDecoder or TransDecoder.LongOrfs and TransDecoder.Predict) in PATH.\n" 
}

# Pfam settings
my $pfam_search;
my $pfam_cpus = $free_cpus;

# Allowed models 
my %models = (NG => 1, LWL => 1, LPB => 1, MLWL => 1, MLPB => 1, 
              GY => 1, YN => 1, MYN => 1, MS => 1, MA => 1);

# Default settings
my %ks;
my $lower_ks;
my $upper_ks;
my $model = "YN";
my $threshold = 300;
my $min_contig_length = 300;
my @transcriptomes;
GetOptions(
#	"ks|k:s"                 => \&get_ks_range,
	"help|h"                  => \&help,
	"model|m=s"               => \$model,
	"output|o=s"              => \$output_dir,
	"pfam-cpus|p=i"           => \$pfam_cpus,
	"pfam-search|s=s"         => \$pfam_search,
	"transcriptomes|t:s{2,2}" => \@transcriptomes,
);

# Check that hmmscan can be located if user wants pfam search
my $hmmscan = check_path_for_exec("hmmscan") if ($pfam_search);

# Error check input
die "Unknown arguments specified.\n".&help if (@ARGV);
die "You need to specify two fasta files containing transcriptomes.\n".&usage if (scalar(@transcriptomes) != 2);
foreach my $transcriptome (@transcriptomes) {
	die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);
}

# Translate each transcriptome into best frame
my @pids;
my @gene_pairs;

foreach my $index (0 .. $#transcriptomes) {
#	$lower_ks = $ks{"$index-LOWER"};
#	$upper_ks = $ks{"$index-UPPER"};

	# Needed?
#	$lower_ks = (sort {$a <=> $b} keys %ks)[0] if (!defined($lower_ks)); 
#	$upper_ks = (sort {$a <=> $b} keys %ks)[1] if (!defined($upper_ks)); 
#
#	die "You must specify a Ks range to investigate.\n".&usage if (!defined($lower_ks) || !defined($upper_ks));

	my $transcriptome = $transcriptomes[$index];
	#(my $base_name = $transcriptomes[$index]) =~ s/(.*)\..*/$1/;
	(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\..*/$2/;

	#push(@gene_pairs, $base_name.".ks_$lower_ks-$upper_ks.fasta");
	push(@gene_pairs, $base_name.".pairs.fasta");

	my $pid = fork();
	if ($pid == 0) {

		my $proteome = $transcriptome.".transdecoder.pep";

		if (!-e $proteome) {
			print "Running TransDecoder on '$transcriptome'...\n";
			my $return = system("$transdecoder -t $transcriptome >/dev/null");
			die "Error occurred while running TransDecoder: '$return'.\n" if ($return);
			print "Finished TransDecoder for '$transcriptome.\n";
		}
		else {
			print "Using TransDecoder output in file '$proteome'...\n";
		}

		my %align = parse_fasta("$transcriptome.transdecoder.cds");

		# Perform self-blat
		my $return = system("$blat $proteome $proteome -prot -out=pslx -noHead $proteome.pslx") if (!-e "$proteome.pslx");
		die "Error running self-blat: '$return'.\n" if ($return);

		# Parse blat output 
		my (%queries, %matches);
		open(my $blat_output, "<", "$proteome.pslx");
		while (my $line = <$blat_output>) {
			chomp($line);

			# split the line on tabs and extract the information we want
			my @line = split("\t", $line);
			#my ($query_name, $match_name, $query_align, $match_align) = 
			#	($line[9], $line[13], $line[21], $line[22]);
			my ($match, $mismatch, $query_name, $match_name, $query_align, $match_align) = 
				($line[0], $line[1], $line[9], $line[13], $line[21], $line[22]);

			next if ("$query_name" eq "$match_name");

			# exclude self match
			my $match_length = $match + $mismatch;
			#if ($query_name ne $match_name) {

				next if ($match_length * 3 < $threshold);

				my @query_align =  split(",", $query_align);
				my @match_align =  split(",", $match_align);
				my $query_align_prot = join("", @query_align);

				my $query_nuc_align = $align{$query_name};
				my $match_nuc_align = $align{$match_name};

				# Reverse translate protein matches to corresponding nucleotide sequences
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

				my @names;
				push(@names, $query_name);
				push(@names, $match_name);
				@names = sort { "$a" cmp "$b" } @names;

				#if (length($query_align) >= $threshold && length($match_align) >= $threshold) {
					#&& length($query_align) / $query_len >= $hit_frac && length($match_align) / $match_len >= $hit_frac) {

				#	# Check that the match hasn't already been extracted in reverse order
					if (!exists($queries{"$match_name-$query_name"})) {

						my $name = "q_$query_name"."_t_$match_name";
						#my $name = "q_$names[0]"."_t_$names[1]";
						my $pair = {'QUERY_ALIGN' => $query_align,
									'QUERY_ALIGN_PROT' => $query_align_prot,
									'MATCH_ALIGN' => $match_align,
									'LENGTH' => length($query_align)};

						# Check if there is already a match between these two sequences
						# if there is a match, the longer one will be output
						my $current_hit = $matches{$name};
						#if (exists($matches{$name})) {
						if (defined($current_hit)) {
							#my $current_length = $matches{$name}->{'LENGTH'};
							my $current_length = $current_hit->{'LENGTH'};
							if ($current_length <= length($query_align)) {
								$matches{$name} = $pair;
							}
						}
						else {
							$matches{$name} = $pair;
						}
						$queries{"$query_name-$match_name"}++;
					}
				#}
			#}
		}
		close($blat_output);

		my $num_matches = scalar(keys %matches);
		die "No homologs meeting the nucleotide length threshold found.\n" if (!$num_matches);
		print "$num_matches pair(s) met the requirements.\n";

		# Write .atx input file for KaKs_Calculator
		open(my $blat_atx, ">", "$base_name.atx");
		foreach my $key (sort { $a cmp $b} keys %matches) {
			my $pair = $matches{$key};

			(my $query_name = $key) =~ s/q_(.*)_t_.*/$1/;
			(my $match_name = $key) =~ s/q_.*_t_(.*)/$1/;

			print {$blat_atx} ">$key\n";
			print {$blat_atx} "$pair->{QUERY_ALIGN}\n";
			print {$blat_atx} "$pair->{MATCH_ALIGN}\n\n";
		}
		close($blat_atx);

		# use KaKs_Calculator to determine Ks between homologus pairs
		$return = system("$kaks_calc -i $base_name.atx -o $base_name.kaks -m $model >/dev/null") if (!-e "$base_name.kaks");
		die "Error running KaKs_Calculator: '$return'.\n" if ($return);

		# Select gene pairs with Ks between $lower_ks and $upper_ks
		my %genes;

		open(my $kaks_output, "<", $base_name.".kaks");
		#open(my $gene_subset, ">", $base_name.".ks_$lower_ks-$upper_ks.fasta");
		open(my $gene_subset, ">", $base_name.".pairs.fasta");

		while (my $line = <$kaks_output>) {
			chomp($line);

			next if ($line =~ /^Sequence/);

			my @line = split(/\s+/, $line);
			my $ks = $line[3];
			$ks = 0 if ($ks eq "NA");

	#		if ($ks >= $lower_ks && $ks <= $upper_ks) {
				(my $pair_name = $line[0]) =~ s/^>//;

				if ($pair_name =~ /q_(.*)_t_(.*)/) {
					my ($query_name, $match_name) = ($1, $2);
#					print $query_name, " ", $match_name,"\n";

					#my $found_match = 0;
#					foreach my $key (keys %matches) {
#						if ($key eq $pair_name) {
#							my $pair = $matches{$key};
#							print {$gene_subset} ">$pair_name\n";
#							print {$gene_subset} $pair->{QUERY_ALIGN_PROT},"\n";
#							$found_match++;
#							last;
#						}
#					}
					my $pair = $matches{$pair_name};
					print {$gene_subset} ">$pair_name\n";
					print {$gene_subset} $pair->{QUERY_ALIGN_PROT},"\n";
					#die "Could not locate match: '$pair_name'\n" if (!$found_match);
				}
	#		}
		}
		close($kaks_output);
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
#unlink(@pids);

# need to make it so gene names in proteinortho output have q_t_ format for easier parsing

# run proteinortho using subset of genes of base transcriptome against all other transcriptomes 
system("$protein_ortho @gene_pairs -clean -project=wgd-plot") if (!-e "wgd-plot.proteinortho");
#system("$protein_ortho @gene_pairs -clean");

# parse proteinortho output and create a new multisequence alignment for each hit with proteins common to all transcriptomes
my %families;
my $count = 0;
my @species_order;
open(my $ortho_output, "<", "wgd-plot.proteinortho");
while (my $line = <$ortho_output>) {
	chomp($line);

	my @line = split("\t", $line);
	my ($num_species, $num_genes) = ($line[0], $line[1]);

	if ($. == 1) {
		@species_order = @line[3 .. $#line]; 
		next;
	}

	# only want matches which include all species
	#if ($num_species == scalar(@transcriptomes)) {
	if ($num_species == scalar(@transcriptomes) && $num_genes == scalar(@transcriptomes)) {

		my @contigs = @line[3 .. $#line];	

		#my $num_genes = 0;
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
	#(my $file = $species_order[$index]) =~ s/_k.*?-.*?\.fasta$/.cds/;

	# Nucleotide-based
	#(my $file = $species_order[$index]) =~ s/\.ks_.*?-.*?\.fasta$/.fasta.transdecoder.cds/;
	#(my $file = $species_order[$index]) =~ s/\.ks_.*?-.*?\.fasta$/.fasta.transdecoder.pep/;
	(my $file = $species_order[$index]) =~ s/\.pairs\.fasta$/.fasta.transdecoder.pep/;
	$transcriptomes[$index] = $file;
}

mkdir($output_dir) if (!-e $output_dir);
unlink("against.txt", "support.txt");

my @unlink;
$SIG{CHLD} = 'IGNORE';
foreach my $family (sort { $a <=> $b } keys %families) {
#foreach my $family (92 .. 92) {
	my $members = $families{$family};

	#$max_forks = 1;
	until(okay_to_run()) {};

	my $pid = fork();
	if ($pid == 0) {
		setpgrp();
		analyze_family({'ID' => $family, 'MEMBERS' => $members});
		exit(0);
	}
	else {
		push(@pids, $pid);
	}
	#die;
}

# Create the plot using R
#chdir($output_dir);


sub analyze_family {
	my $args = shift;

	my $id = $args->{ID};
	my @members = @{$args->{MEMBERS}};

	my %names;
	my %sequences;
	my %nuc_sequences;
	my %original_names;
	my $total_members = 0;
	foreach my $index (0 .. $#members) {
		my $member = $members[$index];
		#(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\.fasta\.transdecoder\.cds/$2/;
		(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\.fa(sta)?\.transdecoder\.pep/$2/;
		(my $cds_file = $transcriptomes[$index]) =~ s/\.pep$/.cds/;
		#(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*?)\..*/$2/;
		foreach my $index2 (0 .. scalar(@{$member}) - 1) {
			my $contig = @{$member}[$index2];
			my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $contig});
			my $nuc_sequence = get_seq_from_fasta({'FILE' => $cds_file, 'TAXON' => $contig});
			#my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index].".transdecoder.cds", 'TAXON' => $contig});
			$sequences{$base_name."_$index2"} = $sequence;
			$nuc_sequences{$base_name."_$index2"} = $nuc_sequence;
			$names{$base_name}++;
			$original_names{$base_name."_$index2"} = $contig;
			$total_members++;
		}
	}

	my $family_raw = "$output_dir/$id-raw.fasta";
	my $family_aligned = "$output_dir/$id-aligned.nex";
	my $family_reduced = "$output_dir/$id-reduced.fasta";

	#push(@unlink, $family_raw, $family_reduced, $family_raw.".nhr", $family_raw.".nin", $family_raw.".nsq", $family_raw.".out");
	push(@unlink, $family_raw, $family_reduced, $family_raw.".phr", $family_raw.".pin", $family_raw.".psq", $family_raw.".out");

	$SIG{INT} = sub { unlink(@unlink); exit(0) };

	open(my $raw, ">", $family_raw);
	foreach my $taxon (sort { $a cmp $b } keys %sequences) {
		print {$raw} ">$taxon\n";
		print {$raw} "$sequences{$taxon}\n";
	}
	close($raw);

	#print "[$id] Reducing sequence to homologous sites.\n";
	system("$makeblastdb -in $family_raw -input_type fasta -dbtype prot >/dev/null") || die;
	#system("$blastp -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend' -evalue 1e-05 >/dev/null") || die;
	system("$blastp -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstart send' -evalue 1e-05 >/dev/null") || die;
	#system("$blastp -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend' -evalue 1e-10 >/dev/null") || die;

	#reorient_contigs($id, \%sequences);
	#%sequences = parse_fasta($family_raw);

	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

#	# remove matches in reverse direction
#	foreach my $index (reverse(0 .. $#lines)) {
#		my @line = split("\t", $lines[$index]);
#		my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
#		if ($s_strand eq "minus") {
#			splice(@lines, $index, 1);	
#		}
#	}

	my %counts;
	#my %strands;
	my %matches;
	my %quartet;
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];

		my @line = split("\t", $line);
		#my ($query, $match, $q_start, $q_end, $s_strand) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		#my ($query, $match, $q_start, $q_end) = ($line[0], $line[1], $line[2], $line[3]);
		my ($query, $match, $q_start, $q_end, $s_start, $s_end) = ($line[0], $line[1], $line[2], $line[3], $line[4], $line[5]);
		next if ($query eq $match);

		$counts{$query}++;

		my $next_line = " \t \t";
		if ($index + 1 <= $#lines) {
			$next_line = $lines[$index + 1];
		}
		my @next_line = split("\t", $next_line);
		my ($next_query, $next_match) = ($next_line[0], $next_line[1]);

		$q_start--; $q_end--; # blast uses 1-based indexing and we don't want that
		$s_start--; $s_end--;

		my $q_range = "$q_start-$q_end";
		my $s_range = "$s_start-$s_end";

		if (!exists($matches{"$query-$match"})) {
			$matches{"$query-$match"} = {$q_range => $s_range};
		}
		else {

			my $current_match_q_range = (keys %{$matches{"$query-$match"}})[0];
			my $current_match_s_range = (values %{$matches{"$query-$match"}})[0];

			my $s_intersection = intersection($current_match_s_range, $s_range);
			if ($s_intersection eq "0-0") {
				my $q_intersection = intersection($current_match_q_range, $q_range);

				if ($q_intersection eq "0-0") {
					$q_range = union($q_range, $current_match_q_range);
					$s_range = union($s_range, $current_match_s_range);
					$matches{"$query-$match"} = {$q_range => $s_range};
				}
			}

		}

		if ($query.$match ne $next_query.$next_match) {
			if (!exists($quartet{$query})) {
				$quartet{$query} = (keys %{$matches{"$query-$match"}})[0];
			}
			else {
				$quartet{$query} = intersection((keys %{$matches{"$query-$match"}})[0], $quartet{$query});	
			}
		}
	}
	if (scalar(keys %quartet) < 4) {
	#if (scalar(keys %quartet) < scalar(@transcriptomes) * 2) {
		unlink(@unlink);
		die "[$id] Could not identify homologous sites for a species in this family.\n";
	}

	foreach my $taxon (keys %counts) {
		my $count = $counts{$taxon};
		if ($count < $total_members - 1) {
			die "[$id] Some taxa did not share homology.\n";
			unlink(@unlink);
			exit(0);
		}
	}

	# Output unaligned sequences of family
	open(my $out, ">", $family_reduced);
	foreach my $contig (sort { $a cmp $b } keys %quartet) {

		my $seq = $sequences{$contig};
		my $range = $quartet{$contig};

		print {$out} ">$contig\n";

		my $reduced_length = 0;
		foreach my $segment (split(",", $range)) {
			my $sequence = get_partial_seq($segment, $seq);
			print {$out} "$sequence\n";

			$reduced_length += length($sequence);
		}

		# Nucleotides
		#if ($reduced_length < $min_contig_length) {

		# AA's
		if ($reduced_length < $min_contig_length / 3) {
			unlink(@unlink);
			die "[$id] A sequence did not meet the minimum length requirements ($reduced_length).\n";
		}
		print {$out} "\n";
	}
	close($out);

	print "[$id] Aligning reduced sites with muscle.\n";
	system("$muscle -in $family_reduced -out $family_aligned >/dev/null 2>&1") || die;

	my %family_aligned = parse_fasta($family_aligned);

	foreach my $index1 (0 .. scalar(keys %family_aligned) - 2) {
		my $taxon1 = (keys %family_aligned)[$index1];
		foreach my $index2 ($index1 + 1 .. scalar(keys %family_aligned) - 1) {
			my $taxon2 = (keys %family_aligned)[$index2];
			my $pid = get_pid($family_aligned{$taxon1}, $family_aligned{$taxon2});

			die "[$id] Pairwise similarity between sequences was too low ($pid).\n" if ($pid < $pid_cutoff);
		}
	}

	my %family_aligned_nuc;
	foreach my $taxon (keys %family_aligned) {
		$family_aligned_nuc{$taxon} = reverse_translate({"DNA" => $nuc_sequences{$taxon}, "PROT" => $family_aligned{$taxon}});
	}

	#write_phylip({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned});
	write_phylip({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned_nuc});

	chdir("$output_dir");
	system("raxmlHPC -f a -m GTRGAMMA -p 123123 -x 123123 -#100 -s ../$family_aligned -n $id >/dev/null");
	#system("raxmlHPC -f a -m PROTGAMMALG -p 123123 -x 123123 -#100 -s ../$family_aligned -n $id >/dev/null");
	#chdir("..");

	#open(my $raxml_out, "<", "families/RAxML_bipartitions.$id");
	open(my $raxml_out, "<", "RAxML_bipartitions.$id");
	chomp(my $tree = <$raxml_out>);
	close($raxml_out);

	my @names = keys %names;
	my $taxon_1 = $names[0];
	my $taxon_2 = $names[1];

	# Remove branch lengths
	(my $tree_no_bls = $tree) =~ s/:\d+(\.\d+)?//g;

	my $out_file_name;
	if ($tree_no_bls =~ /.*\((\S+?)_\d,(\S+?)_\d\)(\d+)/) {
		my $taxon_1 = $1;
		my $taxon_2 = $2;
		my $bootstrap = $3;

		print "[$id] tree: $tree_no_bls\n";
		print "[$id] $taxon_1 and $taxon_2 are in a clade with $bootstrap% support.\n";

		die "[$id] RAxML tree does not have high enough support.\n" if ($bootstrap < $cutoff);

		if ($taxon_1 eq $taxon_2) {
			$out_file_name = "against.txt";
		}
		else {
			$out_file_name = "support.txt";
		}
	}

	foreach my $name (@names) {

		my %contigs;
		foreach my $contig (keys %original_names) {
			if ($contig =~ /\Q$name\E_\d+/) {
				#$contigs{$contig} = $original_names{$contig};
				$contigs{$original_names{$contig}} = $contig;
			}
		}

		open(my $topo_output, ">>", "$name.csv");
		flock($topo_output, LOCK_EX) || die "Could not lock '$topo_output': $!.\n";
		seek($topo_output, 0, SEEK_END) || die "Could not seek '$topo_output': $!.\n";

		open(my $kaks_output, "<", "../$name.kaks");
		while (my $line = <$kaks_output>) {
			chomp($line);

			next if ($line =~ /^Sequence/);

			my @line = split(/\s+/, $line);

			if ($line[0] =~ /q_(.*)_t_(.*)/) {
				if (exists($contigs{$1}) && exists($contigs{$2})) {
					my $ks = $line[3];
					$ks = 0 if ($ks eq "NA");
					
					if ($out_file_name =~ /support/) {
						print {$topo_output} "$ks,1\n";
					}
					else {
						print {$topo_output} "$ks,0\n";
					}
					last;
				}
			}

		}
		close($kaks_output);
		flock($topo_output, LOCK_UN) || die "Could not unlock '$topo_output': $!.\n";
		close($topo_output);

#		open(my $atx_file, ">", "$id-$name.atx");
#		print {$atx_file} ">$name\n";
#		#foreach my $sequence (keys %sequences) {
#		foreach my $sequence (keys %family_aligned_nuc) {
#			if ($sequence =~ /^\Q$name\E_\d+$/) {
#				print {$atx_file} "$family_aligned_nuc{$sequence}\n";
#			}
#		}
#		print {$atx_file} "\n";
#		close($atx_file);
#
#		# use KaKs_Calculator to determine Ks between homologus pairs
#		my $return = system("$kaks_calc -i $id-$name.atx -o $id-$name.kaks -m $model >/dev/null");
#		#my $return = system("$kaks_calc -i $id-$name.atx -o $id-$name.kaks -m $model");
#		#die "Error running KaKs_Calculator: '$return'.\n" if ($return);
#
#		open(my $topo_output, ">>", "$name.csv");
#		flock($topo_output, LOCK_EX) || die "Could not lock '$topo_output': $!.\n";
#		seek($topo_output, 0, SEEK_END) || die "Could not seek '$topo_output': $!.\n";
#
#		#open(my $topo_output, ">", "$id-$name.csv");
#		open(my $kaks_output, "<", "$id-$name.kaks");
#		while (my $line = <$kaks_output>) {
#			chomp($line);
#
#			next if ($line =~ /^Sequence/);
#
#			my @line = split(/\s+/, $line);
#			my $ks = $line[3];
#			$ks = 0 if ($ks eq "NA");
#
#			if ($out_file_name =~ /support/) {
#				print {$topo_output} "$ks,1\n";
#			}
#			else {
#				print {$topo_output} "$ks,0\n";
#			}
#		}
#		close($kaks_output);
#
#		flock($topo_output, LOCK_UN) || die "Could not unlock '$topo_output': $!.\n";
#		close($topo_output);
	}

	open(my $out_file, ">>", $out_file_name);
	flock($out_file, LOCK_EX) || die "Could not lock '$out_file': $!.\n";
	seek($out_file, 0, SEEK_END) || die "Could not seek '$out_file': $!.\n";

	print {$out_file} "$id: $tree\n";

	flock($out_file, LOCK_UN) || die "Could not unlock '$out_file': $!.\n";
	close($out_file);

	#unlink(@unlink);
}

sub get_pid {
	my ($seq1, $seq2) = @_;

	my $same = 0;
	my $total = 0;
	for (my $i = length($seq1); $i > 0; $i--) {
		my $char1 = chop($seq1);
		my $char2 = chop($seq2);

		if ($char1 ne "-" && $char2 ne "-") {
			if ($char1 eq $char2) {
				$same++;
			}
			$total++;
		}
	}

	return ($same / $total) * 100;
}

sub reorient_contigs { 
	my ($id, $sequences) = (@_);

	my $family_raw = "$output_dir/$id-raw.fasta";

	my %sequences = %{$sequences};

	system("$makeblastdb -in $family_raw -input_type fasta -dbtype nucl >/dev/null") || die;
	#system("$blastn -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstrand' -evalue 1e-05 >/dev/null") || die;
	#system("$blastn -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstrand' -evalue 1e-02 >/dev/null") || die;

	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

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
			$strands{$match}++;
		}
		$hits{"$query-$match"}++;
	}
	foreach my $taxon (keys %quartet) {
		my $count = $quartet{$taxon};
		if ($count < scalar(@transcriptomes) - 1) {
			#unlink(@unlink);
			die "[$id] Some taxa did not share homology.\n";
			exit(0);
		}
	}

	if (keys %strands > 0) {
		print "[$id] Reorienting contig(s).\n";
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

sub get_ks_range {
	my ($arg, $ks_range) = @_;

	my @ks_ranges = split(",", $ks_range);
	foreach my $index (0 .. $#ks_ranges) {
		my $ks_range = $ks_ranges[$index];
		my @ks_range = split("-", $ks_range);
		@ks_range = sort { $a <=> $b } @ks_range;

		$ks{"$index-LOWER"} = $ks_range[0];
		$ks{"$index-UPPER"} = $ks_range[1];

		die "Lower Ks bound '$ks_range[0]' is not a number.\n" if ($ks_range[0] !~ /^(\d*\.)?\d+$/);
		die "Upper Ks bound '$ks_range[1]' is not a number.\n" if ($ks_range[1] !~ /^(\d*\.)?\d+$/);
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
			#$taxon =~ s/-/_/g;
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

sub get_free_cpus {

	my $os_name = $^O;

	# Returns a two-member array containing CPU usage observed by top,
	# top is run twice as its first output is usually inaccurate
	my @percent_free_cpu;
	if ($os_name eq "darwin") {
		# Mac OS
		chomp(@percent_free_cpu = `top -i 1 -l 2 | grep "CPU usage"`);
	}
	else {
		# Linux
		chomp(@percent_free_cpu = `top -b -n2 -d0.05 | grep "Cpu(s)"`);
	}

	my $percent_free_cpu = pop(@percent_free_cpu);

	if ($os_name eq "darwin") {
		# Mac OS
		$percent_free_cpu =~ s/.*?(\d+\.\d+)%\s+id.*/$1/;
	}
	else {
		# linux 
		$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
		$percent_free_cpu = eval($percent_free_cpu);
	}

	my $total_cpus;
	if ($os_name eq "darwin") {
		# Mac OS
		$total_cpus = `sysctl -n hw.ncpu`;
	}
	else {
		# linux
		$total_cpus = `grep --count 'cpu' /proc/stat` - 1;
	}

	my $free_cpus = int($total_cpus * $percent_free_cpu / 100) + 1;

	if ($free_cpus == 0 || $free_cpus !~ /^\d+$/) {
		$free_cpus = 1; # assume that at least one cpu can be used
	}
	
	return $free_cpus;
}

sub help {
	return "";
}

sub usage {
	return "";
}

