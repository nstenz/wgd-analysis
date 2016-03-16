#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use POSIX qw(:sys_wait_h);
use Time::HiRes qw(usleep);
use Fcntl qw(:flock SEEK_END);
use File::Path qw(remove_tree);

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
my $pid_cutoff = 80; # minimum amino acid percent identity in alignments

# Where we will output results
my $output_dir = "wgd-test-".int(time());

# check that dependencies are present in user's PATH
my $blat = check_path_for_exec("blat");
my $blastp = check_path_for_exec("blastp");
my $muscle = check_path_for_exec("muscle");
my $raxml = check_path_for_exec("raxmlHPC");
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
my $threshold = 300; # TODO: more specific name
my $min_contig_length = 300;
my @transcriptomes;
GetOptions(
	"help|h"                  => \&help,
	"ks|k=s"                  => \&get_ks_range,
	"model|m=s"               => \$model,
	"output|o=s"              => \$output_dir,
	"pfam-cpus|p=i"           => \$pfam_cpus,
	"pfam-search|s=s"         => \$pfam_search,
	"transcriptomes|t=s{2,2}" => \@transcriptomes,
);

# Check that hmmscan can be located if user wants pfam search
my $hmmscan = check_path_for_exec("hmmscan") if ($pfam_search);

# Error check input
die "Unknown arguments specified.\n".&help if (@ARGV);
die "You need to specify two fasta files containing transcriptomes.\n".&usage if (scalar(@transcriptomes) != 2);
die "The model '$model' does not exist or is not usable by this script.\n".&usage if (!exists($models{$model}));
die "Could not locate hmmscan binary file: '$pfam_search'.\n" if (defined($pfam_search) && !-e $pfam_search);

# Check that specified input exists
foreach my $transcriptome (@transcriptomes) {
	die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);

	$transcriptome = abs_path($transcriptome);	
}

# Create the output directory
mkdir($output_dir) if (!-e $output_dir);

# Move to output directory, remove result files
chdir($output_dir);
unlink("against.txt", "support.txt");

# Translate transcriptomes into most likely proteome
my @pids;
my @gene_pairs;
$SIG{INT} = sub { foreach my $pid (@pids) { kill -9, $pid }; exit(0); };
foreach my $index (0 .. $#transcriptomes) {
	my $transcriptome = $transcriptomes[$index];

	# Fetch upper and lower Ks user specified for this file
	$lower_ks = $ks{"$index-LOWER"};
	$upper_ks = $ks{"$index-UPPER"};

	# Wait until a CPU is free
	until(okay_to_run(\@pids)) {};

	my $pid = fork();
	if ($pid == 0) {
		setpgrp();
		extract_protein_pairs($transcriptome, $lower_ks, $upper_ks);
		exit(0);
	}
	else {
		# Store output filename
		(my $transcriptome_name = $transcriptome) =~ s/(.*\/)?(.*)/$2/;
		#my $output_file = $transcriptome_name.".ks_$lower_ks-$upper_ks.fasta";
		my $output_file = $transcriptome_name.".ks_$lower_ks-$upper_ks";
		push(@gene_pairs, $output_file);

		# Store child pid
		push(@pids, $pid);
	}
}

# Wait until everything is complete
foreach my $pid (@pids) {
	waitpid($pid, 0);
}

# Run proteinortho using subset of genes of base transcriptome against all other transcriptomes 
my $return = system("$protein_ortho @gene_pairs -clean -project=wgd-test") if (!-e "wgd-test.proteinortho");
die "Error running ProteinOrtho: '$return'.\n" if ($return);

# Parse proteinortho output and create a new multisequence alignment for each hit with proteins common to all transcriptomes
my %families;
my $count = 0;
my @species_order;
open(my $ortho_output, "<", "wgd-test.proteinortho");
while (my $line = <$ortho_output>) {
	chomp($line);

	my @line = split("\t", $line);
	my ($num_species, $num_genes) = ($line[0], $line[1]);

	# Use first line to reorder transcriptomes
	if ($. == 1) {
		@species_order = @line[3 .. $#line]; 
		next;
	}

	# Only want matches which include all species
	if ($num_species == scalar(@transcriptomes) && $num_genes == scalar(@transcriptomes)) {

		# Extract clustered contig names
		my @contigs = @line[3 .. $#line];	

		# Separate contig names 
		foreach my $index (0 .. $#contigs) {

			# Parse out the two separate sequences present in the name
			my $contigs = $contigs[$index];
			if ($contigs =~ /q_(.*)_t_(.*)/) {
				my ($contig1, $contig2) = ($1, $2);

				$contigs[$index] = [$contig1, $contig2];
			}
		}
		$families{$count} = \@contigs;
		$count++;

		print "$count: ";
		foreach my $contig (@contigs) {
			print "@{$contig} ";
		}
		print "\n";
	}
}
close($ortho_output);

print "$count quartets found.\n";

# Reorient transcriptome array to match protein ortho output
foreach my $index (0 .. $#species_order) {

	# Modify to TransDecoder output
	print $transcriptomes[$index],"\n";
	print $species_order[$index],"\n";
	#(my $file = $species_order[$index]) =~ s/\.ks_.*?-.*?\.fasta$/.fasta.transdecoder.pep/;
	(my $file = $species_order[$index]) =~ s/\.fasta\.ks_.*?-.*?$/.fasta.transdecoder.pep/;
	print $file,"\n";
	$transcriptomes[$index] = $file;
}

# Extract alignable sequence for each quartet, align, then run RAxML
undef(@pids);
#$SIG{CHLD} = 'IGNORE';
foreach my $family (sort { $a <=> $b } keys %families) {
	my $members = $families{$family};

	# Wait until a CPU is free
	until(okay_to_run(\@pids)) {};

	my $pid = fork();
	if ($pid == 0) {
		setpgrp();
		analyze_family({'ID' => $family, 'MEMBERS' => $members});
		exit(0);
	}
	else {
		push(@pids, $pid);
	}
}

# Wait until each job is done
foreach my $pid (@pids) {
	waitpid($pid, 0);
}

sub extract_protein_pairs {
	my ($transcriptome, $lower_ks, $upper_ks) = @_;

	(my $base_name = $transcriptome) =~ s/(.*\/)?(.*)\..*/$2/;
	(my $transcriptome_name = $transcriptome) =~ s/(.*\/)?(.*)/$2/;
	#my $output_file = $transcriptome_name.".ks_$lower_ks-$upper_ks.fasta";
	my $output_file = $transcriptome_name.".ks_$lower_ks-$upper_ks";

	# TransDecoder output filename
	my $proteome = $transcriptome_name.".transdecoder.pep";

	# Check if TransDecoder output already exists, run if it doesn't
	if (!-e $proteome) {
		print "Running TransDecoder on '$transcriptome'...\n";
		#run_transdecoder($transcriptome);
		print "Finished TransDecoder for '$transcriptome.\n";
	}
	else {
		print "Using TransDecoder output in file '$proteome'...\n";
	}

	# Parse TransDecoder CDS
	my %cds = parse_fasta("$transcriptome_name.transdecoder.cds");

	# Perform self-blat at protein level
	my $return = system("$blat '$proteome' '$proteome' -prot -out=pslx -noHead $proteome.pslx") if (!-e "$proteome.pslx");
	die "Error running self-blat: '$return'.\n" if ($return);

	# Parse blat output 
	my (%queries, %matches);
	open(my $blat_output, "<", "$proteome.pslx");
	while (my $line = <$blat_output>) {
		chomp($line);

		# Split the line on tabs and extract the information we want
		my @line = split("\t", $line);
		my ($match, $mismatch, $query_name, $match_name, $query_align, $match_align) = 
			($line[0], $line[1], $line[9], $line[13], $line[21], $line[22]);

		# Exclude self matches
		next if ("$query_name" eq "$match_name");

		# Calculate full length of alignment
		my $match_length = $match + $mismatch;

		# Skip hit if it is too short (3 * match_length = length nucleotides);
		next if ($match_length * 3 < $threshold);

		# Extract query and match protein sequences
		my @query_align =  split(",", $query_align);
		my @match_align =  split(",", $match_align);
		my $query_align_prot = join("", @query_align);

		# Extract query and match nucleotides sequences
		my $query_nuc_align = $cds{$query_name};
		my $match_nuc_align = $cds{$match_name};

		# Reverse translate protein matches to corresponding nucleotide sequences
		my $trans_query_align;
		foreach my $align (@query_align) {
			$trans_query_align .= reverse_translate({"DNA" => $query_nuc_align, "PROT" => $align});
		}

		my $trans_match_align;
		foreach my $align (@match_align) {
			$trans_match_align .= reverse_translate({"DNA" => $match_nuc_align, "PROT" => $align});
		}

		# Set aligned sequences to their corresponding reverse translated sequence
		$query_align = $trans_query_align;
		$match_align = $trans_match_align;

		# Check that the match hasn't already been extracted in reverse order
		if (!exists($queries{"$match_name-$query_name"})) {

			my $name = "q_$query_name"."_t_$match_name";
			my $pair = {'QUERY_ALIGN' => $query_align,
						'QUERY_ALIGN_PROT' => $query_align_prot,
						'MATCH_ALIGN' => $match_align,
						'LENGTH' => length($query_align)};

			# Check if there is already a match between these two sequences
			# if there is a match, the longer one will be output
			my $current_hit = $matches{$name};
			if (defined($current_hit)) {
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
	}
	close($blat_output);

	# Calculate how many hits we found
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

	# Use KaKs_Calculator to determine Ks between homologus pairs
	$return = system("$kaks_calc -i '$base_name.atx' -o '$base_name.kaks' -m $model >/dev/null") if (!-e "$base_name.kaks");
	die "Error running KaKs_Calculator: '$return'.\n" if ($return);

	# Select gene pairs with Ks between $lower_ks and $upper_ks
	my %genes;

	open(my $kaks_output, "<", $base_name.".kaks");
	open(my $gene_subset, ">", $output_file);

	# Parse KaKs_Calculator output
	while (my $line = <$kaks_output>) {
		chomp($line);

		# Skip header line
		next if ($line =~ /^Sequence/);

		# Extract Ks for each pair
		my @line = split(/\s+/, $line);
		my $ks = $line[3];
		$ks = 0 if ($ks eq "NA");

		# Check that Ks is within specified range
		if ($ks >= $lower_ks && $ks <= $upper_ks) {
			(my $pair_name = $line[0]) =~ s/^>//;

			if ($pair_name =~ /q_.*_t_.*/) {
				my $pair = $matches{$pair_name};
				print {$gene_subset} ">$pair_name\n";
				print {$gene_subset} $pair->{QUERY_ALIGN_PROT},"\n";
			}
		}
	}
	close($kaks_output);
	close($gene_subset);

	return;
}

sub analyze_family {
	my $args = shift;

	# Extract method arguments
	my $id = $args->{ID};
	my @members = @{$args->{MEMBERS}};

	# Extract each members' alignment data
	my %names;
	my %sequences;
	my %nuc_sequences;
	my $total_members = 0;
	foreach my $index (0 .. $#members) {
		my $member = $members[$index];

		(my $transcriptome_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)/$2/;
		(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\.fa(sta)?\.transdecoder\.pep/$2/;

		(my $cds_file = $transcriptomes[$index]) =~ s/\.pep$/.cds/;
		foreach my $index2 (0 .. scalar(@{$member}) - 1) {
			my $contig = @{$member}[$index2];
			my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $contig});
			my $nuc_sequence = get_seq_from_fasta({'FILE' => $cds_file, 'TAXON' => $contig});

			# Store extracted sequences
			$sequences{$base_name."_$index2"} = $sequence;
			$nuc_sequences{$base_name."_$index2"} = $nuc_sequence;
			$names{$base_name}++;
			$total_members++;
		}
	}

#	# Set output filenames
#	my $family_raw = "$output_dir/$id-raw.fasta";
#	my $family_aligned = "$output_dir/$id-aligned.nex";
#	my $family_reduced = "$output_dir/$id-reduced.fasta";
	my $family_raw = "$id-raw.fasta";
	my $family_aligned = "$id-aligned.nex";
	my $family_reduced = "$id-reduced.fasta";

	# For cleaning up temp files in case of an interruption
	my @unlink;
	push(@unlink, $family_raw, $family_reduced, $family_raw.".phr", $family_raw.".pin", $family_raw.".psq", $family_raw.".out");
	$SIG{INT} = sub { unlink(@unlink); exit(0) };

	# Write the raw sequences to file
	open(my $raw, ">", $family_raw);
	foreach my $taxon (sort { $a cmp $b } keys %sequences) {
		print {$raw} ">$taxon\n";
		print {$raw} "$sequences{$taxon}\n";
	}
	close($raw);

	# Reduce family to homologous sites
	print "[$id] Reducing sequence to homologous sites.\n";

	# Run makeblastdb on the raw sequences
	my $return = system("$makeblastdb -in '$family_raw' -input_type fasta -dbtype prot >/dev/null");
	die "Error creating blastdb: '$return'.\n" if ($return);

	# Run blastp on the raw sequences
	$return = system("$blastp -db '$family_raw' -query '$family_raw' -out='$family_raw.out' -outfmt '6 qseqid sseqid qstart qend sstart send' -evalue 1e-05 >/dev/null");
	die "Error running blastp: '$return'.\n" if ($return);

	# Read in blast output
	open(my $blast, "<", "$family_raw.out");
	chomp(my @lines = <$blast>);
	close($blast);

	# Parse blast output
	my %counts;
	my %matches;
	my %quartet;
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];

		# Split line on tabs
		my @line = split("\t", $line);
		my ($query, $match, $q_start, $q_end, $s_start, $s_end) = ($line[0], $line[1], $line[2], $line[3], $line[4], $line[5]);

		# Skip self-matches
		next if ($query eq $match);

		$counts{$query}++;

		# For checking if the next line is from the same species
		my $next_line = " \t \t";
		if ($index + 1 <= $#lines) {
			$next_line = $lines[$index + 1];
		}
		my @next_line = split("\t", $next_line);
		my ($next_query, $next_match) = ($next_line[0], $next_line[1]);

		# Convert from Blast's 1-based indexing
		$q_start--; $q_end--;
		$s_start--; $s_end--;

		# Ranges of sites included in alignment
		my $q_range = "$q_start-$q_end";
		my $s_range = "$s_start-$s_end";

		# Check if a hit already exists between the two species
		if (!exists($matches{"$query-$match"})) {
			$matches{"$query-$match"} = {$q_range => $s_range};
		}
		else {

			# Since a hit already exists, combine it with the newer one
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

	# Check that we had a hit from each member of the family
	if (scalar(keys %quartet) < 4) {
		unlink(@unlink);
		die "[$id] Could not identify homologous sites for a species in this family.\n";
	}

	# Check that each member shares sequence with other members
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

		# Sequences are currently amino acids so divide length by 3
		if ($reduced_length < $min_contig_length / 3) {
			unlink(@unlink);
			die "[$id] A sequence did not meet the minimum length requirements ($reduced_length).\n";
		}
		print {$out} "\n";
	}
	close($out);

	# Align with Muscle
	print "[$id] Aligning reduced sites with muscle.\n";
	$return = system("$muscle -in '$family_reduced' -out '$family_aligned' >/dev/null 2>&1");
	die "Error running muscle: '$return'.\n" if ($return);

	# Parse Muscle alignment
	my %family_aligned = parse_fasta($family_aligned);

	# Check pairwise identity between all members
	foreach my $index1 (0 .. scalar(keys %family_aligned) - 2) {
		my $taxon1 = (keys %family_aligned)[$index1];

		foreach my $index2 ($index1 + 1 .. scalar(keys %family_aligned) - 1) {
			my $taxon2 = (keys %family_aligned)[$index2];
			my $pid = get_pid($family_aligned{$taxon1}, $family_aligned{$taxon2});

			#die "[$id] Pairwise similarity between sequences was too low ($pid).\n" if ($pid < $pid_cutoff);
			
			# Stop analysis if pid between a pair of sequences is too
			if ($pid < $pid_cutoff) {
				die "[$id] Pairwise similarity between sequences was too low ($pid).\n";
			}
		}
	}

	# Reverse translate from amino acid to nucleotide
	my %family_aligned_nuc;
	foreach my $taxon (keys %family_aligned) {
		$family_aligned_nuc{$taxon} = reverse_translate({"DNA" => $nuc_sequences{$taxon}, "PROT" => $family_aligned{$taxon}});
	}

	# Write alignment to file
	write_phylip({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned_nuc});

	# Run RAxML
	#chdir($output_dir);
	#$return = system("$raxml -f a -m GTRGAMMA -p 123123 -x 123123 -#100 -s '../$family_aligned' -n $id >/dev/null");
	$return = system("$raxml -f a -m GTRGAMMA -p 123123 -x 123123 -#100 -s '$family_aligned' -n $id >/dev/null");
	die "Error running raxml: '$return'.\n" if ($return);

	# Parse RAxML output
	open(my $raxml_out, "<", "RAxML_bipartitions.$id");
	chomp(my $tree = <$raxml_out>);
	close($raxml_out);

	# Parse taxon names
	my @names = keys %names;
	my $taxon_1 = $names[0];
	my $taxon_2 = $names[1];

	# Remove branch lengths
	(my $tree_no_bls = $tree) =~ s/:\d+(\.\d+)?//g;

	# Determine whether topology is against or supports WGD
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
	else {
		print "REGEX ERROR\n";
	}

	# Add family to corresponding output file
	open(my $out_file, ">>", $out_file_name);
	flock($out_file, LOCK_EX) || die "Could not lock '$out_file': $!.\n";
	seek($out_file, 0, SEEK_END) || die "Could not seek '$out_file': $!.\n";

	print {$out_file} "$id: $tree\n";

	flock($out_file, LOCK_UN) || die "Could not unlock '$out_file': $!.\n";
	close($out_file);

	unlink(@unlink);
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

	# Split on commas
	my @ks_ranges = split(",", $ks_range);
	foreach my $index (0 .. $#ks_ranges) {
		my $ks_range = $ks_ranges[$index];
		my @ks_range = split("-", $ks_range);
		@ks_range = sort { $a <=> $b } @ks_range;

		$ks{"$index-LOWER"} = $ks_range[0];
		$ks{"$index-UPPER"} = $ks_range[1];

		die "Error parsing Ks from '$ks_range'.\n" if (!defined($ks{"$index-LOWER"}) || !defined($ks{"$index-UPPER"}));

		# Check that parsed Ks is a number
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
			die "Invalid character in taxon name: '$taxon'.\n" if (!$taxon =~ /-/);
		}
		else {
			$align{$taxon} .= $line;
		}
	}
	return %align;
}

sub check_path_for_exec {
	my ($exec, $continue) = @_;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec && -x $dir.$exec && !-d $dir.$exec);
	}

	if (!defined($exec_path) && $continue) {
		return;
	}
	elsif (defined($exec_path)) {
		return $exec_path;
	}
	else {
		die "Could not find the following executable: '$exec'. This script requires this program in your path.\n";
	}
	return $exec_path;
}

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
		$translation = substr($translation, 0, $gap)."---".substr($translation, $gap, length($translation) - $gap);
		$count++;
	}

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
		$union .= $starts[$index]."-".$ends[$index].",";
	}
	chop($union);

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

	my $seq;
	my $is_desired_seq;

	open(my $align, "<", $file) or die "Could not open '$file': $!.\n";
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
	my $pids = shift;

	# Free up a CPU by sleeping for 10 ms
	usleep(10000);

	my $current_forks = scalar(@{$pids});
	foreach my $index (reverse(0 .. $current_forks - 1)) {
		next if ($index < 0);

		my $pid = @{$pids}[$index];
		my $wait = waitpid($pid, WNOHANG);

		# Successfully reaped child
		if ($wait > 0) {
			$current_forks--;
			#splice(@pids, $index, 1);
			splice(@{$pids}, $index, 1);
		}
	}

	return ($current_forks < $free_cpus);
}

sub run_transdecoder {
	my $transcriptome = shift;

	# TransDecoder working directory name
	my $transdecoder_out_dir = "transdecoder-$transcriptome";

	# Check if the files we need from TransDecoder are present
	if (!-e "$transcriptome.transdecoder.pep" || !-e "$transcriptome.transdecoder.mRNA") {
		print "\nRunning TransDecoder on '$transcriptome'...";

		# Opt for newer version of TransDecoder if present
		if ($transdecoder_orfs) {
			system("$transdecoder_orfs -t '$transcriptome'");

			# Need to run things differently to include pfam
			if (defined($pfam_search)) {
				system("$hmmscan --cpu $pfam_cpus --domtblout pfam.domtblout '$pfam_search' '$transcriptome.transdecoder_dir/longest_orfs.pep'");
				system("$transdecoder_predict -t '$transcriptome' --retain_pfam_hits pfam.domtblout");
			}
			else {
				system("$transdecoder_predict -t '$transcriptome'");
			}
			remove_tree("$transcriptome.transdecoder_dir/");
		}
		# Old version of TransDecoder
		else {
			if (defined($pfam_search)) {
				system("$transdecoder -t '$transcriptome' --workdir '$transdecoder_out_dir' --search_pfam $pfam_search --CPU $pfam_cpus");
			}
			else {
				system("$transdecoder -t '$transcriptome' --workdir '$transdecoder_out_dir'");
			}
			remove_tree($transdecoder_out_dir);
		}

		# Clean up unneeded files
		remove_tree($transdecoder_out_dir);

		print "Completed TransDecoder.\n";
	}
	else {
		print "\nUsing TransDecoder output from previous run.";
	}
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

