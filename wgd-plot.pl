#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use Time::HiRes qw(usleep);
use Fcntl qw(:flock SEEK_END);
use POSIX qw(:sys_wait_h ceil);
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
my $pid_cutoff = 60; # amino acid level
my $boot_cutoff = 70; # bootstrap support for RAxML

# Where we will output results
my $kaks_calc_in_dir = "kaks-in";
my $kaks_calc_out_dir = "kaks-out";
my $output_dir = "wgd-plot-".int(time());

# Check that dependencies are present in user's PATH
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
	die "Could not locate required TransDecoder executables (TransDecoder or TransDecoder.LongOrfs and TransDecoder.Predict) in \$PATH.\n" 
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
#my $threshold = 300;
my $min_contig_length = 300;
my @transcriptomes;
GetOptions(
	"help|h"                  => \&help,
	"model|m=s"               => \$model,
	"output|o=s"              => \$output_dir,
	"n-threads|T=i"           => \$free_cpus,
	"pid-cut|c=i"             => \$pid_cutoff,
	"boot-cut|b=i"            => \$boot_cutoff,
	"pfam-cpus|p=i"           => \$pfam_cpus,
	"pfam-search|s=s"         => \$pfam_search,
	"min-length|l=i"          => \$min_contig_length,
	"transcriptomes|t:s{2,2}" => \@transcriptomes,
);

# Check that hmmscan can be located if user wants pfam search
my $hmmscan = check_path_for_exec("hmmscan") if ($pfam_search);

# Error check input
die "Unknown arguments specified.\n".&help if (@ARGV);
die "You need to specify two fasta files containing transcriptomes.\n".&usage if (scalar(@transcriptomes) != 2);

# Check that specified input exists
foreach my $transcriptome (@transcriptomes) {
	die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);

	$transcriptome = abs_path($transcriptome);	
}

# Check that specified input exists
foreach my $transcriptome (@transcriptomes) {
	die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);

	$transcriptome = abs_path($transcriptome);	
}

# Create the output directories
mkdir($output_dir) if (!-e $output_dir);
mkdir("$output_dir/genes") if (!-e "$output_dir/genes");

# Move to output directory, remove result files
chdir($output_dir);

# Translate each transcriptome into best frame
my @pids;
my @gene_pairs;
$SIG{INT} = sub { foreach my $pid (@pids) { kill -9, $pid }; exit(0); };
foreach my $index (0 .. $#transcriptomes) {
	my $transcriptome = $transcriptomes[$index];

	# Wait until a CPU is free
	until(okay_to_run(\@pids)) {};

	# Perform fork
	my $pid = fork();
	if ($pid == 0) {
		setpgrp();
		extract_protein_pairs($transcriptome);
		exit(0);
	}
	else {
		# Store output filename
		(my $transcriptome_name = $transcriptome) =~ s/(.*\/)?(.*)/$2/;
		my $output_file = $transcriptome_name.".$model.pairs.fasta";
		push(@gene_pairs, $output_file);

		# Store child pid
		push(@pids, $pid);
	}

#	my $pid = fork();
#	if ($pid == 0) {
#
#		my $proteome = $transcriptome.".transdecoder.pep";
#
#		if (!-e $proteome) {
#			print "Running TransDecoder on '$transcriptome'...\n";
#			my $return = system("$transdecoder -t $transcriptome >/dev/null");
#			die "Error occurred while running TransDecoder: '$return'.\n" if ($return);
#			print "Finished TransDecoder for '$transcriptome.\n";
#		}
#		else {
#			print "Using TransDecoder output in file '$proteome'...\n";
#		}
#
#		my %align = parse_fasta("$transcriptome.transdecoder.cds");
#
#		# Perform self-blat
#		my $return = system("$blat $proteome $proteome -prot -out=pslx -noHead $proteome.pslx") if (!-e "$proteome.pslx");
#		die "Error running self-blat: '$return'.\n" if ($return);
#
#		# Parse blat output 
#		my (%queries, %matches);
#		open(my $blat_output, "<", "$proteome.pslx");
#		while (my $line = <$blat_output>) {
#			chomp($line);
#
#			# split the line on tabs and extract the information we want
#			my @line = split("\t", $line);
#			#my ($query_name, $match_name, $query_align, $match_align) = 
#			#	($line[9], $line[13], $line[21], $line[22]);
#			my ($match, $mismatch, $query_name, $match_name, $query_align, $match_align) = 
#				($line[0], $line[1], $line[9], $line[13], $line[21], $line[22]);
#
#			next if ("$query_name" eq "$match_name");
#
#			# exclude self match
#			my $match_length = $match + $mismatch;
#			#if ($query_name ne $match_name) {
#
#				next if ($match_length * 3 < $threshold);
#
#				my @query_align =  split(",", $query_align);
#				my @match_align =  split(",", $match_align);
#				my $query_align_prot = join("", @query_align);
#
#				my $query_nuc_align = $align{$query_name};
#				my $match_nuc_align = $align{$match_name};
#
#				# Reverse translate protein matches to corresponding nucleotide sequences
#				my $trans_query_align;
#				foreach my $align (@query_align) {
#					$trans_query_align .= reverse_translate({"DNA" => $query_nuc_align, "PROT" => $align});
#				}
#
#				my $trans_match_align;
#				foreach my $align (@match_align) {
#					$trans_match_align .= reverse_translate({"DNA" => $match_nuc_align, "PROT" => $align});
#				}
#
#				$query_align = $trans_query_align;
#				$match_align = $trans_match_align;
#
#				my @names;
#				push(@names, $query_name);
#				push(@names, $match_name);
#				@names = sort { "$a" cmp "$b" } @names;
#
#				#if (length($query_align) >= $threshold && length($match_align) >= $threshold) {
#					#&& length($query_align) / $query_len >= $hit_frac && length($match_align) / $match_len >= $hit_frac) {
#
#				#	# Check that the match hasn't already been extracted in reverse order
#					if (!exists($queries{"$match_name-$query_name"})) {
#
#						my $name = "q_$query_name"."_t_$match_name";
#						#my $name = "q_$names[0]"."_t_$names[1]";
#						my $pair = {'QUERY_ALIGN' => $query_align,
#									'QUERY_ALIGN_PROT' => $query_align_prot,
#									'MATCH_ALIGN' => $match_align,
#									'LENGTH' => length($query_align)};
#
#						# Check if there is already a match between these two sequences
#						# if there is a match, the longer one will be output
#						my $current_hit = $matches{$name};
#						#if (exists($matches{$name})) {
#						if (defined($current_hit)) {
#							#my $current_length = $matches{$name}->{'LENGTH'};
#							my $current_length = $current_hit->{'LENGTH'};
#							if ($current_length <= length($query_align)) {
#								$matches{$name} = $pair;
#							}
#						}
#						else {
#							$matches{$name} = $pair;
#						}
#						$queries{"$query_name-$match_name"}++;
#					}
#				#}
#			#}
#		}
#		close($blat_output);
#
#		my $num_matches = scalar(keys %matches);
#		die "No homologs meeting the nucleotide length threshold found.\n" if (!$num_matches);
#		print "$num_matches pair(s) met the requirements.\n";
#
#		# Write .atx input file for KaKs_Calculator
#		open(my $blat_atx, ">", "$base_name.atx");
#		foreach my $key (sort { $a cmp $b} keys %matches) {
#			my $pair = $matches{$key};
#
#			(my $query_name = $key) =~ s/q_(.*)_t_.*/$1/;
#			(my $match_name = $key) =~ s/q_.*_t_(.*)/$1/;
#
#			print {$blat_atx} ">$key\n";
#			print {$blat_atx} "$pair->{QUERY_ALIGN}\n";
#			print {$blat_atx} "$pair->{MATCH_ALIGN}\n\n";
#		}
#		close($blat_atx);
#
#		# use KaKs_Calculator to determine Ks between homologus pairs
#		$return = system("$kaks_calc -i $base_name.atx -o $base_name.kaks -m $model >/dev/null") if (!-e "$base_name.kaks");
#		die "Error running KaKs_Calculator: '$return'.\n" if ($return);
#
#		# Select gene pairs with Ks between $lower_ks and $upper_ks
#		my %genes;
#
#		open(my $kaks_output, "<", $base_name.".kaks");
#		#open(my $gene_subset, ">", $base_name.".ks_$lower_ks-$upper_ks.fasta");
#		open(my $gene_subset, ">", $base_name.".pairs.fasta");
#
#		while (my $line = <$kaks_output>) {
#			chomp($line);
#
#			next if ($line =~ /^Sequence/);
#
#			my @line = split(/\s+/, $line);
#			my $ks = $line[3];
#			$ks = 0 if ($ks eq "NA");
#
#	#		if ($ks >= $lower_ks && $ks <= $upper_ks) {
#				(my $pair_name = $line[0]) =~ s/^>//;
#
#				if ($pair_name =~ /q_(.*)_t_(.*)/) {
#					my ($query_name, $match_name) = ($1, $2);
##					print $query_name, " ", $match_name,"\n";
#
#					#my $found_match = 0;
##					foreach my $key (keys %matches) {
##						if ($key eq $pair_name) {
##							my $pair = $matches{$key};
##							print {$gene_subset} ">$pair_name\n";
##							print {$gene_subset} $pair->{QUERY_ALIGN_PROT},"\n";
##							$found_match++;
##							last;
##						}
##					}
#					my $pair = $matches{$pair_name};
#					print {$gene_subset} ">$pair_name\n";
#					print {$gene_subset} $pair->{QUERY_ALIGN_PROT},"\n";
#					#die "Could not locate match: '$pair_name'\n" if (!$found_match);
#				}
#		}
#		close($kaks_output);
#		close($gene_subset);
#
#		exit(0);
#	}
#	else {
#		push(@pids, $pid);
#	}
}

# Wait until everything is complete
foreach my $pid (@pids) {
	waitpid($pid, 0);
}

# Run proteinortho using subset of genes of base transcriptome against all other transcriptomes 
my $return = system("$protein_ortho @gene_pairs -clean -project=wgd-plot") if (!-e "wgd-plot.proteinortho");
die "Error running ProteinOrtho: '$return'.\n" if ($return);

# Parse proteinortho output and create a new multisequence alignment for each hit with proteins common to all transcriptomes
my %families;
my $count = 0;
my @species_order;
open(my $ortho_output, "<", "wgd-plot.proteinortho");
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

		# Output info about group
		print "$count: ";
		foreach my $contig (@contigs) {
			print "@{$contig} ";
		}
		print "\n";
	}
}
close($ortho_output);

# Reorient transcriptome array to match protein ortho output
foreach my $index (0 .. $#species_order) {

	# Modify to TransDecoder output
	(my $file = $species_order[$index]) =~ s/\.fasta\.$model\.pairs\.fasta$/.fasta.transdecoder.pep/;
	$transcriptomes[$index] = $file;
}

mkdir($output_dir) if (!-e $output_dir);
#unlink("against.txt", "support.txt");

# Analyze families concurrently
my @unlink;
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
	my $transcriptome = shift;

	# Extract name of transcriptome file
	(my $base_name = $transcriptome) =~ s/(.*\/)?(.*)\..*/$2/;
	(my $transcriptome_name = $transcriptome) =~ s/(.*\/)?(.*)/$2/;
	my $output_file = $transcriptome_name.".$model.pairs.fasta";

	# TransDecoder output filenames
	my $cds = "$transcriptome_name.transdecoder.cds";
	my $proteome = $transcriptome_name.".transdecoder.pep";

	# Check if TransDecoder output already exists, run if it doesn't
	if (!-e $proteome) {
		print "Running TransDecoder on '$transcriptome'...\n";
		run_transdecoder($transcriptome);
		print "Finished TransDecoder for '$transcriptome.\n";
	}
	else {
		print "Using TransDecoder output in file '$proteome'...\n";
	}

	# Check whether or not we have KaKs_Calculator output for this model
	if (!-e $output_file) {

		# Run all versus all self-blat
		my %matches = run_self_blat($proteome, $cds);

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
		run_kaks_calc("$base_name.atx", "$base_name.$model.kaks");

		# Calculate how many hits we found
		my $num_matches = scalar(keys %matches);
		die "No homologs meeting the nucleotide length threshold found.\n" if (!$num_matches);
		print "$num_matches pair(s) met the requirements.\n";

		# Select gene pairs with Ks between $lower_ks and $upper_ks
		my %genes;

		# Parse and filter KaKs_Calculator output
		open(my $gene_subset, ">", $output_file);
		open(my $kaks_output, "<", "$base_name.$model.kaks");
		while (my $line = <$kaks_output>) {
			chomp($line);

			# Skip header line
			next if ($line =~ /^Sequence/);

			# Extract Ks for each pair
			my @line = split(/\s+/, $line);
			my $ks = $line[3];
			$ks = 0 if ($ks eq "NA");

			# Output representative sequence for pair
			(my $pair_name = $line[0]) =~ s/^>//;
			if ($pair_name =~ /q_.*_t_.*/) {
				my $pair = $matches{$pair_name};
				print {$gene_subset} ">$pair_name\n";
				print {$gene_subset} $pair->{QUERY_ALIGN_PROT},"\n";
			}
		}
		close($kaks_output);
		close($gene_subset);
	}

	return;
}

sub run_self_blat {
	my ($proteome, $cds) = @_;

	# Load CDS sequences
	my %cds = parse_fasta($cds);

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
		#next if ($match_length * 3 < $threshold);
		next if ($match_length * 3 < $min_contig_length);

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

	return %matches;
}

sub run_kaks_calc {
	my ($full_atx, $output) = @_;

	# Run KaKs_Calculator to get Ks values
	print "Running KaKs_Calculator...\n";

	# Determine how to split the input file
	chomp(my $total_hits = `wc -l $full_atx`);
	$total_hits =~ s/(^\s*\d+).*/$1/;
	$total_hits = $total_hits / 4;

	# Calculate how many lines each fork gets
	my $lines_per_thread = ceil($total_hits / $free_cpus) * 4;

	# Create the directories that will hold the split files and their KaKs_Calculator output
	mkdir("$kaks_calc_in_dir") || die "Could not create directory '$kaks_calc_in_dir': $!.\n";
	mkdir("$kaks_calc_out_dir") || die "Could not create directory '$kaks_calc_out_dir': $!.\n";

	# Split the files
	#run_cmd("split -l $lines_per_thread $paralogs_seqs $kaks_calc_in_dir/$paralogs_seqs.");
	my $return = system("split -l $lines_per_thread '$full_atx' '$kaks_calc_in_dir/$full_atx.'");
	die "Error splitting '$full_atx': '$return'.\n" if ($return);

	# Run each instance of KaKs_Calculator
	my @pids;
	foreach my $file (glob("$kaks_calc_in_dir/$full_atx.*")) {

		(my $file_id = $file) =~ s/.*\.(\S+$)/$1/;

		# Fork
		my $pid = fork();
		if ($pid == 0) {
			# Child 
			$return = system("$kaks_calc -i '$file' -o '$kaks_calc_out_dir/$output.$file_id' -m $model >/dev/null");

			exit(0);
		}
		else {
			# Parent
			push(@pids, $pid);
		}
	}

	# Don't let forks become zombies
	foreach my $pid (@pids) {
		waitpid($pid, 0);
	}

	# Combine all the output files generated by the multiple forks
	$return = system("cat $kaks_calc_out_dir/$output.* > $output");
	die "Error concatenating KaKs output files.\n" if ($return);

	# Clean up
	remove_tree($kaks_calc_in_dir);
	remove_tree($kaks_calc_out_dir);

	print "Completed KaKs_Calculator.\n";
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
	my %original_names;
	my $total_members = 0;
	foreach my $index (0 .. $#members) {
		my $member = $members[$index];

		# Get name info for this member
		(my $base_name = $transcriptomes[$index]) =~ s/(.*\/)?(.*)\.fa(sta)?\.transdecoder\.pep/$2/;
		(my $cds_file = $transcriptomes[$index]) =~ s/\.pep$/.cds/;

		foreach my $index2 (0 .. scalar(@{$member}) - 1) {

			# Get sequence for this member
			my $contig = @{$member}[$index2];
			my $sequence = get_seq_from_fasta({'FILE' => $transcriptomes[$index], 'TAXON' => $contig});
			my $nuc_sequence = get_seq_from_fasta({'FILE' => $cds_file, 'TAXON' => $contig});

			# Store extracted sequences
			$sequences{$base_name."_$index2"} = $sequence;
			$nuc_sequences{$base_name."_$index2"} = $nuc_sequence;
			$names{$base_name}++;
			$original_names{$base_name."_$index2"} = $contig;
			$total_members++;
		}
	}

	# Move into gene directory
	chdir("genes");

	# Set output filenames
	my $family_raw = "$output_dir/$id-raw.fasta";
	my $family_aligned = "$output_dir/$id-aligned.nex";
	my $family_reduced = "$output_dir/$id-reduced.fasta";

	# For cleaning up temp files in case of an interruption
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

	# Run blastp on the raw sequences
	system("$makeblastdb -in $family_raw -input_type fasta -dbtype prot >/dev/null") and die "Error running makeblastdb: $!.\n";
	system("$blastp -db $family_raw -query $family_raw -out=$family_raw.out -outfmt '6 qseqid sseqid qstart qend sstart send' -evalue 1e-05 >/dev/null") and die "Error running blastp: $!.\n";

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

		# Split line on tabs, assign corresponding values
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

			# If the hits do not overlap, get their union
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
	system("$muscle -in $family_reduced -out $family_aligned >/dev/null 2>&1") and die "Error running muscle: $!.\n";

	# Parse Muscle alignment
	my %family_aligned = parse_fasta($family_aligned);

	# Check pairwise identity between all members
	foreach my $index1 (0 .. scalar(keys %family_aligned) - 2) {
		my $taxon1 = (keys %family_aligned)[$index1];
		foreach my $index2 ($index1 + 1 .. scalar(keys %family_aligned) - 1) {
			my $taxon2 = (keys %family_aligned)[$index2];
			my $pid = get_pid($family_aligned{$taxon1}, $family_aligned{$taxon2});

			# Stop analysis if pid between a pair of sequences is too
			if ($pid < $pid_cutoff) {
				unlink(@unlink);
				die "[$id] Pairwise similarity between sequences was too low ($pid).\n";
			}
		}
	}

	# Reverse translate from amino acid to nucleotide
	my %family_aligned_nuc;
	foreach my $taxon (keys %family_aligned) {
		$family_aligned_nuc{$taxon} = reverse_translate({"DNA" => $nuc_sequences{$taxon}, "PROT" => $family_aligned{$taxon}});
	}

	#write_phylip({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned});
	write_phylip({'OUT' => $family_aligned, 'ALIGN' => \%family_aligned_nuc});

	# Run RAxML
	system("raxmlHPC -f a -m GTRGAMMA -p 123123 -x 123123 -#100 -s '$family_aligned' -n $id >/dev/null");

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

		# Check if we met the threshold
		if ($bootstrap < $boot_cutoff) {
			unlink(@unlink);
			die "[$id] RAxML tree does not have high enough support.\n";
		}

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

	# Create CSV table for each species
	foreach my $name (@names) {

		my %contigs;
		foreach my $contig (keys %original_names) {
			if ($contig =~ /\Q$name\E_\d+/) {
				$contigs{$original_names{$contig}} = $contig;
			}
		}

		# Open CSV table for appending, use flock to prevent simultaneous writes
		open(my $topo_output, ">>", "../$name.csv") || die "Could not open '$name.csv': $!.\n";
		flock($topo_output, LOCK_EX) || die "Could not lock '$name.csv': $!.\n";
		seek($topo_output, 0, SEEK_END) || die "Could not seek '$name.csv': $!.\n";

		# Parse KaKs_Calculator output to determine original Ks of species pair
		#open(my $kaks_output, "<", "../$name.kaks");
		open(my $kaks_output, "<", "../$name.$model.kaks");
		while (my $line = <$kaks_output>) {
			chomp($line);

			# Skip header
			next if ($line =~ /^Sequence/);

			my @line = split(/\s+/, $line);

			if ($line[0] =~ /q_(.*)_t_(.*)/) {
				if (exists($contigs{$1}) && exists($contigs{$2})) {
					my $ks = $line[3];
					$ks = 0 if ($ks eq "NA");
					
					# Output Ks as well as its WGD support status
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
	}

#	open(my $out_file, ">>", $out_file_name);
#	flock($out_file, LOCK_EX) || die "Could not lock '$out_file': $!.\n";
#	seek($out_file, 0, SEEK_END) || die "Could not seek '$out_file': $!.\n";
#
#	print {$out_file} "$id: $tree\n";
#
#	flock($out_file, LOCK_UN) || die "Could not unlock '$out_file': $!.\n";
#	close($out_file);

	# Cleanup
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

sub run_transdecoder {
	my $transcriptome = shift;

	# TransDecoder working directory name
	(my $transcriptome_name = $transcriptome) =~ s/(.*\/)?(.*)/$2/;
	my $transdecoder_out_dir = "transdecoder-$transcriptome_name";

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
			remove_tree("$transcriptome_name.transdecoder_dir/");
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
			#$taxon =~ s/-/_/g;
			#$align->{$taxon} .= $line;
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

sub usage {
	return "Usage: perl wgd-plot.pl -t [TRANSCRIPTOME1 TRANSCRIPTOME2]...\n";
}

sub help {
print <<EOF; 
@{[usage()]}
Generate a Ks plot for a given transcriptome in fasta format

  -t, --transcriptomes              file names of at least two transcriptomes (in FASTA format) to use for analyses
  -m, --model                       model used by KaKs_Calculator to determine Ks (default: YN)
                                        valid models: NG, LWL, LPB, MLWL, MLPB, GY, YN, MYN, MS, MA. See KaKs_Calculator doc for details.
  -o, --output                      name of the directory to store output files in (default: "wgd-plot-" + Unix time of script invocation)
  -l, --min-length                  the minimum alignment length of paralogous sequences (default: 300 bp)
  -p, --pid-cut                     minimum pairwise percent identity allowed between two members of a quartet (default: 60)
  -b, --boot-cut                    minimum RAxML bootstrap support required to use a quartet (default: 70)
  --pfam-cpus                       the number of CPUs to let hmmscan using during pfam analysis (default: current number of free CPUs)
  --pfam-search                     full path to pfam binary for usage in pfam search
  -T, --n-threads                   the number of CPUs to use during analysis (default: current number of free CPUs)
  -h, --help                        display this help and exit

Example:
  perl wgd-plot.pl -t trans1.fa trans2.fa                                 create a sliding window plot for trans1.fa trans2.fa showing WGD support as a function of Ks 


Mail bug reports and suggestions to <noah.stenz.github\@gmail.com>
EOF
exit(0);
}
