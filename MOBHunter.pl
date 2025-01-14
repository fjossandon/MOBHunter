#!/usr/bin/perl
=head1 NAME

MobHunter

=head1 AUTHOR

Francisco J. Ossandon <fco.j.ossandon(at)gmail.com>

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
# Version 1.0.0 Jan 2025: First public release

use strict;
use warnings;
use File::Copy;
use File::Path;
use Bio::Root::Version 1.006925;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::SearchIO;
use Cwd;
our $VERSION = 'v1.0.0';


# Check if the working directory is provided as a command-line argument
if (@ARGV != 1) {
    die "Usage: perl $0 <working_directory>\n";
}

# Get the working directory from the command-line argument
my $work_dir = $ARGV[0];

# Change to the provided working directory
chdir $work_dir or die "Could not change to directory '$work_dir': $!\n";


# Programs folders
my $alien_hunter_dir = '/home/bioinfo/mobhunter/MobHunter_v1/MobHunter/Alien_Hunter-1.8';
my $aragorn_dir      = '/home/bioinfo/mobhunter/MobHunter_v1/MobHunter/Aragorn1.2.36';
my $island_path_dir  = '/home/bioinfo/mobhunter/MobHunter_v1/MobHunter/IslandPath-DIMOB_v0.3';
my $paida_dir        = '/home/bioinfo/mobhunter/MobHunter_v1/MobHunter/PAI-DA';
my $trna_scan_se_dir = '/home/bioinfo/mobhunter/MobHunter/tRNAscan-SE_1.3.2';
my $phage_finder_dir = '/home/bioinfo/mobhunter/MobHunter_v1/MobHunter/Phage_Finder_v2.2';
my $phispy_dir       = '/home/bioinfo/mobhunter/MobHunter_v1/MobHunter/phiSpyNov11_v2.3';

# Databases
my $conj_scan_hmm = '/home/bioinfo/mobhunter/MobHunter_v1/MobHunter/CONJscan-T4SSscan/CONJscan.hmm3';
my $tnp_pred_hmm  = '/home/bioinfo/mobhunter/MobHunter_v1/MobHunter/TnpPred_hmm/TnpPred_HMM_Profiles.hmm';

# Optional file 'keywords.txt' for keywords search on CDS "/product" tag from annotation,
# assumed to be located in the same folder than sequences. Format is 1 keyword or phrase per line
# using case insensitive complete match ('m/\b$keyword\b/i').
# Examples: 'integrase', 'hypothetical protein', 'EmrB/QacA'
# NOTE FOR DEVS: For ease of use the 'quotemeta' function is used, if you want to use
#  regular expressions on keywords you must comment it out just before the match is done
my $keyword_file = 'keywords.txt';
my %keywords;
if (-e $keyword_file) {
    open my $fh_KEYWORDS, '<', $keyword_file or die "Could not read '$keyword_file': $!\n";
    while (my $line = <$fh_KEYWORDS>) {
        chomp $line;
        next if ($line eq '');
        # Lines that start with '#' are considered comments and skipped
        next if ($line =~ m/^#/);

        # Clear extra whitespaces just in case
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        $keywords{$line} = '';
    }
    close $fh_KEYWORDS;
}

# Check which Blast package is available
my $blast_cmd   = `rpsblast 2>&1`; # Only binary present in both packages
my $blast_pkg   = ($blast_cmd =~ m#BLAST query/options error:#)  ? 'blast+'
                : ($blast_cmd =~ m#rpsblast\s+\S+\s+arguments:#) ? 'blast'
                :                                                   die "No BLAST program seems available, please install BLAST+ or BLAST to continue!\n";

# Genbank (.gbk) or Fasta (.fna) files are needed, others will be derived from .gbk if available
opendir INPUT_DIR, $work_dir or die "Could not read folder '$work_dir': $!\n";
my @input_files = grep {/\.(?:fa|fna|fasta|faa|ffn|gbk?|ptt)$/i} readdir INPUT_DIR;
closedir INPUT_DIR;

# Group together multiple filetypes of the same sequences
my %input_sequences;
foreach my $input_file (@input_files) {
    my ($header, $ext) = ($input_file =~ m/^(.+)\.(\w+)$/);
    $ext = lc $ext;

      ($ext eq 'fa')    ? ($input_sequences{$header}{fasta_fna} = $input_file)
    : ($ext eq 'fna')   ? ($input_sequences{$header}{fasta_fna} = $input_file)
    : ($ext eq 'fasta') ? ($input_sequences{$header}{fasta_fna} = $input_file)
    : ($ext eq 'faa')   ? ($input_sequences{$header}{fasta_faa} = $input_file)
    : ($ext eq 'ffn')   ? ($input_sequences{$header}{fasta_ffn} = $input_file)
    : ($ext eq 'gb')    ? ($input_sequences{$header}{genbank}   = $input_file)
    : ($ext eq 'gbk')   ? ($input_sequences{$header}{genbank}   = $input_file)
    : ($ext eq 'ptt')   ? ($input_sequences{$header}{table_ptt} = $input_file)
    :                     next;
}
# Remove Genbank files produced by the subroutines that could be already present from the hash keys
foreach my $key (keys %input_sequences) {
    delete $input_sequences{$key}
        if ($key =~ m/(?: AlienHunter  | Aragorn    | Blastn      | CONJscan    | gaps
                       |  ICEberg_ICEs | IslandPath | keywords    | PAI-DA      | Phage_Finder
                       |  PhiSpy       | TnpPred    | tRNA_frags  | tRNAscan-SE | Consensus )/x);
}

# Process each input sequence
my $counter    = 0;
my $total_keys = scalar keys %input_sequences;
foreach my $sequence (sort {$a cmp $b} keys %input_sequences) {
    $counter++;
    my %results;
    my %db_results;
    my $results_counter    = 0;
    my $db_results_counter = 0;
    my $codontable         = 11; # Default value, later updated if possible

    # Prepare input files necessary for the subroutines
    print "\nPreparing '$sequence' for analysis... ($counter / $total_keys)\n";

    # Store description and use taxonomy separated in columns
    my $taxonomy    = join "\t", ('unclassified') x 9;
    my $description = '';

    # PhiSpy will die if an input Genbank file contains 40 CDS features or less,
    # so keep record of the CDS number before trying to run it
    my $genbank_CDS_number = 0;

    # Store possible hits from optional keyword search on annotation to print them later
    my %keyword_hits;

    # If the .gbk is available, generate .fna, .ffn and .faa if any is missing
    if (exists $input_sequences{$sequence}{genbank}) {
        my $file_obj_gbk = Bio::SeqIO->new(-file   => $input_sequences{$sequence}{genbank},
                                           -format => 'genbank');

        while (my $seq_obj = $file_obj_gbk->next_seq) {
            # In NCBI genomes, primary_id is the genome GI
            my $sequence_gi = $seq_obj->primary_id;
            my $locus       = $seq_obj->display_id;
            my $ncbi_locus  = "$locus." . ($seq_obj->version || 1);
            my $seq_length  = $seq_obj->length;
            my @tax_levels  = (eval { $seq_obj->species->classification   }) ? $seq_obj->species->classification : ();
            my $ncbi_spp    = (eval { $seq_obj->species->binomial('FULL') }) ? $seq_obj->species->binomial('FULL')
                            :                                                  '';
            $description    = $seq_obj->desc || $locus; # Copy $locus value if there is not description for the table
            $description    =~ s/\.$//;

            # Prepare Taxonomy columns
            @tax_levels = reverse @tax_levels; # BioPerl stores taxonomy in inverted order with Kingdom in the end, so reverse it
            my $kingdom       = 'unclassified';
            my $phylum        = 'unclassified';
            my $class         = 'unclassified';
            my $order         = 'unclassified';
            my $family        = 'unclassified';
            my $tribe         = 'unclassified';
            my $genus         = 'unclassified';
            my $species_group = 'unclassified';
            my $species       = 'unclassified';

            # Uncommon 2 taxonomy levels
            if (scalar @tax_levels == 2) {
                # Assume that only Kingdom and Species are known
                ($kingdom, $species) = @tax_levels;
            }
            # Uncommon 3 taxonomy levels
            elsif (scalar @tax_levels == 3) {
                # Based in the cases seen so far, just assume that the values are
                # Kingdom, Genus and Species
                ($kingdom, $genus, $species) = @tax_levels;
            }
            # Uncommon 4 taxonomy levels
            elsif (scalar @tax_levels == 4) {
                # Based in the cases seen so far, just assume that the values are
                # Kingdom, Phylum, Genus and Species
                ($kingdom, $phylum, $genus, $species) = @tax_levels;
            }
            # Uncommon 5 taxonomy levels
            elsif (scalar @tax_levels == 5) {
                # Based in the cases seen so far, just assume that the values are
                # Kingdom, Phylum, Class, Genus and Species
                ($kingdom, $phylum, $class, $genus, $species) = @tax_levels;
            }
            # Uncommon 6 taxonomy levels
            elsif (scalar @tax_levels == 6) {
                # Check if Class is missing by testing if the string in the respective column
                # ends in "ales" (typical of Order names)
                if ($tax_levels[2] =~ m/ales$/) {
                    ($kingdom, $phylum, $order, $family, $genus, $species) = @tax_levels;
                }
                # Check if Family is missing by checking if Order is present and the
                # last 2 columnsare genus and species
                elsif (    $tax_levels[3]  =~ m/ales$/ # Order name present in the right position
                       and $tax_levels[-2] =~ m/^(?:Candidatus )?\S+$/ # usually 1 word (excluding Candidatus)
                       and $tax_levels[-1] =~ m/^$tax_levels[-2]/      # typically for genus and species
                    ) {
                    ($kingdom, $phylum, $class, $order, $genus, $species) = @tax_levels;
                }
                # Here are cases where Order is missing, since it was not detected already at Class position
                elsif ($tax_levels[3]  !~ m/ales$/) {
                    ($kingdom, $phylum, $class, $family, $genus, $species) = @tax_levels;
                }
                # All the other are cases where Genus name is missing
                else {
                    ($kingdom, $phylum, $class, $order, $family, $species) = @tax_levels;
                }
            }
            # Most commonly used 7 taxonomy levels
            elsif (scalar @tax_levels == 7) {
                ($kingdom, $phylum, $class, $order, $family, $genus, $species) = @tax_levels;
            }
            # Uncommon 8 taxonomy levels
            elsif (scalar @tax_levels == 8) {
                # For 8 taxonomy levels, try to determinate if the extra level is 'tribe'
                # or 'species group'. For this check if the last 2 levels are genus and species
                if (    $tax_levels[-2] =~ m/^(?:Candidatus )?\S+$/ # usually 1 word (excluding Candidatus)
                    and $tax_levels[-1] =~ m/^$tax_levels[-2]/      # typically for genus and species
                    ) {
                    ($kingdom, $phylum, $class, $order, $family, $tribe, $genus, $species) = @tax_levels;
                }
                else {
                    ($kingdom, $phylum, $class, $order, $family, $genus, $species_group, $species) = @tax_levels;
                }
            }
            # Uncommon 9 taxonomy levels, examples from Rickettsia, Mycobacterium and Agrobacterium
            elsif (scalar @tax_levels == 9) {
                ($kingdom, $phylum, $class, $order, $family, $tribe, $genus, $species_group, $species) = @tax_levels;
            }

            # Join the all the taxonomy columns
            $taxonomy = join "\t", $kingdom, $phylum, $class, $order, $family, $tribe, $genus, $species_group, $species;

            # Get all CDS features and count them
            my @feat_CDS = grep {$_->primary_tag eq 'CDS'} $seq_obj->get_SeqFeatures;
            $genbank_CDS_number = scalar @feat_CDS;

            # Check .fna
            if (not exists $input_sequences{$sequence}{fasta_fna} or -z $input_sequences{$sequence}{fasta_fna}) {
                my $sequence_file = "$sequence.fna";
                my $file_obj_fna  = Bio::SeqIO->new(-file   => ">$sequence_file",
                                                    -format => 'fasta',
                                                    -width  => 70); # NCBI value

                my $display_id = $locus;
                if ($sequence_gi =~ m/^\d+$/) {
                    $display_id = ($ncbi_locus =~ m/^N[A-Z]_/) ? "gi|$sequence_gi|ref|$ncbi_locus|"
                                   :                             "gi|$sequence_gi|gb|$ncbi_locus|";
                }
                my $seq_obj_fna = Bio::Seq->new (-display_id => $display_id,
                                                 -desc       => $description,
                                                 -seq        => $seq_obj->seq);
                $file_obj_fna->write_seq($seq_obj_fna);
                close $file_obj_fna->_fh;

                $input_sequences{$sequence}{fasta_fna} = $sequence_file;
            }

            # Check .ptt .ffn and .faa
            my $feats_checked = 0;
            if (   not exists $input_sequences{$sequence}{fasta_ffn} or -z $input_sequences{$sequence}{fasta_ffn}
                or not exists $input_sequences{$sequence}{fasta_faa} or -z $input_sequences{$sequence}{fasta_faa}
                or not exists $input_sequences{$sequence}{table_ptt} or -z $input_sequences{$sequence}{table_ptt}
                or scalar keys %keywords > 0
                ) {
                $feats_checked = 1;
                # Generate .ffn and .faa together
                my ($ffn_file, $file_obj_ffn);
                my ($faa_file, $file_obj_faa);
                if (   not exists $input_sequences{$sequence}{fasta_ffn} or -z $input_sequences{$sequence}{fasta_ffn}
                    or not exists $input_sequences{$sequence}{fasta_faa} or -z $input_sequences{$sequence}{fasta_faa}
                    ) {
                    $ffn_file = "$sequence.ffn";
                    $faa_file = "$sequence.faa";

                    $file_obj_ffn = Bio::SeqIO->new(-file   => ">$ffn_file",
                                                    -format => 'fasta',
                                                    -width  => 70); # NCBI value
                    $file_obj_faa = Bio::SeqIO->new(-file   => ">$faa_file",
                                                    -format => 'fasta',
                                                    -width  => 70); # NCBI value
                }

                # Generate .ptt
                my $table_file;
                my %table_entries;
                my $table_counter = 0;
                my $num_genes     = 0;
                if (not exists $input_sequences{$sequence}{table_ptt} or -z $input_sequences{$sequence}{table_ptt}) {
                    $table_file = "$sequence.ptt";
                }

                foreach my $feat (@feat_CDS) {
                    # Feature location
                    my $feat_loc_obj = $feat->location;
                    my $feat_strand  = $feat_loc_obj->strand;
                    my $ptt_strand   = ($feat_strand == -1) ? '-' : '+';

                    my ($feat_start, $feat_end, $ffn_coords);
                    if ( $feat_loc_obj->isa('Bio::Location::SplitLocationI') ) {
                        my @sub_locs = $feat_loc_obj->sub_Location;
                        $feat_start = $sub_locs[0]->start;
                        $feat_end   = $sub_locs[-1]->end;

                        my @loc_strings;
                        if ($feat_strand == -1) {
                            @sub_locs = reverse @sub_locs;
                            foreach my $sub_loc (@sub_locs) {
                                my $sub_start = $sub_loc->start;
                                my $sub_end   = $sub_loc->end;
                                push @loc_strings, "c$sub_end-$sub_start";
                            }
                        }
                        else {
                            foreach my $sub_loc (@sub_locs) {
                                my $sub_start = $sub_loc->start;
                                my $sub_end   = $sub_loc->end;
                                push @loc_strings, "$sub_start-$sub_end";
                            }
                        }
                        $ffn_coords = ":" . join ",", @loc_strings;
                    }
                    else {
                        $feat_start = $feat_loc_obj->start;
                        $feat_end   = $feat_loc_obj->end;
                        $ffn_coords = ($feat_strand == -1) ? ":c$feat_end-$feat_start" : ":$feat_start-$feat_end";
                    }

                    if ($feat_start > $feat_end) {
                        # Reduce $genbank_CDS_number by 1 because this gene will by skipped
                        # by PhiSpy script, since its gene coordinates limitations don't work
                        # well with split by origin genes coordinates (e.g. join(43500..44459,1..3192)
                        $genbank_CDS_number--;
                    }

                    # Tags
                    my $locus_tag    = $feat->has_tag('locus_tag')    ? ($feat->get_tag_values('locus_tag'))[0]
                                     :                                   '';
                    my $gene_tag     = $feat->has_tag('gene')         ? ($feat->get_tag_values('gene'))[0]
                                     :                                   '';
                    my $product      = $feat->has_tag('product')      ? ($feat->get_tag_values('product'))[0]
                                     :                                   '';
                    my $codon_start  = $feat->has_tag('codon_start')  ? ($feat->get_tag_values('codon_start'))[0]
                                     :                                   1;
                    my $pseudo       = $feat->has_tag('pseudo')       ?  1
                                     :                                   0;
                    my $transl_table = $feat->has_tag('transl_table') ? ($feat->get_tag_values('transl_table'))[0]
                                     :                                   11;
                    $codontable      = $transl_table;
                    my $prot_id      = '';

                    # For genomes like Pandoraea pnomenusa 3kgm (NC_022904), that have CDS with
                    # 2 protein ids (e.g. "REF_PRJNA226227:U875_00030" and "YP_008834357.1")
                    if ($feat->has_tag('protein_id')) {
                        my @prot_ids = $feat->get_tag_values('protein_id');
                        my $best_id  = '';
                        if (scalar @prot_ids > 1) {
                            foreach my $id (@prot_ids) {
                                if (   $id =~ m/^[A-Z]{2}_\d+\.\d+$/ # RefSeq format
                                    or $id =~ m/^[A-Z]+\d+\.\d+$/    # Genbank format
                                    ) {
                                    $best_id = $id;
                                    last;
                                }
                            }
                            $prot_id = ($best_id ne '') ? $best_id : $prot_ids[0];
                        }
                        else {
                            $prot_id = $prot_ids[0];
                        }
                    }
                    else {
                        $prot_id   = $feat->has_tag('locus_tag')  ? ($feat->get_tag_values('locus_tag'))[0]
                                   : $feat->has_tag('gene')       ? ($feat->get_tag_values('gene'))[0]
                                   :                                 '';
                    }

                    $product =~ s/(?<=\S)-\s/-/;
                    $prot_id =~ s/\s+/_/g;

                    # Perform optional keyword search now before continuing
                    if (scalar keys %keywords > 0) {
                        # Check if product match with 1 or more keywords
                        my $match = 0;
                        foreach my $keyword (keys %keywords) {
                            $keyword = quotemeta $keyword;
                            $match++ if ($product =~ m/\b$keyword\b/i);
                        }
                        if ($match > 0) {
                            my $hit_loc_str   = $feat_loc_obj->to_FTstring;
                            my $hit_loc_index = $feat_start + (1 - $feat_end / $seq_length);
                            my $feature       = Bio::SeqFeature::Generic->new(-primary_tag => 'misc_feature',
                                                                              -location    => $feat_loc_obj,
                                                                              -tag => { inference => 'keyword_match',
                                                                                        note      => "ID: $prot_id; Product: $product",
                                                                                        colour    => '200 200 100',
                                                                                       },
                                                                              );
                            # Store feature
                            $keyword_hits{$hit_loc_str}{feature} = $feature;
                            $keyword_hits{$hit_loc_str}{loc_idx} = $hit_loc_index;
                        }
                    }

                    # Skip pseudogenes without protein ID
                    if ($prot_id eq '' and $pseudo == 1) {
                        next;
                    }
                    elsif ($prot_id eq '') {
                        # Since it have no protein ID, create a custom ID
                        # for it to not cause trouble downstream
                        $prot_id = $feat_start . "_" . $feat_end;
                    }

                    # Nucleotide and aminoacid sequence
                    my $feat_seq = $feat->spliced_seq(-nosort => 1);
                    my $nt_seq   = $feat_seq->seq;
                    my $aa_seq   = $feat_seq->translate(-codontable_id => $transl_table,
                                                        -frame         => $codon_start -1,
                                                        )->seq;
                    if ($feat_seq->seq =~ m/^[CTG]TG/) {
                        $aa_seq =~ s/^\w/M/;
                    }
                    # Delete stop codon before measuring peptide length
                    $aa_seq =~ s/\*+$//;
                    my $aa_len = length $aa_seq;

                    # A GI number is needed in one of the .ptt table columns, and is also
                    # needed in the display_id field of .ffn and .faa sequences to comply
                    # with IslandPath-DIMOB code (.ffn use Genome GI, .faa use Protein GI)
                    my $gi_id = '';
                    my @db_xref = $feat->has_tag('db_xref') ? $feat->get_tag_values('db_xref') : ();
                    # For NCBI annotations with GI
                    if (@db_xref) {
                        foreach my $xref (@db_xref) {
                            if ($xref =~ m/^GI:(\d+)/) {
                                $gi_id = $1;
                                last;
                            }
                        }
                    }
                    # For ORFminer annotations without GI
                    if ($gi_id eq '') {
                        if ($locus_tag =~ m/^aorf_(\d+)/) {
                            $gi_id = $1;
                            if ($locus_tag =~ m/\d+_(\d+)/) {
                                $gi_id .= $1;
                            }
                        }
                        else {
                            $gi_id = "$feat_start$feat_end";
                        }
                    }

                    # Complete .ffn and .faa process for the feature
                    if (defined $ffn_file and defined $faa_file) {
                        my $prod = '';
                        if (   (   $product eq 'uncharacterized protein'
                                or $product eq 'hypothetical protein')
                            and $locus_tag ne ''
                            ) {
                            $prod = "$product $locus_tag";
                        }
                        else {
                            $prod = $product;
                        }

                        # Assemble complete display_ids
                        my $display_id_faa  = '';
                        my $display_id_ffn  = '';
                        my $description_ffn = '';
                        my $description_faa = ($ncbi_spp ne '') ? "$prod [$ncbi_spp]" : $prod;

                        if ($prot_id =~ m/^[A-Z]{2}_\d+\.\d+$/) { # RefSeq format
                            $display_id_faa  = "gi|$gi_id|ref|$prot_id|";
                            $display_id_ffn  = "gi|$sequence_gi|ref|$ncbi_locus|";
                            $description_ffn = $description;
                        }
                        elsif ($prot_id =~ m/^[A-Z]+\d+\.\d+$/) { # Genbank format
                            $display_id_faa = "gi|$gi_id|gb|$prot_id|";
                            $display_id_ffn = "gi|$sequence_gi|gb|$ncbi_locus|";
                            $description_ffn = $description;
                        }
                        else {
                            $display_id_faa = "gi|$gi_id|$prot_id|";
                            $display_id_ffn = "gi|$gi_id|$prot_id|";
                            $description_ffn = ($prod ne '') ? $prod : $description;
                        }

                        # Display_id for .ffn files must incorporate coordinates
                        $display_id_ffn .= $ffn_coords;

                        # Fasta output
                        my $nt_seq_obj = Bio::Seq->new (-display_id => $display_id_ffn,
                                                        -desc       => $description_ffn,
                                                        -seq        => $nt_seq);
                        my $aa_seq_obj = Bio::Seq->new (-display_id => $display_id_faa,
                                                        -desc       => $description_faa,
                                                        -seq        => $aa_seq);
                        $file_obj_ffn->write_seq($nt_seq_obj);
                        $file_obj_faa->write_seq($aa_seq_obj) if ($pseudo == 0);
                    }

                    # Complete .ptt process for the feature
                    if (defined $table_file and $pseudo == 0) {
                        $num_genes++;

                        $locus_tag = '-' if $locus_tag eq '';
                        $gene_tag  = '-' if $gene_tag  eq '';
                        $product   = '-' if $product   eq '';
                        # Note: COG code is not established because its absent from the usual annotation
                        my $entry = join("\t", "$feat_start..$feat_end", $ptt_strand, $aa_len, $gi_id,
                                                $gene_tag, $locus_tag, '-', '-', $product
                                         ) . "\n";

                        my $entry_loc_index = $feat_start + (1 - $feat_end / $seq_length);
                        $table_counter++;
                        $table_entries{$table_counter}{loc_idx} = $entry_loc_index;
                        $table_entries{$table_counter}{entry}   = $entry;
                    }
                }

                # Finish .ffn and .faa process
                if (defined $ffn_file and defined $faa_file) {
                    close $file_obj_ffn->_fh;
                    close $file_obj_faa->_fh;

                    # Only register the sequence files if they are not empty.
                    # Delete them if the Genbank file didn't have any CDS
                    if (-z $ffn_file or -z $faa_file) {
                        unlink $ffn_file or die "Could not delete file '$ffn_file': $!\n";
                        unlink $faa_file or die "Could not delete file '$faa_file': $!\n";
                    }
                    else {
                        $input_sequences{$sequence}{fasta_ffn} = $ffn_file;
                        $input_sequences{$sequence}{fasta_faa} = $faa_file;
                    }
                }

                # Finish .ptt process
                if (defined $table_file) {
                    # Phage Finder splits first line by space and use $a[1] =~ /^\D/ in "check_infofile"
                    # to certify that the file is a Genbank ".ptt" file, so if there is a second word
                    # and it start with a number, add an extra undescore to prevent failing that match
                    my $temp_desc = $description;
                    my @words     = split /\s+/, $temp_desc;
                    if (scalar @words > 1 and $words[1] =~ /^\d/) {
                        $words[1]  = "_$words[1]";
                        $temp_desc = join " ", @words;
                    }

                    # Prepare the PTT table
                    open my $fh_PTT_TABLE, '>', $table_file or die "Could not write file '$table_file': $!\n";
                    print $fh_PTT_TABLE "$temp_desc - 1..$seq_length\n";
                    print $fh_PTT_TABLE "$num_genes proteins\n";
                    print $fh_PTT_TABLE "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";

                    # Print sorted by Start
                    foreach my $key (sort{$table_entries{$a}{loc_idx} <=> $table_entries{$b}{loc_idx}} keys %table_entries) {
                        my $entry = $table_entries{$key}{entry};
                        print $fh_PTT_TABLE $entry;
                    }
                    close $fh_PTT_TABLE;

                    $input_sequences{$sequence}{table_ptt} = $table_file;
                }
            }
            # If features where not processed because .faa, .ffn and .ptt files already exists,
            # make a quick check to adjust $genbank_CDS_number for PhiSpy
            if ($feats_checked == 0) {
                foreach my $feat (@feat_CDS) {
                    my $feat_start  = $feat->start;
                    my $feat_end    = $feat->end;
                    if ($feat_start > $feat_end) {
                        # Reduce $genbank_CDS_number by 1 because this gene will by skipped
                        # by PhiSpy script, since its gene coordinates limitations don't work
                        # well with split by origin genes coordinates (e.g. join(43500..44459,1..3192)
                        $genbank_CDS_number--;
                    }

                    # Recover transl_table value
                    my $transl_table = $feat->has_tag('transl_table') ? ($feat->get_tag_values('transl_table'))[0]
                                     :                                   11;
                    $codontable      = $transl_table;
                }
            }

            # Only 1 sequence per file allowed
            last;
        }
        close $file_obj_gbk->_fh;
    }

    if (not exists $input_sequences{$sequence}{fasta_fna}) {
        print "  No .fna fasta file found, skipping it...\n";
        next;
    }

    # Prepare the fasta sequence object to pass it to the subs
    my $sequence_obj;
    my $sequence_length = 0;
    my $fa_obj = Bio::SeqIO->new(-file   => $input_sequences{$sequence}{fasta_fna},
                                 -format => 'fasta');
    my $seq_counter = 0;
    while (my $seq_obj = $fa_obj->next_seq) {
        $seq_counter++;
        my $display_id   = $seq_obj->display_id;
        my $description  = $seq_obj->desc;
        my $sequence     = $seq_obj->seq;
        $sequence_length = $seq_obj->length;

        if ($display_id =~ m/^gi\|(\d+)\|[a-z]+\|([\.\w]+)\|/) {
            my $gi_id      = $1;
            my $ncbi_locus = $2;
            $sequence_obj  = Bio::Seq::RichSeq->new(-primary_id => $gi_id,
                                                    -desc       => $description,
                                                    -seq        => $sequence);

            my ($locus, $version) = ($ncbi_locus =~ m/^(\w+)\.(\d+)$/);
            if (defined $locus) {
                $sequence_obj->accession_number($locus);
                $sequence_obj->display_id($locus);
                $sequence_obj->seq_version($version);
            }
            else {
                $sequence_obj->display_id($ncbi_locus);
            }
        }
        else {
            $sequence_obj = Bio::Seq->new(-display_id => $display_id,
                                          -desc       => $description,
                                          -seq        => $sequence);
        }
    }
    close $fa_obj->_fh;

    if ($seq_counter == 0) {
        die "No sequences found in the fasta file '$input_sequences{$sequence}{fasta_fna}'!\n";
    }
    elsif ($seq_counter > 1) {
        die "Only 1 sequence allowed becase some programs (like Alien Hunter or PAI-DA)\n"
          . "can't process more than 1 correctly or will only process 1!\n";
    }

    # Result table collector for the sequence
    my $result_table_file = "$work_dir/$sequence.result_table.txt";
    my @result_table;
    my $table_header = join("\t", 'Description', 'Locus', 'Kingdom', 'Phylum', 'Class', 'Order',
                                  'Family', 'Tribe', 'Genus', 'Species group', 'Species',
                                  'Start', 'End', 'Coordinates', 'Predictor', 'Score or Note') . "\n";
    push @result_table, $table_header;

    # If a keyword search was performed, generate a Genbank file that map the hits
    my $keyword_result = "$sequence.keywords.gbk";
    if (scalar keys %keywords > 0) {
        print "Reporting optional keyword search in '$sequence'...\n";

        # Print found keyword hits
        if (scalar keys %keyword_hits > 0) {
            my $out_obj = $sequence_obj->clone;

            # Add features sorted by a location index
            foreach my $loc (sort {
                                   ($keyword_hits{$a}{loc_idx})
                                   <=>
                                   ($keyword_hits{$b}{loc_idx})
                                   }
                             keys %keyword_hits
                ) {
                my $feature = $keyword_hits{$loc}{feature};
                $out_obj->add_SeqFeature($feature);
            }

            # Print object with results in a Genbank file
            my $gbk_obj = Bio::SeqIO->new(-file   => ">$keyword_result",
                                          -format => 'genbank');
            $gbk_obj->write_seq($out_obj);
            close $gbk_obj->_fh;
        }
    }
    if (-e $keyword_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $keyword_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;
                my $note    = $feat->has_tag('note')  ? ($feat->get_tag_values('note'))[0]
                            :                            '';

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'keywords', $note) . "\n";
                push @result_table, $line;
            }
        }
        close $result_obj->_fh;
    }

    # If sequence have Ns (unknown nucleotide), generate a Genbank file that maps them
    my $gaps_result  = "$sequence.gaps.gbk";
    print "Looking for gaps in '$sequence'...\n";
    if (not -e $gaps_result or -z $gaps_result) {
        my $check_seq = $sequence_obj->seq;
        if ($check_seq =~ m/N/i) {
            my %gaps;

            # First check if the entry sequence is a pseudochromosome (multiple contigs fused by
            # a common joining sequence) by searching the joining sequence. If exists, keep
            # their location in the anotation using the 'assembly_gap' primary tag
            my $seq_start    = 1;
            my $seq_end      = length $check_seq;
            my $join_seq     = 'NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN';

            while ($check_seq =~ m/(?:$join_seq){1,}/poig) {
                my $end_pos   = pos($check_seq);
                my $start_pos = $end_pos - length( ${^MATCH} ) + 1;

                my $gap_loc_string = "$start_pos..$end_pos";
                my $gap_loc_index  = $start_pos + (1 - $end_pos / $seq_end);

                # Build a custom feature GAP using black color
                my $feature = Bio::SeqFeature::Generic->new (
                                -start       => $start_pos,
                                -end         => $end_pos,
                                -primary_tag => 'assembly_gap',
                                -tag => {estimated_length => 'unknown',
                                         gap_type         => 'between scaffolds',
                                         colour           => '0 0 0'},
                                );
                # Store feature
                $gaps{$gap_loc_string}{feature} = $feature;
                $gaps{$gap_loc_string}{domain}  = 'ASSEMBLY_GAP';
                $gaps{$gap_loc_string}{loc_idx} = $gap_loc_index;
            }
            my @assembly_locs = keys %gaps;

            # Then check for Ns strings and label them as internal contig gaps if their length
            # is above a minimum cutoff (to reduce noise from single/very short Ns)
            while ($check_seq =~ m/(?:N){1,}/pig) {
                my $end_pos   = pos($check_seq);
                my $match_len = length( ${^MATCH} );
                my $start_pos = $end_pos - $match_len + 1;
                next if $match_len < 10;

                # Check if the gap found belongs to an already registered assembly gap,
                # and skip it in that case
                my $match_found = 0;
                foreach my $assembly_loc (@assembly_locs) {
                    my ($assembly_start, $assembly_end) = ($assembly_loc =~ m#(\d+)\.\.(\d+)#);
                    if (   ($assembly_start <= $start_pos and $start_pos <= $assembly_end)
                        or ($assembly_start <= $end_pos   and $end_pos   <= $assembly_end)
                        ) {
                        $match_found = 1;
                        last;
                    }
                }
                next if ($match_found == 1);

                my $i_gap_loc_string = "$start_pos..$end_pos";
                my $i_gap_loc_index  = $start_pos + (1 - $end_pos / $seq_end);

                # Construir un feature del gap a medida usando el tag primario "gap"
                my $estimated_length = ($match_len == 100) ? 'unknown' : $match_len;
                my $feature = Bio::SeqFeature::Generic->new (
                                -start       => $start_pos,
                                -end         => $end_pos,
                                -primary_tag => 'gap',
                                -tag => {estimated_length => $estimated_length},
                                );

                # Store feature
                $gaps{$i_gap_loc_string}{feature} = $feature;
                $gaps{$i_gap_loc_string}{domain}  = 'GAP';
                $gaps{$i_gap_loc_string}{loc_idx} = $i_gap_loc_index;
            }

            # Print found gaps
            if (scalar keys %gaps > 0) {
                my $out_obj = $sequence_obj->clone;

                # Add features sorted by a location index
                foreach my $loc (sort {
                                       ($gaps{$a}{loc_idx})
                                       <=>
                                       ($gaps{$b}{loc_idx})
                                       }
                                 keys %gaps
                    ) {
                    my $feature = $gaps{$loc}{feature};
                    $out_obj->add_SeqFeature($feature);
                }

                # Print object with results in a Genbank file
                my $gbk_obj = Bio::SeqIO->new(-file   => ">$gaps_result",
                                              -format => 'genbank');
                $gbk_obj->write_seq($out_obj);
                close $gbk_obj->_fh;
            }
        }
    }
    if (-e $gaps_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $gaps_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                my $primary_tag = $feat->primary_tag;
                my $loc_str     = $feat->location->to_FTstring;
                my $start       = $feat->start;
                my $end         = $feat->end;
                my $length      = $feat->has_tag('estimated_length') ? ($feat->get_tag_values('estimated_length'))[0]
                                :                                       '';
                my $length_str  = ($length ne 'unknown') ? "$length nts" : $length;
                my $gap_type = ($primary_tag eq 'gap') ? 'gap of Ns' : $primary_tag;

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'gaps', "$gap_type; length: $length_str") . "\n";
                push @result_table, $line;

                $results_counter++;
                my $loc_index = $start + (1 - $end / $sequence_length);
                $results{$results_counter}{program} = 'gaps';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = $length;
            }
        }
        close $result_obj->_fh;
    }

    # Prediction of Horizontal Gene Transfer (HGT)
    print "Analizing '$sequence' through Alien Hunter...\n";
    my $alien_result = "$sequence.AlienHunter.gbk";
    if (not -e $alien_result or -z $alien_result) {
        alien_hunter($alien_hunter_dir, $input_sequences{$sequence}{fasta_fna}, $sequence_obj);
    }
    if (-e $alien_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $alien_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                next if $feat->primary_tag ne 'misc_feature';

                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;
                my $score   = $feat->has_tag('score') ? ($feat->get_tag_values('score'))[0]
                            :                            '';
                my $note    = $feat->has_tag('note')  ? ($feat->get_tag_values('note'))[0]
                            :                            '';
                my ($threshold) = ($note =~ m/threshold: ([\.\d]+)/);
                my $ratio = $score / $threshold;
                $ratio    = sprintf "%.3f", $ratio;

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'Alien Hunter', "$ratio ($score)") . "\n";
                push @result_table, $line;

                $results_counter++;
                my $loc_index = $start + (1 - $end / $sequence_length);
                $results{$results_counter}{program} = 'AlienHunter';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = $ratio;
                if ($note =~ m/rRNA operon/) {
                    $results{$results_counter}{rRNA_operon} = '';
                }
            }
        }
        close $result_obj->_fh;
    }

    # Prediction of tRNA, mtRNA, and tmRNA genes
    print "Analizing '$sequence' through Aragorn...\n";
    my %aragorn_data; # To help discard redundant data from tRNAscan-SE for island consensus analysis
    my $aragorn_result = "$sequence.Aragorn.gbk";
    if (not -e $aragorn_result or -z $aragorn_result) {
        aragorn($aragorn_dir, $input_sequences{$sequence}{fasta_fna}, $sequence_obj, $codontable);
    }
    if (-e $aragorn_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $aragorn_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                next if $feat->primary_tag eq 'CDS';

                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;
                my $product = $feat->has_tag('product') ? ($feat->get_tag_values('product'))[0]
                            :                              '';

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'Aragorn', $product) . "\n";
                push @result_table, $line;

                my $loc_index = $start + (1 - $end / $sequence_length);
                $aragorn_data{$loc_str} = '';

                $results_counter++;
                $results{$results_counter}{program} = 'Aragorn';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = $product;
            }
        }
        close $result_obj->_fh;
    }

    # Prediction of tRNA, mtRNA, and tmRNA gene FRAGMENTS.
    # NOTE: An important difference from above is that only
    # non-intron predictions are used to search for fragments
    print "Analizing '$sequence' through Aragorn and Blastn to look for tRNA fragments...\n";
    my $aragorn_frags_result = "$sequence.tRNA_frags.gbk";
    if (not -e $aragorn_frags_result or -z $aragorn_frags_result) {
        aragorn_frags($aragorn_dir, $input_sequences{$sequence}{fasta_fna}, $sequence_obj);
    }
    if (-e $aragorn_frags_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $aragorn_frags_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                next if $feat->primary_tag eq 'CDS';

                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;
                my $product = $feat->has_tag('product') ? ($feat->get_tag_values('product'))[0]
                            :                              '';

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'tRNA_frags', $product) . "\n";
                push @result_table, $line;

                my $note = $feat->has_tag('note') ? ($feat->get_tag_values('note'))[0]
                         :                           '';
                my ($original_loc) = ($note =~ m/Location: (\S+)/);

                $results_counter++;
                my $loc_index = $start + (1 - $end / $sequence_length);
                $results{$results_counter}{program} = 'tRNA_frags';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = $note;
                if ($original_loc ne $loc_str) {
                    $results{$results_counter}{original_loc} = $original_loc;
                }
            }
        }
        close $result_obj->_fh;
    }

    # Genomic island prediction through IslandPath-DIMOB
    my $island_path_result = "$sequence.IslandPath.gbk";
    if (    exists $input_sequences{$sequence}{fasta_fna}
        and exists $input_sequences{$sequence}{fasta_faa}
        and exists $input_sequences{$sequence}{fasta_ffn}
        and exists $input_sequences{$sequence}{table_ptt}
        ) {
        # IMPORTANT: It needs annotation!!! For this program, the subroutine
        # will expect that $sequence_file is the '.fna' with the genome info
        # (that will be needed when generating the Genbank output), and it will
        # look for the presence of '.faa', '.ffn' and '.ptt' files with the same name
        # (that will be used by the IslandPath program).
        # IslandPath command example: "./dimob.pl your_file.faa your_file.ffn your_file.ptt"
        print "Analizing '$sequence' through IslandPath-DIMOB...\n";
        if (not -e $island_path_result or -z $island_path_result) {
            island_path($island_path_dir, $input_sequences{$sequence}{fasta_fna}, $sequence_obj);
        }
    }
    if (-e $island_path_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $island_path_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'IslandPath', 'N/A') . "\n";
                push @result_table, $line;

                $results_counter++;
                my $loc_index = $start + (1 - $end / $sequence_length);
                $results{$results_counter}{program} = 'IslandPath';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = 'N/A';
            }
        }
        close $result_obj->_fh;
    }

    # Genomic island prediction through PAI-DA
    my $paida_result = "$sequence.PAI-DA.gbk";
    if (exists $input_sequences{$sequence}{genbank}){
        # IMPORTANT: It needs annotation!!! In this case an annotated Genbank file
        # as input for the compositional analysis of the coding sequences
        # PAI-DA command example:
        # "./compscan.pl -i sequence.gbk -o seq.dat"
        # "./compscore.pl -i seq.dat -o seq.score"
        print "Analizing '$sequence' through PAI-DA...\n";
        if (not -e $paida_result or -z $paida_result) {
            paida($paida_dir, $input_sequences{$sequence}{genbank}, $sequence_obj);
        }
    }
    if (-e $paida_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $paida_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;
                my $score   = $feat->has_tag('score') ? ($feat->get_tag_values('score'))[0]
                            :                            '';

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'PAI-DA', $score) . "\n";
                push @result_table, $line;

                $results_counter++;
                my $loc_index = $start + (1 - $end / $sequence_length);
                $results{$results_counter}{program} = 'PAIDA';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = $score;
            }
        }
        close $result_obj->_fh;
    }

    # Prediction of tRNA
    print "Analizing '$sequence' through tRNAscan-SE...\n";
    my $trna_scan_result = "$sequence.tRNAscan-SE.gbk";
    if (not -e $trna_scan_result or -z $trna_scan_result) {
        trna_scan($trna_scan_se_dir, $input_sequences{$sequence}{fasta_fna}, $sequence_obj);
    }
    if (-e $trna_scan_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $trna_scan_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
          FEAT:
            foreach my $feat (@features) {
                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;
                my $product = $feat->has_tag('product') ? ($feat->get_tag_values('product'))[0]
                            :                              '';

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'tRNAscan-SE', $product) . "\n";
                push @result_table, $line;

                # Discard tRNAs already covered by Aragorn
                if (exists $aragorn_data{$loc_str}) {
                    # Exact match with known tRNA
                    next FEAT;
                }
                else {
                    my @known_trna_locs = keys %aragorn_data;
                    foreach my $loc (@known_trna_locs) {
                        my ($loc_start, $loc_end) = ($loc =~ m/(\d+)\.\.(\d+)/);

                        if (   ($start >= $loc_start  and  $start <= $loc_end)
                            or ($end   >= $loc_start  and  $end   <= $loc_end)
                            ) {
                            # Overlapping match with known tRNA
                            next FEAT;
                        }
                    }
                }

                $results_counter++;
                my $loc_index = $start + (1 - $end / $sequence_length);
                $results{$results_counter}{program} = 'tRNAscan';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = $product;
            }
        }
        close $result_obj->_fh;
    }

    # Prediction of prophage regions through Phage Finder
    my $phage_finder_result = "$sequence.Phage_Finder.gbk";
    if (    exists $input_sequences{$sequence}{fasta_fna}
        and exists $input_sequences{$sequence}{fasta_faa}
        and exists $input_sequences{$sequence}{table_ptt}
        ) {
        print "Analizing '$sequence' through Phage Finder...\n";
        if (not -e $phage_finder_result or -z $phage_finder_result) {
            phage_finder($phage_finder_dir, $input_sequences{$sequence}{fasta_fna}, $sequence_obj, $blast_pkg, $codontable);
        }
    }
    if (-e $phage_finder_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $phage_finder_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;
                my $note    = $feat->has_tag('note') ? ($feat->get_tag_values('note'))[0]
                            :                           '';
                my $species = ($note =~ m/no species hit/) ? 'unknown species'
                            :                                ($note =~ m/DB hit is ([^\(]+) \(/)[0];

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'Phage Finder', $species) . "\n";
                push @result_table, $line;

                $results_counter++;
                my $loc_index = $start + (1 - $end / $sequence_length);
                $results{$results_counter}{program} = 'Phage_Finder';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = $note;
                if ($note =~  m/att[LR]/) {
                    $results{$results_counter}{attach} = $loc_str;
                }
            }
        }
        close $result_obj->_fh;
    }

    # Prediction of prophage regions through PhiSpy (Python program!)
    my $phispy_result  = "$sequence.PhiSpy.gbk";
    my $python_version = `python -V 2>&1`;
    my $python_ready   = ($python_version =~ m/^Python 2\./) ? 1 : 0;
    if (    $python_ready == 1
        and exists $input_sequences{$sequence}{genbank}
        and $genbank_CDS_number > 40
        ) {
        print "Analizing '$sequence' through PhiSpy...\n";
        if (not -e $phispy_result or -z $phispy_result) {
            phispy($phispy_dir, $input_sequences{$sequence}{genbank}, $sequence_obj);
        }
    }
    if (-e $phispy_result) {
        my $result_obj = Bio::SeqIO->new(-file   => $phispy_result,
                                         -format => 'genbank');
        while (my $seq_obj = $result_obj->next_seq) {
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                my $loc_str = $feat->location->to_FTstring;
                my $start   = $feat->start;
                my $end     = $feat->end;
                my $note    = $feat->has_tag('note')  ? ($feat->get_tag_values('note'))[0]
                            :                            '';

                my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                      $loc_str, 'PhiSpy', $note) . "\n";
                push @result_table, $line;

                $results_counter++;
                my $loc_index = $start + (1 - $end / $sequence_length);
                $results{$results_counter}{program} = 'PhiSpy';
                $results{$results_counter}{loc_str} = $loc_str;
                $results{$results_counter}{loc_idx} = $loc_index;
                $results{$results_counter}{details} = $note;
            }
        }
        close $result_obj->_fh;
    }

    # Database sequence hits and parsing
    my %databases = (CONJscan      => $conj_scan_hmm,
                     TnpPred       => $tnp_pred_hmm);
    foreach my $database (sort {$a cmp $b} keys %databases) {
        my $query_file = $input_sequences{$sequence}{fasta_faa};
        next if (not defined $query_file);

        print "Searching '$sequence' for $database database sequences...\n";
        my $database_file   = $databases{$database};
        my $database_result = "$sequence.$database.gbk";
        if (not -e $database_result or -z $database_result) {
            database_search($database, $database_file, $query_file, $sequence_obj);
        }

        if (-e $database_result) {
            my @conjscan_hits;

            my $result_obj = Bio::SeqIO->new(-file   => $database_result,
                                             -format => 'genbank');
            while (my $seq_obj = $result_obj->next_seq) {
                my @features = $seq_obj->get_SeqFeatures;
                foreach my $feat (@features) {
                    my $loc_str = $feat->location->to_FTstring;
                    my $start   = $feat->start;
                    my $end     = $feat->end;
                    my $note    = $feat->has_tag('note') ? ($feat->get_tag_values('note'))[0]
                                :                           '';

                    my $line = join("\t", $description, $sequence, $taxonomy, $start, $end,
                                          $loc_str, $database, $note) . "\n";
                    push @result_table, $line;

                    my $loc_index = $start + (1 - $end / $sequence_length);
                    $db_results_counter++;
                    $db_results{$db_results_counter}{database} = $database;
                    $db_results{$db_results_counter}{loc_str}  = $loc_str;
                    $db_results{$db_results_counter}{loc_idx}  = $loc_index;
                    $db_results{$db_results_counter}{details}  = $note;

                    # For CONJscan hits, try to group T4SS genes into single regions that will later
                    # treated as ICE "islands" in the consensus island analysis. Is necessary to
                    # gather all hits before doing it.
                    if ($database eq 'CONJscan') {
                        my ($conjscan_hit) = ($note =~ m/; Hit: (.+);$/);
                        push @conjscan_hits, "$start\t$end\t$conjscan_hit";
                    }
                }
            }
            close $result_obj->_fh;

            # Analyze CONJscan hits to create a region an expand using a 30000 nt
            # distance window (based on examples found in ICEberg like 'CTn9343').
            # This value can reunite conjugated genes separated from long lengths,
            # but it also can fuse 2 different neighboring ICEs, so for distances
            # longer than 15000 nt test if the found hit belongs to other cluster.
            # This assumes that genes are already sorted by start.
            my %conjscan_groups;
            my $conjscan_counter = 1;

            # Keep track of the 3 key genes to decide on hits separated by very long distances
            my $relax_gene    = 0;
            my $coupling_gene = 0;
            my $atpase_gene   = 0;

            for (my $i = 0; $i < scalar @conjscan_hits; $i++) {
                my $curr_hit = $conjscan_hits[$i];
                my ($hit_start, $hit_end, $conjscan_id) = split /\t/, $curr_hit;
                my $hit_distance = 0;
                if (scalar keys %conjscan_groups != 0) {
                    $hit_distance = $hit_start - $conjscan_groups{$conjscan_counter}{end};
                }

                # First cluster
                if (scalar keys %conjscan_groups == 0) {
                    $conjscan_groups{$conjscan_counter}{start} = $hit_start;
                    $conjscan_groups{$conjscan_counter}{end}   = $hit_end;
                    $conjscan_groups{$conjscan_counter}{genes}{$conjscan_id} = '';
                    $relax_gene++    if $conjscan_id =~ m/MOB\w+/i;
                    $coupling_gene++ if $conjscan_id =~ m/(?:T4CP|TcpA)/i;
                    $atpase_gene++   if $conjscan_id =~ m/(?:VirB4|TraU)/i;
                }
                # There can be some variation on the separation of T4SS genes,
                # so use a maximum distance of 15000
                elsif ($hit_distance <= 15000) {
                    $conjscan_groups{$conjscan_counter}{end} = $hit_end if ($hit_end > $conjscan_groups{$conjscan_counter}{end});
                    $conjscan_groups{$conjscan_counter}{genes}{$conjscan_id} = '';
                    $relax_gene++    if $conjscan_id =~ m/MOB\w+/i;
                    $coupling_gene++ if $conjscan_id =~ m/(?:T4CP|TcpA)/i;
                    $atpase_gene++   if $conjscan_id =~ m/(?:VirB4|TraU)/i;
                }
                # For very long distances be more cautious and check if the hit belongs
                # to the current cluster or a new cluster (e.g. conjugation genes from ICEVchBan9
                # in CP001485 genome are internally separated by ~21000 nts at one point)
                elsif ($hit_distance > 15000 and $hit_distance <= 30000) {
                    my %next_cluster;

                    # Keep track of the 3 key genes
                    my $next_relax_gene    = 0;
                    my $next_coupling_gene = 0;
                    my $next_atpase_gene   = 0;

                    # Guess current gene cluster type from its mating pair formation (mfc) genes
                    my %curr_mfp_types;
                    foreach my $hit (keys %{ $conjscan_groups{$conjscan_counter}{genes} }) {
                        my $type = ($hit =~ m/^(\w|FA|FATA)_/)[0] || '';
                        next if $type eq '';

                        $curr_mfp_types{$type} = 0 if not exists $curr_mfp_types{$type};
                        $curr_mfp_types{$type}++;
                    }
                    my $curr_highest_value = (scalar keys %curr_mfp_types > 0) ? (sort {$b <=> $a} values %curr_mfp_types)[0]
                                           :                                      0;
                    # Use of Hash for cases when there is a tie for the top represented types
                    my %curr_top_types;
                    foreach my $type (keys %curr_mfp_types) {
                        $curr_top_types{$type} = '' if ($curr_mfp_types{$type} == $curr_highest_value);
                    }

                    for (my $j = $i; $j < scalar @conjscan_hits; $j++) {
                        my $check_hit = $conjscan_hits[$j];
                        my ($check_start, $check_end, $check_id) = split /\t/, $check_hit;
                        my $check_distance = 0;
                        if (scalar keys %next_cluster != 0) {
                            $check_distance = $check_start - $next_cluster{end};
                        }

                        if (scalar keys %next_cluster == 0) {
                            $next_cluster{start} = $check_start;
                            $next_cluster{end}   = $check_end;
                            $next_cluster{genes}{$check_id} = '';
                            $next_relax_gene++    if $check_id =~ m/MOB\w+/i;
                            $next_coupling_gene++ if $check_id =~ m/(?:T4CP|TcpA)/i;
                            $next_atpase_gene++   if $check_id =~ m/(?:VirB4|TraU)/i;
                        }
                        # Conservative distances for this check, longer distances will be re-evaluated in the next cycle
                        elsif ($check_distance <= 15000) {
                            $next_cluster{end} = $check_end if ($check_end > $next_cluster{end});
                            $next_cluster{genes}{$check_id} = '';
                            $next_relax_gene++    if $check_id =~ m/MOB\w+/i;
                            $next_coupling_gene++ if $check_id =~ m/(?:T4CP|TcpA)/i;
                            $next_atpase_gene++   if $check_id =~ m/(?:VirB4|TraU)/i;
                        }
                        else {
                            last;
                        }
                    }

                    # Guess nexxt gene cluster type from its mating pair formation (mfc) genes
                    my %next_mfp_types;
                    foreach my $hit (keys %{ $next_cluster{genes} }) {
                        my $type = ($hit =~ m/^(\w|FA|FATA)_/)[0] || '';
                        next if $type eq '';

                        $next_mfp_types{$type} = 0 if not exists $next_mfp_types{$type};
                        $next_mfp_types{$type}++;
                    }
                    my $next_highest_value = (scalar keys %next_mfp_types > 0) ? (sort {$b <=> $a} values %next_mfp_types)[0]
                                           :                                      0;
                    # Use of Hash for cases when there is a tie for the top represented types
                    my %next_top_types;
                    foreach my $type (keys %next_mfp_types) {
                        $next_top_types{$type} = '' if ($next_mfp_types{$type} == $next_highest_value);
                    }

                    # Check if gene cluster types match
                    my $type_match = 0;
                    foreach my $type (keys %next_top_types) {
                        $type_match = 1 if (exists $curr_top_types{$type});
                    }

                    # A new cluster must be formed if the original hit found belonged to a new
                    # cluster (at least 3 members), there are repeated genes, or it repeats one
                    # of the key genes (e.g. both clusters have the relaxase)
                    my $repated_key  = ($relax_gene    > 0 and $next_relax_gene    > 0) ? 1
                                     : ($coupling_gene > 0 and $next_coupling_gene > 0) ? 1
                                     : ($atpase_gene   > 0 and $next_atpase_gene   > 0) ? 1
                                     :                                                    0;
                    my $repeated_acc = 0;
                    foreach my $key (keys %{ $next_cluster{genes} } ) {
                        $repeated_acc++ if exists $conjscan_groups{$conjscan_counter}{genes}{$key};
                    }

                    # Expand current cluster when its small and there are no repeated elements found
                    # (e.g. new ICE found in CP001392 genome)
                    if (    scalar keys %{ $conjscan_groups{$conjscan_counter}{genes} } <= 2
                        and $repated_key == 0 and $repeated_acc == 0
                        ) {
                        $conjscan_groups{$conjscan_counter}{end}   = $hit_end if ($hit_end > $conjscan_groups{$conjscan_counter}{end});
                        $conjscan_groups{$conjscan_counter}{genes}{$conjscan_id} = '';
                        $relax_gene++    if $conjscan_id =~ m/MOB\w+/i;
                        $coupling_gene++ if $conjscan_id =~ m/(?:T4CP|TcpA)/i;
                        $atpase_gene++   if $conjscan_id =~ m/(?:VirB4|TraU)/i;
                    }
                    # New cluster formation
                    elsif (   ($type_match == 0 and scalar keys %curr_top_types > 0 and scalar keys %next_top_types > 0)
                           or  $repated_key > 0
                           or $repeated_acc > 0
                        ) {
                        $conjscan_counter++;
                        $conjscan_groups{$conjscan_counter}{start} = $hit_start;
                        $conjscan_groups{$conjscan_counter}{end}   = $hit_end;
                        $conjscan_groups{$conjscan_counter}{genes}{$conjscan_id} = '';
                        $relax_gene    = 0;
                        $coupling_gene = 0;
                        $atpase_gene   = 0;
                        $relax_gene++    if $conjscan_id =~ m/MOB\w+/i;
                        $coupling_gene++ if $conjscan_id =~ m/(?:T4CP|TcpA)/i;
                        $atpase_gene++   if $conjscan_id =~ m/(?:VirB4|TraU)/i;
                    }
                    # Current cluster expansion
                    else {
                        $conjscan_groups{$conjscan_counter}{end}   = $hit_end if ($hit_end > $conjscan_groups{$conjscan_counter}{end});
                        $conjscan_groups{$conjscan_counter}{genes}{$conjscan_id} = '';
                        $relax_gene++    if $conjscan_id =~ m/MOB\w+/i;
                        $coupling_gene++ if $conjscan_id =~ m/(?:T4CP|TcpA)/i;
                        $atpase_gene++   if $conjscan_id =~ m/(?:VirB4|TraU)/i;
                    }
                }
                # New cluster formation
                else {
                    $conjscan_counter++;
                    $conjscan_groups{$conjscan_counter}{start} = $hit_start;
                    $conjscan_groups{$conjscan_counter}{end}   = $hit_end;
                    $conjscan_groups{$conjscan_counter}{genes}{$conjscan_id} = '';
                    $relax_gene    = 0;
                    $coupling_gene = 0;
                    $atpase_gene   = 0;
                    $relax_gene++    if $conjscan_id =~ m/MOB\w+/i;
                    $coupling_gene++ if $conjscan_id =~ m/(?:T4CP|TcpA)/i;
                    $atpase_gene++   if $conjscan_id =~ m/(?:VirB4|TraU)/i;
                }
            }

            # Copy %conjscan_groups data in %results
            foreach my $key (sort {$a <=> $b} keys %conjscan_groups) {
                my @group_gene_keys = keys %{ $conjscan_groups{$key}{genes} };
                if (scalar @group_gene_keys > 1) {
                    my $group_start  = $conjscan_groups{$key}{start};
                    my $group_end    = $conjscan_groups{$key}{end};

                    # From Guglielmini 2011/2014 papers (CONJscan-related). MPF (mating-pair formation)
                    # Relaxases (mobilization): MOBB, MOBC, MOBF, MOBH, MOBP, MOBQ, MOBT, MOBV
                    # Coupling protein T4CP: VirD4 and TcpA
                    # Major ATPases: VirB4 and TraU
                    # MPFT (vir system of the Ti plasmid): VirB 1-11
                    # MPFF (vir system of the F plasmid): Tra LEKBVRCIWUCNEFHGD, Trb C
                    # MPFI (vir system of the R64 Incl plasmid): Tra EIJKLMNOPQRSTUVWXY, Trb AB
                    # MPFG (vir system of the ICE - ICEHIN1056): Tfc 1-24
                    # MPFB (vir system of the Bacteroidetes): Tra EFHIJKLMNOPQ
                    # MPFC (vir system of the Cyanobacteria): Alr 7204-7213
                    # MPF_FA (vir system of the Firmicutes and Actinobacteria): Orf 13-24
                    # MPF_FATA (vir system of the FA + Tenericutes and Archaea): trs CDEFGJLK, prg BUCFGHIJKL, CD 410-424, gbs 1346-1369
                    # NOTE: MPF_FATA above have 4 subclasses, but they are very similar to each other, so there can be mixed hits in an island
                    my $relaxases      = 0; # 10 HMM models in CONJscan
                    my $coupling       = 0; # 3 HMM models in CONJscan
                    my $atpases        = 0; # 2 HMM models in CONJscan
                    my $mpf_T_genes    = 0; # 9 HMM models in CONJscan (not counting VirB4)
                    my $mpf_F_genes    = 0; # 12 HMM models in CONJscan
                    my $mpf_I_genes    = 0; # 16 HMM models in CONJscan
                    my $mpf_G_genes    = 0; # 18 HMM models in CONJscan
                    my $mpf_B_genes    = 0; # 12 HMM models in CONJscan
                    my $mpf_C_genes    = 0; # 8 HMM models in CONJscan
                    my $mpf_FA_genes   = 0; # 7 HMM models in CONJscan
                    my $mpf_FATA_genes = 0; # 27 HMM models in CONJscan between the 4 subfamilies
                    foreach my $gene (@group_gene_keys) {
                        $relaxases++      if $gene =~ m/MOB\w+/i;
                        $coupling++       if $gene =~ m/(?:T4CP|TcpA)/i;
                        $atpases++        if $gene =~ m/(?:VirB4|TraU)/i;
                        $mpf_T_genes++    if $gene =~ m/T_virB[0-3,5-9]+/i; # VirB4 not included
                        $mpf_F_genes++    if $gene =~ m/(?:F_tra\w|F_trb\w)/i;
                        $mpf_I_genes++    if $gene =~ m/(?:I_tra\w|I_trb\w)/i;
                        $mpf_G_genes++    if $gene =~ m/G_tfc\d+/i;
                        $mpf_B_genes++    if $gene =~ m/B_tra\w/i;
                        $mpf_C_genes++    if $gene =~ m/C_alr\d+/i;
                        $mpf_FA_genes++   if $gene =~ m/FA_orf\w+/i;
                        $mpf_FATA_genes++ if $gene =~ m/(?:FATA_cd\w+|FATA_gbs\d+|FATA_prg\w+|FATA_trs\w)/i;
                    }

                    # Guglielmini classification: VirB4-T4CP-MOBP: ICE | VirB4 (+/- T4CP): Putative protein-exporting T4SS | MOBP: IME
                    if (   $atpases        >= 1
                        or $relaxases      >= 1
                        or $mpf_T_genes    >= 2
                        or $mpf_F_genes    >= 2
                        or $mpf_B_genes    >= 2
                        or $mpf_C_genes    >= 2
                        or $mpf_FA_genes   >= 2
                        or $mpf_FATA_genes >= 2
                        or scalar @group_gene_keys >= 3
                        ) {
                        my $loc_index     = $group_start + (1 - $group_end / $sequence_length);
                        my $group_loc_str = "$group_start..$group_end";
                        $results_counter++;
                        $results{$results_counter}{program} = 'CONJscan';
                        $results{$results_counter}{loc_str} = $group_loc_str;
                        $results{$results_counter}{loc_idx} = $loc_index;
                        $results{$results_counter}{details} = 'N/A';
                    }
                }
            }
        }
    }

    # Results merging
    islands_analysis($sequence, $sequence_obj, \%results, \%db_results);

    # Print table with results from all programs
    open my $TABLE, '>', $result_table_file or die "Could not write file '$result_table_file': $!\n";
    if (scalar @result_table == 1) {
        # If there is nothing else than the header, print a No Results line
        my $nothing_line = join("\t", $description, $sequence, $taxonomy, ('No results') x 5) . "\n";
        push @result_table, $nothing_line;
    }
    print $TABLE @result_table;
    close $TABLE;
}

print "\nFinished!\n";
exit;


# Algorithm that decides if an island candidate is conserved or rejected
sub finalize_island {
    my ($islands, $island_buffer, $db_data, $island_counter) = @_;

    # For islands finalizing because of truncation,
    # remove hits that could belong to the next island
    my %removed_members;
    my @removed_members_data;
    if (exists $$island_buffer{truncation}) {
        my %clean_members;
        my @clean_members_data;
        foreach my $check_data ( @{ $$island_buffer{members_data} } ) {
            my $check_program = $$check_data{program};
            my $check_loc_str = $$check_data{loc_str};
            my $current_start = (exists $$island_buffer{hard_start}) ? $$island_buffer{hard_start} : $$island_buffer{start};
            my $current_end   = (exists $$island_buffer{hard_end})   ? $$island_buffer{hard_end}   : $$island_buffer{end};

            my ($check_start, $check_end) = ($check_loc_str =~ m/(\d+)\.\.(\d+)/);
            if (   ($check_start >= $current_start and $check_start <= $current_end)
                or ($check_end   >= $current_start and $check_end   <= $current_end)
                or ($check_start <  $current_start and $check_end   >  $current_end)
                ) {
                # For results that have only one end inside of the island, check if its significant or just a bit.
                # Exclude tRNAs/tmRNAS from this check though, since they are small and important for DNA insertion.
                if ( $check_program !~ m/(?:Aragorn|tRNAscan)/
                    and (   ($check_end   > $current_end   and $check_start >= $current_start and $check_start <= $current_end )
                         or ($check_start < $current_start and $check_end   >= $current_start and $check_end   <= $current_end ) )
                      ) {
                    my $inside_length  = (    $check_start >= $current_start
                                          and $check_start <= $current_end) ? ($current_end - $check_start + 1)
                                       :                                      ($check_end - $current_start + 1);
                    my $inside_percent = $inside_length / ($check_end - $check_start + 1) * 100;
                    if ($inside_percent < 50) {
                        if ($check_program eq 'Phage_Finder') {
                            if (exists $removed_members{Phage_Finder} and exists $$check_data{attach}) {
                                $removed_members{Phage_Finder} = $$check_data{attach};
                            }
                            elsif (not exists $removed_members{Phage_Finder}) {
                                $removed_members{Phage_Finder} = (exists $$check_data{attach}) ? $$check_data{attach}
                                                               : '';
                            }
                        }
                        elsif ($check_program eq 'tRNA_frags') {
                            if (exists $removed_members{tRNA_frags} and exists $$check_data{original_loc}) {
                                $removed_members{tRNA_frags} = $$check_data{loc_str};
                            }
                            elsif (not exists $removed_members{tRNA_frags}) {
                                $removed_members{tRNA_frags} = (exists $$check_data{original_loc}) ? $$check_data{loc_str}
                                                             : '';
                            }
                        }
                        else {
                            $removed_members{$check_program} = '';
                        }
                        push @removed_members_data, $check_data;
                        next;
                    }
                }

                if ($check_program eq 'Phage_Finder') {
                    if (exists $clean_members{Phage_Finder} and exists $$check_data{attach}) {
                        $clean_members{Phage_Finder} = $$check_data{attach};
                    }
                    elsif (not exists $clean_members{Phage_Finder}) {
                        $clean_members{Phage_Finder} = (exists $$check_data{attach}) ? $$check_data{attach}
                                                     : '';
                    }
                }
                elsif ($check_program eq 'tRNA_frags') {
                    if (exists $clean_members{tRNA_frags} and exists $$check_data{original_loc}) {
                        $clean_members{tRNA_frags} = $$check_data{loc_str};
                    }
                    elsif (not exists $clean_members{tRNA_frags}) {
                        $clean_members{tRNA_frags} = (exists $$check_data{original_loc}) ? $$check_data{loc_str}
                                                   : '';
                    }
                }
                else {
                    $clean_members{$check_program} = '';
                }
                push @clean_members_data, $check_data;
            }
            # Hit completely out of the new range
            else {
                if ($check_program eq 'Phage_Finder') {
                    if (exists $removed_members{Phage_Finder} and exists $$check_data{attach}) {
                        $removed_members{Phage_Finder} = $$check_data{attach};
                    }
                    elsif (not exists $removed_members{Phage_Finder}) {
                        $removed_members{Phage_Finder} = (exists $$check_data{attach}) ? $$check_data{attach}
                                                       : '';
                    }
                }
                elsif ($check_program eq 'tRNA_frags') {
                    if (exists $removed_members{tRNA_frags} and exists $$check_data{original_loc}) {
                        $removed_members{tRNA_frags} = $$check_data{loc_str};
                    }
                    elsif (not exists $removed_members{tRNA_frags}) {
                        $removed_members{tRNA_frags} = (exists $$check_data{original_loc}) ? $$check_data{loc_str}
                                                     : '';
                    }
                }
                else {
                    $removed_members{$check_program} = '';
                }
                push @removed_members_data, $check_data;
            }
        }
        $$island_buffer{members} = \%clean_members;
        @{ $$island_buffer{members_data} } = @clean_members_data;
    }

    # Check if the island contains an important gap
    # that could give an island false positive
    my $island_gap_percent = 0;
    my $discard_island     = 0;
    if (exists $$island_buffer{members}{gaps}) {
        $island_gap_percent = gap_percent($$island_buffer{members_data}, $$island_buffer{start}, $$island_buffer{end});

        # Discard islands that contains an unknown length gap
        # because the effect on detection cannot be evaluated,
        # also discard if the length is more than 25% of the island
        if (   $island_gap_percent eq 'unknown'
            or $island_gap_percent >= 25
            ) {
            $discard_island = 1;
        }
    }

    if ($discard_island == 0) {
        # Decide to store the island seed depending on its members
        my $added_flag = 0;

        # For AlienHunter results that are not supported
        # by previous programs (or only by PAIDA)
        if (    (   exists $$island_buffer{members}{AlienHunter}
                 or exists $$island_buffer{members}{IslandPath}
                 or exists $$island_buffer{members}{PAIDA}
                 )
            and not exists $$island_buffer{members}{Phage_Finder}
            and not exists $$island_buffer{members}{PhiSpy}
            and not exists $$island_buffer{members}{CONJscan}
            ) {
            # Rescue results that are supported by a tRNA fragment region,
            # or have a good score without the involvement of an overlapping rRNA operon
            my $highest_ah_ratio    = 0;
            my $highest_paida_score = 0;
            my $rRNA_operon         = 0;
            my $ah_hits             = 0;
            my $paida_hits          = 0;
            my $tRNA_hits           = 0;
            foreach my $ah_data ( @{ $$island_buffer{members_data} } ) {
                if ($$ah_data{program} eq 'AlienHunter') {
                    $ah_hits++;
                    my $ratio = $$ah_data{details};
                    $highest_ah_ratio = $ratio if ($ratio > $highest_ah_ratio);
                    $rRNA_operon   = 1      if (exists $$ah_data{rRNA_operon});
                }
                elsif ($$ah_data{program} eq 'Aragorn') {
                    $tRNA_hits++;
                }
                elsif ($$ah_data{program} eq 'PAIDA') {
                    $paida_hits++;
                    my $paida_score = $$ah_data{details};
                    $highest_paida_score = $paida_score if ($paida_score > $highest_paida_score);
                }
            }

            # If there are 2 or more Alien Hunter or PAI-DA hits,
            # it means that its a long island, so in that case
            # revert the possible rRNA interference
            if ($ah_hits > 1 or $paida_hits > 1) {
                $rRNA_operon = 0;
            }

            # Check if the score is good enough depending on who support the island.
            # Note: IslandPath hits without support of AlienHunter or PAI-DA will not pass.
            my $good_score = 0;
            if (# AlienHunter high cutoff
                   (        exists $$island_buffer{members}{AlienHunter} and $highest_ah_ratio    >= 2.8
                    )
                # PAI-DA high cutoff
                or (        exists $$island_buffer{members}{PAIDA}       and $highest_paida_score >= 3.2
                    )
                # AlienHunter lower cutoff with IslandPath support
                or (        exists $$island_buffer{members}{AlienHunter} and $highest_ah_ratio    >= 1.7
                    and     exists $$island_buffer{members}{IslandPath}
                    and not exists $$island_buffer{members}{PAIDA}
                    )
                # PAI-DA lower cutoff with IslandPath support
                or (        exists $$island_buffer{members}{PAIDA}       and $highest_paida_score >= 2.8
                    and     exists $$island_buffer{members}{IslandPath}
                    and not exists $$island_buffer{members}{AlienHunter}
                    )
                # AlienHunter and PAI-DA combined cutoffs, with and without IslandPath support
                or (    exists $$island_buffer{members}{AlienHunter}
                    and exists $$island_buffer{members}{PAIDA}
                    and (   (    not exists $$island_buffer{members}{IslandPath}
                             and $highest_ah_ratio    >= 1.5 and $highest_paida_score >= 2.6
                             )
                         or (    exists $$island_buffer{members}{IslandPath}
                             and $highest_ah_ratio    >= 1.3 and $highest_paida_score >= 2.5
                             )
                         )
                    )
                ) {
                $good_score = 1;
            }

            # Only store if its supported by a tRNA fragment region,
            # or if it have a high score ratio without a predicted
            # rRNA operon overlapping
            if (   (    exists $$island_buffer{members}{tRNA_frags}
                    and exists $$island_buffer{hard_start}
                    )
                or (    $rRNA_operon == 0 and $tRNA_hits <= 6
                    and $good_score  == 1
                    )
                ) {
                $added_flag = 1;
            }
        }
        elsif (   exists $$island_buffer{members}{Phage_Finder}
               or exists $$island_buffer{members}{PhiSpy}
               or exists $$island_buffer{members}{CONJscan}
            ) {
            $added_flag = 1;
        }

        # If the island was added, calculate its score
        if ($added_flag == 1) {
            $island_counter++;

            $$islands{$island_counter}{start}             = $$island_buffer{start};
            $$islands{$island_counter}{end}               = $$island_buffer{end};
            @{ $$islands{$island_counter}{members_data} } = @{ $$island_buffer{members_data} };
            foreach my $key (keys %{ $$island_buffer{members} })  {
                $$islands{$island_counter}{members}{$key} = '';
            }
            if (exists $$island_buffer{hard_start}) {
                $$islands{$island_counter}{hard_start}    = $$island_buffer{hard_start};
                $$islands{$island_counter}{hard_end}      = $$island_buffer{hard_end};
            }
            $$islands{$island_counter}{truncation}        = '' if exists $$island_buffer{truncation};
            $$islands{$island_counter}{conjscan_attach}   = '' if exists $$island_buffer{conjscan_attach};

            # PAI-DA results give unrefined 5 kilobases windows (1..5000, 10001..15000, etc),
            # so if PAI-DA is the first result of the island, use the first coordinate of the next
            # program prediction if its located within the first 5kb window of the PAI-DA result
            # (except when there are "hard limits" set on the island)
            if (    scalar @{ $$islands{$island_counter}{members_data} } > 1
                and $$islands{$island_counter}{members_data}[0]{program} eq 'PAIDA'
                and not exists $$islands{$island_counter}{hard_start}
                ) {
                my ($paida_start) = ($$islands{$island_counter}{members_data}[0]{loc_str} =~ m/(\d+)\.\./);
                my ($next_start)  = ($$islands{$island_counter}{members_data}[1]{loc_str} =~ m/(\d+)\.\./);

                # When the first element added is a CONJscan gene group, it leaves no trace in "members_data"
                # and PAI-DA get the first array slot, but PAI-DA start will be different from island start.
                # So only continue if both start points match
                if (    $paida_start == $$islands{$island_counter}{start}
                    and $next_start  <= $paida_start + 4999
                    ) {
                    $$islands{$island_counter}{start} = $next_start;
                }
            }

            # Now add any extra gene from %db_data that are located within the island,
            # and collect information from CONJscan hits for classification
            my $added_genes = 0;
            my $relaxases   = 0; # 10 HMM models in CONJscan
            my $coupling    = 0; # 3 HMM models in CONJscan
            my $atpases     = 0; # 2 HMM models in CONJscan
            my %mfp_types;
            my %mfp_hits;

            foreach my $number (sort {
                                      $$db_data{$a}{loc_idx}
                                      <=>
                                      $$db_data{$b}{loc_idx}
                                      }
                                keys %$db_data
                ) {
                my $db_loc_idx = $$db_data{$number}{loc_idx};
                next if (sprintf("%d",$db_loc_idx) < $$islands{$island_counter}{start});

                my $database   = $$db_data{$number}{database};
                my $db_loc_str = $$db_data{$number}{loc_str};
                my $db_details = $$db_data{$number}{details};
                my ($db_start, $db_end) = ($db_loc_str =~ m/(\d+)\.\.(\d+)/);
                my ($db_hit)   = ($db_details =~ m/; Hit: (.+);$/);
                last if ($db_start > $$islands{$island_counter}{end});

                if (   ($db_start >= $$islands{$island_counter}{start} and $db_start <= $$islands{$island_counter}{end})
                    or ($db_end   >= $$islands{$island_counter}{start} and $db_end   <= $$islands{$island_counter}{end})
                    or ($db_start <  $$islands{$island_counter}{start} and $db_end   >  $$islands{$island_counter}{end})
                    ) {
                    $added_genes++;

                    # Update islands coordinates if no "hard limits" are set
                    if (not exists $$islands{$island_counter}{hard_start}) {
                        $$islands{$island_counter}{start} = $db_start if ($db_start < $$islands{$island_counter}{start});
                        $$islands{$island_counter}{end}   = $db_end   if ($db_end   > $$islands{$island_counter}{end});
                    }

                    my %temp_data = (program => $database,
                                     loc_str => $db_loc_str,
                                     loc_idx => $db_loc_idx,
                                     details => $db_details);
                    $$islands{$island_counter}{members}{$database} = '';
                    push @{ $$islands{$island_counter}{members_data} }, \%temp_data;

                    # For ICE/IME classification
                    if ($database eq 'CONJscan') {
                        $relaxases++ if $db_hit =~ m/MOB\w+/i;
                        $coupling++  if $db_hit =~ m/(?:T4CP|TcpA)/i;
                        $atpases++   if $db_hit =~ m/(?:VirB4|TraU)/i;

                        my $type = ($db_hit =~ m/^(\w|FA|FATA)_/)[0] || '';
                        if ($type ne '') {
                            $mfp_types{$type} = 0 if not exists $mfp_types{$type};
                            $mfp_types{$type}++;
                            $mfp_hits{$type}{$db_hit} = '';
                        }
                    }
                }
            }
            if ($added_genes > 0) {
                my $key = 0;
                my %temp_data   = map { ($key++, $_) } @{ $$islands{$island_counter}{members_data} };
                my @sorted_keys = sort { $temp_data{$a}{loc_idx}
                                         <=>
                                         $temp_data{$b}{loc_idx}
                                        }
                                  keys %temp_data;
                my @sorted_values = map { $temp_data{$_} } @sorted_keys;
                @{ $$islands{$island_counter}{members_data} } = @sorted_values;

                # CONJscan classification
                my $highest_value = (scalar keys %mfp_types > 0) ? (sort {$b <=> $a} values %mfp_types)[0]
                                  :                                 0;
                # Use of Hash for cases when there is a tie for the top represented types
                my %top_types;
                my $conjscan_type = '';
                foreach my $type (keys %mfp_types) {
                    $top_types{$type} = '' if ($mfp_types{$type} == $highest_value);
                }
                if (scalar keys %top_types > 0) {
                    if (scalar keys %top_types == 1) {
                        $conjscan_type   = (keys %top_types)[0];
                    }
                    else {
                        $conjscan_type = join '-', sort {$a cmp $b} keys %top_types;
                    }
                }

                # Guglielmini classification: VirB4-T4CP-MOBP: ICE | VirB4 (+/- T4CP): Putative protein-exporting T4SS | MOBP: IME
                my $possible_ice = (    $atpases   > 0
                                    and $coupling  > 0
                                    and $relaxases > 0
                                    )                   ? 'ICE'
                                 : (    $relaxases > 0
                                    )                   ? 'IME'
                                 : (    $atpases    > 0
                                    and $relaxases == 0
                                    )                   ? 'Putative protein-exporting T4SS'
                                 :                        '';
                if ($possible_ice ne '') {
                    $$islands{$island_counter}{ice_type} = ($conjscan_type ne '') ? "$possible_ice MPF_$conjscan_type family" : $possible_ice;
                }
            }

            # Calculate the Score
            my $island_score = island_score($$islands{$island_counter}{members_data}, $island_gap_percent);

            # Revert island addition if it have gaps and the resulting score falls below 3, or if the score
            # is less than 2 (e.g. 1 CONJscan hit separated from the other clustered genes by truncation)
            if (   (    exists $$islands{$island_counter}{members}{gaps}
                    and $island_score < 3
                    )
                or $island_score < 2
                ) {
                delete $$islands{$island_counter};
                # Revert island counter addition
                $island_counter--;
            }
            else {
                $island_score = sprintf("%.2f", $island_score);
                $$islands{$island_counter}{score} = $island_score;
            }
        }
    }

    return ($island_counter, \%removed_members, \@removed_members_data);
}

# Island score
sub island_score {
    my ($members_data, $island_gap_percent) = @_;

    my $score = 0;
    # Score per program
    foreach my $prediction (@{ $members_data }) {
        if ($$prediction{program} eq 'AlienHunter') {
            $score += 1 + $$prediction{details}; # Ratio
        }
        elsif ($$prediction{program} eq 'Aragorn') {
            $score += 1;
        }
        elsif ($$prediction{program} eq 'CONJscan') {
            $score += 1;
        }
        elsif ($$prediction{program} eq 'IslandPath') {
            $score += 4;
        }
        elsif ($$prediction{program} eq 'PAIDA') {
            $score += 1 + ($$prediction{details} * 2); # Score
        }
        elsif ($$prediction{program} eq 'Phage_Finder') {
            $score += 6;
        }
        elsif ($$prediction{program} eq 'PhiSpy') {
            $score += 5;
        }
        elsif ($$prediction{program} eq 'TnpPred') {
            $score += 1;
        }
        elsif ($$prediction{program} eq 'tRNA_frags') {
            # Give a best score to predictions where a possible
            # range could be established with the similar tRNA
            if (exists $$prediction{original_loc}) {
                $score += 3;
            }
            else {
                $score += 1;
            }
        }
        elsif ($$prediction{program} eq 'tRNAscan') {
            $score += 1;
        }
    }
    # Penalization by total gap percent
    $score -= ($island_gap_percent / 4); # -1 per every 4% gap

    return $score;
}

# Island gap percent
sub gap_percent {
    my ($members_data, $start, $end) = @_;
    my $total_gap_len = 0;
    foreach my $gap_data ( @{ $members_data } ) {
        next if $$gap_data{program} ne 'gaps';

        my $gap_length = $$gap_data{details};
        if ($gap_length eq 'unknown') {
            $total_gap_len = 'unknown';
            last;
        }
        else {
            $total_gap_len += $gap_length;
        }
    }

    my $gap_percent;
    if ($total_gap_len eq 'unknown') {
        $gap_percent = 'unknown'
    }
    else {
        $gap_percent = $total_gap_len / ($end - $start + 1) * 100;
    }
    return $gap_percent;
}

# Consensus analysis on islands found
sub islands_analysis {
    my ($header, $seq_obj, $data, $db_data) = @_;
    my $desc = $seq_obj->desc;
    my %data = %$data;

    my %islands;
    my %island_buffer;
    my $island_counter  = 0;
    my @consensus_table = ( "Description\tLocus\tLocation\tPossible Type\tScore\tPrograms\tPrediction Elements\n" );
    my $consensus_gbk   = "$header.Consensus.gbk";
    my $consensus_txt   = "$header.Consensus_table.txt";

    # Check the results sorted by location index
    my $tRNA_frag_size    = 0;
    my $phage_finder_size = 0;
    foreach my $num (sort {
                           ($data{$a}{loc_idx})
                           <=>
                           ($data{$b}{loc_idx})
                           }
                     keys %data
        ) {
        my $program = $data{$num}{program};
        my $loc_str = $data{$num}{loc_str};
        my $details = $data{$num}{details};

        my ($start, $end) = ($loc_str =~ m/(\d+)\.\.(\d+)/);
        next if (not defined $start); # For ocasional single nucleotide predictions (Alien Hunter)

        # To store data from hits removed from an island by truncation
        my ($removed_members, $removed_members_data);

        # If the buffer have an island seed, check if the result
        # can be added to the current or a new island must be initiated
        if (scalar keys %island_buffer != 0) {
            my $truncate_island         = 0;
            my $expand_island           = 0;
            my $island_overlap_percent  = 0;
            my $margin_difference_ratio = 0;
            my $size_difference_ratio   = 0;

            # If it overlaps with the previous island, choose if add it to the previous island
            # or let it continue to be merged with the current island
            if (exists $island_buffer{trunc_start} and $start >= $islands{$island_counter}{start} and $start <= $islands{$island_counter}{end}) {
                my $prev_island_overlap_percent;
                if ($end <= $islands{$island_counter}{end}) {
                    $prev_island_overlap_percent = 100;
                }
                else {
                    $prev_island_overlap_percent = ($islands{$island_counter}{end} - $start + 1) / ($end - $start + 1) * 100;
                }

                if ($prev_island_overlap_percent >= 50) {
                    if ($program eq 'Phage_Finder') {
                        if (exists $islands{$island_counter}{members}{Phage_Finder} and exists $data{$num}{attach}) {
                            $islands{$island_counter}{members}{Phage_Finder} = $data{$num}{attach};
                        }
                        elsif (not exists $islands{$island_counter}{members}{Phage_Finder}) {
                            $islands{$island_counter}{members}{Phage_Finder} = (exists $data{$num}{attach}) ? $data{$num}{attach}
                                                                             : '';
                        }
                    }
                    elsif ($program eq 'tRNA_frags') {
                        if (exists $islands{$island_counter}{members}{tRNA_frags} and exists $data{$num}{original_loc}) {
                            $islands{$island_counter}{members}{tRNA_frags} = $data{$num}{loc_str};
                        }
                        elsif (not exists $islands{$island_counter}{members}{tRNA_frags}) {
                            $islands{$island_counter}{members}{tRNA_frags} = (exists $data{$num}{original_loc}) ? $data{$num}{loc_str}
                                                                           : '';
                        }
                    }
                    else {
                        $islands{$island_counter}{members}{$program} = '';
                    }
                    push @{ $islands{$island_counter}{members_data} }, $data{$num};

                    # Sort again the members with the new addition
                    my $key = 0;
                    my %temp_data   = map { ($key++, $_) } @{ $islands{$island_counter}{members_data} };
                    my @sorted_keys = sort { $temp_data{$a}{loc_idx}
                                             <=>
                                             $temp_data{$b}{loc_idx}
                                            }
                                      keys %temp_data;
                    my @sorted_values = map { $temp_data{$_} } @sorted_keys;
                    @{ $islands{$island_counter}{members_data} } = @sorted_values;

                    ## Pass the data to decide if the island will be kept or rejected
                    #($island_counter, $removed_members, $removed_members_data) = finalize_island(\%islands, \%island_buffer, $db_data, $island_counter);

                    # Check again if the island contains an important gap that could give
                    # an island false positive (which could be hit just added)
                    my $island_gap_percent = 0;
                    if (exists $island_buffer{members}{gaps}) {
                        my $island_gap_percent = gap_percent($island_buffer{members_data}, $island_buffer{start}, $island_buffer{end});

                        # Discard islands that contains an unknown length gap
                        # because the effect on detection cannot be evaluated,
                        # also discard if the length is more than 25% of the island
                        if (   $island_gap_percent eq 'unknown'
                            or $island_gap_percent >= 25
                            ) {
                            delete $islands{$island_counter};
                            $island_counter--;
                            next;
                        }
                    }

                    # Recalculate the Score
                    my $island_score = island_score($islands{$island_counter}{members_data}, $island_gap_percent);

                    # Revert island addition if it have gaps and the resulting score falls below 3, or if the score
                    # is less than 2 (e.g. 1 CONJscan hit separated from the other clustered genes by truncation)
                    if (   (    exists $islands{$island_counter}{members}{gaps}
                            and $island_score < 3
                            )
                        or $island_score < 2
                        ) {
                        delete $islands{$island_counter};
                        $island_counter--;
                    }
                    else {
                        $island_score = sprintf("%.2f", $island_score);
                        $islands{$island_counter}{score} = $island_score;
                    }
                    next;
                }
            }

            # If the result start coordinate is located within the island region, add it and expand
            # the island region. Make a special treatment for Phage Finder with att data and tRNA_frags;
            # since they establish "hard limits", check if the island in current form is not lost by
            # a small overlap with the new limits.
            if (   ($start >= $island_buffer{start}    and $start <= $island_buffer{end})
                # For islands with truncated start, include those hits that were not taken by previous island
                or (exists $island_buffer{trunc_start} and $start < $island_buffer{start})
                ) {
                if ($end <= $island_buffer{end}) {
                    $island_overlap_percent = 100;
                }
                else {
                    $island_overlap_percent = ($island_buffer{end} - $start + 1) / ($end - $start + 1) * 100;
                }
                $margin_difference_ratio = ($start - $island_buffer{start} + 1) / ($end - $start + 1);
                $size_difference_ratio   = ($island_buffer{end} - $island_buffer{start} + 1) / ($end - $start + 1);

                if (   (    $program eq 'Phage_Finder' and $details =~ m/att[LR]/
                        and exists $island_buffer{members}{PhiSpy}
                        )
                    or (    (   ($program eq 'Phage_Finder' and $details =~ m/att[LR]/)
                             or ($program eq 'tRNA_frags'   and exists $data{$num}{original_loc})
                             # ICEs overlapping putative Phages
                             or ($program eq 'CONJscan'     and exists $island_buffer{hard_start})
                             # Avoids unreasonable fusion of 2 PhiSpy predictions
                             or ($program eq 'PhiSpy'       and exists $island_buffer{members}{PhiSpy})
                             # Avoids excesive lost of a PhiSpy prediction by tRNA_frags when overlap is small
                             or ($program eq 'PhiSpy'       and not exists $island_buffer{members}{PhiSpy} and exists $island_buffer{hard_start})
                             )
                        and (# Overlap with small islands
                                ($size_difference_ratio <= 0.5 and $island_overlap_percent >= 20)
                             # Overlap with medium-big islands
                             or ($size_difference_ratio  > 0.5 and $island_overlap_percent >= 40 and $margin_difference_ratio <= 1.5)
                             )
                        )
                    or (    (   ($program eq 'tRNA_frags' and     exists $data{$num}{original_loc}       and exists $island_buffer{hard_start})
                             or ($program eq 'CONJscan'   and     exists $island_buffer{hard_start})
                             or ($program eq 'PhiSpy'     and     exists $island_buffer{members}{PhiSpy})
                             or ($program eq 'PhiSpy'     and not exists $island_buffer{members}{PhiSpy} and exists $island_buffer{hard_start})
                             )
                        # When overlap is significant, add to the island regardless of size/margin conditions above
                        and $island_overlap_percent >= 60
                    )
                    or (   (   ($program eq 'Phage_Finder' and $details !~ m/att[LR]/)
                            or ($program eq 'tRNA_frags'   and not exists $data{$num}{original_loc})
                            or ($program eq 'CONJscan'     and not exists $island_buffer{hard_start})
                            or ($program eq 'PhiSpy'       and not exists $island_buffer{members}{PhiSpy} and not exists $island_buffer{hard_start})
                            or  $program !~ m/(?:CONJscan|Phage_Finder|PhiSpy|tRNA_frags)/
                            )
                        and $island_overlap_percent >  0
                        )
                ) {
                    $expand_island = 1;
                }
                else {
                    $truncate_island = 1;
                }
            }
            # Following conditions allows isolated T4SS gene groups to attach to adjacent islands
            elsif (   (   $program eq 'CONJscan'
                       or (exists $island_buffer{members}{CONJscan} and keys %{ $island_buffer{members} } == 1)
                       )
                   and ($start > $island_buffer{end} and $start <= $island_buffer{end} + 1000)
                ) {
                $expand_island = 1;
                $island_buffer{conjscan_attach} = '';
            }

            if ($expand_island == 1) {
                my $changed_start = 0;
                if ($program eq 'tRNA_frags' and exists $data{$num}{original_loc}) {
                    $tRNA_frag_size = $end - $start + 1;
                }
                elsif ($program eq 'Phage_Finder') {
                    $phage_finder_size = $end - $start + 1;
                }

                # Look for the presence of "hard linits" and update if necessary
                if (   (    not exists $island_buffer{hard_start}
                        and (   (    $program eq 'tRNA_frags'
                                 and exists $data{$num}{original_loc}
                                 and $tRNA_frag_size >= $phage_finder_size)
                             or (    $program eq 'Phage_Finder'
                                 and $details =~ m/att[LR]/
                                 and not exists $island_buffer{members}{Phage_Finder}
                                 and (    not exists $island_buffer{members}{PhiSpy}
                                      or (    exists $island_buffer{members}{PhiSpy}
                                          and (# Overlap with small islands
                                                  ($size_difference_ratio <= 0.5 and $island_overlap_percent >= 20)
                                               # Overlap with medium-big islands
                                               or ($size_difference_ratio  > 0.5 and $island_overlap_percent >= 40 and $margin_difference_ratio <= 1.5)
                                               )
                                          )
                                      )
                                 )
                             )
                        )
                    # If limits were set by tRNA_frags, allow an update
                    # if Phage Finder have new data
                    or (    exists $island_buffer{hard_start}
                        and $island_buffer{hard_program} eq 'tRNA_frags'
                        and (    $program eq 'Phage_Finder'
                             and $details =~ m/att[LR]/)
                        )
                    ) {
                    $island_buffer{hard_program} = $program;
                    $island_buffer{hard_start}   = $start;
                    $island_buffer{hard_end}     = $end;
                    $island_buffer{hard_length}  = ($end - $start + 1);
                    if ($island_buffer{start} != $start) {
                        $changed_start = 1;
                    }
                    $island_buffer{start}        = $start;
                    $island_buffer{end}          = $end;
                }
                # Found cases where the tRNA_frags truncated good Phage_Finder
                # predictions because it lacked the att[LR] data, so allow
                # Phage_Finder to undo "hard limits" in some cases
                elsif (    exists $island_buffer{hard_start}
                       and $island_buffer{hard_program} eq 'tRNA_frags'
                       and (    $program eq 'Phage_Finder'
                            and $details !~ m/att[LR]/
                            and $island_buffer{hard_length} < $phage_finder_size)
                    ) {
                    delete $island_buffer{hard_program};
                    delete $island_buffer{hard_start};
                    delete $island_buffer{hard_end};
                    delete $island_buffer{hard_length};

                    # Backtrack $start to the first start found and look for the furthest $end
                    # among the island seed members, since they could have been altered by tRNA_frags
                    (my $new_start) = ($island_buffer{members_data}[0]{loc_str} =~ m/(\d+)\.\./);
                    if ($island_buffer{start} != $new_start) {
                        $changed_start = 1;
                    }
                    my $new_end = 0;
                    foreach my $member ( @{ $island_buffer{members_data} } ) {
                        (my $member_end) = $$member{loc_str} =~ m/\.\.(\d+)/;
                        $new_end = $member_end if ($member_end > $new_end);
                    }

                    $island_buffer{start} = $new_start;
                    $island_buffer{end}   = $new_end;
                }

                # If the island Start changed, remove the results
                # that could be now out of range
                if ($changed_start == 1) {
                    my %clean_members;
                    my @clean_members_data;
                    foreach my $check_data ( @{ $island_buffer{members_data} } ) {
                        my $check_program = $$check_data{program};
                        my $check_loc_str = $$check_data{loc_str};
                        my $current_start = (exists $island_buffer{hard_start}) ? $island_buffer{hard_start} : $island_buffer{start};
                        my $current_end   = (exists $island_buffer{hard_end})   ? $island_buffer{hard_end}   : $island_buffer{end};

                        my ($check_start, $check_end) = ($check_loc_str =~ m/(\d+)\.\.(\d+)/);
                        if (   ($check_start >= $current_start and $check_start <= $current_end)
                            or ($check_end   >= $current_start and $check_end   <= $current_end)
                            or ($check_start <  $current_start and $check_end   >  $current_end)
                            ) {
                            # For results that have only one end inside of the island, check if its significant or just a bit.
                            # Exclude tRNAs/tmRNAS from this check though, since they are small and important for DNA insertion.
                            if ( $check_program !~ m/(?:Aragorn|tRNAscan)/
                                and (   ($check_end   > $current_end   and $check_start >= $current_start and $check_start <= $current_end )
                                     or ($check_start < $current_start and $check_end   >= $current_start and $check_end   <= $current_end ) )
                                  ) {
                                my $inside_length  = (    $check_start >= $current_start
                                                      and $check_start <= $current_end) ? ($current_end - $check_start + 1)
                                                   :                                      ($check_end - $current_start + 1);
                                my $inside_percent = $inside_length / ($check_end - $check_start + 1) * 100;
                                next if $inside_percent < 40;
                            }

                            if ($check_program eq 'Phage_Finder') {
                                if (exists $clean_members{Phage_Finder} and exists $$check_data{attach}) {
                                    $clean_members{Phage_Finder} = $$check_data{attach};
                                }
                                elsif (not exists $clean_members{Phage_Finder}) {
                                    $clean_members{Phage_Finder} = (exists $$check_data{attach}) ? $$check_data{attach}
                                                                 : '';
                                }
                            }
                            elsif ($check_program eq 'tRNA_frags') {
                                if (exists $clean_members{tRNA_frags} and exists $$check_data{original_loc}) {
                                    $clean_members{tRNA_frags} = $$check_data{loc_str};
                                }
                                elsif (not exists $clean_members{tRNA_frags}) {
                                    $clean_members{tRNA_frags} = (exists $$check_data{original_loc}) ? $$check_data{loc_str}
                                                               : '';
                                }
                            }
                            else {
                                $clean_members{$check_program} = '';
                            }
                            push @clean_members_data, $check_data;
                        }
                    }
                    $island_buffer{members} = \%clean_members;
                    @{ $island_buffer{members_data} } = @clean_members_data;
                }

                # Add the new data to the buffer, but for CONJscan gene groups don't add 'members_data'
                # since the individual CONJscan genes will be added later on the %db_data check
                if ($program eq 'Phage_Finder') {
                    if (exists $island_buffer{members}{Phage_Finder} and exists $data{$num}{attach}) {
                        $island_buffer{members}{Phage_Finder} = $data{$num}{attach};
                    }
                    elsif (not exists $island_buffer{members}{Phage_Finder}) {
                        $island_buffer{members}{Phage_Finder} = (exists $data{$num}{attach}) ? $data{$num}{attach}
                                                              : '';
                    }
                }
                elsif ($program eq 'tRNA_frags') {
                    if (exists $island_buffer{members}{tRNA_frags} and exists $data{$num}{original_loc}) {
                        $island_buffer{members}{tRNA_frags} = $data{$num}{loc_str};
                    }
                    elsif (not exists $island_buffer{members}{tRNA_frags}) {
                        $island_buffer{members}{tRNA_frags} = (exists $data{$num}{original_loc}) ? $data{$num}{loc_str}
                                                            : '';
                    }
                }
                else {
                    $island_buffer{members}{$program} = '';
                }

                if ($program ne 'CONJscan') {
                    push @{ $island_buffer{members_data} }, $data{$num};
                }
                if ($end > $island_buffer{end} and not exists $island_buffer{hard_end}) {
                    $island_buffer{end} = $end;
                }
            }
            # If there is no overlap, close the current island and clean the buffer
            else {
                # If the island was truncated as consequence of a tRNA_frag region or
                # Phage Finder prediction with att sites present in the current or next island,
                # set the end depending on who have the "hard limits"
                if ($truncate_island == 1) {
                    # Special scenario: a CONJscan region is getting split by truncation caused by
                    # "hard limits" on current island, but in this case there are both a Phage Finder
                    # prediction with att sites and a tRNA_frag region with different sizes; in that
                    # case try to use the region that leave it out to avoid splitting the CONJscan region
                    if (    $program eq 'CONJscan'
                        and exists $island_buffer{hard_start}
                        and exists $island_buffer{members}{tRNA_frags}   and $island_buffer{members}{tRNA_frags}   ne ''
                        and exists $island_buffer{members}{Phage_Finder} and $island_buffer{members}{Phage_Finder} ne ''
                        ) {
                        my ($tRNA_frags_end)   = ($island_buffer{members}{tRNA_frags}   =~ m/\.\.(\d+)/);
                        my ($Phage_Finder_end) = ($island_buffer{members}{Phage_Finder} =~ m/\.\.(\d+)/);

                        my $chosen = '';
                        if ($start > $tRNA_frags_end and $start <= $Phage_Finder_end) {
                            $chosen = 'tRNA_frags';
                        }
                        elsif ($start > $Phage_Finder_end and $start <= $tRNA_frags_end) {
                            $chosen = 'Phage_Finder';
                        }

                        # If somebody was selected, modify the buffer to use the new hard coordinates,
                        # if there is no good candidate then just truncate next island start
                        if ($chosen ne '') {
                            my ($new_start, $new_end) = ($island_buffer{members}{$chosen} =~ m/(\d+)\.\.(\d+)/);
                            $start = $new_end + 1;

                            $island_buffer{hard_program} = $chosen;
                            $island_buffer{start}        = $new_start;
                            $island_buffer{hard_start}   = $new_start;
                            $island_buffer{end}          = $new_end;
                            $island_buffer{hard_end}     = $new_end;
                            $island_buffer{hard_length}  = ($new_end - $new_start + 1);
                            $island_buffer{truncation}   = '';
                        }
                        else {
                            $start = $island_buffer{hard_end} + 1;
                            $data{$num}{truncation} = '';
                        }
                    }
                    # Only current island, so truncate next island start
                    elsif (    exists $island_buffer{hard_start}
                        and not (   $program eq 'tRNA_frags'   and exists $data{$num}{original_loc}
                                 or $program eq 'Phage_Finder' and $details =~ m/att[LR]/)
                        ) {
                        $start = $island_buffer{hard_end} + 1;
                        $data{$num}{truncation} = '';
                    }
                    # Only next island, so truncate current island end
                    elsif (    not exists $island_buffer{hard_start}
                           and (   $program eq 'tRNA_frags'   and exists $data{$num}{original_loc}
                                or $program eq 'Phage_Finder' and $details =~ m/att[LR]/)
                        ) {
                        $island_buffer{end} = $start - 1;
                        $island_buffer{truncation} = '';
                    }
                    # Both or none have (so no clear preference), truncate current island end
                    else {
                        $island_buffer{end} = $start - 1;
                        $island_buffer{truncation} = '';
                    }
                }

                # Pass the data to decide if the island will be kept or rejected,
                # and recover possible removed hits from truncation.
                # Sometimes the first island causes truncation but later is rejected,
                # making the second island unnecessarily lose information, so monitor
                # the acceptance/rejection to revert if needed.
                my $counter_monitor = $island_counter;
                ($island_counter, $removed_members, $removed_members_data) = finalize_island(\%islands, \%island_buffer, $db_data, $island_counter);

                if (exists $data{$num}{truncation} and $counter_monitor == $island_counter) {
                    delete $data{$num}{truncation};
                    ($start) = ($data{$num}{loc_str} =~ m/(\d+)\.\./);
                }

                # Clear buffer for next island
                %island_buffer     = ();
                $tRNA_frag_size    = 0;
                $phage_finder_size = 0;
            }
        }

        # If buffer is empty, initialize the new possible island
        # and go to the next result on the list
        if (scalar keys %island_buffer == 0) {
            # First include data that could have been removed by
            # previous island and therefore is located before current hit
            my $added_hits = 0;
            if (scalar keys %$removed_members > 0) {
                $added_hits++;
                $island_buffer{members}           = $removed_members;
                @{ $island_buffer{members_data} } = @$removed_members_data;
            }

            # Now continue with current hit
            $island_buffer{start}             = $start;
            $island_buffer{end}               = $end;
            if ($program eq 'Phage_Finder') {
                if (exists $island_buffer{members}{Phage_Finder} and exists $data{$num}{attach}) {
                    $island_buffer{members}{Phage_Finder} = $data{$num}{attach};
                }
                elsif (not exists $island_buffer{members}{Phage_Finder}) {
                    $island_buffer{members}{Phage_Finder} = (exists $data{$num}{attach}) ? $data{$num}{attach}
                                                          : '';
                }
            }
            elsif ($program eq 'tRNA_frags') {
                if (exists $island_buffer{members}{tRNA_frags} and exists $data{$num}{original_loc}) {
                    $island_buffer{members}{tRNA_frags} = $data{$num}{loc_str};
                }
                elsif (not exists $island_buffer{members}{tRNA_frags}) {
                    $island_buffer{members}{tRNA_frags} = (exists $data{$num}{original_loc}) ? $data{$num}{loc_str}
                                                        : '';
                }
            }
            else {
                $island_buffer{members}{$program} = '';
            }

            # For island seeds starting from truncation
            if (exists $data{$num}{truncation}) {
                delete $data{$num}{truncation};
                $island_buffer{truncation} = '';
                # Next key purpose is to signal next hits that they must
                # check if they fall in the previous island or this one
                $island_buffer{trunc_start} = '';
            }

            # For CONJscan gene groups don't add 'members_data' since the individual
            # CONJscan genes will be added later on the %db_data check
            if ($program eq 'CONJscan' and $added_hits == 0) {
                @{ $island_buffer{members_data} } = ();
            }
            elsif ($program ne 'CONJscan' and $added_hits == 0) {
                @{ $island_buffer{members_data} } = ( $data{$num} );
            }
            elsif ($program ne 'CONJscan' and $added_hits  > 0) {
                push @{ $island_buffer{members_data} }, $data{$num};
            }

            # If the result is a tRNA fragment with region coords (associated with
            # its predicted tRNA source), or a Phage Finder result with established
            # insertion points (attL/R), set "hard limits" to the island
            if (   (    $program eq 'tRNA_frags'
                    and exists $data{$num}{original_loc})
                or (    $program eq 'Phage_Finder'
                    and $details =~ m/att[LR]/)
                ) {
                $island_buffer{hard_program} = $program;
                $island_buffer{hard_start}   = $start;
                $island_buffer{hard_end}     = $end;
                $island_buffer{hard_length}  = ($end - $start + 1);
            }

            if ($program eq 'tRNA_frags' and exists $data{$num}{original_loc}) {
                $tRNA_frag_size = $end - $start + 1;
            }
            elsif ($program eq 'Phage_Finder') {
                $phage_finder_size = $end - $start + 1;
            }

            next;
        }
    }
    # Process the last buffer
    if (scalar keys %island_buffer != 0) {
        # Pass the data to decide if the island will be kept or rejected
        ($island_counter) = finalize_island(\%islands, \%island_buffer, $db_data, $island_counter);

        %island_buffer     = ();
        $tRNA_frag_size    = 0;
        $phage_finder_size = 0;
    }

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;
    foreach my $key (sort {$a <=> $b} keys %islands) {
        my $island_start = $islands{$key}{start};
        my $island_end   = $islands{$key}{end};
        my $island_score = $islands{$key}{score};
        my @island_data  = @{ $islands{$key}{members_data} };

        my $note = "Consensus Score: $island_score;";
        my $prediction_counter = 0;
        my $possible_type      = 'Genomic island';
        foreach my $prediction (@island_data) {
            $prediction_counter++;
            if ($$prediction{program} eq 'AlienHunter') {
                $note .= " $prediction_counter) Alien Hunter: $$prediction{loc_str}, Score Ratio: $$prediction{details};";
            }
            elsif ($$prediction{program} eq 'Aragorn') {
                $note .= " $prediction_counter) Aragorn: $$prediction{loc_str}, Product: $$prediction{details};";
            }
            elsif ($$prediction{program} eq 'CONJscan') {
                my $details        = $$prediction{details};
                my ($conjscan_hit) = ($details =~ m/; Hit: (.+);$/);
                $note .= " $prediction_counter) CONJscan: $$prediction{loc_str}, Product: $conjscan_hit;";
            }
            elsif ($$prediction{program} eq 'gaps') {
                $note .= " $prediction_counter) Ns or Assembly gap: $$prediction{loc_str}, Length: $$prediction{details} nts;";
            }
            elsif ($$prediction{program} eq 'IslandPath') {
                $note .= " $prediction_counter) IslandPath-DIMOB: $$prediction{loc_str};";
            }
            elsif ($$prediction{program} eq 'PAIDA') {
                $note .= " $prediction_counter) PAI-DA: $$prediction{loc_str}, Score: $$prediction{details};";
            }
            elsif ($$prediction{program} eq 'Phage_Finder') {
                $possible_type = 'Prophage';
                my $details    = $$prediction{details};
                $details       =~ s/;/,/g;
                $note .= " $prediction_counter) Phage Finder: $$prediction{loc_str}, Note: $details;";
            }
            elsif ($$prediction{program} eq 'PhiSpy') {
                my $details    = $$prediction{details};
                $details       =~ s/;/,/g;
                $note .= " $prediction_counter) PhiSpy: $$prediction{loc_str}, Note: $details;";
            }
            elsif ($$prediction{program} eq 'TnpPred') {
                my $details       = $$prediction{details};
                my ($tnppred_hit) = ($details =~ m/; Hit: (.+);$/);
                $note .= " $prediction_counter) TnpPred: $$prediction{loc_str}, Product: $tnppred_hit;";
            }
            elsif ($$prediction{program} eq 'tRNA_frags') {
                my ($type, $details) = ($$prediction{details} =~ m/(tm?RNA) fragment: ([^\.]+)/);
                if (exists $$prediction{original_loc}) {
                    $note .= " $prediction_counter) $type fragment: $$prediction{original_loc}, From: $details, Predicted Range: $$prediction{loc_str};";
                }
                else {
                    $note .= " $prediction_counter) $type fragment: $$prediction{loc_str}, From: $details;";
                }
            }
            elsif ($$prediction{program} eq 'tRNAscan') {
                $note .= " $prediction_counter) tRNAscan-SE: $$prediction{loc_str}, Product: $$prediction{details};";
            }
        }

        # Re-classification by CONJscan data
        if (exists $islands{$key}{ice_type}) {
            if ($possible_type eq 'Genomic island') {
                $possible_type = $islands{$key}{ice_type};
            }
            else {
                $possible_type = "$possible_type - $islands{$key}{ice_type}";
            }
        }


        # Add violet color grading to feature depending on their score,
        # currently >= 30 is the score for maximum saturation
        my $value = ($island_score - 3) / 30;
        $value    = 1 if ($value > 1);
        $value    = 0 if ($value < 0);

        my $red   = 255;
        my $green = int(255 - ($value * 255));
        my $blue  = 255;

        # Build feature with the information collected
        my $feature = Bio::SeqFeature::Generic
            ->new( -primary_tag => 'misc_feature',
                   -start       => $island_start,
                   -end         => $island_end,
                   -tag         => { inference => 'prediction consensus',
                                     colour    => "$red $green $blue",
                                     note      => "$possible_type; $note" },
        );
        $out_obj->add_SeqFeature($feature);

        my $programs = join ", ", sort{$a cmp $b} keys %{ $islands{$key}{members} };
        push @consensus_table, "$desc\t$header\t$island_start..$island_end\t$possible_type\t$island_score\t$programs\t$note\n";

        print "CHECK TRUNCATION HERE!!!\n"          if (exists $islands{$key}{truncation});
        print "CHECK new CONJscan attach HERE!!!\n" if (exists $islands{$key}{conjscan_attach});
        print "Coords: $island_start..$island_end; Score: $island_score; Possible Type: $possible_type\n";
    }
    my $gbk_obj = Bio::SeqIO->new(-file   => ">$consensus_gbk",
                                  -format => 'genbank');
    $gbk_obj->write_seq($out_obj);

    open my $CONS_TABLE, '>', $consensus_txt or die "Could not write file '$consensus_txt': $!\n";
    if (scalar @consensus_table == 1) {
        # If there is nothing else than the header, print a No Results line
        my $nothing_line = join("\t", $desc, $header, ('No results') x 5) . "\n";
        push @consensus_table, $nothing_line;
    }
    print $CONS_TABLE @consensus_table;
    close $CONS_TABLE;

    print "Islands predicted: $island_counter\n";
    return;
}

# Prediction of tRNA, mtRNA, and tmRNA gene FRAGMENTS
sub database_search {
    my ($database, $database_file, $file, $seq_obj) = @_;
    if (not -e $file) {
        warn "Could not find the provided query file for $database: '$file'\n";
        return;
    }
    my $seq_length = $seq_obj->length;

    # Accessory file
    my ($out_header) = ($file =~ m/^(.+)\.(?:faa|fna)$/i);
    my $file_ptt     = "$out_header.ptt";
    if ($file =~ m/\.faa$/i and not -e $file_ptt) {
        warn "Could not find needed file for $database: '$file_ptt'\n";
        return;
    }

    my %genes_loc;
    if ($file =~ m/\.faa$/i) {
        open my $PTT_TABLE, '<', $file_ptt or die "Could not read file '$file_ptt': $!\n";
        while (my $line = <$PTT_TABLE>) {
            chomp $line;
            next if ($line !~ m/^\d+\.\.\d+\t/);
            next if ($line eq '');

            my @columns  = split /\t/, $line;
            my $location = ($columns[1] eq '-') ? "complement($columns[0])" : $columns[0];
            my ($gene_start, $gene_end) = ($location =~ m/(\d+)\.\.(\d+)/);
            my $loc_idx = $gene_start + (1 - $gene_end / $seq_length);
            $genes_loc{$columns[3]}{$location} = $loc_idx;
        }
        close $PTT_TABLE;
    }

    # Run Blast using the predicted tRNAs as query and the original sequence as target
    my $temp_file = "$out_header.$database.txt";
    my $format    = '';
    my $null      = ($^O =~ m/mswin/i) ? 'NUL' : '/dev/null';
    if (   $database eq 'CONJscan'
        or $database eq 'TnpPred'
        ) {
        $format = 'hmmer';
        if (not -e $temp_file or -z $temp_file) {
            system("hmmscan --cpu 20 -o \"$temp_file\" \"$database_file\" \"$file\" ");
        }
    }

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Feature collector
    my %out_feats;

    # Process output
    my $search_obj = Bio::SearchIO->new(-file   => "<$temp_file",
                                        -format => $format);
    my $number = 0;
    my $flag   = 0;
  RESULT:
    while (my $result = $search_obj->next_result) {
        # faa display_id examples:
        # >gi|16082468|ref|NP_393478.1|           # RefSeq format (ptt ID: 16082468)
        # >gi|339832456|gb|EGQ60373.1|            # GenBank format (ptt ID: 339832456)
        # >gi|00001|lcl|Contig33|aorf_00001|M3|   # ORFminer custom format (ptt ID: 00001)
        # >gi|12092751209412|ref|1209275_1209412| # Old format conversion for custom seqs without any ID (ptt ID: 12092751209412)
        # >gi|5501083|550_1083|                   # New format conversion for custom seqs without any ID (ptt ID: 5501083)
        my $qry_name   = $result->query_name;
        my $qry_desc   = $result->query_description;
        my $qry_length = $result->query_length;

        # Hmmseach format 3.1b1 added a final "[ok]" line after the last "//" line.
        # This makes current BioPerl version to get confused and generate an empty result,
        # so skip result if it have empty $qry_name and $qry_desc
        next if ($qry_name eq '' and $qry_desc eq '');

        # Note: For protein searches there may be more than 1 location in %genes_loc for
        # the same ID (like WP_xxxx ID proteins), but for nucleotide searches it will be empty
        my $gene_loc    = '';
        my $gene_start  = 0;
        my $gene_end    = 0;
        my $gene_strand = 0;
        if (scalar keys %genes_loc > 0) {
            # Considering that there can be non-unique protein ids, always assume that
            # the current hit is the one with the smaller Start and remove it from %genes_loc
            my ($ptt_id) = ($qry_name =~ m/^gi\|(\d+)\|/);
            ($gene_loc)  = sort {$genes_loc{$ptt_id}{$a} <=> $genes_loc{$ptt_id}{$b}} keys %{ $genes_loc{$ptt_id} };
            $gene_strand = ($gene_loc =~ m/complement/) ? -1 : 1;

            ($gene_start, $gene_end) = ($gene_loc =~ m/(\d+)\.\.(\d+)/);
            delete $genes_loc{$ptt_id}{$gene_loc};

            # Skip rare cut-by-origin genes (e.g. 1934779..110)
            next RESULT if ($gene_start > $gene_end);
        }

        while (my $hit = $result->next_hit) {
            my $hit_name = $hit->name;
            my $hit_desc = $hit->description;

            while (my $hsp = $hit->next_hsp) {
                $number++;
                my $score      = ($format eq 'hmmer') ? $hsp->score : $hsp->bits;
                my $evalue     = $hsp->evalue;
                my $qry_start  = $hsp->start('query');
                my $qry_end    = $hsp->end('query');
                my $qry_strand = $hsp->strand('query');
                my $hit_start  = $hsp->start('hit');
                my $hit_end    = $hsp->end('hit');
                my $hit_strand = $hsp->strand('hit');
                my $identity   = $hsp->percent_identity;
                my $similarity = $hsp->frac_conserved * 100;
                my $align_len  = $qry_end - $qry_start + 1;

                my ($qry_cov, $hit_cov);
                if ($format eq 'blast') {
                    $qry_cov    = $hsp->length('query') / $qry_length  * 100;
                    $hit_cov    = $hsp->length('hit')   / $hit->length * 100;
                }
                elsif ($format eq 'hmmer') {
                    $qry_cov    = ($qry_end - $qry_start + 1) / $qry_length * 100;
                    if ($hit->length != 0) {
                        $hit_cov = $hsp->length('hit') / $hit->length * 100;
                    }
                    else {
                        $hit_cov = $hsp->length('hit') / $hsp->hit->end * 100;
                    }
                }
                $identity   = sprintf("%.1f", $identity);
                $similarity = sprintf("%.1f", $similarity);
                $qry_cov    = sprintf("%.1f", $qry_cov);
                $hit_cov    = sprintf("%.1f", $hit_cov);

                if (    $evalue     <= 1e-10
                    and $qry_cov    >= 20
                    and $hit_cov    >= 20
                    and $similarity >= 40
                    ) {
                    # Set hit location
                    my $loc_idx = '';
                    my $loc_str = '';
                    if ($gene_loc ne '') {
                        # Calculate genomic HSP location by converting
                        # and adjusting aminoacids to nucleotides number
                        my ($hsp_start, $hsp_end);
                        if ($gene_strand == 1) {
                            $hsp_start = $gene_start + ( $qry_start * 3 ) - 3;
                            $hsp_end   = $gene_start + ( $qry_end   * 3 ) - 1;
                        }
                        elsif ($gene_strand == -1) {
                            $hsp_start = $gene_end   - ( $qry_end   * 3 ) + 1;
                            $hsp_end   = $gene_end   - ( $qry_start * 3 ) + 3;
                        }
                        $loc_idx = $hsp_start + (1 - $hsp_end / $seq_length);
                        $loc_str = ($gene_strand == -1) ? "complement($hsp_start..$hsp_end)" : "$hsp_start..$hsp_end";
                    }
                    else {
                        $loc_idx = $qry_start + (1 - $qry_end / $seq_length);
                        $loc_str = ($qry_strand == -1) ? "complement($qry_start..$qry_end)" : "$qry_start..$qry_end";
                    }

                    # Set product and description for $note
                    my $description = ($database =~ m/(?:CONJscan|TnpPred)/) ? $hit_name : "$hit_name $hit_desc";
                    my $product     = ($database eq 'TnpPred') ? "$hit_name transposase"
                                    :                            "$hit_name protein";

                    # BLAST/HMMER statistics
                    my $comparison_stats = ($file =~ m/\.faa$/i) ? "Score: $score; E-value: $evalue; Identity: $identity%; Similarity: $similarity%; Qry Coverage: $qry_cov%; Hit Coverage: $hit_cov%"
                                         :                         "Score: $score; E-value: $evalue; Identity: $identity%; Hit Coverage: $hit_cov%; Length: $align_len";

                    my $note    = ($file =~ m/\.faa$/i) ? "$database hit; Orig Annot: $qry_desc; $comparison_stats; Hit: $description;"
                                :                         "$database hit; $comparison_stats; Hit: $description;";
                    my $colour  = ($database eq 'CONJscan')      ? '150 150 150'
                                : ($database eq 'TnpPred')       ? '255 150 0'
                                :                                  '255 255 255';
                    my $loc_obj = Bio::Factory::FTLocationFactory->from_string($loc_str);
                    my $feature = Bio::SeqFeature::Generic->new(-location    => $loc_obj,
                                                                -primary_tag => 'misc_feature',
                                                                -tag => {note    => $note,
                                                                         product => $product,
                                                                         colour  => $colour,},
                                                                );

                    # Store the finished feature
                    $out_feats{$number}{feature} = $feature;
                    $out_feats{$number}{loc_idx} = $loc_idx;

                    my $feat_start = $feature->start;
                    my $feat_end   = $feature->end;
                }

                # For protein searches, only keep the top hit for each query
                next RESULT if ($file =~ m/\.faa$/i);
            }
        }
    }
    close $search_obj->_fh;

    # Add features sorted by a location index
    foreach my $id_num (sort {
                              ($out_feats{$a}{loc_idx})
                              <=>
                              ($out_feats{$b}{loc_idx})
                              }
                        keys %out_feats
        ) {
        my $feature = $out_feats{$id_num}{feature};
        $out_obj->add_SeqFeature($feature);
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.$database.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;
    }

    # Now that all data is in the Genbank file, delete some intermediary files
    unlink $temp_file or die "Could not delete '$temp_file': $!\n";
    return;
}

# Prediction of prophage regions by PhiSpy
sub phispy {
    my ($folder, $file, $seq_obj) = @_;
    if (not -e $folder) {
        warn "Could not find the provided folder for PhiSpy: '$folder'\n";
        return;
    }
    if (not -e $file) {
        warn "Could not find the provided file for PhiSpy: '$file'\n";
        return;
    }

    my $desc       = $seq_obj->desc;
    my $seq_length = $seq_obj->length;

    # Temporary data folders that will be created by PhiSpy on execution,
    # needed to be set before adding the full path to $file
    my $input_data_folder  = $file . "_input";
    my $output_data_folder = $file . "_output";
    $input_data_folder  =~ s/[\.\s]/_/g;
    $output_data_folder =~ s/[\.\s]/_/g;

    # Keep previous working directory to return to it after program execution
    my $old_folder  = getcwd;
    if ($file !~ m'[\\/]') {
        $file = "$old_folder/$file";
    }
    # Change to program directory
    chdir $folder or die "Could not change directory to '$folder': $!\n";

    # Define output
    my $out_header = ($file =~ m/^(.+)\.gbk?$/i) ? $1 : $file;
    my $out_file   = "$out_header.PhiSpy.prophage_tbl.txt";

    # Training sets list available from PhiSpy, choose one if description matches
    # or leave the default generic set (0)
    my $training_set = 0;
    my %training_genomes = (
        'Bacillus halodurans'        => 1,
        'Bacillus subtilis'          => 2,
        'Bifidobacterium longum'     => 3,
        'Brucella melitensis'        => 4,
        'Caulobacter crescentus'     => 5,
        'Clostridium perfringens'    => 6,
        'Clostridium tetani'         => 7,
        'Deinococcus radiodurans'    => 8,
        'Escherichia coli'           => 9,
        'Haemophilus influenzae'     => 10,
        'Lactococcus lactis'         => 11,
        'Listeria innocua'           => 12,
        'Listeria monocytogenes'     => 13,
        'Mesorhizobium loti'         => 14,
        'Mycobacterium tuberculosis' => 15,
        'Neisseria meningitidis'     => 16,
        'Pasteurella multocida'      => 17,
        'Pseudomonas aeruginosa'     => 18,
        'Pseudomonas putida'         => 19,
        'Ralstonia solanacearum'     => 20,
        'Salmonella enterica'        => 21,
        'Shewanella oneidensis'      => 22,
        'Shigella flexneri'          => 23,
        'Staphylococcus aureus'      => 24,
        'Streptococcus pyogenes'     => 25,
        'Streptomyces coelicolor'    => 26,
        'Vibrio cholerae'            => 27,
        'Xanthomonas axonopodis'     => 28,
        'Xylella fastidiosa'         => 29,
        'Yersinia pestis'            => 30,
    );
    foreach my $set (keys %training_genomes) {
        if ($desc =~ m/^$set/) {
            $training_set = $training_genomes{$set};
            last;
        }
    }

    # Make sure that $input_data_folder don't exists to avoid errors
    # with the genbank_to_seed.py script
    if (-d "$folder/$input_data_folder") {
        File::Path->remove_tree("$folder/$input_data_folder");
    }

    # First run the script for converting the Genbank file into the format used by PhiSpy,
    # and the run PhiSpy on the output folder
    my $null = ($^O =~ m/mswin/i) ? 'NUL' : '/dev/null';
    system("python genbank_to_seed.py \"$file\" \"$input_data_folder\" > $null");
    system("python phiSpy.py -i \"$input_data_folder\" -o \"$output_data_folder\" -t $training_set > $null");

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Feature collector
    my %out_feats;

    # Read result table
    my $results_file = "$folder/$output_data_folder/prophage_tbl.txt";
    open my $RESULTS, '<', $results_file or die "Could not read file '$results_file': $!\n";
    my $feat_counter = 0;
    while (my $line = <$RESULTS>) {
        chomp $line;
        next if $line =~ m/^fig_no/;
        next if $line eq '';

        my ($internal_id, $function, $contig,         $start,       $end,
            $position,    $rank,     $initial_status, $probability, $final_status,
            $attL_start,  $attL_end, $attR_start,     $attR_end,    $attL_seq,
            $attR_seq
            ) = split /\t/, $line;
        next if $final_status == 0;

        $feat_counter++;

        # Feature Note
        my $note     = "prophage";
        # This variable will be appended to $note and will
        # contain att data if it spans the whole region
        my $att_note = '';

        # Process attL/R sequences
        if (defined $attL_seq and defined $attR_seq) {
            $attL_seq  = lc $attL_seq;
            $attR_seq  = lc $attR_seq;
            $att_note .= "; attL location: $attL_start..$attL_end, attR location: $attR_start..$attR_end. ";

            # Revert one sequence if attL and attR have opposite direction (+/- or -/+),
            # or both if both are in the reverse strand
            if ($attL_start < $attL_end and $attR_start > $attR_end) {    # +/-
                $attR_seq = reverse $attR_seq;
                $attR_seq =~ tr/atcg/tagc/;
            }
            elsif ($attL_start > $attL_end and $attR_start < $attR_end) { # -/+
                $attL_seq = reverse $attL_seq;
                $attL_seq =~ tr/atcg/tagc/;
            }
            elsif ($attL_start > $attL_end and $attR_start > $attR_end) { # -/-
                $attR_seq = reverse $attR_seq;
                $attL_seq = reverse $attL_seq;
                $attR_seq =~ tr/atcg/tagc/;
                $attL_seq =~ tr/atcg/tagc/;
            }

            my $att_string = '';
            if ($attL_seq eq $attR_seq) {
                $att_string = "attL/R=$attL_seq.";
            }
            else {
                if (length $attL_seq == length $attR_seq) {
                    my @attL_nts = split //, $attL_seq;
                    my @attR_nts = split //, $attR_seq;
                    my $consensus_seq = '';

                    for (my $i = 0; $i < scalar @attL_nts; $i++) {
                        $consensus_seq .= ($attL_nts[$i] eq $attR_nts[$i]) ? $attL_nts[$i]
                                        :                                  "[$attL_nts[$i]$attR_nts[$i]]";
                    }
                    $att_string = "attL/R=$consensus_seq.";
                }
                else {
                    $att_string = "attL=$attL_seq; attR=$attR_seq.";
                }
            }

            # Beware of strings longer than the Genbank format (57 characters)
            if (length $att_string > 57) {
                my @segments = ($att_string =~ m/(.{1,57})/g);
                $att_string = "@segments";
            }

            # Add attachment information to $att_note
            $att_note .= $att_string;
        }

        if ($start > $end) {
            my $new_start = $end;
            my $new_end   = $start;
            $start = $new_start;
            $end   = $new_end;
        }

        # Phage region coordinates
        my $phage_start = 0;
        my $phage_end   = 0;
        if (defined $attL_start) {
            $phage_start = ($attL_start < $attL_end) ? $attL_start : $attL_end;
            $phage_end   = ($attR_start < $attR_end) ? $attR_end   : $attR_start;
        }
        else {
            $phage_start = $start;
            $phage_end   = $end;
        }

        # Keep track of how many CDS are involved
        my $CDS_counter = 1;

        # Fast forward to find phage end coordinate, but keep
        # an end copy in case an att data is discarded
        my $furthest_end = $end;
        while (my $next_line = <$RESULTS>) {
            chomp $next_line;

            my @columns = split /\t/, $next_line;

            # Keep record of the End while final_status is still 1.
            # $columns[3] => start, $columns[4] => end
            if ($columns[9] == 1) {
                $CDS_counter++;
                $end          = ($columns[3] < $columns[4]) ? $columns[4] : $columns[3];
                $furthest_end = $end if ($end > $furthest_end);

                # Only update $phage_end if there is no known attR
                $phage_end = $end if (not defined $attR_start and $end > $phage_end);

                # Remove att data if it does not span the whole region
                $att_note = '' if (scalar @columns == 10);
            }
            else {
                last;
            }
        }
        $note .= $att_note;
        # If the region started with att data but it didn't last until the end,
        # discard att data and use the region genes start/end
        if (defined $attR_start and $att_note eq '') {
            $phage_start = $start;
            $phage_end   = $furthest_end;
        }

        my $loc_idx = $phage_start + (1 - $phage_end / $seq_length);

        # Build feature with the information collected
        my $feature = Bio::SeqFeature::Generic
            ->new( -primary_tag => 'misc_feature',
                   -start       => $phage_start,
                   -end         => $phage_end,
                   -tag         => { inference => 'profile:PhiSpy:2.3',
                                     note      => $note,
                                     colour    => '200 100 100',},
        );

        $out_feats{$feat_counter}{loc_idx}     = $loc_idx;
        $out_feats{$feat_counter}{feature}     = $feature;
        $out_feats{$feat_counter}{CDS_counter} = $CDS_counter;
        # Keep $start/$end without attL/attR modification for summary file validation
        $out_feats{$feat_counter}{start}       = $start;
        $out_feats{$feat_counter}{end}         = $furthest_end;
    }
    close $RESULTS;

    # Confirm if the features data is consistent with the summary table prophage.tbl
    my $summary_counter = 0;
    my $summary_file    = "$folder/$output_data_folder/prophage.tbl";
    open my $SUMMARY, '<', $summary_file or die "Could not read file '$summary_file': $!\n";
    while (my $line = <$SUMMARY>) {
        chomp $line;
        next if $line eq '';

        $summary_counter++;
        my ($internal_id, $location) = split /\t/, $line;

        # Note: PhiSpy internally checks which coordinate is bigger to always print Start < End
        my ($start, $end) = ($location =~ m/_(\d+)_(\d+)$/);

        my $phage_start = $out_feats{$summary_counter}{start};
        my $phage_end   = $out_feats{$summary_counter}{end};

        if ($phage_start != $start or $phage_end != $end) {
            die "  ERROR!!! Result from prophage.tbl ($start..$end) does not match data from table ($phage_start..$phage_end)! Please check!\n";
        }
    }
    close $SUMMARY;

    # Add features sorted by a location index
    foreach my $key ( sort { $out_feats{$a}{loc_idx} <=> $out_feats{$b}{loc_idx} }
                      keys %out_feats
        ) {
        # Don't keep single CDS results, which are a direct result of
        # PhiSpy refining step, where a candidate phage area can be reduced
        # and may leave a fragment of the original unrefined area
        my $CDS_counter = $out_feats{$key}{CDS_counter};
        next if $CDS_counter == 1;

        my $out_feat = $out_feats{$key}{feature};
        $out_obj->add_SeqFeature($out_feat);
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.PhiSpy.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;

        # Keep a copy of the result table
        copy($results_file, $out_file);
    }

    # Now that all data is in the Genbank file, delete the temporary folders
    File::Path->remove_tree("$folder/$input_data_folder");
    File::Path->remove_tree("$folder/$output_data_folder");

    # Return to previous working folder
    chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
    return;
}

# Prediction of prophage regions by Phage Finder
sub phage_finder {
    my ($folder, $file, $seq_obj, $blast_pkg, $codontable) = @_;
    if (not -e $folder) {
        warn "Could not find the provided folder for Phage_Finder: '$folder'\n";
        return;
    }
    if (not -e $file) {
        warn "Could not find the provided file for Phage_Finder: '$file'\n";
        return;
    }

    my $seq_length = $seq_obj->length;

    # Keep previous working directory to return to it after program execution
    my $old_folder = getcwd;
    if ($file !~ m'[\\/]') {
        $file = "$old_folder/$file";
    }

    # Accessory files needed for Phage_Finder
    my $out_header = ($file =~ m/^(.+)\.(?:fa|fna|fasta)$/i) ? $1 : $file;
    my $file_faa   = "$out_header.faa";
    my $file_ptt   = "$out_header.ptt";
    if (not -e $file_faa) {
        warn "Could not find needed file for Phage_Finder: '$file_faa'\n";
        return;
    }
    if (not -e $file_ptt) {
        warn "Could not find needed file for Phage_Finder: '$file_ptt'\n";
        return;
    }

    # Define temporary output files (they must be located in
    # the same folder and not use paths, only simple filenames)
    my $file_header = ($file =~ m#( [^\\/]+ ) \.\w+ $#ix)[0];
    my $temp_folder = "$old_folder/$file_header-Phage_Finder_process";
    my $hmm_file    = "temp_file.hmm_out.txt";
    my $blast_file  = "temp_file.blast_out.txt";
    my $tRNA_file   = "temp_file.tRNAscan_out.txt";
    my $tmRNA_file  = "temp_file.aragorn_out.txt";

    # Phage Finder use some very long filenames on temporary files, so when running
    # on Windows (which have a very narrow limit compared to other OSes) try to
    # estimate if those filenames can stay below the maximum local path length
    # (1+2+256+1 or [drive][:][path][null] = 260)
    if ($^O =~ m/mswin/i) {
        # Example (297): D:\very_long_temp_folder\
        #                CLST_final_pseudochr_dir\
        #                CLST_final_pseudochr.4.phtest.10000.5000.1e-006.phage.tRNA-Pro-2.r2.10.tail.seq.nsq
        my $long_dir        = "$temp_folder/$file_header" . "_dir/";
        my $max_path_length = length($long_dir) + length($file_header) + 70;
        my $length_diff     = $max_path_length - 260;

        if ($max_path_length >= 260) {
            die "Current '$file_header' path is so long that is possible that Phage Finder temporary files "
              . "can't be created because of Windows 260 chars length path limit, please reduce path length "
              . "by at least $length_diff chars by using shorter file/folder names or moving the sequence "
              . "data to a less nested folder!\n";
        }
    }

    # Run accessory programs to generate input files...
    my $null = ($^O =~ m/mswin/i) ? 'NUL' : '/dev/null';
    if (not -e $temp_folder) {
        mkdir $temp_folder or die "Could not create temp folder '$temp_folder': $!\n";
    }

    # Copy the original sequence files that will be used by Phage Finder to the temporary folder
    my ($local_file)     = ($file     =~ m#([^\\/]+)$#);
    my ($local_file_ptt) = ($file_ptt =~ m#([^\\/]+)$#);
    copy($file, "$temp_folder/$local_file")         if (not -e "$temp_folder/$local_file");
    copy($file_ptt, "$temp_folder/$local_file_ptt") if (not -e "$temp_folder/$local_file_ptt");

    # Produce HMM3 results. Mimic HMM3_searches.sh
    if (not -e "$temp_folder/$hmm_file") {
        my $hmmsearch   = `hmmsearch -h`;
        if ($hmmsearch !~ m/# HMMER 3/) {
            warn "Phage_Finder only works with data from HMMER3!\n";
            return;
        }

        my $hmm_file_db = "$temp_folder/HMM3_db.hmm";
        open my $HMM_DB, '>', $hmm_file_db or die "Could not write file '$hmm_file_db': $!\n";

        open my $HMM_LIST, '<', "$folder/hmm3.lst" or die "Could not read file '$folder/hmm3.lst': $!\n";
        my @hmm_list = map { (my $tmp = $_) =~ s/[\s\n]//g; $tmp; } <$HMM_LIST>;
        close $HMM_LIST;

        foreach my $hmm (@hmm_list) {
            my $hmm_path = "$folder/PHAGE_HMM3s_dir/$hmm.HMM";
            if (not -e $hmm_path) {
                warn "  WARNING!!! Could not find file '$hmm_path' for the HMM3 DB!\n";
                next;
            }
            open my $HMM_IN, '<', $hmm_path or die "Could not read file '$hmm_path': $!\n";
            while (my $line = <$HMM_IN>) {
                print $HMM_DB $line;
            }
            close $HMM_IN;
        }
        close $HMM_DB;

        # Run HMMer3
        system("hmmsearch --cpu 20 -o \"$temp_folder/$hmm_file\" \"$hmm_file_db\" \"$file_faa\" ");

        # Remove temporary database
        unlink $hmm_file_db or die "Could not delete file '$hmm_file_db': $!\n";
    }

    # Produce BLAST results (choose between blastall and blastp)
    if (not -e "$temp_folder/$blast_file") {
        if ($blast_pkg eq 'blast+') {
            system(  "blastp -db \"$folder/DB/phage_10_02_07_release.db_new\" -outfmt 6 "
                   . "-query \"$file_faa\" -out \"$temp_folder/$blast_file\" -evalue 0.001 "
                   . "-num_alignments 4 -num_threads 20 -seg no 2> $null" );
        }
        else {
            system(  "blastall -p blastp -d \"$folder/DB/phage_10_02_07_release.db\" -m 8 "
                   . "-i \"$file_faa\" -o \"$temp_folder/$blast_file\" -e 0.001 "
                   . "-v 4 -b 4 -a 20 -F F 2> $null" );
        }
    }
    # Keep a list of the proteins that gave phage hits, so when Phage Finder fails to report
    # a candidate species its possible to discriminate if its an internal bug or there are no hits
    my %phage_hits;
    open my $PHAGE_BLAST, '<', "$temp_folder/$blast_file" or die "Could not read file '$blast_file': $!\n";
    while (my $line = <$PHAGE_BLAST>) {
        chomp $line;
        next if ($line eq '');

        my ($gi_id) = ($line =~ m/^gi\|(\d+)\|/);
        $phage_hits{$gi_id}{hit} = '';
    }
    close $PHAGE_BLAST;
    if (scalar keys %phage_hits > 0) {
        my $seq_counter = 0;
        open my $SEQS, '<', $file_faa or die "Could not read file '$file_faa': $!\n";
        while (my $line = <$SEQS>) {
            if ($line =~ m/^>gi\|(\d+)\|/) {
                my $gi_id = $1;
                $seq_counter++;
                $phage_hits{$gi_id}{pos} = $seq_counter;
            }
        }
        close $SEQS;
    }

    # Produce tRNAscan-SE results
    if (not -e "$temp_folder/$tRNA_file") {
        chdir $trna_scan_se_dir or die "Could not change directory to '$trna_scan_se_dir': $!\n";
        system(  "perl /home/bioinfo/mobhunter/MobHunter_v1/MobHunter/tRNAscan-SE_1.3.2/tRNAscan-SE.pl "
               . "--quiet --forceow --general -o \"$temp_folder/$tRNA_file\" \"$file\" ");
    }

    # Produce Aragorn tmRNA results
    if (not -e "$temp_folder/$tmRNA_file") {
        chdir $aragorn_dir or die "Could not change directory to '$aragorn_dir': $!\n";

        # Linux needs './' before the executable binary name
        my $aragorn_exe    = ($^O =~ m/mswin/i) ? 'aragorn_alpha3' : './aragorn_alpha3';

        # Get version, since option -wunix is added in v1.2.37
        my ($version)      = (`$aragorn_exe -h` =~ m/ARAGORN v([\.\d]+) /);
        my ($v1, $v2, $v3) = ($version =~ m/\d+/g);
        my $ver_number     = ($v1*100000) + ($v2*1000) + $v3;
        # NOTE: -wunix option causes massive slowdown for files with a
        # big ammount of sequences. Only use for files with 1 to few sequences.
        my $newline_fix    = ($^O =~ m/mswin/i and $ver_number >= 102037) ? '-wunix' : '';

        # If the -wunix option is unavailable, is very important to
        # create temporary files with native line-endings to avoid Aragorn's
        # coordinates bug
        my $input_file = "$temp_folder/$local_file";
        if ($newline_fix eq '') {
            # Create temporary file
            my $temp_file   = "$temp_folder/$local_file.copy.fna";

            open my $IN,  '<', $input_file or die "Could not read file '$input_file': $!\n";
            open my $OUT, '>', $temp_file   or die "Could not write file '$temp_file': $!\n";

            while (my $line = <$IN>) {
                chomp $line;
                print $OUT "$line\n";
            }
            close $IN;
            close $OUT;
            $input_file = $temp_file;
        }
        system("$aragorn_exe $newline_fix -m -o \"$temp_folder/$tmRNA_file\" \"$input_file\" ");

        # If a temporary file was used, delete it now
        if ($input_file eq "$temp_folder/$local_file.copy.fna") {
            unlink $input_file or die "Could not delete temporary file '$input_file': $!\n";
        }
    }

    # Change to Phage_Finder script folder (so it can compile correctly)
    chdir "$folder/bin" or die "Could not change to '$folder/bin': $!\n";

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Phage Finder will work only if the Blast output or PTT file are not empty
    if (-s "$temp_folder/$blast_file" and -s "$temp_folder/$local_file_ptt") {
        # Run program with all the data and remove some leftovers
        system(  "perl Phage_Finder_v2.1.pl -b \"$temp_folder\" -t \"$blast_file\" "
               . "-i \"$local_file_ptt\" -m \"$hmm_file\" -r \"$tRNA_file\" -n \"$tmRNA_file\" "
               . "-A \"$local_file\" -c $codontable -S > $null" );

        if (-e "$folder/bin/error.log") {
            unlink "$folder/bin/error.log" or die "Could not delete '$folder/bin/error.log': $!\n";
        }
        if (-e "$old_folder/error.log") {
            unlink "$old_folder/error.log" or die "Could not delete '$old_folder/error.log': $!\n";
        }
        if (-e "$folder/bin/formatdb.log") {
            unlink "$folder/bin/formatdb.log" or die "Could not delete '$folder/bin/formatdb.log': $!\n";
        }

        # Feature collector
        my %out_feats;

        # Phage_Finder internal database, read to convert later
        # genome IDs into species names
        my %species;
        my $db_file = "$folder/DB/phage_10_02_07_release.crib";
        open my $DB, '<', $db_file or die "Could not read file '$db_file': $!\n";
        while (my $line = <$DB>) {
            chomp $line;
            next if $line eq '';

            my @data = split /\t/, $line;
            $species{$data[0]} = $data[2];
        }
        close $DB;

        # Phage_Finder have a bug where the predicted prophage region
        # coordinates are sometimes wrong when the begin or end gene is located
        # in the reverse strand, so read "phage_finder_info.txt" to get
        # the original genes coordinates and apply a fix later if needed
        my %gene_coords;
        my $info_file = "$temp_folder/phage_finder_info.txt";
        open my $INFO, '<', $info_file or die "Could not read file '$info_file': $!\n";
        while (my $line = <$INFO>) {
            chomp $line;
            next if $line eq '';

            my ($assemble_id, $size, $acc, $a_pos, $b_pos, $product) = split /\t/, $line;
            # The coordinates are reversed for genes in the reverse strand
            if ($a_pos < $b_pos) {
                $gene_coords{$acc}{end_5}  = $a_pos;
                $gene_coords{$acc}{end_3}  = $b_pos;
                $gene_coords{$acc}{strand} = '+';
            }
            else {
                $gene_coords{$acc}{end_5} = $b_pos;
                $gene_coords{$acc}{end_3} = $a_pos;
                $gene_coords{$acc}{strand} = '-';
            }
        }
        close $INFO;

        # Read result table
        my $results_file = "$temp_folder/PFPR_tab.txt";
        open my $RESULTS, '<', $results_file or die "Could not read file '$results_file': $!\n";
        my $feat_counter = 0;
        while (my $line = <$RESULTS>) {
            chomp $line;
            next if $line =~ m/^#asmbl_id/;
            next if $line eq '';

            $feat_counter++;
            my ($assemble_id,   $genome_size,     $genome_gc,       $begin_region,
                $end_region,    $size_region,     $label,           $type,
                $prime5_att,    $prime3_att,      $target,          $region_gc,
                $best_db_match, $begin_gene,      $end_gene,        $integrase_HMMs,
                $core_HMMs,     $over_noise_HMMs, $lytic_HMMs,      $tail_HMMs,
                $Mu_HMMs,       $region_orient,   $dist_int_to_att, $num_genes,
                $ser_recombs
                ) = split /\t/, $line;

            # Check for bugs in the region coordinates. Note: When begin/end gene is a tRNA,
            # it will not appear in phage_finder_info.txt, so no record in %gene_coords
            if (exists $gene_coords{$begin_gene}) {
                if ($prime5_att eq 'N.D.' and $begin_region == $gene_coords{$begin_gene}{end_3}) {
                    $begin_region = $gene_coords{$begin_gene}{end_5};
                }
            }
            if (exists $gene_coords{$end_gene}) {
                if ($prime3_att eq 'N.D.' and $end_region == $gene_coords{$end_gene}{end_5}) {
                    $end_region = $gene_coords{$end_gene}{end_3};
                }
            }
            # For reverse hits, sometimes Phage Finder has printed "1..10" and others "10..1",
            # so sort them first just in case
            my $region_start = ($begin_region > $end_region) ? $end_region   : $begin_region;
            my $region_end   = ($begin_region > $end_region) ? $begin_region : $end_region;
            my $location     = ($region_orient eq '-') ? "complement($region_start..$region_end)"
                             :                           "$region_start..$region_end";
            my $loc_idx      = $begin_region + (1 - $end_region / $seq_length);

            # Build tags
            my $species  = $species{$best_db_match};
            if (not defined $species) {
                # Check if there was a blast hit for the any of the involved genes
                # but it was not reported by Phage Finder because a bug in its code,
                # if nothing happens then there is really no hit
                my @genes     = sort {$phage_hits{$a}{pos} <=> $phage_hits{$b}{pos}} keys %phage_hits;
                my $evaluate  = 0;
                foreach my $gene (@genes) {
                    $evaluate = 1 if ($gene eq $begin_gene);
                    if ($evaluate == 1) {
                        if (exists $phage_hits{$gene}{hit}) {
                            die "  There was a species blast hit in Phage_Finder DB that didn't appear in results at $location! Check bug!\n";
                        }
                    }
                    $evaluate = 0 if ($gene eq $end_gene);
                }
                $species  = 'unknown species';
            }
            my $note     = ($species eq 'unknown species') ? "$type; $label; no species hit in Phage_Finder DB"
                         :                                   "$type; $label; best Phage_Finder DB hit is $species ($best_db_match)";
            if ($target ne 'N.D.') {
                $note .= "; Insertion into ";
                $note .= ($target =~ m/^\d+$/) ? "PID $target. " : "$target. ";

                my $att_string = '';
                if ($prime5_att eq $prime3_att) {
                    $att_string = "attL/R=" . lc($prime5_att) . '.';
                }
                else {
                    if (length $prime5_att == length $prime3_att) {
                        my @att5_nts = split //, lc $prime5_att;
                        my @att3_nts = split //, lc $prime3_att;
                        my $consensus_seq = '';

                        for (my $i = 0; $i < scalar @att5_nts; $i++) {
                            $consensus_seq .= ($att5_nts[$i] eq $att3_nts[$i]) ? $att5_nts[$i]
                                            :                                  "[$att5_nts[$i]$att3_nts[$i]]";
                        }
                        $att_string = "attL/R=$consensus_seq.";
                    }
                    else {
                        $att_string = "attL=" . lc($prime5_att) . "; attR=" . lc($prime3_att). '.';
                    }
                }

                # Beware of strings longer than the Genbank format (57 characters)
                if (length $att_string > 57) {
                    my @segments = ($att_string =~ m/(.{1,57})/g);
                    $att_string = "@segments";
                }

                # Add attachment information to $note
                $note .= $att_string;
            }

            # Build feature with the information collected
            my $feat_loc_obj = Bio::Factory::FTLocationFactory->from_string($location);
            my $feature = Bio::SeqFeature::Generic
                ->new(-primary_tag => 'misc_feature',
                      -location    => $feat_loc_obj,
                      -tag => { inference           => 'profile:Phage_Finder:2.1',
                                note                => $note,
                                colour              => '100 200 100',
                               },
                     );
            $out_feats{$feat_counter}{loc_idx} = $loc_idx;
            $out_feats{$feat_counter}{feature} = $feature;
            $out_feats{$feat_counter}{note}    = $note;
        }
        close $RESULTS;

        # Add features sorted by a location index. NOTE: There is at least 1 case in FN433596
        # where an unresolved bug in an island results in 2 overlaped predictions,
        # so check adjacent predictions to reject one of the overlaped regions
        my @sorted_keys = sort { $out_feats{$a}{loc_idx} <=> $out_feats{$b}{loc_idx} }
                          keys %out_feats;
        for (my $i = 0; $i < scalar @sorted_keys; $i++) {
            my $key        = $sorted_keys[$i];
            my $out_feat   = $out_feats{$key}{feature};
            my $curr_start = $out_feat->start;
            my $curr_end   = $out_feat->end;

            if ($i < $#sorted_keys) {
                my $next_key   = $sorted_keys[$i+1];
                my $next_feat  = $out_feats{$next_key}{feature};
                my $next_start = $next_feat->start;
                my $next_end   = $next_feat->end;

                # Check if the prediction overlaps with the next one
                if ($curr_end >= $next_start and $curr_end <= $next_end) {
                    my $curr_note = $out_feats{$key}{note};
                    my $next_note = $out_feats{$next_key}{note};
                    my $curr_len  = $curr_end - $curr_start + 1;
                    my $next_len  = $next_end - $next_start + 1;

                    # If one prediction have att data, remove the other;
                    # if none or both have att data, keep the bigger one
                    if ($curr_note !~ m/attL/ and $next_note =~ m/attL/) {
                        $out_feats{$key}{remove} = '';
                    }
                    elsif ($curr_note =~ m/attL/ and $next_note !~ m/attL/) {
                        $out_feats{$next_key}{remove} = '';
                    }
                    else {
                        if ($curr_len <= $next_len) {
                            $out_feats{$key}{remove} = '';
                        }
                        else {
                            $out_feats{$next_key}{remove} = '';
                        }
                    }
                }
            }

            # Add feature if it has not been marked for removal
            if (not exists $out_feats{$key}{remove}) {
                $out_obj->add_SeqFeature($out_feat);
            }
        }
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.Phage_Finder.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;
    }

    # Now that all data is in the Genbank file, delete the temporary folder
    File::Path->remove_tree($temp_folder);

    # Return to previous working folder
    chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
    return;
}

# Prediction of tRNA genes
sub trna_scan {
    my ($folder, $file, $seq_obj) = @_;
    if (not -e $folder) {
        warn "Could not find the provided folder for tRNAscan-SE: '$folder'\n";
        return;
    }
    if (not -e $file) {
        warn "Could not find the provided file for tRNAscan-SE: '$file'\n";
        return;
    }

    my $seq_length = $seq_obj->length;

    # Keep previous working directory to return to it after program execution
    my $old_folder = getcwd;
    if ($file !~ m'[\\/]') {
        $file = "$old_folder/$file";
    }
    # Change to program directory
    chdir $folder or die "Could not change directory to '$folder': $!\n";

    # Define output
    my $out_header = ($file =~ m/^(.+)\.(?:fa|fna|fasta)$/i) ? $1 : $file;
    my $out_file   = "$out_header.struct.txt";

    # Run program
    system( "perl /home/bioinfo/mobhunter/MobHunter_v1/MobHunter/tRNAscan-SE_1.3.2/tRNAscan-SE.pl "
           . "--quiet --forceow --general --struct \"$out_file\" "
           . "-o \"$out_header.txt\" \"$file\" " );

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Feature collector
    my %out_feats;

    # Process output
    open my $RESULT, '<', $out_file or die "Could not read file '$out_file': $!\n";
    my $number = 0;

    # NOTE: tRNAscan-SE does not make predictions over circular sequences
    # (with sequences passing through the origin: Start-SEQ_END -> 1-End),
    # also possible introns are forced to locate outside of the anticodon,
    # so coordinates processing is simpler than Aragorn
    while (my $line = <$RESULT>) {
        chomp $line;
        # Result start
        if ($line =~ m/^\S+.trna(\d+) \s+ \( (\d+) - (\d+) \) \s+ Length:/x) {
            $number   = $1;
            my $a_pos = $2;
            my $b_pos = $3;

            if (exists $out_feats{$number}) {
                die "Do not process more than 1 query sequence per run!\n";
            }

            # Prepare coordinates
            my $start  = 0;
            my $end    = 0;
            my $strand = 0;
            if ($a_pos > $b_pos) {
                $start  = $b_pos;
                $end    = $a_pos;
                $strand = -1;
            }
            else {
                $start  = $a_pos;
                $end    = $b_pos;
                $strand = 1;
            }
            my $location = ($strand == -1) ? "complement($start..$end)"
                         :                   "$start..$end";
            my $loc_idx  = $start + (1 - $end / $seq_length);

            # Prepare the hash keys
            $out_feats{$number}{coords}        = $location;
            $out_feats{$number}{loc_idx}       = $loc_idx;
            $out_feats{$number}{tRNA}          = '';
            $out_feats{$number}{score}         = 0;
            $out_feats{$number}{anticodon}     = '';
            $out_feats{$number}{anti_coords  } = '';
            $out_feats{$number}{intron_coords} = '';
            $out_feats{$number}{intron_seq}    = '';
            $out_feats{$number}{pseudogene}    = 0;
        }
        # tRNA data
        elsif ($line =~ m/^Type: (\S+)	Anticodon: (\S+) at \d+-\d+ \((\d+)-(\d+)\)	Score: ([\.\d]+)/) {
            my $tRNA                       = $1;
            $out_feats{$number}{anticodon} = $2;
            my $anti_a_pos                 = $3;
            my $anti_b_pos                 = $4;
            $out_feats{$number}{score}     = $5;

            # For Selenocysteine tRNA
            $tRNA =~ s#SeC (?: \( [pe] \) )?#Sec#x; #SeC(p), SeC(e), SeC
            $out_feats{$number}{tRNA} = $tRNA;

            # Prepare coordinates
            my $anti_codon_start  = 0;
            my $anti_codon_end    = 0;
            my $anti_codon_strand = 0;
            if ($anti_a_pos > $anti_b_pos) {
                $anti_codon_start  = $anti_b_pos;
                $anti_codon_end    = $anti_a_pos;
                $anti_codon_strand = -1;
            }
            else {
                $anti_codon_start  = $anti_a_pos;
                $anti_codon_end    = $anti_b_pos;
                $anti_codon_strand = 1;
            }
            my $location = ($anti_codon_strand == -1) ? "complement($anti_codon_start..$anti_codon_end)"
                         :                              "$anti_codon_start..$anti_codon_end";
            $out_feats{$number}{anti_coords} = $location;
        }
        # Intron data
        elsif ($line =~ m/^Possible intron: (\d+)-(\d+) \((\d+)-(\d+)\)/) {
            my $intron_rel_a_pos = $1;
            my $intron_rel_b_pos = $2;
            my $intron_abs_a_pos = $3;
            my $intron_abs_b_pos = $4;

            # Prepare coordinates
            my $intron_start  = 0;
            my $intron_end    = 0;
            my $intron_strand = 0;
            if ($intron_abs_a_pos > $intron_abs_b_pos) {
                $intron_start  = $intron_abs_b_pos;
                $intron_end    = $intron_abs_a_pos;
                $intron_strand = -1;
            }
            else {
                $intron_start  = $intron_abs_a_pos;
                $intron_end    = $intron_abs_b_pos;
                $intron_strand = 1;
            }
            my $location = ($intron_strand == -1) ? "complement($intron_start..$intron_end)"
                         :                          "$intron_start..$intron_end";
            $out_feats{$number}{intron_coords} = $location;
            # Store relative positions to use at seq line
            $out_feats{$number}{intron_a_pos} = $intron_rel_a_pos;
            $out_feats{$number}{intron_b_pos} = $intron_rel_b_pos;
        }
        elsif ($line =~ m/^Possible pseudogene:  HMM Sc=/) {
            $out_feats{$number}{pseudogene} = 1;
        }
        # Only check seq line if there is an intron to get its sequence
        elsif ($out_feats{$number}{intron_coords} ne '' and $line =~ m/^Seq:\s+(\S+)/) {
            my $tRNA_sequence = $1;
            my $intron_rel_a_pos = $out_feats{$number}{intron_a_pos};
            my $intron_rel_b_pos = $out_feats{$number}{intron_b_pos};
            my $intron_seq = substr $tRNA_sequence, $intron_rel_a_pos - 1,
                                                    $intron_rel_b_pos - $intron_rel_a_pos + 1;

            # Introns should be all lowercase
            if ($intron_seq =~ m/[A-Z]/) {
                die "Given intron coordinates $intron_rel_a_pos-$intron_rel_b_pos does not "
                  . "produce an all lowercase sequence: $intron_seq\n";
            }
            $out_feats{$number}{intron_seq} = $intron_seq;
        }
    }
    close $RESULT;

    # Add features sorted by a location index
    foreach my $id_num (sort {
                              ($out_feats{$a}{loc_idx})
                              <=>
                              ($out_feats{$b}{loc_idx})
                              }
                        keys %out_feats
        ) {
        my $tRNA_type     = "tRNA-$out_feats{$id_num}{tRNA}";
        my $coords        = $out_feats{$id_num}{coords};
        my $score         = $out_feats{$id_num}{score};
        my $anticodon     = $out_feats{$id_num}{anticodon};
        my $anti_coords   = $out_feats{$id_num}{anti_coords};
        my $intron_coords = $out_feats{$id_num}{intron_coords};
        my $intron_seq    = $out_feats{$id_num}{intron_seq};
        my $pseudogene    = $out_feats{$id_num}{pseudogene};

        # Tags to build the feature object
        my %tags = (inference => 'profile:tRNAscan-SE:1.3.1');

        my $product     = ($tRNA_type eq 'tRNA-Undet') ? 'tRNA-OTHER' : $tRNA_type;
        my ($aminoacid) = ($product   =~ m/^tRNA-(\w+)$/);

        my $dna_codon = reverse $anticodon;
        $dna_codon    =~ tr/ATCG/TAGC/;

        my $note = "codon recognized: ";
        $note   .= ($tRNA_type eq 'tRNA-Undet') ? 'Undet' : $dna_codon;

        # Add /anticodon tag only where the information is available
        if ($tRNA_type ne 'tRNA-Undet') {
            # Check if given anticodon coordinates are correct
            my $anti_loc_obj     = Bio::Factory::FTLocationFactory->from_string($anti_coords);
            my $check_anti_codon = $seq_obj->subseq($anti_loc_obj);

            if (uc $check_anti_codon ne uc $anticodon) {
                die "Problems with anticodon position for result number $id_num!\n";
            }

            if ($intron_coords ne '') {
                # Check if given intron coordinates are correct
                my $anti_loc_obj = Bio::Factory::FTLocationFactory->from_string($intron_coords);
                my $check_intron = $seq_obj->subseq($anti_loc_obj);

                if (uc $check_intron ne uc $intron_seq) {
                    die "Problems with intron position for result number $id_num!\n";
                }

                # If an intron was found in tRNA, add the info to $note
                $note .= "; 1 intron: $intron_coords";

            }
            $anticodon = lc $anticodon;
            $tags{anticodon} = "(pos:$anti_coords,aa:$aminoacid,seq:$anticodon)";
        }

        # Add collected information to %tags
        $tags{product} = $product;
        $tags{note}    = "$note; Score: $score";
        # Note: '_no_value' tells BioPerl to print the tag without any value
        $tags{pseudo}  = '_no_value' if ($pseudogene == 1);

        # Build and add feature for the rest of the cases with the information collected in %tags
        my $feat_loc_obj = Bio::Factory::FTLocationFactory->from_string($coords);
        my $feature = Bio::SeqFeature::Generic->new(-primary_tag => 'tRNA',
                                                    -location    => $feat_loc_obj,
                                                    -tag         => \%tags,
                                                    );
        $out_obj->add_SeqFeature($feature);
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.tRNAscan-SE.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;
    }

    # Now that all data is in the Genbank file, delete some intermediary files
    unlink  $out_file        or die "Could not delete '$out_file': $!\n";
    unlink "$out_header.txt" or die "Could not delete '$out_header.txt': $!\n";

    # Return to previous working folder
    chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
    return;
}

# Prediction of genomic islands from a Genbank file.
# NOTE: PAI-DA will only check CDS features which are
# more than 300 bp long
sub paida {
    my ($folder, $file, $seq_obj) = @_;
    if (not -e $folder) {
        warn "Could not find the provided folder for IslandPath-DIMOB: '$folder'\n";
        return;
    }
    if (not -e $file) {
        warn "Could not find the provided file for IslandPath-DIMOB: '$file'\n";
        return;
    }

    my $seq_length = $seq_obj->length;

    # Keep previous working directory to return to it after program execution
    my $old_folder = getcwd;
    if ($file !~ m'[\\/]') {
        $file = "$old_folder/$file";
    }
    # Change to program directory
    chdir $folder or die "Could not change directory to '$folder': $!\n";

    # Define output
    my $out_header = ($file =~ m/^(.+)\.gbk?$/i) ? $1 : $file;
    my $out_file1  = "$out_header.txt";
    my $out_file2  = "$out_header.PAI-DA.score";

    # Program scan parameters
    my $window = '20000'; # Default: 20000
    my $step   = '5000';  # Default:  5000

    # Run the program Scan and Score scripts
    my $null = ($^O =~ m/mswin/i) ? 'NUL' : '/dev/null';
    system("perl compscan.pl -w $window -s $step -i \"$file\" -o \"$out_file1\" > $null");
    system("perl compscore.pl -i \"$out_file1\" -o \"$out_file2\" ");

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Feature collector
    my %out_feats;

    # The PAI-DA paper use a score cutoff of 3.9, but after trying to replicate
    # E. Coli O157:H7 EDL933 results with an updated annotation, the same predicted
    # islands needed a score cutoff of 3.6
    my $score_cutoff = 2.3;

    # Is necessary to combine adjacent regions who are above the score cutoff into
    # single feature that span over them and average their score
    my $save_state  = 0;
    my $regions     = 0;
    my $start       = 0;
    my $total_score = 0;
    my $max_score   = 0;

    # Get genomic island predictions from the Score file
    my $feat_counter = 0;
    open my $RESULT, '<', $out_file2 or die "Could not read file '$out_file2': $!\n";
    while (my $line = <$RESULT>) {
        chomp $line;
        next if $line =~ m/^#/;
        next if $line eq '';

        # The output table have many columns, but only capture
        # the location (first) and score (last)
        $line =~ s/ //g; # remove whitespaces but leave tabs
        my @columns = split /\t/, $line;
        my $score   = $columns[-1];

        if ($save_state == 1 and $score < $score_cutoff) {
            # A genomic island have just finished (and the present region
            # is below cutoff), so add the finished island to the results
            $feat_counter++;
            my $end     = $start + $regions*$step - 1;
            my $loc_idx = $start + (1 - $end / $seq_length);

            my $feat_score = $total_score / $regions;
            $feat_score    = sprintf "%.6f", $feat_score;
            $max_score     = $feat_score if ($feat_score > $max_score);

            # Build feature with the information collected in %tags
            my $feature = Bio::SeqFeature::Generic
                ->new( -primary_tag => 'misc_feature',
                       -start       => $start,
                       -end         => $end,
                       -tag         => { inference => 'profile:PAI-DA:1.1',
                                         score     => $feat_score,
                                         note      => 'genomic island' },
            );

            $out_feats{$feat_counter}{loc_idx} = $loc_idx;
            $out_feats{$feat_counter}{feature} = $feature;
            $out_feats{$feat_counter}{score}   = $feat_score;

            # Reset variables
            $save_state  = 0;
            $regions     = 0;
            $start       = 0;
            $total_score = 0;
        }
        elsif ($save_state == 1) {
            # A genomic island in growing state
            $regions++;
            $total_score += $score;
        }
        elsif ($save_state == 0 and $score < $score_cutoff) {
            # Discarded region, no genomic island in the previous region
            next;
        }
        else {
            # A new genomic island detected
            $save_state  = 1;
            $regions     = 1;
            $start       = $columns[0];
            $total_score = $score;
        }
    }
    close $RESULT;

    # Check if there was a final genomic island pending to be added
    if ($save_state == 1) {
        # A genomic island have just finished (and the present region
        # is below cutoff), so add the finished island to the results
        $feat_counter++;
        my $end     = $start + $regions*$step - 1;
        my $loc_idx = $start + (1 - $end / $seq_length);

        my $feat_score = $total_score / $regions;
        $feat_score    = sprintf "%.6f", $feat_score;
        $max_score     = $feat_score if ($feat_score > $max_score);

        # Build feature with the information collected in %tags
        my $feature = Bio::SeqFeature::Generic
            ->new( -primary_tag => 'misc_feature',
                   -start       => $start,
                   -end         => $end,
                   -tag         => { inference => 'profile:PAI-DA:1.1',
                                     score     => $feat_score,
                                     note      => 'genomic island', },
        );

        $out_feats{$feat_counter}{loc_idx} = $loc_idx;
        $out_feats{$feat_counter}{feature} = $feature;
        $out_feats{$feat_counter}{score}   = $feat_score;

        # Reset variables
        $save_state  = 0;
        $regions     = 0;
        $start       = 0;
        $total_score = 0;
    }

    # Add features sorted by a location index
    foreach my $key ( sort { $out_feats{$a}{loc_idx} <=> $out_feats{$b}{loc_idx} }
                      keys %out_feats
        ) {
        my $loc_idx    = $out_feats{$key}{loc_idx};
        my $out_feat   = $out_feats{$key}{feature};
        my $feat_score = $out_feats{$key}{score};

        # Add blue color grading to feature depending of their score
        my $value = ($feat_score - $score_cutoff) / ($max_score - $score_cutoff);
        my $red   = int(255 - ($value * 255));
        my $green = int(255 - ($value * 255));
        my $blue  = 255;
        $out_feat->add_tag_value('colour', "$red $green $blue");

        $out_obj->add_SeqFeature($out_feat);
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.PAI-DA.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;
    }
    else {
        # Delete score file if there are no significant results
        unlink $out_file2 or die "Could not delete '$out_file2': $!\n";
    }

    # Now that all data is in the Genbank file, delete some intermediary files
    unlink $out_file1 or die "Could not delete '$out_file1': $!\n";

    # Return to previous working folder
    chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
    return;
}

# Prediction of genomic islands
sub island_path {
    my ($folder, $file, $seq_obj) = @_;
    if (not -e $folder) {
        warn "Could not find the provided folder for IslandPath-DIMOB: '$folder'\n";
        return;
    }
    if (not -e $file) {
        warn "Could not find the provided file for IslandPath-DIMOB: '$file'\n";
        return;
    }

    my $seq_length = $seq_obj->length;

    # Keep previous working directory to return to it after program execution
    my $old_folder = getcwd;
    if ($file !~ m'[\\/]') {
        $file = "$old_folder/$file";
    }
    # Change to program directory
    chdir $folder or die "Could not change directory to '$folder': $!\n";

    # Accessory files needed for IslandPath-DIMOB
    my $out_header = ($file =~ m/^(.+)\.(?:fa|fna|fasta)$/i) ? $1 : $file;
    my $file_faa   = "$out_header.faa";
    my $file_ffn   = "$out_header.ffn";
    my $file_ptt   = "$out_header.ptt";
    if (not -e $file_faa) {
        warn "Could not find needed file for IslandPath-DIMOB: '$file_faa'\n";
        return;
    }
    if (not -e $file_ffn) {
        warn "Could not find needed file for IslandPath-DIMOB: '$file_ffn'\n";
        return;
    }
    if (not -e $file_ptt) {
        warn "Could not find needed file for IslandPath-DIMOB: '$file_ptt'\n";
        return;
    }

    # Run program and capture standard output
    my $output = `perl dimob.pl "$file_faa" "$file_ffn" "$file_ptt"`;

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Feature collector
    my %out_feats;

    # Get genomic island predictions from program output
    my $feat_counter = 0;
    my @results = ($output =~ m/(\d+\t\d+)/g);
    foreach my $result (@results) {
        $feat_counter++;

        my ($start, $end) = ($result =~ m/^(\d+)\t(\d+)$/);
        my $loc_idx = $start + (1 - $end / $seq_length);

        # Build feature with the information collected in %tags
        my $feature = Bio::SeqFeature::Generic
            ->new( -primary_tag => 'misc_feature',
                   -start       => $start,
                   -end         => $end,
                   -tag         => { inference => 'profile:IslandPath-DIMOB:0.2',
                                     note      => 'genomic island',
                                     colour    => '0 255 0',},
        );

        $out_feats{$feat_counter}{loc_idx} = $loc_idx;
        $out_feats{$feat_counter}{feature} = $feature;
    }

    # Add features sorted by a location index
    foreach my $key ( sort { $out_feats{$a}{loc_idx} <=> $out_feats{$b}{loc_idx} }
                      keys %out_feats
        ) {
        my $out_feat = $out_feats{$key}{feature};
        $out_obj->add_SeqFeature($out_feat);
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.IslandPath.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;
    }

    # Now that all data is in the Genbank file, delete some intermediary files
    unlink "$file_faa.mob" or die "Could not delete '$file_faa.mob': $!\n";

    # Return to previous working folder
    chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
    return;
}

# Prediction of tRNA, mtRNA, and tmRNA gene FRAGMENTS
sub aragorn_frags {
    my ($folder, $file, $seq_obj) = @_;
    if (not -e $folder) {
        warn "Could not find the provided folder for Aragorn: '$folder'\n";
        return;
    }
    if (not -e $file) {
        warn "Could not find the provided file for tRNA fragments search: '$file'\n";
        return;
    }

    my $seq_length = $seq_obj->length;

    # Keep previous working directory to return to it after program execution
    my $old_folder = getcwd;
    if ($file !~ m'[\\/]') {
        $file = "$old_folder/$file";
    }
    # Change to program directory
    chdir $folder or die "Could not change directory to '$folder': $!\n";

    # Define output
    my ($out_header, $file_type) = ($file =~ m/^(.+)\.([^\.]+)$/i);
    if ($file_type =~ m/^txt$/i) {
        # The temporary files uses ".txt" extension, so use the full
        # filename in this case
        $out_header = $file;
    }
    my $out_file = "$out_header.txt";

    # Get tRNA/tmRNA sequences directly from the Aragorn subroutine results,
    # since running Aragorn again could include results that were already discarded
    my $aragorn_file = "$out_header.Aragorn.gbk";
    # If there are no results, abort and return
    if (not -e $aragorn_file) {
        # Return to previous working folder
        chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
        return;
    }
    my $file_in_obj = Bio::SeqIO->new(-file   => "<$aragorn_file",
                                      -format => 'genbank');
    my $file_out_obj = Bio::SeqIO->new(-file   => ">$out_file",
                                       -format => 'fasta',
                                       -width  => 50);
    while (my $seq_obj = $file_in_obj->next_seq) {
        my $tmRNA_start = 0;
        my @features = $seq_obj->get_SeqFeatures;
        foreach my $feat (@features) {
            my $primary_tag = $feat->primary_tag;
            next if $primary_tag eq 'CDS';

            # Description
            my $start   = $feat->start;
            my $end     = $feat->end;
            my $strand  = $feat->strand;
            my $loc_str = ($strand == -1) ? "c[$start,$end]" : "[$start,$end]";
            # Skip rare split-by-origin tRNA predictions
            # (e.g. complement(join(2298937..2298983,1..29)))
            next if ($start > $end);

            # Display_id
            my $anticodon = '';
            if ($feat->has_tag('anticodon')) {
                my $tag = ($feat->get_tag_values('anticodon'))[0];
                ($anticodon) = ($tag =~ m/,\s*seq\s*:([\w\s]+)/);
                $anticodon =~ s/\s+//g;
            }

            my $product = $feat->has_tag('product') ? ($feat->get_tag_values('product'))[0]
                        :                              '';
            if ($product eq 'tRNA-OTHER') {
                my $note = $feat->has_tag('note') ? ($feat->get_tag_values('note'))[0]
                         :                              '';
                my ($new_product) = ($note =~ m/;\s+(.+)$/);
                $new_product =~ s/\s+//g;
                $product = $new_product;
            }
            elsif (     $tmRNA_start == 0
                and ($product eq 'tmRNA acceptor piece' or $product eq 'tmRNA coding piece')
                ) {
                $tmRNA_start  = $start;
                next;
            }
            elsif (     $tmRNA_start != 0
                   and ($product eq 'tmRNA acceptor piece' or $product eq 'tmRNA coding piece')
                ) {
                $start       = $tmRNA_start;
                $tmRNA_start = 0; # Reset
                $loc_str     = ($strand == -1) ? "c[$start,$end]" : "[$start,$end]";
                $loc_str     = "(Permuted) $loc_str";
            }

            my $display_id = ($primary_tag eq 'tmRNA') ? $primary_tag
                           :                             "$product($anticodon)";

            # Sequence
            my $feat_seq = $seq_obj->subseq($start, $end);
            $feat_seq    = lc $feat_seq;
            if ($strand == -1) {
                $feat_seq = reverse $feat_seq;
                $feat_seq =~ tr/atcg/tagc/;
            }
            my $fa_obj = Bio::Seq->new(-display_id => $display_id,
                                       -desc       => $loc_str,
                                       -seq        => $feat_seq);
            $file_out_obj->write_seq($fa_obj);
        }
    }
    close $file_in_obj->_fh;
    close $file_out_obj->_fh;

    # Run Blast using the predicted tRNAs as query and the original sequence as target
    my $blast_out_file = "$out_header.blast.txt";
    my $null           = ($^O =~ m/mswin/i) ? 'NUL' : '/dev/null';
    system("blastn -task dc-megablast -query \"$out_file\" -subject \"$file\" -out \"$blast_out_file\" 2> $null");

    # Read first the tRNA prediction output to store their locations
    # and discard later any fragment that overlaps them
    my @predicted_locations;
    open my $ARAGORN_FILE, '<', $out_file or die "Could not read file '$out_file': $!\n";
    while (my $line = <$ARAGORN_FILE>) {
        if ($line =~ m/^>/) {
            chomp $line;
            my ($location) = ($line =~ m/(c?\[\d+,\d+\])$/);
            push @predicted_locations, $location;
        }
    }
    close $ARAGORN_FILE;

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Feature collector
    my %out_feats;

    # For cases where 2 or more fragments points to the same tRNA
    my %tangled_tRNAs_loc;

    # Process output
    my $blast_obj = Bio::SearchIO->new(-file   => "<$out_header.blast.txt",
                                      -format => 'blast');
    my $number = 0;
    while (my $result = $blast_obj->next_result) {
        # Aragorn Fasta header examples:
        # >tRNA-Val(cac) [226619,226693]
        # >tmRNA [254349,254713]
        my $qry_name   = $result->query_name;
        my $qry_desc   = $result->query_description;
        my $qry_length = $result->query_length;

        while (my $hit = $result->next_hit) {
            while (my $hsp = $hit->next_hsp) {
                $number++;
                my $identity    = $hsp->percent_identity;
                my $hsp_qry_len = $hsp->length('query');

                my $qry_cov = $hsp_qry_len
                            / $qry_length
                            * 100;
                $qry_cov  = sprintf("%.1f", $qry_cov);
                $identity = sprintf("%.1f", $identity);

                if (    $hsp_qry_len > 20
                    and $qry_cov     < 100
                    and $identity   >= 91
                    ) {
                    my $hit_start   = $hsp->start('hit');
                    my $hit_end     = $hsp->end('hit');

                    # First discard any hit where the query and hit strands mismatch,
                    # since strands should be the same for tRNA fragments produced by insertions
                    my ($tRNA_strand, $tRNA_start, $tRNA_end) = ($qry_desc =~ m#(c?)\[(\d+),(\d+)\]#);
                    $tRNA_strand    = ($tRNA_strand eq 'c') ? -1 : 1;
                    my $hit_strand  = $hsp->strand('hit');
                    next if ($hit_strand != $tRNA_strand);

                    # Then discard any hit that overlaps a known tRNA/tmRNA predicted location
                    my $overlap = 0;
                    foreach my $location (@predicted_locations) {
                        my ($loc_start, $loc_end) = ($location =~ m/\[(\d+),(\d+)\]/);
                        if (   ($hit_start >= $loc_start and $hit_start <= $loc_end)
                            or ($hit_end   >= $loc_start and $hit_end <= $loc_end)
                            ) {
                            $overlap = 1;
                            last;
                        }
                    }
                    next if $overlap == 1;

                    # Continue filling data for those hits who passed the previous filter
                    my $loc_idx     = $hit_start + (1 - $hit_end / $seq_length);
                    my $loc_str     = ($hit_strand == -1) ? "complement($hit_start..$hit_end)" : "$hit_start..$hit_end";
                    my $frag_type   = ($qry_name =~ m/tmRNA/) ? 'tmRNA fragment' : 'tRNA fragment';
                    my $description = "$frag_type: $qry_name $qry_desc. Length: $hsp_qry_len nt. "
                                    . "Qry Coverage: $qry_cov%. Identity: $identity%. Location: $loc_str";

                    # Create region coordinates that goes between
                    # the tRNA fragment and the possible tRNA source
                    my $region_start  = ($tRNA_start < $hit_start) ? $tRNA_start : $hit_start;
                    my $region_end    = ($hit_end < $tRNA_end)     ? $tRNA_end   : $hit_end;
                    my $region_len    = $region_end - $region_start + 1;
                    my $region_coords = '';
                    my $loc_index     = '';
                    # Trust the proposed region when the region length is not too short
                    # or too large (best length cutoffs still being refined)
                    if ($region_len >= 4_000 and $region_len <= 300_000) {
                        $region_coords = ($hit_strand != $tRNA_strand) ? "$region_start..$region_end"
                                       : ($hit_strand == -1)           ? "complement($region_start..$region_end)"
                                       :                                 "$region_start..$region_end";
                        $loc_index     = $region_start + (1 - $region_end / $seq_length);
                    }
                    # For the other cases use only the frag coordinates
                    else {
                        $region_coords = $loc_str;
                        $loc_index     = $loc_idx;
                    }

                    my $loc_obj = Bio::Factory::FTLocationFactory->from_string($region_coords);
                    my $feature = Bio::SeqFeature::Generic->new(-location    => $loc_obj,
                                                                -primary_tag => 'misc_feature',
                                                                -tag => {
                                                                         inference => 'profile:Aragorn:1.2.36',
                                                                         note      => $description,
                                                                         product   => "$qry_name fragment",
                                                                         colour    => '255 255 0',},
                                                                );

                    # Store the finished feature
                    $out_feats{$number}{feature} = $feature;
                    $out_feats{$number}{loc_idx} = $loc_index;

                    if ($region_coords ne $loc_str) {
                        $tangled_tRNAs_loc{$qry_desc}{$number}{details} = $description;
                    }
                }
            }
        }
    }
    close $blast_obj->_fh;

    # Only let the best fragment hit to use region coordinates
    # with their predicted source tRNA
    foreach my $tRNA_loc (%tangled_tRNAs_loc) {
        my @fragments = keys %{ $tangled_tRNAs_loc{$tRNA_loc} };
        next if scalar @fragments < 2;

        my $best_hit      = '';
        my $best_location = '';
        my $best_identity = 0;
        my $best_coverage = 0;
        my $best_length   = 0;
        my $best_distance = 0;
        my %non_best_keys;
        foreach my $frag (@fragments) {
            my ($target_loc, $len, $cov, $ident, $loc)
                = ( $tangled_tRNAs_loc{$tRNA_loc}{$frag}{details}
                    =~ m/(\S+)\. Length: (\d+) nt. Qry Coverage: ([\.\d]+)%. Identity: ([\.\d]+)%. Location: (\S+)/);
            my ($target_start, $target_end) = ($target_loc =~ m/(\d+),(\d+)/);
            my ($frag_start,   $frag_end)   = ($loc        =~ m/(\d+)\.\.(\d+)/);
            my $dist = ($frag_start < $target_start) ? ($target_start - $frag_end) : ($frag_start - $target_end);

            if (   (    $ident >  $best_identity)
                or (    $ident == $best_identity
                    and $cov    > $best_coverage)
                or (    $ident == $best_identity
                    and $cov   == $best_coverage
                    and $len    > $best_length)
                or (    $ident == $best_identity
                    and $cov   == $best_coverage
                    and $len   == $best_length
                    and $dist   > $best_distance)
                ) {
                $non_best_keys{$best_hit} = $best_location if ($best_hit ne '');

                $best_hit      = $frag;
                $best_location = $loc;
                $best_identity = $ident;
                $best_coverage = $cov;
                $best_length   = $len;
                $best_distance = $dist;
            }
            else {
                $non_best_keys{$frag} = $loc;
            }
        }

        # Revert the region coordinates to their original for non best hits
        foreach my $non_best (keys %non_best_keys) {
            my $original_loc = $non_best_keys{$non_best};
            my $temp_feature = $out_feats{$non_best}{feature};

            my ($orig_start, $orig_end) = ($original_loc =~ m/(\d+)\.\.(\d+)/);
            my $orig_loc_index          = $orig_start + (1 - $orig_end / $seq_length);

            my $loc_obj = Bio::Factory::FTLocationFactory->from_string($original_loc);
            $temp_feature->location($loc_obj);
            $out_feats{$non_best}{feature} = $temp_feature;
            $out_feats{$non_best}{loc_idx} = $orig_loc_index;
        }
    }

    # Add features sorted by a location index
    foreach my $id_num (sort {
                              ($out_feats{$a}{loc_idx})
                              <=>
                              ($out_feats{$b}{loc_idx})
                              }
                        keys %out_feats
        ) {
        my $feature = $out_feats{$id_num}{feature};
        $out_obj->add_SeqFeature($feature);
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.tRNA_frags.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;
    }

    # Now that all data is in the Genbank file, delete some intermediary files
    unlink $out_file       or die "Could not delete '$out_file': $!\n";
    unlink $blast_out_file or die "Could not delete '$blast_out_file': $!\n";

    # Return to previous working folder
    chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
    return;
}

# Prediction of tRNA, mtRNA, and tmRNA genes
sub aragorn {
    my ($folder, $file, $seq_obj, $codontable) = @_;
    if (not -e $folder) {
        warn "Could not find the provided folder for Aragorn: '$folder'\n";
        return;
    }
    if (not -e $file) {
        warn "Could not find the provided file for Aragorn: '$file'\n";
        return;
    }

    my $desc       = $seq_obj->desc;
    my $sequence   = $seq_obj->seq;
    my $seq_length = $seq_obj->length;

    # Keep previous working directory to return to it after program execution
    my $old_folder = getcwd;
    if ($file !~ m'[\\/]') {
        $file = "$old_folder/$file";
    }
    # Change to program directory
    chdir $folder or die "Could not change directory to '$folder': $!\n";

    # Define output
    my ($out_header, $file_type) = ($file =~ m/^(.+)\.([^\.]+)$/i);
    if ($file_type =~ m/^txt$/i) {
        # The temporary file uses ".txt" extension, so use the full
        # filename in this case
        $out_header = $file;
    }
    my $out_file   = "$out_header.txt";

    # Linux needs './' before the executable binary name
    my $aragorn_exe    = ($^O =~ m/mswin/i) ? 'aragorn_alpha3' : './aragorn_alpha3';

    # Get version, since option -wunix is added in v1.2.37
    my ($version)      = (`$aragorn_exe -h` =~ m/ARAGORN v([\.\d]+) /);
    my ($v1, $v2, $v3) = ($version =~ m/\d+/g);
    my $ver_number     = ($v1*100000) + ($v2*1000) + $v3;
    # NOTE: -wunix option causes massive slowdown for files with a
    # big ammount of sequences. Only use for files with 1 to few sequences.
    my $newline_fix    = ($^O =~ m/mswin/i and $ver_number >= 102037) ? '-wunix' : '';

    # When processing the file with Aragorn v1.2.36 or older, its important
    # for Genbank files to contain the DEFINITION line (this was fixed in v1.2.37),
    # or Aragorn will not recognize it and will return an empty result (some older
    # ORFminer results can miss that line); if its missing it must be added to
    # the file. Also, if the -wunix option is unavailable, is very important to
    # create temporary files with native line-endings to avoid Aragorn's
    # coordinates bug, so the report of false positives caused by overlapping
    # between tRNAs and CDS genes is trustable.
    my $definition_line = 0;
    if ($file_type =~ m/gbk?$/i) {
        open my $IN,  '<', $file or die "Could not read file '$file': $!\n";
        while (my $line = <$IN>) {
            if ($line =~ m/^DEFINITION  /) {
                $definition_line = 1;
                last;
            }
        }

        # Create a temporary file only if DEFINITION is missing
        # or -wunix option is unavailable
        if ($definition_line == 0 or $newline_fix eq '') {
            # Reset filehandle position to file start
            seek $IN, 0, 0;

            # Create temporary file
            my $temp_file = "$out_header.copy.gbk";
            open my $OUT, '>', $temp_file or die "Could not write file '$temp_file': $!\n";

            my $header = '';
            my $flag   = 'header';
            while (my $line = <$IN>) {
                chomp $line;
                if ($flag eq 'header') {
                    if (   $line =~ m/^ACCESSION   /
                        or $line =~ m/^FEATURES    /
                        or $line =~ m/^ORIGIN      /
                        ) {
                        $flag = 'flush';
                        $header .= "DEFINITION  \n" if ($header !~ m/DEFINITION  /);
                    }
                    $header .= "$line\n";
                }
                else {
                    if ($flag eq 'flush') {
                        print $OUT $header;
                        $flag = 'done';
                    }
                    print $OUT "$line\n";
                }
            }
            close $OUT;
            $file = $temp_file;

            # Since the temporary file have native
            # line-endings, -wunix is not needed
            $newline_fix = '';
        }
        close $IN;
    }

    # Run program. e: Show prediction Score, br: Show secondary structure
    # of the tRNA gene, rp: Flag possible pseudogenes
    system("$aragorn_exe $newline_fix -e -br -rp -o \"$out_file\" \"$file\" ");

    # If a temporary file was used, delete it now
    if ($file eq "$out_header.copy.gbk") {
        unlink $file or die "Could not delete temporary file '$file': $!\n";
    }

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Feature collector
    my %out_feats;

    # Process output
    open my $RESULT, '<', $out_file or die "Could not read file '$out_file': $!\n";
    my $number = 0;
    while (my $line = <$RESULT>) {
        # NOTE: Aragorn print the secondary structure for tRNA but not for tmRNA
        chomp $line;

        # Result start
        if ($line =~ m/^(\d+)\.$/) {
            $number = $1;

            if (exists $out_feats{$number}) {
                die "Do not process more than 1 query sequence per run!\n";
            }

            # Prepare the hash keys
            $out_feats{$number}{tRNA}           = '';
            $out_feats{$number}{coords}         = '';
            $out_feats{$number}{loc_idx}        = 0;
            $out_feats{$number}{score}          = 0;
            $out_feats{$number}{primary_seq}    = '';
            $out_feats{$number}{sec_struct}     = '';
            $out_feats{$number}{tag_nt_seq}     = '';
            $out_feats{$number}{tag_aa_seq}     = '';
            $out_feats{$number}{intron}         = 0;
            $out_feats{$number}{pseudo}         = 0;
            $out_feats{$number}{possible_false} = 0;
        }
        # Possible tRNA/tmRNA Pseudogene
        elsif ($line =~ m/^\s*Possible Pseudogene/) {
            $out_feats{$number}{pseudo} = 1;
        }
        # tRNA identity line. Example: tRNA-?(Ala|Cys)(gc)
        elsif ($line =~ m/^\s{4}(tRNA-\S+)/) {
            my $tRNA = $1;
            # For Selenocysteine tRNA
            $tRNA =~ s/SeC/Sec/;

            $out_feats{$number}{tRNA} = $tRNA;
        }
        # tmRNA identity line
        elsif ($line =~ m/^tmRNA Sequence in /) {
            $out_feats{$number}{tRNA} = 'tmRNA';
        }
        # tRNA/tmRNA coordinates
        elsif ($line =~ m/^\s*(?:Sequence|Location) (c?)\[(\d+),(\d+)\]$/) {
            my ($compl, $start, $end) = ($1, $2, $3);

            my $location = ($compl eq 'c') ? "complement($start..$end)" : "$start..$end";
            my $loc_idx  = $start + (1 - $end / $seq_length);
            $out_feats{$number}{coords}  = $location;
            $out_feats{$number}{loc_idx} = $loc_idx;
        }
        # tRNA/tmRNA prediction score
        elsif ($line =~ m/^\s*Score = ([\d\.]+)/) {
            $out_feats{$number}{score} = $1;
        }
        # Possible false positive tRNA/tmRNA (intron or possible pseudogene and overlaps with another gene)
        elsif ($line =~ m/^\s*(?:Annotation false positive|Not annotated)$/) {
            # Only consider as false positive if there is a "false positive" or
            # "Not annotated" flag (tRNA/tmRNA is not present in the original
            # annotation) associated with the overlapping of another gene,
            # AND if its a pseudogene, possess an intron or have a weak prediction
            # Score (that will be checked later)
            $line  = <$RESULT>;
            # There can be multiple overlaps, so use while() to pass through them
            while ($line =~ m/^\s*Overlap with /) {
                $out_feats{$number}{possible_false} = 1;
                $line  = <$RESULT>;
            }
        }
        # tRNA primary sequence and secondary structure
        elsif ($line =~ m/^Primary sequence for /) {
            # Pass to primary sequence line
            $line  = <$RESULT>;
            chomp $line;
            $out_feats{$number}{primary_seq} = $line;

            # Pass to secondary structure line
            $line  = <$RESULT>;
            chomp $line;
            $out_feats{$number}{sec_struct}  = $line;
        }
        # tRNA intron existence
        elsif ($line =~ m/^Intron from /) {
            $out_feats{$number}{intron} = 1;
        }
        # tmRNA permuted 2-pieces
        elsif ($line =~ m/^Permuted$/) {
            $out_feats{$number}{tRNA} .= '-permuted';
        }
        # tmRNA primary sequence OR tRNA intron
        elsif ($line =~ m/^1   .   10    .   20    .   30/) {
            # Skip the sequence if its a tRNA intron
            next if $out_feats{$number}{intron} == 1;

            # Pass to first primary sequence line
            $line  = <$RESULT>;
            chomp $line;

            # Append the sequence until an empty line
            my $primary_seq = $line;
            while ($line = <$RESULT>) {
                chomp $line;
                last if $line eq '';

                $primary_seq .= $line;
            }
            $out_feats{$number}{primary_seq} = $primary_seq;
        }
        # tmRNA domains location (optional, starting from v1.2.37),
        # used later to double-check permuted acceptor and coding pieces
        elsif ($line =~ m/^(\d)' tRNA domain at \[(\d+,\d+)\]$/) {
            my ($dom, $loc_str) = ($1, $2);
            $out_feats{$number}{"tmRNA $dom"} = $loc_str;
        }
        # tmRNA resume consensus sequence check
        elsif ($line =~ m/^Resume consensus sequence/) {
            my $cons_seq   = '';
            my $cons_start = 0;
            my $cons_end   = 0;
            # v1.2.36
            if ($line =~ m/^Resume consensus sequence \(at (\d+)\): (\w+)/) {
                $cons_start = $1;
                $cons_seq   = $2;
            }
            # starting from v1.2.37
            elsif ($line =~ m/^Resume consensus sequence at \[(\d+),(\d+)\]: (\w+)/) {
                $cons_start = $1;
                $cons_end   = $2;
                $cons_seq   = $3;
            }
            else {
                die "Could not retrieve consensus sequence location from description for result $number!!!\n";
            }

            my $qry_string = $out_feats{$number}{primary_seq} || '';
            if ($qry_string ne '') {
                my $match_found = 0;
                while ($qry_string =~ m/$cons_seq/pig) {
                    $match_found++;
                    my $length    = length( ${^MATCH} );
                    my $end_pos   = pos($qry_string);
                    my $start_pos = $end_pos - $length + 1;
                    if ($cons_start != $start_pos) {
                        warn "Consensus sequence start conflict at result $number: Aragorn=$cons_start, found=$start_pos\n";
                    }
                    if ($cons_end != 0 and $cons_end != $end_pos) {
                        warn "Consensus sequence end conflict at result $number: Aragorn=$cons_end, found=$end_pos\n";
                    }
                }
                if ($match_found == 0) {
                    die "Could not find consensus sequence in tmRNA sequence for result $number!!!\n";
                }
            }
            else {
                die "No tmRNA sequence found for result $number!!!\n";
            }
        }
        # tmRNA tag peptide
        elsif ($line =~ m/^Tag peptide /) {
            # Preserve location line for later checkup
            my $tag_location_line = $line;

            # Pass to tag nucletide sequence line
            $line  = <$RESULT>;
            chomp $line;
            my ($tag_nt_seq) = ($line =~ m/Tag sequence: (\S+)/);
            $tag_nt_seq =~ s/-//g;

            # Pass to tag peptide sequence line (the one with 1-letter abbreviations)
            $line  = <$RESULT>; # 3-letter peptide sequence line
            $line  = <$RESULT>; # 1-letter peptide sequence line
            chomp $line;
            my ($tag_aa_seq) = ($line =~ m/Tag peptide:  (\S+)/);

            # Sometimes there are multiple stop codons in the tag peptide,
            # so only leave the first stop codon in $tag_nt_seq and then
            # remove stops from $tag_aa_seq
            if ($tag_aa_seq =~ m/\*{2,}$/) {
                my $stop_count = $tag_aa_seq =~ tr/\*//;
                my $extra_nts  = ($stop_count - 1) * 3;

                $tag_nt_seq =~ s/.{$extra_nts}$//;
            }
            $tag_aa_seq =~ s/\*//g;

            $out_feats{$number}{tag_nt_seq} = $tag_nt_seq;
            $out_feats{$number}{tag_aa_seq} = $tag_aa_seq;

            # Tag peptide location check
            my $tag_start = 0;
            my $tag_end   = 0;
            # v1.2.36
            if ($tag_location_line =~ m/^Tag peptide \(at (\d+)\)/) {
                $tag_start = $1;
            }
            # starting from v1.2.37
            elsif ($tag_location_line =~ m/^Tag peptide at \[(\d+),(\d+)\]/) {
                $tag_start = $1;
                $tag_end   = $2;
            }
            else {
                die "Could not retrieve consensus sequence location from description for result $number!!!\n";
            }

            my $qry_string = $out_feats{$number}{primary_seq} || '';
            if ($qry_string ne '') {
                my $match_found = 0;
                while ($qry_string =~ m/$tag_nt_seq/pig) {
                    $match_found++;
                    my $length    = length( ${^MATCH} );
                    my $end_pos   = pos($qry_string);
                    my $start_pos = $end_pos - $length + 1;
                    if ($tag_start != $start_pos) {
                        warn "Tag peptide start conflict at result $number: Aragorn=$tag_start, found=$start_pos\n";
                    }
                    if ($tag_end != 0 and $tag_end != $end_pos) {
                        warn "Tag peptide end conflict at result $number: Aragorn=$tag_end, found=$end_pos\n";
                    }
                }
                if ($match_found == 0) {
                    die "Could not find tag peptide in tmRNA sequence for result $number!!!\n";
                }
            }
            else {
                die "No tmRNA sequence found for result $number!!!\n";
            }
        }
    }
    close $RESULT;

    # Add features sorted by a location index
    foreach my $id_num (sort {
                              ($out_feats{$a}{loc_idx})
                              <=>
                              ($out_feats{$b}{loc_idx})
                              }
                        keys %out_feats
        ) {
        # Skip if hit is a possible false positive, and possess
        # introns in its sequence or is flagged as pseudogene
        if (    $out_feats{$id_num}{possible_false} == 1
            and (   $out_feats{$id_num}{sec_struct} =~ m/i+/
                 or $out_feats{$id_num}{pseudo}     == 1
                 # Score >= 110 cutoff suggested by NCBI annotation staff
                 or $out_feats{$id_num}{score}       < 110)
            ) {
            next;
        }

        # Hit features
        my $tRNA        = $out_feats{$id_num}{tRNA};
        my $coords      = $out_feats{$id_num}{coords};
        my $score       = $out_feats{$id_num}{score};
        my $primary_seq = $out_feats{$id_num}{primary_seq};
        my $sec_struct  = $out_feats{$id_num}{sec_struct};
        my $tag_nt_seq  = $out_feats{$id_num}{tag_nt_seq};
        my $tag_aa_seq  = $out_feats{$id_num}{tag_aa_seq};
        my $pseudo      = $out_feats{$id_num}{pseudo};
        my $primary_tag = ($tRNA =~ m/tmRNA/) ? 'tmRNA' : 'tRNA';

        if (    $tRNA        eq ''
            or  $coords      eq ''
            or  $score       == 0
            or  $primary_seq eq ''
            or ($sec_struct  eq '' and $primary_tag eq 'tRNA')
            or ($tag_nt_seq  eq '' and $primary_tag eq 'tmRNA')
            or ($tag_aa_seq  eq '' and $primary_tag eq 'tmRNA')
            ) {
            die "  ERROR!!! Some data is missing after parsing result number $id_num:\n"
              . "* tRNA         = $tRNA\n"
              . "* coords       = $coords\n"
              . "* score        = $score\n"
              . "* primaary_seq = $primary_seq\n"
              . "* sec_struct   = $sec_struct\n"
              . "* tag_nt_seq   = $tag_nt_seq\n"
              . "* tag_aa_seq   = $tag_aa_seq\n"
        }

        # Bug reported in v1.2.36 and fixed in v1.2.37 through option -wunix:
        # While the printed coordinates are fine in Linux, they are displaced
        # in Windows when the file use linux-line-ending (\n) instead of windows
        # line-ending (\r\n), which causes cumulative coordinate errors for multiple
        # query sequences. When processing only 1 query sequence, Aragorn almost
        # always MISSES THE FIRST NT, but in at least 1 case it adds in the sequence
        # CHARACTERS FROM THE >HEADER LINE!!! Because of this, in Windows is imperative
        # to compare the subsequence obtained directly from the original sequence
        # using the result coordinates, with the result printed sequence and check
        # they match (an add the missing nucleotide to $primary_seq if necessary).
        # Note: can't use "Bio::Factory::FTLocationFactory->from_string()" here,
        # it won't catch divided-by-origin predictions (e.g. c[2298936,28])
        my ($start) = ($coords =~ m/(\d+)\.\./);
        my ($end)   = ($coords =~ m/\.\.(\d+)/);
        my $strand  = ($coords =~ m/complement/) ? -1 : 1;

        if ($^O =~ m/mswin/i) {
            # Sequence given by the coordinates
            my $coord_seq = '';
            if ($end < $start) {
                # This happens when Aragorn prediction location is split
                # by origin in a circular sequence (e.g. c[2298936,28])
                $coord_seq  = substr $sequence, $start-1, $seq_length-$start+1;
                $coord_seq .= substr $sequence, 0, $end;
            }
            else {
                $coord_seq = substr $sequence, $start-1, $end-$start+1;
            }
            $coord_seq = uc $coord_seq;
            if ($strand == -1) {
                $coord_seq = reverse $coord_seq;
                $coord_seq =~ tr/ATCG/TAGC/;
            }

            # Aragorn transforms ambiguous nucleotides in '.',
            # so replace any non-ATCG for '.' before comparing
            $coord_seq =~ s/[^ATCG]/\./g;

            if (uc $coord_seq ne uc $primary_seq) {
                # Die if there is an error despite
                # using the new -wunix option
                if ($newline_fix eq '-wunix') {
                    die "  ERROR!!! Result number $id_num coords $coords are wrong despite using -wunix!\n";
                }

                # Adjust coordinates
                $start++;
                $end++;

                # Proceed to check if adjustment was successful
                my $first_nt = ''; # Will be needed in split-by-origin cases
                if ($end < $start) {
                    # This happens when Aragorn prediction location is split
                    # by origin in a circular sequence (e.g. c[2298936,28])
                    $coord_seq  = substr $sequence, $start-1, $seq_length-$start+1;
                    $coord_seq .= substr $sequence, 0, $end;

                    # Recover the missing nucleotide to add it to $primary_seq if the
                    # range involves the first nucleotide of the original sequence
                    $first_nt = lc substr $sequence, 0, 1;
                }
                else {
                    $coord_seq = substr $sequence, $start-1, $end-$start+1;
                }
                $coord_seq = uc $coord_seq;

                if ($strand == -1) {
                    $coord_seq = reverse $coord_seq;
                    $coord_seq =~ tr/ATCG/TAGC/;
                    $first_nt  =~ tr/atcg/tagc/;
                }

                # Aragorn transforms ambiguous nucleotides in '.',
                # so replace any non-ATCG for '.' before comparing
                $first_nt  =~ s/[^atcg]/\./g;
                $coord_seq =~ s/[^ATCG]/\./g;

                # Add the missing nucleotide to $primary_seq in split-by-origin
                # cases, and add an empty space to $sec_struct to correct
                # anticodon location (but beware that $sec_struct is shorter
                # than $primary_seq)
                if ($end < $start) {
                    my $insertion_pos = ($strand == -1) ? ($end-1) : ($seq_length-$start+1);

                    # Correct $primary_seq
                    substr $primary_seq, $insertion_pos, 0, $first_nt;

                    # Correct $sec_struct
                    if ($insertion_pos <= length $sec_struct) {
                        # When the missing nucleotide maps into the anticodon,
                        # there are 2 options. If the insertion position or the previous
                        # position is used by an intron that is splitting the anticodon,
                        # then the missing nucleotide is part of an intron and an 'i'
                        # is inserted; if there is no intron there, then add the empty space
                        # just before the anticodon to avoid braking it
                        my $pos = 0;
                        my $len = 0;
                        while ($sec_struct =~ m/(A+i*A+)/g ) {
                            $len = length $1;
                            $pos = pos $sec_struct;
                        }

                        # If insertion occurs outside the anticodon
                        if (   $insertion_pos >= $pos
                            or $insertion_pos <= $pos - $len
                            ) {
                            substr $sec_struct, $insertion_pos, 0, ' ';
                        }
                        else {
                            # If insertion occurs in the middle of the anticodon
                            my $insert_neighborhood = substr $sec_struct, $insertion_pos-1, 2;

                            # Possibilities: Ai, iA, ii
                              ($insert_neighborhood =~ m/i/) ? (substr $sec_struct, $insertion_pos, 0, 'i')
                            :                                  (substr $sec_struct, $pos-$len,      0, ' ');
                        }
                    }
                }

                if (uc $coord_seq ne uc $primary_seq) {
                    # If there is still no match we have the worst scenario,
                    # Aragorn probably did add characters from the >header line
                    # to the sequence. I don't like it, but now is necessary
                    # to SEARCH the $primary_seq in the whole query $sequence.
                    # Since there is always a possibility that there were 2 identicals
                    # sequences in different locations, ONLY KEEP the result that
                    # appears in an close +/- range given by $desc length
                    my @guess_offset;
                    my $test_seq     = $sequence;
                    my $desc_length  = length $desc;

                    # If it were a case of split-by-origin, just drop it because the prediction
                    # would include non-sequence characters and would be too altered to make it valid
                    if ($end < $start) {
                        die "Result number $id_num coords $coords are wrong and could not be adjusted properly!\n";
                    }

                    # Proceed with the location rescue attempt
                    if ($strand == -1) {
                        $test_seq = reverse $test_seq;
                        $test_seq =~ tr/ATCG/TAGC/;
                        while ($test_seq =~ m/$primary_seq/gi) {
                            my $pos = pos $test_seq;
                            # Transform $pos from reverse coord to forward coord
                            $pos = abs ($pos - $seq_length) + 1;

                            if (    $pos >= $start-$desc_length-10
                                and $pos <= $start+$desc_length+10
                                ) {
                                my $offset = $pos - $start;
                                push @guess_offset, $offset;
                            }
                        }
                    }
                    else {
                        while ($test_seq =~ m/$primary_seq/gi) {
                            my $pos = pos $test_seq;

                            if (    $pos >= $end-$desc_length-10
                                and $pos <= $end+$desc_length+10
                                ) {
                                my $offset = $pos - $end;
                                push @guess_offset, $offset;
                            }
                        }
                    }

                    # Continue only if there is a match in only 1 location
                    if (scalar @guess_offset == 0 or scalar @guess_offset >= 2) {
                        die "Result number $id_num coords $coords are wrong and could not be adjusted properly!\n";
                    }

                    # Now make the final check which hopefully will pass
                    $start += $guess_offset[0];
                    $end   += $guess_offset[0];

                    my $guess_seq = substr $sequence, $start-1, $end-$start+1;
                    $guess_seq = uc $guess_seq;

                    if ($strand == -1) {
                        $guess_seq = reverse $guess_seq;
                        $guess_seq =~ tr/ATCG/TAGC/;
                    }

                    # Aragorn transforms ambiguous nucleotides in '.',
                    # so replace any non-ATCG for '.' before comparing
                    $guess_seq =~ s/[^ATCG]/\./g;

                    if (uc $guess_seq ne uc $primary_seq) {
                        die "Result number $id_num coords $coords are wrong and could not be adjusted properly!\n";
                    }
                }
            }
        }

        # For use after next IF when defining tmRNA
        # tag peptide location in its own CDS
        my $tag_peptide_start = 0;
        my $tag_peptide_end   = 0;

        # Tags to build the feature object
        my %tags = (inference => 'profile:Aragorn:1.2.36');
        if ($primary_tag eq 'tRNA') {
            # Ex: tRNA-?(Ala|Cys)(gc), tRNA-Lys(ctt), tRNA-???(.ag)
            my ($tRNA_type, $anticodon) = ($tRNA =~ m/^(tRNA-\S+)
                                                      \(
                                                        ( [\.\w]+ ) # '.' is used for ambiguous nts
                                                       \)/x);
            my $product     = ($tRNA_type =~ m/^tRNA-\w+$/)   ? $tRNA_type : 'tRNA-OTHER';
            my ($aminoacid) = ($product   =~ m/^tRNA-(\w+)$/);

            my $dna_codon = uc reverse $anticodon;
            $dna_codon    =~ tr/ATCG/TAGC/;

            my $note = "codon recognized: $dna_codon";
            $note   .= "; $tRNA_type" if $product eq 'tRNA-OTHER';

            # Calculate anticodon position
            my $anti_codon_pos   = 0;
            my $anti_codon_len   = 0; # Anticodon may be 2, 3 or 4 nts

            # Beware of introns splitting apart the anticodon
            my $intron_pos   = 0;
            my $intron_len   = 0;
            my $intron_split = 0; # Flag for anticodon split by intron

            while ($sec_struct =~ m/(A+(i*)A+)/g ) {
                my $alignment   = $1;
                my $intron      = $2;
                $anti_codon_len = length $alignment;
                $anti_codon_pos = pos $sec_struct;

                # Intron in the middle of anticodon
                if ($intron ne '') {
                    $intron_len   = length $intron;
                    $intron_pos   = index($sec_struct, $intron) + $intron_len ;
                    $intron_split = 1;
                }
            }
            if ($anti_codon_pos == 0) {
                die "Problems detecting anticodon position in the primary seq for result number $id_num!\n";
            }

            # If there is no intron in the anticodon, check for its presence
            # in other places of the tRNA structure
            if ($intron_len == 0 and $sec_struct =~ m/(i+)/g) {
                $intron_len = length $1;
                $intron_pos = pos $sec_struct;
            }

            # Calculate coordinates for anticodon and intron (if exists)
            my $anti_codon_start   = 0;
            my $anti_codon_end     = 0;
            my $intron_start       = 0;
            my $intron_end         = 0;
            # For split anticodon locations
            my $anti_codon_inner_a = 0;
            my $anti_codon_inner_b = 0;

            if ($intron_len != 0) {
                # Intron coordinates
                if ($strand == -1) {
                    $intron_start = $end - $intron_pos + 1;
                    $intron_end   = $end - $intron_pos + $intron_len;
                }
                else {
                    $intron_start = $start + $intron_pos - $intron_len;
                    $intron_end   = $start + $intron_pos - 1;
                }

                # Inner coordinates for anticodons split by an intron
                if ($intron_split == 1) {
                    $anti_codon_inner_a = $intron_start - 1;
                    $anti_codon_inner_b = $intron_end   + 1;
                }
            }

            # Anticodon coordinates
            if ($strand == -1) {
                $anti_codon_start = $end - $anti_codon_pos + 1;
                $anti_codon_end   = $end - $anti_codon_pos + $anti_codon_len;
            }
            else {
                $anti_codon_start = $start + $anti_codon_pos - $anti_codon_len;
                $anti_codon_end   = $start + $anti_codon_pos - 1;
            }

            # Adjustments for split-by-origin locations,
            # where start or end could be out of bounds
            $anti_codon_start += $seq_length if ($anti_codon_start < 1);
            $anti_codon_start -= $seq_length if ($anti_codon_start > $seq_length);
            $anti_codon_end   += $seq_length if ($anti_codon_end   < 1);
            $anti_codon_end   -= $seq_length if ($anti_codon_end   > $seq_length);
            if ($intron_len != 0) {
                $intron_start       += $seq_length if ($intron_start < 1);
                $intron_start       -= $seq_length if ($intron_start > $seq_length);
                $intron_end         += $seq_length if ($intron_end   < 1);
                $intron_end         -= $seq_length if ($intron_end   > $seq_length);
            }
            if ($intron_split == 1) {
                $anti_codon_inner_a += $seq_length if ($anti_codon_inner_a < 1);
                $anti_codon_inner_a -= $seq_length if ($anti_codon_inner_a > $seq_length);
                $anti_codon_inner_b += $seq_length if ($anti_codon_inner_b < 1);
                $anti_codon_inner_b -= $seq_length if ($anti_codon_inner_b > $seq_length);
            }

            # Set of the intron coodinate string
            my $intron_coords = '';
            if ($intron_len != 0) {
                # Intron passing through sequence origin (End..1)
                if ($intron_end < $intron_start) {
                    my $a_region
                        = ($intron_start == $seq_length) ? "$intron_start"
                        :                                  "$intron_start..$seq_length";
                    my $b_region
                        = ($intron_end   == 1)           ? "1"
                        :                                  "1..$intron_end";

                    $intron_coords = ($strand == -1) ? "complement(join($a_region,$b_region))"
                                   :                   "join($a_region,$b_region)";
                }
                # Normal intron
                else {
                    $intron_coords = ($strand == -1) ? "complement($intron_start..$intron_end)"
                                   :                   "$intron_start..$intron_end";
                }
            }

            # Set of the anticodon coodinate string
            my $anti_codon_coords = '';
            # /anticodon=(pos:34..36,aa:Phe,seq:aaa)
            # /anticodon=(pos:join(5,495..496),aa:Leu,seq:taa)
            # /anticodon=(pos:complement(4156..4158),aa:Gln,seq:ttg)
            if ($anti_codon_end < $anti_codon_start) {
                # Split anticodon passing through sequence origin (End..1)
                if ($intron_split == 1) {
                    my $a_region
                        = ($anti_codon_start == $seq_length)         ? "$anti_codon_start"
                        : ($anti_codon_start == $anti_codon_inner_a) ? "$anti_codon_start"
                        :                                              "$anti_codon_start..$anti_codon_inner_a";
                    my $b_region
                        = ($anti_codon_end   == 1)                   ? "1"
                        : ($anti_codon_end   == $anti_codon_inner_b) ? "$anti_codon_end"
                        :                                              "$anti_codon_inner_b..$anti_codon_end";

                    $anti_codon_coords = ($strand == -1) ? "complement(join($a_region,$b_region))"
                                       :                   "join($a_region,$b_region)";
                }
                # Non-split anticodon passing through sequence origin (End..1)
                else {
                    my $a_region
                        = ($anti_codon_start == $seq_length) ? "$anti_codon_start"
                        :                                      "$anti_codon_start..$seq_length";
                    my $b_region
                        = ($anti_codon_end   == 1)           ? "1"
                        :                                      "1..$anti_codon_end";

                    $anti_codon_coords = ($strand == -1) ? "complement(join($a_region,$b_region))"
                                       :                   "join($a_region,$b_region)";
                }
            }
            else {
                # Normal split anticodon
                if ($intron_split == 1) {
                    my $a_region = ($anti_codon_start == $anti_codon_inner_a) ? "$anti_codon_start"
                                 :                                              "$anti_codon_start..$anti_codon_inner_a";
                    my $b_region = ($anti_codon_end   == $anti_codon_inner_b) ? "$anti_codon_end"
                                 :                                              "$anti_codon_inner_b..$anti_codon_end";

                    $anti_codon_coords = ($strand == -1) ? "complement(join($a_region,$b_region))"
                                       :                   "join($a_region,$b_region)";
                }
                # Normal non-split anticodon
                else {
                    $anti_codon_coords = ($strand == -1) ? "complement($anti_codon_start..$anti_codon_end)"
                                       :                   "$anti_codon_start..$anti_codon_end";
                }
            }

            # Check if calculated coordinates are correct. If the coords pass
            # through the origin or there is an intron in the anticodon,
            # concatenate the different segments
            my $check_anti_codon = '';
            if ($anti_codon_end < $anti_codon_start) {
                my $a_region = '';
                my $b_region = '';
                if ($intron_split == 1) {
                    $a_region = substr $sequence,
                                       $anti_codon_start - 1,
                                       $anti_codon_inner_a - $anti_codon_start + 1;
                    $b_region = substr $sequence,
                                       $anti_codon_inner_b - 1,
                                       $anti_codon_end - $anti_codon_inner_b + 1
                }
                else {
                    $a_region = substr $sequence,
                                       $anti_codon_start - 1,
                                       $seq_length - $anti_codon_start + 1;
                    $b_region = substr $sequence,
                                       0,
                                       $anti_codon_end;
                }
                $check_anti_codon = $a_region . $b_region;
            }
            else {
                if ($intron_split == 1) {
                    my $a_region = substr $sequence,
                                          $anti_codon_start - 1,
                                          $anti_codon_inner_a - $anti_codon_start + 1;
                    my $b_region = substr $sequence,
                                          $anti_codon_inner_b - 1,
                                          $anti_codon_end - $anti_codon_inner_b + 1;
                    $check_anti_codon = $a_region . $b_region;
                }
                else {
                    $check_anti_codon = substr $sequence, $anti_codon_start-1, $anti_codon_len;
                }
            }
            $check_anti_codon = uc $check_anti_codon;
            if ($strand == -1) {
                $check_anti_codon = reverse $check_anti_codon;
                $check_anti_codon =~ tr/ATCG/TAGC/;
            }

            # Aragorn transforms ambiguous nucleotides in '.',
            # so replace any non-ATCG for '.' before comparing
            $check_anti_codon =~ s/[^ATCG]/\./g;

            if (uc $check_anti_codon ne uc $anticodon) {
                die "Problems calculating anticodon position for result number $id_num!\n";
            }

            # If an intron was found in tRNA, add the info to $note
            $note .= "; 1 intron: $intron_coords" if $intron_len > 0;

            # Add collected information to %tags
            $tags{product}   = $product;
            $tags{note}      = $note;
            $tags{anticodon} = "(pos:$anti_codon_coords,aa:$aminoacid,seq:$anticodon)";
            $tags{pseudo}    = '_no_value' if ($pseudo == 1);

        }
        elsif ($primary_tag eq 'tmRNA') {
            # Calculate tag_peptide position
            my $tag_peptide_pos = 0;
            my $tag_peptide_len = length $tag_nt_seq;
            while ($primary_seq =~ m/$tag_nt_seq/g) {
                $tag_peptide_pos = pos $primary_seq;
            }
            if ($tag_peptide_pos == 0) {
                die "Problems detecting tag peptide position in the primary seq for result number $id_num!\n";
            }

            if ($strand == -1) {
                $tag_peptide_start = $end - $tag_peptide_pos + 1;
                $tag_peptide_end   = $tag_peptide_start + $tag_peptide_len - 1;
            }
            else {
                $tag_peptide_start = $start + $tag_peptide_pos - $tag_peptide_len;
                $tag_peptide_end   = $tag_peptide_start + $tag_peptide_len - 1;
            }
            # Adjustments for split-by-origin locations,
            # where start or end could be out of bounds
            $tag_peptide_start += $seq_length if ($tag_peptide_start < 1);
            $tag_peptide_start -= $seq_length if ($tag_peptide_start > $seq_length);
            $tag_peptide_end   += $seq_length if ($tag_peptide_end   < 1);
            $tag_peptide_end   -= $seq_length if ($tag_peptide_end   > $seq_length);

            # Set of the coodinates string
            my $coord_string = '';
            if ($tag_peptide_end < $tag_peptide_start) {
                # Example: c[2298936,28]
                $coord_string = ($strand == -1) ? "complement(join($tag_peptide_start..$seq_length,1..$tag_peptide_end))"
                              :                   "join($tag_peptide_start..$seq_length,1..$tag_peptide_end)";
            }
            else {
                $coord_string = ($strand == -1) ? "complement($tag_peptide_start..$tag_peptide_end)"
                              :                   "$tag_peptide_start..$tag_peptide_end";
            }

            # Check if calculated coordinates are correct
            my $check_tag_nt_seq = '';
            if ($tag_peptide_end < $tag_peptide_start) {
                # Concatenate the sequence with a few nucleotides from
                # the beginning to emulate the passing through origin,
                # and get the sequence from there
                my $temp_seq = $sequence . substr $sequence, 0, $tag_peptide_len;
                $check_tag_nt_seq = uc substr $temp_seq, $tag_peptide_start-1, $tag_peptide_len;
            }
            else {
                $check_tag_nt_seq = uc substr $sequence, $tag_peptide_start-1, $tag_peptide_len;
            }
            if ($strand == -1) {
                $check_tag_nt_seq = reverse $check_tag_nt_seq;
                $check_tag_nt_seq =~ tr/ATCG/TAGC/;
            }

            # Aragorn transforms ambiguous nucleotides in '.',
            # so replace any non-ATCG for '.' before comparing
            $check_tag_nt_seq =~ s/[^ATCG]/\./g;

            if (uc $check_tag_nt_seq ne uc $tag_nt_seq) {
                die "Problems calculating tag peptide position for result number $id_num!\n";
            }

            # Add collected information to %tags
            $tags{gene}        = 'ssrA';
            $tags{note}        = 'tmRNA';
            $tags{product}     = 'tmRNA, 10Sa RNA';
            $tags{tag_peptide} = $coord_string;
            $tags{pseudo}      = '_no_value' if ($pseudo == 1);
        }

        # Build the feature location object
        my $feat_coord_str = '';
        if ($end < $start) {
            # Example: c[2298936,28]
            $feat_coord_str = ($strand == -1) ? "complement(join($start..$seq_length,1..$end))"
                            :                   "join($start..$seq_length,1..$end)";
        }
        else {
            $feat_coord_str = ($strand == -1) ? "complement($start..$end)"
                            :                   "$start..$end";
        }
        my $feat_loc_obj = Bio::Factory::FTLocationFactory->from_string($feat_coord_str);

        # Build and add features for special two-piece tmRNA cases
        if ($tRNA eq 'tmRNA-permuted') {
            # Calculate separate coordinates for acceptor and coding pieces
            my ($acceptor_seq, $intervening_seq, $coding_seq)
                = ($primary_seq =~ m/^( [A-Z\.]+ [a-z\.]+ )  # Acceptor piece
                                      ( [A-Z\.]+ )           # Intervening Segment (IVS)
                                      ( [a-z\.]+ [A-Z\.]+
                                        [a-z\.]+ [A-Z\.]+ )$ # Coding piece
                                     /x);

            if ($acceptor_seq eq '' or $intervening_seq eq '' or $coding_seq eq '') {
                die "  ERROR!!! Can't find relevant sequences for permuted tmRNA on result number $id_num:\n"
                  . "* tRNA         = $tRNA\n"
                  . "* Primary Seq  = $primary_seq\n";
            }

            # Double-check pieces locations if available
            if (exists $out_feats{$id_num}{'tmRNA 5'}) {
                my ($dom_5_start, $dom_5_end) = ($out_feats{$id_num}{'tmRNA 5'} =~ m/(\d+),(\d+)/);
                my ($dom_3_start, $dom_3_end) = ($out_feats{$id_num}{'tmRNA 3'} =~ m/(\d+),(\d+)/);

                if ($dom_3_end != length $acceptor_seq) {
                    die "  ERROR!!! Mismatch between permuted tmRNA sequence and "
                      . "reported 3' tRNA domain location [$dom_3_start,$dom_3_end]\n"
                }
                if ($dom_5_start != 1+length "$acceptor_seq$intervening_seq") {
                    die "  ERROR!!! Mismatch between permuted tmRNA sequence and "
                      . "reported 5' tRNA domain location [$dom_5_start,$dom_5_end]\n"
                }
            }

            my $acceptor_start = 0;
            my $acceptor_end   = 0;
            my $coding_start   = 0;
            my $coding_end     = 0;
            if ($strand == -1) {
                $coding_start   = $start;
                $coding_end     = $start      - 1 + length $coding_seq;
                $acceptor_start = $coding_end + 1 + length $intervening_seq;
                $acceptor_end   = $end;
            }
            else {
                $acceptor_start = $start;
                $acceptor_end   = $start        - 1 + length $acceptor_seq;
                $coding_start   = $acceptor_end + 1 + length $intervening_seq;
                $coding_end     = $end;
            }
            # Adjustments for split-by-origin locations,
            # where inners start or end could be out of bounds
            $acceptor_end   += $seq_length if ($acceptor_end   < 1);
            $acceptor_end   -= $seq_length if ($acceptor_end   > $seq_length);
            $coding_start   += $seq_length if ($coding_start   < 1);
            $coding_start   -= $seq_length if ($coding_start   > $seq_length);

            # Build the acceptor location object
            my $acc_coord_str = '';
            if ($acceptor_end < $acceptor_start) {
                # Example: c[2298936,28]
                $acc_coord_str = ($strand == -1) ? "complement(join($acceptor_start..$seq_length,1..$acceptor_end))"
                               :                   "join($acceptor_start..$seq_length,1..$acceptor_end)";
            }
            else {
                $acc_coord_str = ($strand == -1) ? "complement($acceptor_start..$acceptor_end)"
                               :                   "$acceptor_start..$acceptor_end";
            }
            my $acc_loc_obj = Bio::Factory::FTLocationFactory->from_string($acc_coord_str);

            # Build the acceptor location object
            my $cod_coord_str = '';
            if ($coding_end < $coding_start) {
                # Example: c[2298936,28]
                $cod_coord_str = ($strand == -1) ? "complement(join($coding_start..$seq_length,1..$coding_end))"
                               :                   "join($coding_start..$seq_length,1..$coding_end)";
            }
            else {
                $cod_coord_str = ($strand == -1) ? "complement($coding_start..$coding_end)"
                               :                   "$coding_start..$coding_end";
            }
            my $cod_loc_obj = Bio::Factory::FTLocationFactory->from_string($cod_coord_str);

            # Build features
            my $acc_feature = Bio::SeqFeature::Generic->new(-primary_tag => $primary_tag,
                                                            -location    => $acc_loc_obj,
                                                            -tag         => { gene      => 'ssrA',
                                                                              inference => 'profile:Aragorn:1.2.36',
                                                                              note      => 'tmRNA acceptor piece',
                                                                              product   => 'tmRNA acceptor piece',
                                                                             },
                                                            );
            my $cod_feature = Bio::SeqFeature::Generic->new(-primary_tag => $primary_tag,
                                                            -location    => $cod_loc_obj,
                                                            -tag         => { gene        => 'ssrA',
                                                                              inference   => 'profile:Aragorn:1.2.36',
                                                                              note        => 'tmRNA coding piece',
                                                                              product     => 'tmRNA coding piece',
                                                                              tag_peptide => $tags{tag_peptide},
                                                                             },
                                                            );
            # If detected as pseudo
            if ($pseudo == 1) {
                $acc_feature->add_tag_value('pseudo', '_no_value');
                $cod_feature->add_tag_value('pseudo', '_no_value');
            }

            # Add features by location order
            if ($strand == -1) {
                $out_obj->add_SeqFeature($cod_feature);
                $out_obj->add_SeqFeature($acc_feature);
            }
            else {
                $out_obj->add_SeqFeature($acc_feature);
                $out_obj->add_SeqFeature($cod_feature);
            }
        }
        # Build and add feature for the rest of the cases with the information collected in %tags
        else {
            my $feature = Bio::SeqFeature::Generic->new(-primary_tag => "$primary_tag",
                                                        -location    => $feat_loc_obj,
                                                        -tag         => \%tags,
                                                        );
            $out_obj->add_SeqFeature($feature);
        }

        # For tmRNA only, per INSDC recommendation
        # (http://www.insdc.org/documents/feature_table.html),
        # the amino acid sequence corresponding to the /tag_peptide
        # is additionally annotated by describing a
        # 5' partial CDS feature (e.g. CDS    <90..122)
        if ($primary_tag eq 'tmRNA') {
            my %cds_tags = (gene         => 'ssrA',
                            inference    => 'profile:Aragorn:1.2.36',
                            note         => 'tmRNA tag peptide',
                            codon_start  => 1,
                            transl_table => $codontable,
                            product      => 'tmRNA tag peptide',
                            translation  => $tag_aa_seq,
            );
            $cds_tags{pseudo} = '_no_value' if ($pseudo == 1);

            # Build location object
            my $tag_coord_str = '';
            if ($tag_peptide_end < $tag_peptide_start) {
                # Example: c[2298936,28]
                $tag_coord_str = ($strand == -1) ? "complement(join($tag_peptide_start..$seq_length,1..>$tag_peptide_end))"
                               :                   "join(<$tag_peptide_start..$seq_length,1..$tag_peptide_end)";
            }
            else {
                $tag_coord_str = ($strand == -1) ? "complement($tag_peptide_start..>$tag_peptide_end)"
                               :                   "<$tag_peptide_start..$tag_peptide_end";
            }
            my $tag_loc_obj = Bio::Factory::FTLocationFactory->from_string($tag_coord_str);

            # Build and add feature
            my $tag_pep_feat = Bio::SeqFeature::Generic->new(-primary_tag => 'CDS',
                                                             -location    => $tag_loc_obj,
                                                             -tag         => \%cds_tags,
                                                             );
            $out_obj->add_SeqFeature($tag_pep_feat);
        }
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.Aragorn.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;
    }

    # Now that all data is in the Genbank file, delete some intermediary files
    unlink $out_file or die "Could not delete '$out_file': $!\n";

    # Return to previous working folder
    chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
    return;
}

# Prediction of Horizontal Gene Transfer (HGT)
sub alien_hunter {
    my ($folder, $file, $seq_obj) = @_;
    if (not -e $folder) {
        warn "Could not find the provided folder for Alien Hunter: '$folder'\n";
        return;
    }
    if (not -e $file) {
        warn "Could not find the provided file for Alien Hunter: '$file'\n";
        return;
    }

    my $display_id = $seq_obj->display_id;
    my $seq_length = $seq_obj->length;

    # Prepare java environment needed for program optimization -c
    my $env_separator = ($^O =~ m/mswin/i) ? ';' : ':';
    my $cmd_separator = ($^O =~ m/mswin/i) ? '&' : ';';
    my @java_class_paths =("$folder",
                           "$folder/biojava/biojava-1.4_new.jar",
                           "$folder/biojava/bytecode-0.92.jar",
                           "$folder/biojava/commons-cli.jar",
                           "$folder/biojava/commons-collections-2.1.jar",
                           "$folder/biojava/commons-dbcp-1.1.jar",
                           "$folder/biojava/commons-pool-1.1.jar",
                           );
    my $java_class_string = join $env_separator, @java_class_paths;
    $java_class_string    = ($^O =~ m/mswin/i) ? ("set CLASSPATH=$java_class_string")
                          : ("export CLASSPATH=$java_class_string");

    # Keep previous working directory to return to it after program execution
    my $old_folder = getcwd;
    if ($file !~ m'[\\/]') {
        $file = "$old_folder/$file";
    }
    # Change to program directory
    chdir $folder or die "Could not change directory to '$folder': $!\n";

    # Define output
    my $out_header = ($file =~ m/^(.+)\.(?:fa|fna|fasta)$/i) ? $1 : $file;
    $out_header   .= '.AlienHunter';
    my $out_file   = $out_header;

    # Run program, setting java CLASSPATH in the same instance
    my $output = `$java_class_string $cmd_separator perl alien_hunter.pl -f "$file" -o "$out_file" -c`;

    # Prepare the Genbank object for feature additions
    my $out_obj = $seq_obj->clone;

    # Feature collector
    my @out_feats;

    # Use program output to get possible rRNA predictions
    my @possible_rRNAs = ($output =~ m/possible rRNA signature: (\S+)/g);
    foreach my $rRNA (@possible_rRNAs) {
        my $rRNA_loc_obj = Bio::Factory::FTLocationFactory->from_string($rRNA);
        my $feature = Bio::SeqFeature::Generic->new(-primary_tag => 'rRNA',
                                                    -location    => $rRNA_loc_obj,
                                                    -tag => { inference => 'profile:Alien_Hunter:1.7',
                                                              product   => '16S ribosomal RNA signature'
                                                             },
                                                    );
        push @out_feats, $feature;
    }

    # Only process results if the optimized "$out_file.opt" result file is present,
    # since Alien Hunter will not produce results for too small sequences (<20 kb),
    # insufficient data or too many ambiguous bases (>30%)
    if (-e "$out_file.opt") {
        # Alien Hunter produce an EMBL file WITHOUT an ID line,
        # which leads to BioPerl refuse of parsing it,
        # so add the missing line
        my $embl_file_old = "$out_file.opt";
        my $embl_file_new = "$out_file.embl";
        open my $OLD, '<', $embl_file_old or die "Could not read file '$embl_file_old': $!\n";
        open my $NEW, '>', $embl_file_new or die "Could not write file '$embl_file_new': $!\n";

        print $NEW "ID   $display_id; SV 1; linear; genomic DNA; STD; PRO; $seq_length BP.\n";
        while (my $line = <$OLD>) {
            print $NEW $line;
        }
        close $OLD;
        close $NEW;

        # Add genomic islands predictions
        my $embl_obj = Bio::SeqIO->new(-file   => $embl_file_new,
                                       -format => 'embl',);
        while (my $seq_obj = $embl_obj->next_seq) {
            # Alien Hunter only have "misc_feature",
            # so no need to filter features
            my @features = $seq_obj->get_SeqFeatures;
            foreach my $feat (@features) {
                # Warn if there is something else beside the misc_features
                if ( (my $primary_tag = $feat->primary_tag) ne 'misc_feature') {
                    warn "Alien Hunter output '$embl_file_old' have unexpected feature '$primary_tag'!\n"
                }
                my $feat_start = $feat->location->start;
                my $feat_end   = $feat->location->end;
                # Skip single nucleotide features (an apparent artifact of
                # the optimization over small genome-border predictions)
                next if ($feat_start == $feat_end);

                my $colour     = $feat->has_tag('colour') ? ($feat->get_tag_values('colour'))[0]
                               :                             '';
                my $note       = $feat->has_tag('note')   ? ($feat->get_tag_values('note'))[0]
                               :                             '';
                my $score      = $feat->has_tag('score')  ? ($feat->get_tag_values('score'))[0]
                               :                             '';
                my $new_note   = "score: $score; $note";
                my $feature = Bio::SeqFeature::Generic->new(-primary_tag => 'misc_feature',
                                                            -start       => $feat_start,
                                                            -end         => $feat_end,
                                                            -tag => { inference => 'profile:Alien_Hunter:1.7',
                                                                      score     => $score,
                                                                      note      => $new_note,
                                                                      colour    => $colour,
                                                                     },
                                                            );
                push @out_feats, $feature;
            }
        }
        close $embl_obj->_fh;

        # Add features sorted by a location index
        foreach my $out_feat (sort {
                                    ($a->start + (1 - $a->end / $seq_length))
                                    <=>
                                    ($b->start + (1 - $b->end / $seq_length))
                                    }
                              @out_feats
            ) {
            # Discard feature if its composed of more than 20% Ns because of
            # its inherent alteration to the sequence properties
            my $feat_start = $out_feat->location->start;
            my $feat_end   = $out_feat->location->end;
            my $check_seq  = $out_obj->subseq($feat_start, $feat_end);
            my $seq_length = length $check_seq;
            my $N_nts      = $check_seq =~ tr/nN//;
            my $N_percent  = $N_nts / $seq_length * 100;
            next if ($N_percent >= 20);

            $out_obj->add_SeqFeature($out_feat);
        }

        # Now that all data is in the Genbank file, delete some intermediary files.
        # We only keep "$out_file.sco"
        unlink  $out_file           or die "Could not delete '$out_file': $!\n";
        unlink "$out_file.plot"     or die "Could not delete '$out_file.plot': $!\n";
        unlink "$out_file.opt.plot" or die "Could not delete '$out_file.opt.plot': $!\n";
        unlink  $embl_file_old      or die "Could not delete '$embl_file_old': $!\n";
        unlink  $embl_file_new      or die "Could not delete '$embl_file_new': $!\n";
    }
    else {
        # Clean possible temporary files
        unlink  $out_file       or die "Could not delete '$out_file': $!\n"      if (-e  $out_file);
        unlink "$out_file.sco"  or die "Could not delete '$out_file.sco': $!\n"  if (-e "$out_file.sco");
        unlink "$out_file.plot" or die "Could not delete '$out_file.plot': $!\n" if (-e "$out_file.plot");
    }

    # If something was found, print object with results in a Genbank file
    if ($out_obj->get_SeqFeatures) {
        my $gbk_obj = Bio::SeqIO->new(-file   => ">$out_header.gbk",
                                      -format => 'genbank');
        $gbk_obj->write_seq($out_obj);
        $gbk_obj->close;
    }

    # Return to previous working folder
    chdir $old_folder or die "Could not change directory back to '$old_folder': $!\n";
    return;
}
