#!/usr/bin/perl
use strict;
use warnings;
use English qw(-no_match_vars);

use Bio::Seq;
use Bio::SeqIO;


# Check if the working directory is provided as a command-line argument
if (@ARGV != 1) {
    die "Usage: perl $PROGRAM_NAME <working_directory>\n";
}

# Get the working directory from the command-line argument
my $work_dir = $ARGV[0];


# Reference Genbank file
my $locus         = 'NC_009089';
my $reference_gbk = "$locus.gbk";

# "mge_genome_comparison" output file with the consensus regions of the reference
# after sequence alignments with other species from the same genus
my $combined_gbk = "$work_dir/$locus.Blastn_megablast.1_combined.ident_95.gbk";

# my $prediction_gbk = "$work_dir/$locus.Consensus.gbk";
my $max_support = 1;

chdir $work_dir or die "Could not change folder to '$work_dir': $ERRNO\n";

# Look for the Genbank files of the Consensus (MOBHunter output), and the other
# tools results for the analysis
opendir my $input_dir, $work_dir
  or die "Could not read folder '$work_dir': $ERRNO\n";
my @gbk_files = grep {
    m/(?: AlienHunter | IslandPath | PAI-DA    | Phage_Finder | PhiSpy
       |  GIHunter    | GIPSy      | GIST_EGID | Consensus)\.gbk$/x
} readdir $input_dir;
closedir $input_dir;

# Make MOBHunter file the first of the list, and the rest in alphabetical order
my ($mobhunter_file) = grep { /Consensus/ } @gbk_files;
@gbk_files = grep { !/Consensus/ } sort { $a cmp $b } @gbk_files;
unshift @gbk_files, $mobhunter_file;

# Read reference annotation and store genes
print "Reading reference '$reference_gbk'...\n";
my %cds_data;
my $cds_number  = 0;
my $ref_gbk_obj = Bio::SeqIO->new(
    -file   => "<$reference_gbk",
    -format => 'genbank',
);
while ( my $seq_obj = $ref_gbk_obj->next_seq ) {
    my $genome_len = $seq_obj->length;
    my @feats = grep { $_->primary_tag eq 'CDS' } $seq_obj->get_SeqFeatures;

    foreach my $feat_obj (@feats) {
        $cds_number++;
        my $location  = $feat_obj->location->to_FTstring;
        my $cds_start = $feat_obj->start;
        my $cds_end   = $feat_obj->end;
        my $loc_index = $cds_start + ( 1 - $cds_end / $genome_len );

        $cds_data{$location}{start}   = $cds_start;
        $cds_data{$location}{end}     = $cds_end;
        $cds_data{$location}{loc_idx} = $loc_index;
    }
}
close $ref_gbk_obj->_fh;

# Read genome comparison combined file to associate genes to conserved regions
print "Reading combined file '$combined_gbk'...\n";
my $comb_gbk_obj = Bio::SeqIO->new(
    -file   => "<$combined_gbk",
    -format => 'genbank',
);
while ( my $seq_obj = $comb_gbk_obj->next_seq ) {
    my @feats =
      grep { $_->primary_tag eq 'misc_feature' } $seq_obj->get_SeqFeatures;
    foreach my $feat_obj (@feats) {
        my $feat_start = $feat_obj->start;
        my $feat_end   = $feat_obj->end;
        my $note =
          $feat_obj->has_tag('note')
          ? ( $feat_obj->get_tag_values('note') )[0]
          : '';

        # "Or 1" is used when looking individual genomes instead of combined
        my ($support) = ( $note =~ m/^(\d+) supporting genomes;/ )[0] || 1;
        foreach my $cds_loc (
            sort { $cds_data{$a}{loc_idx} <=> $cds_data{$b}{loc_idx} }
            keys %cds_data
          )
        {
            my $cds_start = $cds_data{$cds_loc}{start};
            my $cds_end   = $cds_data{$cds_loc}{end};
            next if $cds_end < $feat_start;
            last if $cds_start > $feat_end;

            if (   ( $cds_start >= $feat_start and $cds_start <= $feat_end )
                or ( $cds_end >= $feat_start  and $cds_end <= $feat_end )
                or ( $cds_start < $feat_start and $cds_end > $feat_end ) )
            {
                my $inside_length =
                    ( $cds_start >= $feat_start and $cds_end <= $feat_end )
                  ? ( $cds_end - $cds_start + 1 )      # Complete smaller
                  : ( $cds_start < $feat_start and $cds_end > $feat_end )
                  ? ( $feat_end - $feat_start + 1 )    # Complete bigger
                  : ( $cds_start >= $feat_start and $cds_start <= $feat_end )
                  ? ( $feat_end - $cds_start + 1 )     # 5' inside
                  : ( $cds_end - $feat_start + 1 );    # 3' inside
                my $cds_inside_percent =
                  $inside_length / ( $cds_end - $cds_start + 1 ) * 100;
                if ( $cds_inside_percent >= 50 ) {
                    $cds_data{$cds_loc}{conserved} = $support;
                }
            }
        }
    }
}
close $comb_gbk_obj->_fh;

# Read prediction file to gather island data
print "\nProgram\tSupport level\tTrue Pos\tFalse Pos\tTrue Neg\tFalse Neg\t"
  . "Predictions\tSensitivity\tSpecificity\tPrecision\tAccuracy\t"
  . "Matthews correlation coefficient\n";
foreach my $prediction_gbk (@gbk_files) {
    my ($tool_name) = ( $prediction_gbk =~ m/([^.]+)\.gbk$/ );
    if ( $tool_name eq 'Consensus' ) {
        $tool_name = 'MOBHunter';
    }
    my %island_data;
    my $pred_gbk_obj = Bio::SeqIO->new(
        -file   => "<$prediction_gbk",
        -format => 'genbank',
    );
    while ( my $seq_obj = $pred_gbk_obj->next_seq ) {
        my $genome_len = $seq_obj->length;
        my @feats =
          grep { $_->primary_tag eq 'misc_feature' } $seq_obj->get_SeqFeatures;

        foreach my $feat_obj (@feats) {
            my $location  = $feat_obj->location->to_FTstring;
            my $cds_start = $feat_obj->start;
            my $cds_end   = $feat_obj->end;
            my $loc_index = $cds_start + ( 1 - $cds_end / $genome_len );

            $island_data{$location}{start}   = $cds_start;
            $island_data{$location}{end}     = $cds_end;
            $island_data{$location}{loc_idx} = $loc_index;
        }
    }
    close $pred_gbk_obj->_fh;

  # Compare island data with gene data to measure True/False Positives/Negatives
    my %data;
    foreach my $cds_loc (
        sort { $cds_data{$a}{loc_idx} <=> $cds_data{$b}{loc_idx} }
        keys %cds_data
      )
    {
        my $cds_start = $cds_data{$cds_loc}{start};
        my $cds_end   = $cds_data{$cds_loc}{end};
        my $support =
          ( exists $cds_data{$cds_loc}{conserved} )
          ? $cds_data{$cds_loc}{conserved}
          : 0;

        my $match = 0;
        foreach my $island_loc (
            sort { $island_data{$a}{loc_idx} <=> $island_data{$b}{loc_idx} }
            keys %island_data
          )
        {
            my $island_start = $island_data{$island_loc}{start};
            my $island_end   = $island_data{$island_loc}{end};
            next if $island_end < $cds_start;
            last if $island_start > $cds_end;

            if (   ( $island_start >= $cds_start and $island_start <= $cds_end )
                or ( $island_end >= $cds_start  and $island_end <= $cds_end )
                or ( $island_start < $cds_start and $island_end > $cds_end ) )
            {
                my $inside_length =
                    ( $island_start >= $cds_start and $island_end <= $cds_end )
                  ? ( $island_end - $island_start + 1 )    # Complete smaller
                  : ( $island_start < $cds_start and $island_end > $cds_end )
                  ? ( $cds_end - $cds_start + 1 )          # Complete bigger
                  : (     $island_start >= $cds_start
                      and $island_start <= $cds_end )
                  ? ( $cds_end - $island_start + 1 )       # 5' inside
                  : ( $island_end - $cds_start + 1 );      # 3' inside
                my $cds_inside_percent =
                  $inside_length / ( $cds_end - $cds_start + 1 ) * 100;
                if ( $cds_inside_percent >= 50 ) {
                    $match = 1;
                    last;
                }
            }
        }

        # Update %data
        if ( $match == 1 ) {
            for my $cycle ( 1 .. $max_support ) {
                if ( $cycle > $support ) {
                    $data{$cycle}{TP}++;
                }
                else {
                    $data{$cycle}{FP}++;
                }
            }
        }
        else {
            for my $cycle ( 1 .. $max_support ) {
                if ( $cycle > $support ) {
                    $data{$cycle}{FN}++;
                }
                else {
                    $data{$cycle}{TN}++;
                }
            }
        }
    }

    # Calculate sensitivity and specificity by support level
    foreach my $supp_level ( sort { $a <=> $b } keys %data ) {
        my $true_pos    = $data{$supp_level}{TP} || 0;
        my $false_pos   = $data{$supp_level}{FP} || 0;
        my $true_neg    = $data{$supp_level}{TN} || 0;
        my $false_neg   = $data{$supp_level}{FN} || 0;
        my $predictions = $true_pos + $false_pos;
        my $total_data  = $true_pos + $false_pos + $true_neg + $false_neg;
        if ( $total_data != $cds_number ) {
            die
              "\nERROR: Counted genes mismatch!!! $total_data vs $cds_number\n";
        }

        my $sensitivity =
          ( $true_pos + $false_neg == 0 )
          ? 0
          : ( $true_pos / ( $true_pos + $false_neg ) * 100 );
        my $specificity =
          ( $true_neg + $false_pos == 0 )
          ? 0
          : ( $true_neg / ( $true_neg + $false_pos ) * 100 );
        my $precision =
          ( $true_pos + $false_pos == 0 )
          ? 0
          : ( $true_pos / ( $true_pos + $false_pos ) * 100 );
        my $accuracy = ( $true_pos + $true_neg ) / $total_data * 100;

        # Matthews correlation coefficient
        my $denominator =
            ( $true_pos + $false_pos == 0 ) ? 1
          : ( $true_pos + $false_neg == 0 ) ? 1
          : ( $true_neg + $false_pos == 0 ) ? 1
          : ( $true_neg + $false_neg == 0 ) ? 1
          : ( ( $true_pos + $false_pos ) *
              ( $true_pos + $false_neg ) *
              ( $true_neg + $false_pos ) *
              ( $true_neg + $false_neg ) );
        my $matthews_cc =
          ( ( $true_pos * $true_neg ) - ( $false_pos * $false_neg ) ) /
          $denominator**0.5;

        $sensitivity = sprintf( "%.3f", $sensitivity );
        $specificity = sprintf( "%.3f", $specificity );
        $accuracy    = sprintf( "%.3f", $accuracy );
        $precision   = sprintf( "%.3f", $precision );
        $matthews_cc = sprintf( "%.3f", $matthews_cc );

        my $data_line = join( "\t",
            $tool_name,   $supp_level, $true_pos,    $false_pos,
            $true_neg,    $false_neg,  $predictions, $sensitivity,
            $specificity, $precision,  $accuracy,    $matthews_cc,
        );
        print "$data_line\n";
    }
}
print "\nFinished!\n";
exit;
