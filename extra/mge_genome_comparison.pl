#!/usr/bin/perl
use strict;
use warnings;
use English qw(-no_match_vars);

use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;


# Check if the working directory is provided as a command-line argument
if (@ARGV != 1) {
    die "Usage: perl $PROGRAM_NAME <working_directory>\n";
}

# Get the working directory from the command-line argument
my $work_dir = $ARGV[0];

# This script assumes all sequences in the folder are going to be compared to
# the selected reference, so specify which will be the reference, while all the
# other sequences will be mapped to it
# Minimum sequence identity percentage to filter blasts alignments
my $identity_cutoff = 95;

# Minimum number of genomes aligned to support a certain combined region
my $support_cutoff = 1;
my $reference_id   = 'NC_009089';
chdir $work_dir or die "Could not change to folder '$work_dir': $ERRNO\n";

# Check which Blast package is available
my $blast_cmd = `rpsblast 2>&1`;    # Only binary present in both packages
my $blast_pkg =
    ( $blast_cmd =~ m#BLAST query/options error:# )  ? 'blast+'
  : ( $blast_cmd =~ m#rpsblast\s+\S+\s+arguments:# ) ? 'blast'
  : die "No BLAST program seems available, "
  . "please install BLAST+ or BLAST to continue!\n";

# For blast+ declare Blastn search type, so a user can change between:
# megablast (default), dc-megablast, blastn and blastn-short.
# This don't apply to the old blast package.
my $blastn_type = 'megablast';

# Genbank (.gbk) or Fasta (.fna) files are needed
opendir my $input_dir, $work_dir
  or die "Could not read folder '$work_dir': $ERRNO\n";
my @input_files = grep { /\.(?:fa|fna|fasta|gbk?)$/i } readdir $input_dir;
closedir $input_dir;

# Group together multiple filetypes of the same sequences
my %input_sequences;
foreach my $input_file (@input_files) {
    my ( $header, $ext ) = ( $input_file =~ m/^(.+)\.(\w+)$/ );
    $ext = lc $ext;

    ( $ext eq 'fa' ) ? ( $input_sequences{$header}{fasta_fna} = $input_file )
      : ( $ext eq 'fna' )
      ? ( $input_sequences{$header}{fasta_fna} = $input_file )
      : ( $ext eq 'fasta' )
      ? ( $input_sequences{$header}{fasta_fna} = $input_file )
      : ( $ext eq 'gb' )  ? ( $input_sequences{$header}{genbank} = $input_file )
      : ( $ext eq 'gbk' ) ? ( $input_sequences{$header}{genbank} = $input_file )
      :                     next;
}

# Remove Genbank files produced by the subroutines that could be already present from the hash keys
foreach my $key ( keys %input_sequences ) {
    delete $input_sequences{$key}
      if (
        $key =~
        m/(?: AlienHunter  | Aragorn    | Blastn      | CONJscan    | gaps
                       |  ICEberg_ICEs | IslandPath | keywords    | PAI-DA      | Phage_Finder
                       |  PhiSpy       | TnpPred    | tRNA_frags  | tRNAscan-SE | Consensus )/x
      );
}

# Make sure that $reference_id sequence is the first to be processed
my @sequences =
  grep { !/$reference_id/ } sort { $a cmp $b } keys %input_sequences;
unshift @sequences, $reference_id;

my %ref_genome_combined;
my $ref_genome_length = 0;
my $ref_genome_seq    = 0;
my $ref_genome_desc   = 0;
foreach my $sequence (@sequences) {
    my $compare_files_header =
      ( $blast_pkg eq 'blast+' )
      ? "$reference_id.Blastn_$blastn_type"
      : "$reference_id.Blastn";

    # Check .fna
    if (
        (
            not exists $input_sequences{$sequence}{fasta_fna}
            or -z $input_sequences{$sequence}{fasta_fna}
        )
        and exists $input_sequences{$sequence}{genbank}
      )
    {
        my $file_obj_gbk = Bio::SeqIO->new(
            -file   => "<$input_sequences{$sequence}{genbank}",
            -format => 'genbank',
        );

        while ( my $seq_obj = $file_obj_gbk->next_seq ) {

            # In NCBI genomes, primary_id is the genome GI
            my $sequence_gi = $seq_obj->primary_id;
            my $locus       = $seq_obj->display_id;
            my $ncbi_locus  = "$locus." . ( $seq_obj->version || 1 );
            my $description = $seq_obj->desc
              || $locus
              ;    # Copy $locus value if there is not description for the table
            $description =~ s/\.$//;

            my $display_id = $locus;
            if ( $sequence_gi =~ m/^\d+$/ ) {
                $display_id =
                  ( $ncbi_locus =~ m/^N[A-Z]_/ )
                  ? "gi|$sequence_gi|ref|$ncbi_locus|"
                  : "gi|$sequence_gi|gb|$ncbi_locus|";
            }
            my $seq_obj_fna = Bio::Seq->new(
                -display_id => $display_id,
                -desc       => $description,
                -seq        => $seq_obj->seq,
            );

            my $sequence_file = "$sequence.fna";
            my $file_obj_fna  = Bio::SeqIO->new(
                -file   => ">$sequence_file",
                -format => 'fasta',
                -width  => 70,
            );    # NCBI value
            $file_obj_fna->write_seq($seq_obj_fna);
            close $file_obj_fna->_fh;

            $input_sequences{$sequence}{fasta_fna} = $sequence_file;
        }
        close $file_obj_gbk->_fh;
    }

    if ( not exists $input_sequences{$sequence}{fasta_fna} ) {
        if ( $sequence eq $reference_id ) {
            die "Missing reference '$reference_id' fasta file!!! Aborting...\n";
        }
        else {
            print "Missing '$sequence' fasta file, skipping it...\n";
            next;
        }
    }
    elsif ( $sequence eq $reference_id ) {
        my $file_obj_fna = Bio::SeqIO->new(
            -file   => "<$input_sequences{$sequence}{fasta_fna}",
            -format => 'fasta',
        );
        while ( my $seq_obj = $file_obj_fna->next_seq ) {
            $ref_genome_length = $seq_obj->length;
            $ref_genome_seq    = $seq_obj->seq;
            $ref_genome_desc   = $seq_obj->desc;
        }
        close $file_obj_fna->_fh;

        print
"Genome comparison for $ref_genome_desc ($reference_id) using $identity_cutoff% identity cutoff:\n";
    }

# After this point, only non-reference sequence will continue to the blast execution
    next if ( $sequence eq $reference_id );

    # Compare the genomes
    my $extra_msg =
      ( $blast_pkg eq 'blast+' ) ? " using \"-task $blastn_type\"" : '';
    print "  Aligning '$sequence' to reference '$reference_id'$extra_msg...";

    my $reference_fasta = $input_sequences{$reference_id}{fasta_fna};
    my $query_fasta     = $input_sequences{$sequence}{fasta_fna};
    my $blast_table     = "$compare_files_header.$sequence" . "_table.txt";
    if ( not -e $blast_table ) {
        if ( $blast_pkg eq 'blast+' ) {
            system(
"blastn -query \"$query_fasta\" -subject \"$reference_fasta\" -out \"$blast_table\" -dust no -outfmt 7 -task $blastn_type"
            );
        }
        elsif ( $blast_pkg eq 'blast' ) {
            system(
"bl2seq -p blastn -i \"$query_fasta\" -j \"$reference_fasta\" -o \"$blast_table\" -F f -D 1"
            );
        }
    }

    # Read result table and store coordinates results
    my %align_coords;
    my %locus
      ; # Must know if there are multiple sequences IDs in the file (like contigs)
    open my $BLAST, '<', $blast_table
      or die "Could not read file '$blast_table': $ERRNO\n";
    while ( my $line = <$BLAST> ) {
        chomp $line;
        next if ( $line eq '' );
        next if ( $line =~ m/^#/ );

        my @columns = split /\t/, $line;
        $locus{ $columns[0] } = '';

        # Only store results with at least $identity_cutoff identity
        if ( $columns[2] >= $identity_cutoff ) {

            # Subject Start/End columns, store using start as lowest value
            if ( $columns[8] < $columns[9] ) {
                $align_coords{"$columns[8]..$columns[9]"} = $columns[0];
            }
            else {
                $align_coords{"$columns[9]..$columns[8]"} = $columns[0];
            }
        }
    }
    close $BLAST;
    my $num_locus = scalar keys %locus;

    # Join individual results to find the big regions
    my %joined_coords;
    foreach my $align_coord (
        sort {
                 ( $a =~ m/^(\d+)/ )[0] <=> ( $b =~ m/^(\d+)/ )[0]
              || ( $a =~ m/(\d+)$/ )[0] <=> ( $b =~ m/(\d+)$/ )[0]
        }
        keys %align_coords
      )
    {
        my ( $align_start, $align_end ) =
          ( $align_coord =~ m/^(\d+)\.\.(\d+)$/ );
        my $overlap = 0;
        foreach my $join_coord (
            sort { ( $a =~ m/^(\d+)/ )[0] <=> ( $b =~ m/^(\d+)/ )[0] }
            keys %joined_coords
          )
        {
            my ( $join_start, $join_end ) =
              ( $join_coord =~ m/^(\d+)\.\.(\d+)$/ );

       # Note: The +1/-1 in the conditions are to fuse together adjacent regions
            if (
                (
                        $join_start >= $align_start
                    and $join_start <= $align_end + 1
                )
                or ( $join_end >= $align_start - 1 and $join_end <= $align_end )
                or ( $join_start < $align_start and $join_end > $align_end )
              )
            {
                $overlap = 1;
                my $desc = $joined_coords{$join_coord};
                my $new_start =
                  ( $align_start < $join_start ) ? $align_start : $join_start;
                my $new_end =
                  ( $align_end > $join_end ) ? $align_end : $join_end;
                delete $joined_coords{$join_coord};

                $joined_coords{"$new_start..$new_end"} =
                  ( $num_locus > 1 )
                  ? "$desc; $align_coords{$align_coord}: $align_coord"
                  : "$desc; $align_coord";
                last;
            }
        }
        if ( $overlap == 0 ) {
            $joined_coords{$align_coord} =
              ( $num_locus > 1 )
              ? "$align_coords{$align_coord}: $align_coord"
              : "$align_coord";
        }
    }

    # Create Genbank file with the joined regions
    my $new_locus = $reference_id;
    $new_locus =~ s/\s+/_/;    # Just in case
    my $out_obj = Bio::Seq->new(
        -display_id => $new_locus,
        -desc       =>
"Genome coverage of reference '$reference_id' from '$sequence' alignments'",
        -seq => $ref_genome_seq,
    );

    my $positive = 0;
    foreach my $join_coord (
        sort { ( $a =~ m/^(\d+)/ )[0] <=> ( $b =~ m/^(\d+)/ )[0] }
        keys %joined_coords
      )
    {
        my $loc_obj = Bio::Factory::FTLocationFactory->from_string($join_coord);
        my $desc    = $joined_coords{$join_coord};
        my $note =
          ( $blast_pkg eq 'blast+' )
          ? "Blastn type: $blastn_type; $desc"
          : "Blastn; $desc";
        my $feature = Bio::SeqFeature::Generic->new(
            -primary_tag => 'misc_feature',
            -location    => $loc_obj,
            -tag         => {
                note   => $note,
                colour => '255 200 200',
            },
        );
        $out_obj->add_SeqFeature($feature);

        # Cumulative genome coverage
        my ( $join_start, $join_end ) = ( $join_coord =~ m/^(\d+)\.\.(\d+)$/ );
        $positive += ( $join_end - $join_start + 1 );

        # Combined multiple sequences regions
        my $overlap = 0;
        foreach my $combined_coord (
            sort { ( $a =~ m/^(\d+)/ )[0] <=> ( $b =~ m/^(\d+)/ )[0] }
            keys %ref_genome_combined
          )
        {
            my ( $combined_start, $combined_end ) =
              ( $combined_coord =~ m/^(\d+)\.\.(\d+)$/ );
            next if ( $combined_end < $join_start );
            last if ( $combined_start > $join_end );

            if (
                (
                        $combined_start >= $join_start
                    and $combined_start <= $join_end
                )
                or
                ( $combined_end >= $join_start and $combined_end <= $join_end )
                or
                ( $combined_start < $join_start and $combined_end > $join_end )
              )
            {
                $overlap = 1;
                my %sequences = %{ $ref_genome_combined{$combined_coord} };
                delete $ref_genome_combined{$combined_coord};

       # All possible overlap conditions. The original region will be fragmented
       # depending on the configuration of the overlap
                if (    $join_start < $combined_start
                    and $join_end < $combined_end )
                {    # 1  a  2  b
                    my $aux_coord1 = $combined_start - 1;
                    my $aux_coord2 = $join_end + 1;
                    $ref_genome_combined{"$join_start..$aux_coord1"}{$sequence}
                      = '';
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$combined_start..$join_end"}{$key}
                          = '';
                    }
                    $ref_genome_combined{"$combined_start..$join_end"}
                      {$sequence} = '';
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$aux_coord2..$combined_end"}{$key}
                          = '';
                    }
                }
                elsif ( $join_start < $combined_start
                    and $join_end == $combined_end )
                {    # 1  a  2b
                    my $aux_coord1 = $combined_start - 1;
                    $ref_genome_combined{"$join_start..$aux_coord1"}{$sequence}
                      = '';
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$combined_start..$join_end"}{$key}
                          = '';
                    }
                    $ref_genome_combined{"$combined_start..$join_end"}
                      {$sequence} = '';
                }
                elsif ( $join_start < $combined_start
                    and $join_end > $combined_end )
                {    # 1  a  b  2
                    my $aux_coord1 = $combined_start - 1;
                    my $aux_coord2 = $combined_end + 1;
                    $ref_genome_combined{"$join_start..$aux_coord1"}{$sequence}
                      = '';
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$combined_start..$combined_end"}
                          {$key} = '';
                    }
                    $ref_genome_combined{"$combined_start..$combined_end"}
                      {$sequence} = '';
                    $ref_genome_combined{"$aux_coord2..$join_end"}{$sequence} =
                      '';
                }
                elsif ( $join_start == $combined_start
                    and $join_end < $combined_end )
                {    # a1  2  b
                    my $aux_coord1 = $join_end + 1;
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$join_start..$join_end"}{$key} =
                          '';
                    }
                    $ref_genome_combined{"$join_start..$join_end"}{$sequence} =
                      '';
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$aux_coord1..$combined_end"}{$key}
                          = '';
                    }
                }
                elsif ( $join_start == $combined_start
                    and $join_end == $combined_end )
                {    # a1  2b
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$join_start..$join_end"}{$key} =
                          '';
                    }
                    $ref_genome_combined{"$join_start..$join_end"}{$sequence} =
                      '';
                }
                elsif ( $join_start == $combined_start
                    and $join_end > $combined_end )
                {    # a1  b  2
                    my $aux_coord1 = $combined_end + 1;
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$join_start..$combined_end"}{$key}
                          = '';
                    }
                    $ref_genome_combined{"$join_start..$combined_end"}
                      {$sequence} = '';
                    $ref_genome_combined{"$aux_coord1..$join_end"}{$sequence} =
                      '';
                }
                elsif ( $join_start > $combined_start
                    and $join_end < $combined_end )
                {    # a  1  2  b
                    my $aux_coord1 = $join_start - 1;
                    my $aux_coord2 = $join_end + 1;
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$combined_start..$aux_coord1"}
                          {$key} = '';
                    }
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$join_start..$join_end"}{$key} =
                          '';
                    }
                    $ref_genome_combined{"$join_start..$join_end"}{$sequence} =
                      '';
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$aux_coord2..$combined_end"}{$key}
                          = '';
                    }
                }
                elsif ( $join_start > $combined_start
                    and $join_end == $combined_end )
                {    # a  1  2b
                    my $aux_coord1 = $join_start - 1;
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$combined_start..$aux_coord1"}
                          {$key} = '';
                    }
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$join_start..$join_end"}{$key} =
                          '';
                    }
                    $ref_genome_combined{"$join_start..$join_end"}{$sequence} =
                      '';
                }
                elsif ( $join_start > $combined_start
                    and $join_end > $combined_end )
                {    # a  1  b  2
                    my $aux_coord1 = $join_start - 1;
                    my $aux_coord2 = $combined_end + 1;
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$combined_start..$aux_coord1"}
                          {$key} = '';
                    }
                    foreach my $key ( keys %sequences ) {
                        $ref_genome_combined{"$join_start..$combined_end"}{$key}
                          = '';
                    }
                    $ref_genome_combined{"$join_start..$combined_end"}
                      {$sequence} = '';
                    $ref_genome_combined{"$aux_coord2..$join_end"}{$sequence} =
                      '';
                }
                else {
                    die
"Something went very wrong with the overlap conditions!\n";
                }
                last;
            }
        }
        if ( $overlap == 0 ) {
            $ref_genome_combined{$join_coord}{$sequence} = '';
        }
    }

    my $out_file_gbk =
      "$compare_files_header.ident_$identity_cutoff.$sequence.gbk";
    my $out_file_obj = Bio::SeqIO->new(
        -file   => ">$out_file_gbk",
        -format => 'genbank',
    );
    $out_file_obj->write_seq($out_obj);
    close $out_file_obj->_fh;

    # Joined regions coverage of the reference genome
    my $percent = $positive / $ref_genome_length * 100;
    $percent = sprintf "%.2f", $percent;
    print "Covered area: $percent%\n";
}

# If there are more than 2 sequences being compared, print combined regions from all sequences
if ( scalar @sequences > 2 ) {
    print "  Combining and fragmenting overlapping features...\n";

    # Get how many times each Start value is repeated, and store
    # the highest repeat value to speed up the next hash self comparison
    my %repeated_keys;
    my $max_repeat = 0;
    foreach my $key ( keys %ref_genome_combined ) {
        my ($key_start) = ( $key =~ m/^(\d+)/ );
        ( exists $repeated_keys{$key_start} )
          ? ( $repeated_keys{$key_start}++ )
          : ( $repeated_keys{$key_start} = 1 );
    }
    foreach my $key ( keys %repeated_keys ) {
        $max_repeat = $repeated_keys{$key}
          if ( $max_repeat < $repeated_keys{$key} );
    }

    # Depuration of %ref_genome_combined, make sure that no features overlap
    my $clean_flag = 0;
    while ( $clean_flag == 0 ) {
        my $overlap_flag = 0;
        my @sorted_keys =
          sort { ( $a =~ m/^(\d+)/ )[0] <=> ( $b =~ m/^(\d+)/ )[0] }
          keys %ref_genome_combined;
      CYCLE:
        for ( my $i = 0 ; $i < scalar @sorted_keys ; $i++ ) {
            my $combined_coord1 = $sorted_keys[$i];

            my ( $combined_start1, $combined_end1 ) =
              ( $combined_coord1 =~ m/^(\d+)\.\.(\d+)$/ );
            my %sequences1 = %{ $ref_genome_combined{$combined_coord1} };

            # No need to read the whole hash, only start to read from
            # $max_repeat positions before the current hash key
            my $position = ( $i - $max_repeat < 0 ) ? 0 : ( $i - $max_repeat );

            for ( my $j = $position ; $j < scalar @sorted_keys ; $j++ ) {
                my $combined_coord2 = $sorted_keys[$j];
                next if ( $combined_coord2 eq $combined_coord1 );

                my ( $combined_start2, $combined_end2 ) =
                  ( $combined_coord2 =~ m/^(\d+)\.\.(\d+)$/ );
                next if ( $combined_end2 < $combined_start1 );
                last if ( $combined_start2 > $combined_end1 );

                if (
                    (
                            $combined_start2 >= $combined_start1
                        and $combined_start2 <= $combined_end1
                    )
                    or (    $combined_end2 >= $combined_start1
                        and $combined_end2 <= $combined_end1 )
                    or (    $combined_start2 < $combined_start1
                        and $combined_end2 > $combined_end1 )
                  )
                {
                    $overlap_flag++;
                    my %sequences2 =
                      %{ $ref_genome_combined{$combined_coord2} };
                    my %sequences_all = ( %sequences1, %sequences2 );
                    delete $ref_genome_combined{$combined_coord1};
                    delete $ref_genome_combined{$combined_coord2};

# Possible overlap conditions, omitted some that cannot happen with current code logic
# (hash compare against itself, sorted by start and fragmentation). The original region
# will be fragmented depending on the configuration of the overlap.
                    if (    $combined_start1 < $combined_start2
                        and $combined_end1 < $combined_end2 )
                    {    # 1  a  2  b
                        my $aux_coord1 = $combined_start2 - 1;
                        my $aux_coord2 = $combined_end1 + 1;
                        foreach my $key ( keys %sequences1 ) {
                            $ref_genome_combined{
                                "$combined_start1..$aux_coord1"}{$key} = '';
                        }
                        foreach my $key ( keys %sequences_all ) {
                            $ref_genome_combined{
                                "$combined_start2..$combined_end1"}{$key} = '';
                        }
                        foreach my $key ( keys %sequences2 ) {
                            $ref_genome_combined{"$aux_coord2..$combined_end2"}
                              {$key} = '';
                        }
                    }
                    elsif ( $combined_start1 < $combined_start2
                        and $combined_end1 == $combined_end2 )
                    {    # 1  a  2b
                        my $aux_coord1 = $combined_start2 - 1;
                        foreach my $key ( keys %sequences1 ) {
                            $ref_genome_combined{
                                "$combined_start1..$aux_coord1"}{$key} = '';
                        }
                        foreach my $key ( keys %sequences_all ) {
                            $ref_genome_combined{
                                "$combined_start2..$combined_end1"}{$key} = '';
                        }
                    }
                    elsif ( $combined_start1 < $combined_start2
                        and $combined_end1 > $combined_end2 )
                    {    # 1  a  b  2
                        my $aux_coord1 = $combined_start2 - 1;
                        my $aux_coord2 = $combined_end2 + 1;
                        foreach my $key ( keys %sequences1 ) {
                            $ref_genome_combined{
                                "$combined_start1..$aux_coord1"}{$key} = '';
                        }
                        foreach my $key ( keys %sequences_all ) {
                            $ref_genome_combined{
                                "$combined_start2..$combined_end2"}{$key} = '';
                        }
                        foreach my $key ( keys %sequences1 ) {
                            $ref_genome_combined{"$aux_coord2..$combined_end1"}
                              {$key} = '';
                        }
                    }
                    elsif ( $combined_start1 == $combined_start2
                        and $combined_end1 < $combined_end2 )
                    {    # a1  2  b
                        my $aux_coord1 = $combined_end1 + 1;
                        foreach my $key ( keys %sequences_all ) {
                            $ref_genome_combined{
                                "$combined_start1..$combined_end1"}{$key} = '';
                        }
                        foreach my $key ( keys %sequences2 ) {
                            $ref_genome_combined{"$aux_coord1..$combined_end2"}
                              {$key} = '';
                        }
                    }
                    elsif ( $combined_start1 == $combined_start2
                        and $combined_end1 > $combined_end2 )
                    {    # a1  b  2
                        my $aux_coord1 = $combined_end2 + 1;
                        foreach my $key ( keys %sequences_all ) {
                            $ref_genome_combined{
                                "$combined_start1..$combined_end2"}{$key} = '';
                        }
                        foreach my $key ( keys %sequences1 ) {
                            $ref_genome_combined{"$aux_coord1..$combined_end1"}
                              {$key} = '';
                        }
                    }
                    else {
                        die
"Something went very wrong with the overlap conditions!\n";
                    }

             # Since this code add and delete keys on the hash that is checking,
             # it must start again after each modification
                    last CYCLE;
                }
            }
        }

        $clean_flag = 1 if ( $overlap_flag == 0 );
    }

    # After %ref_genome_combined cleanup, remove any regions that
    # are supported by less than $support_cutoff genomes
    if ( $support_cutoff > 1 ) {
        foreach my $combined_coord (
            sort { ( $a =~ m/^(\d+)/ )[0] <=> ( $b =~ m/^(\d+)/ )[0] }
            keys %ref_genome_combined
          )
        {
            my %sequences = %{ $ref_genome_combined{$combined_coord} };
            my $number    = scalar keys %sequences;

            if ( $number < $support_cutoff ) {
                delete $ref_genome_combined{$combined_coord};
            }
        }
    }

    # After overlap resolution and filter, fuse any remaining
    # blocks that are adjacent and identical (covered by exactly
    # the same sequences)
    my $ready_flag = 0;
    while ( $ready_flag == 0 ) {
        my $fusion_flag = 0;
        my @sorted_keys =
          sort { ( $a =~ m/^(\d+)/ )[0] <=> ( $b =~ m/^(\d+)/ )[0] }
          keys %ref_genome_combined;
      CYCLE:

        # Since there are no overlaps, compare each block only with the next one
        for ( my $i = 0 ; $i < scalar(@sorted_keys) - 1 ; $i++ ) {
            my $combined_coord1 = $sorted_keys[$i];
            my $combined_coord2 = $sorted_keys[ $i + 1 ];

            my ( $combined_start1, $combined_end1 ) =
              ( $combined_coord1 =~ m/^(\d+)\.\.(\d+)$/ );
            my ( $combined_start2, $combined_end2 ) =
              ( $combined_coord2 =~ m/^(\d+)\.\.(\d+)$/ );

            if ( $combined_end1 + 1 == $combined_start2 ) {
                my %sequences1 = %{ $ref_genome_combined{$combined_coord1} };
                my %sequences2 = %{ $ref_genome_combined{$combined_coord2} };
                my $num_seqs1  = scalar keys %sequences1;
                my $num_seqs2  = scalar keys %sequences2;
                next if ( $num_seqs1 != $num_seqs2 );

                my $matches = 0;
                foreach my $key ( keys %sequences2 ) {
                    $matches++ if ( exists $sequences1{$key} );
                }
                next if ( $matches != $num_seqs1 );

                $fusion_flag++;
                delete $ref_genome_combined{$combined_coord1};
                delete $ref_genome_combined{$combined_coord2};

                %{ $ref_genome_combined{"$combined_start1..$combined_end2"} } =
                  map { $ARG => '' } keys %sequences1;

             # Since this code add and delete keys on the hash that is checking,
             # it must start again after each modification
                last CYCLE;
            }
        }
        $ready_flag = 1 if ( $fusion_flag == 0 );
    }

    print "  Writing all sequences combined file...";
    shift @sequences;    # Remove the reference ID from the top of the list
    my $sequences = join ", ", @sequences;
    my $new_locus = $reference_id;
    $new_locus =~ s/\s+/_/;    # Just in case
    my $out_obj = Bio::Seq->new(
        -display_id => $new_locus,
        -desc => "Genome coverage of reference '$reference_id' from $sequences",
        -seq  => $ref_genome_seq,
    );

    my $combined_length = 0;
    foreach my $combined_coord (
        sort { ( $a =~ m/^(\d+)/ )[0] <=> ( $b =~ m/^(\d+)/ )[0] }
        keys %ref_genome_combined
      )
    {
        my %sequences   = %{ $ref_genome_combined{$combined_coord} };
        my $genomes     = join ", ", sort { $a cmp $b } keys %sequences;
        my $num_genomes = scalar keys %sequences;
        my $loc_obj =
          Bio::Factory::FTLocationFactory->from_string($combined_coord);
        my $method =
          ( $blast_pkg eq 'blast+' ) ? "Blastn type: $blastn_type" : 'Blastn';
        my $note = "$num_genomes supporting genomes; $genomes; $method";

        $combined_length += $loc_obj->end() - $loc_obj->start() + 1;

    # Add color grading to feature depending of the amount of supporting genomes
        my $max_num = scalar keys %input_sequences;
        my $value   = $num_genomes / $max_num;
        my $red     = 255;
        my $green   = int( 255 - ( $value * 255 ) );
        my $blue    = int( 255 - ( $value * 255 ) );
        my $colour  = "$red $green $blue";

        my $feature = Bio::SeqFeature::Generic->new(
            -primary_tag => 'misc_feature',
            -location    => $loc_obj,
            -tag         => {
                note   => $note,
                colour => $colour,
            },
        );
        $out_obj->add_SeqFeature($feature);
    }
    my $out_file_gbk =
      ( $blast_pkg eq 'blast+' )
      ? "$reference_id.Blastn_$blastn_type.$support_cutoff"
      . "_combined.ident_$identity_cutoff.gbk"
      : "$reference_id.Blastn.$support_cutoff"
      . "_combined.ident_$identity_cutoff.gbk";
    my $out_file_obj = Bio::SeqIO->new(
        -file   => ">$out_file_gbk",
        -format => 'genbank',
    );
    $out_file_obj->write_seq($out_obj);
    close $out_file_obj->_fh;

    my $percent = $combined_length / $ref_genome_length * 100;
    $percent = sprintf "%.2f", $percent;
    print "Area covered with at least $support_cutoff sequence(s): $percent%\n";
}

print "\nFinished!\n";
exit;
