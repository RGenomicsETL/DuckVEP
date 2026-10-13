#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::SeqEdit;

# Raw exon-spliced sequences in transcript orientation, including minus-strand RNA.
my @cases = (
    ['plus',  'ATGAAACCCTAA', 5, 5, 'G', 'ATGAGACCCTAA'],
    ['minus', 'TTAGGGTTTCAT', 6, 6, 'A', 'TTAGGATTTCAT'],
    ['insertion', 'ATGAAACCCTAA', 5, 4, 'T', 'ATGATAACCCTAA'],
);
my %vertebrate_mitochondrial = (ATG => 'M', TTC => 'F', TAA => '*', CAA => 'Q', AGG => '*');
sub peptide {
    my ($rna) = @_;
    my $result = '';
    for (my $at = 0; $at + 2 < length($rna); $at += 3) {
        my $codon = substr($rna, $at, 3);
        die "unknown codon $codon\n" unless exists $vertebrate_mitochondrial{$codon};
        my $residue = $vertebrate_mitochondrial{$codon};
        $result .= $residue;
        last if $residue eq '*';
    }
    return $result;
}

for my $case (@cases) {
    my ($name, $raw, $start, $end, $alt, $expected) = @$case;
    my $edit = Bio::EnsEMBL::SeqEdit->new(
        -START => $start, -END => $end, -ALT_SEQ => $alt,
        -CODE => '_rna_edit', -NAME => 'RNA Edit'
    );
    my $edited = $raw;
    $edit->apply_edit(\$edited);
    die "$name: $edited != $expected\n" unless $edited eq $expected;
    print join("\t", $name, $raw, $start, $end, $alt, $edited), "\n";
}

my $raw = 'ATGCCCTAAAGG';
my $programme = Bio::EnsEMBL::SeqEdit->new(
    -START => 4, -END => 5, -ALT_SEQ => 'TT', -CODE => '_rna_edit'
);
my $reference = $raw;
$programme->apply_edit(\$reference);
my $variant = $raw;
substr($variant, 6, 1, 'C');
$programme->apply_edit(\$variant);
die "reference RNA and terminator\n" unless $reference eq 'ATGTTCTAAAGG' && peptide($reference) eq 'MF*';
die "variant RNA and reached terminator\n" unless $variant eq 'ATGTTCCAAAGG' && peptide($variant) eq 'MFQ*';
print join("\t", 'code2-reached', $reference, peptide($reference), $variant, peptide($variant)), "\n";

my $start = 'GTGCCCTAAAGG';
for my $op ([1, 1, 'A'], [4, 5, 'TT']) {
    Bio::EnsEMBL::SeqEdit->new(
        -START => $op->[0], -END => $op->[1], -ALT_SEQ => $op->[2], -CODE => '_rna_edit'
    )->apply_edit(\$start);
}
die "programme-supplied start\n" unless $start eq $reference;
print join("\t", 'code2-declared-start', $start, peptide($start)), "\n";
