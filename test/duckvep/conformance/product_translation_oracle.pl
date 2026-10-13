#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Tools::CodonTable;
use Bio::Root::Version;

my $code = Bio::Tools::CodonTable->new(-id => 1);
my @cases = (
    ['post_stop', 'ATGAAATAATTT', 'ATGAAATAACTT', 1, 'MK*', 'MK*'],
    ['stop_lost', 'ATGAAATAATTTTGA', 'ATGAAACAATTTTGA', 1, 'MK*', 'MKQF*'],
    ['no_stop', 'ATGAAATAATTT', 'ATGAAACAATTT', 1, 'MK*', 'MKQF'],
    ['annotated_ctg', 'CTGAAATAA', 'CTGGAATAA', 1, 'MK*', 'ME*'],
    ['alternative_ctg', 'ATGAAACTGAAATAA', 'ACGAAACTGAAATAA', 7, 'MKLK*', 'MK*'],
);
print "BioPerl $Bio::Root::Version::VERSION\n";
for my $case (@cases) {
    my ($name, $reference, $alternate, $alternate_start, $expected_ref, $expected_alt) = @$case;
    my @products;
    for my $operand ([$reference, 1], [$alternate, $alternate_start]) {
        my ($rna, $start) = @$operand;
        die "undeclared initiator: $name\n" unless $code->is_start_codon(substr($rna, $start - 1, 3));
        my $peptide = 'M';
        for (my $offset = $start + 2; $offset + 3 <= length($rna); $offset += 3) {
            my $aa = $code->translate(substr($rna, $offset, 3));
            $peptide .= $aa;
            last if $aa eq '*';
        }
        push @products, $peptide;
    }
    die "translation mismatch: $name @products\n"
        unless $products[0] eq $expected_ref && $products[1] eq $expected_alt;
    print join("\t", $name, @products), "\n";
}
