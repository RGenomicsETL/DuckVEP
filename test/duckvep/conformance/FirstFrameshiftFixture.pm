package FirstFrameshiftFixture;
use strict;
use warnings;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF;

my $create_transcript = \&Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF::_create_transcript;
{
  no warnings 'redefine';
  *Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF::_create_transcript = sub {
    my $transcript = $create_transcript->(@_);
    if ($transcript && $transcript->stable_id =~ /^TX[1-6]$/) {
      $transcript->add_Attributes(Bio::EnsEMBL::Attribute->new(
        -CODE => 'cds_start_NF', -VALUE => 1));
    }
    return $transcript;
  };
}
1;
