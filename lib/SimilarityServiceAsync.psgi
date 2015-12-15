use Bio::P3::SimilaritySearch::SimilaritySearchImpl;

use strict;

my $impl = Bio::P3::SimilaritySearch::SimilaritySearchImpl->new();

$impl->_sim_service_start();

sub { $impl->_sim_request(@_); };


