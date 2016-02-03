use Bio::P3::SimilaritySearch::SimilaritySearchImpl;
use Plack::Middleware::CrossOrigin;

use strict;

my $impl = Bio::P3::SimilaritySearch::SimilaritySearchImpl->new();

$impl->_sim_service_start();

my $handler = sub { $impl->_sim_request(@_); };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");



