package Bio::P3::SimilaritySearch::SimilaritySearchImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

SimilaritySearch

=head1 DESCRIPTION

API Access to the p3 BLAST service.

Empty for now - initial implementation is a REST / AsyncIO service.

=cut

#BEGIN_HEADER

use Bio::KBase::DeploymentConfig;

use Data::Dumper;

use AnyEvent;
use AnyEvent::HTTP;
use Plack::Request;
use JSON::XS;
use Try::Tiny;

#
# Startup routine for the async similarity search service.
#
sub _sim_service_start
{
    my($self) = @_;

    print STDERR "Sim search service starting\n";
}

#
# Incoming async request.
#
# We support a very simple service.
#
# /search is the endpoint for our main search. We expect a POST message body with
# the request data.
#
sub _sim_request
{
    my($self, $env) = @_;
    my $req = Plack::Request->new($env);
    my $path = $req->path_info;


    if ($path =~ m,^/search$,)
    {
	print STDERR "Have search request for '$path'" . Dumper($req);
	return $self->_handle_search($env, $req);
    }
    else
    {
	return [404, ['Content-Type' => 'text/plain' ], ["Invalid path\n"]];
    }
}

sub _handle_search
{
    my($self, $env, $req) = @_;

    my $fh = $env->{'psgi.input'};

    my $txt;
    {
	local $/;
	undef $/;
	$txt = <$fh>;
    }
    print "Have $txt\n";

    my $ret;

    try {
	my $params = decode_json($txt);

	my $subj_genomes = $params->{subject_genome};
	$subj_genomes = [$subj_genomes] if ($subj_genomes && !ref($subj_genomes));

	my @cmd;
	my $subj_db_type;
	if ($params->{program} eq 'blastp')
	{
	    push(@cmd, 'blastp');
	    $subj_db_type = 'a';
	}
	elsif ($params->{program} eq 'blastn')
	{
	    push(@cmd, 'blastn');
	    $subj_db_type = 'd';
	}
	elsif ($params->{program} eq 'blastx')
	{
	    push(@cmd, 'blastx');
	    $subj_db_type = 'a';
	}
	elsif ($params->{program} eq 'tblastn')
	{
	    push(@cmd, 'tblastn');
	    $subj_db_type = 'd';
	}
	elsif ($params->{program} eq 'tblastx')
	{
	    push(@cmd, 'tblastx');
	    $subj_db_type = 'd';
	}
	else
	{
	    die "Unknown program $params->{program}\n";
	}
	
	my @db_files;
	if (ref($subj_genomes))
	{
	    for my $g (@$subj_genomes)
	    {
		my $f = "$self->{_blast_db_genomes}/$g";
		if ($subj_db_type eq 'a')
		{
		    $f .= "/$g.PATRIC.faa";
		}
		else
		{
		    $f .= "/$g.PATRIC.ffn";
		}
		if (! -f $f)
		{
		    warn "Could not find db file for $g in $f\n";
		}
		else
		{
		    push(@db_files, $f);
		}
	    }
	}
	print STDERR Dumper(\@db_files);
	if (@db_files == 0)
	{
	    die "No database files found\n";
	}
	
	$ret = [200, ['Content-Type' => 'text/plain'], \@db_files];


    } catch {
	print STDERR "Error handling search: $_";
    };

    return $ret if $ret;
    return [500, ['Content-Type' => 'text/plain'], ["Error handling search\n"]];
	   
}
    
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR

    my $cfg = Bio::KBase::DeploymentConfig->new($ENV{KB_SERVICE_NAME} || "SimilaritySearch");

    $self->{_blast_db_genomes} = $cfg->setting('blast-db-genomes');

    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=cut

1;
