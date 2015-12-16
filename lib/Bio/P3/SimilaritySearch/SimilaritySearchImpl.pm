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

use v5.14;
no warnings 'experimental::smartmatch';

use Bio::KBase::DeploymentConfig;

use Data::Dumper;

use IPC::Run;
use AnyEvent;
use AnyEvent::Run;
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

    $self->{_next_job} = "JOB000000";
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

    #
    # If we have query parameters we are going to read our data from the input.
    # Otherwise, we pull from the parsed JSON in the input.

    my $query_fh;
    my $params;
    if ($req->query_string)
    {
	my $q = $req->parameters;
	$params = {};
	print STDERR Dumper($q);
	$params->{subject_genome} = [$q->get_all('subject_genome')];
	for my $k (qw(program output_format max_hits min_coverage evalue))
	{
	    my $v = $q->get($k);
	    $params->{$k} = $v if defined($v);
	}
	$query_fh =  $env->{'psgi.input'}
    }
    else
    {	
	my $fh = $env->{'psgi.input'};

	my $txt;
        {
	    local $/;
	    undef $/;
	    $txt = <$fh>;
	}
	
	try {
	    $params = decode_json($txt);

	} catch {
	    die "Error parsing input JSON: $_";
	};
    }

    my $ret;
	
    try {
	my $subj_genomes = $params->{subject_genome};
	$subj_genomes = [$subj_genomes] if ($subj_genomes && !ref($subj_genomes));
	my @cmd;

	my $subj_db_type;
	if ($params->{program} eq 'blastp')
	{
	    push(@cmd, 'blastp+');
	    $subj_db_type = 'a';
	}
	elsif ($params->{program} eq 'blastn')
	{
	    push(@cmd, 'blastn+');
	    $subj_db_type = 'd';
	}
	elsif ($params->{program} eq 'blastx')
	{
	    push(@cmd, 'blastx+');
	    $subj_db_type = 'a';
	}
	elsif ($params->{program} eq 'tblastn')
	{
	    push(@cmd, 'tblastn+');
	    $subj_db_type = 'd';
	}
	elsif ($params->{program} eq 'tblastx')
	{
	    push(@cmd, 'tblastx+');
	    $subj_db_type = 'd';
	}
	else
	{
	    die "Unknown program $params->{program}\n";
	}

	if (my $v = $params->{evalue})
	{
	    push(@cmd, "-evalue", $v);
	}

	if (my $v = $params->{max_hits})
	{
	    push(@cmd, "-max_hsps", $v);
	}

	if (my $v = $params->{min_coverage})
	{
	    push(@cmd, "-qcov_hsp_perc", $v);
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
	my $db_file;
	my $build_db;
	print STDERR Dumper(\@db_files);
	if (@db_files == 0)
	{
	    die "No database files found\n";
	}
	elsif (@db_files == 1)
	{
	    $db_file = $db_files[0];
	}
	elsif (@db_files > 1)
	{
	    $db_file = File::Temp->new(UNLINK => 0);
	    close($db_file);
	    $build_db = ["blastdb_aliastool+",
			 "-dblist", join(" ", @db_files),
			 "-title", join(" ", @$subj_genomes),
			 "-dbtype", ($subj_db_type eq 'a' ? 'prot' : 'nucl'),
			 "-out", $db_file];
	}   
	
	#
	# Determine output format.
	#
	my $blast_fmt = 13;
	my $content_type = "application/json";
	for ($params->{output_format})
	{
	    when("pairwise")
	    {
		$blast_fmt = 0;
		$content_type = "text/plain";
	    }
	    when ("xml")
	    {
		$blast_fmt = 5;
		$content_type = "text/xml";
	    }
	    when ("tabular")
	    {
		$blast_fmt = 6;
		$content_type = "text/plain";
	    }
	    when ("blast_archive")
	    {
		$blast_fmt = 1;
		$content_type = "text/plain";
	    }
	    when("json")
	    {
		$blast_fmt = 13;
		$content_type = "application/xml";
	    }
	}
    
	#
	# We have done our validation. Set up for async streaming output.
	#

	push(@cmd, "-db", "$db_file");
	push(@cmd, "-outfmt", $blast_fmt);

	$ret = sub {
	    my($responder) = @_;

	    my $job_id = $self->{_next_job}++;

	    my $writer = $responder->([200, ["Content-type" => $content_type]]);

	    my $start_blast = sub {
		my @running = keys %{$self->{_active_jobs}};
		my $n = @running;
		print STDERR "Start @cmd: $n active (@running)\n";

		#
		# Hook in an error handler
		#
		$writer->{handle}->on_error(sub {
		    my($h, $fatal, $error) = @_;
		    print STDERR "Handler for $job_id got error $error\n";
		    delete $self->{_active_jobs}->{$job_id};
		});

		my $h = AnyEvent::Run->new(cmd => \@cmd,
					   on_read => sub {
					       my($h) = @_;
					       $writer->write($h->rbuf);
					       $h->rbuf = '';
					   },
					   on_eof => sub {
					       my($h) = @_;
					       $writer->close();
					       delete $self->{_active_jobs}->{$job_id};
					       my @running = keys %{$self->{_active_jobs}};
					       my $n = @running;
					       print STDERR "Finish @cmd: $n active (@running)\n";
					   },
					   on_error => sub {
					       my($h, $fatal, $msg) = @_;
					       $writer->close();
					       delete $self->{_active_jobs}->{$job_id};
					   });
		
		if ($query_fh)
		{
		    #
		    # Twiggy reads the input into a memory buffer, so we
		    # just dump that into the run handle.
		    my $buf;
		    while (read($query_fh, $buf, 1048576))
		    {
			$h->push_write($buf);
		    }
		    $h->push_shutdown();
		}
		else
		{
		    $h->push_write($params->{query});
		    $h->push_shutdown();
		}
		$self->{_active_jobs}->{$job_id} = $h;
	    };

	    if ($build_db)
	    {
		my $build_h = AnyEvent::Run->new(cmd => $build_db,
						 on_read => sub {
						     my($rh) = @_;
						     print STDERR "Build says: " . $rh->rbuf . "\n";
						     $rh->rbuf = '';
						 },
						 on_eof => sub {
						     my($rh) = @_;
						     print STDERR "Build finished; starting blast\n";
						     delete $self->{_active_builds}->{$job_id};
						     &$start_blast();
						 });
		$self->{_active_builds}->{$job_id} = $build_h;
						 
	    }
	    else
	    {
		&$start_blast();
	    }
	}
    } catch {
	print STDERR "Error handling search: $_";
	$ret = [500, ['Content-Type' => 'text/plain'], ["Error handling search\n"]];
    };

    return $ret;
	   
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
