#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_module_test.t - Are required modules present           |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 12/18/2007                                       |
# UPDATED: 04/14/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test to see if modules required by DAWG-PAWS are         |
#  present.                                                 |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Test::More tests => 20;

diag ("Testing required modules ...\n");

# Need to check on where I am using this
use_ok( 'Carp' );

#-----------------------------+
# Getopt::Long                |
#-----------------------------+
# Used for command line interface
use_ok( 'Getopt::Long' );
# The following is being retired for Getopt::Long 
use_ok( 'Getopt::Std' );

#-----------------------------+
# OS INDEPENDENCE             |
#-----------------------------+
use_ok( 'File::Copy' );        # Copy files without system cmds
use_ok( 'Cwd' );               # Get the current working directory

#-----------------------------+
# MODULES FOR print_help      |
#-----------------------------+
# Modules used for print help subfunction 
use_ok( 'Pod::Select' );       # Print subsections of POD documentation
use_ok( 'Pod::Text' );         # Print POD doc as formatted text file
use_ok( 'IO::Scalar' );
use_ok( 'IO::Pipe' );
use_ok( 'File::Spec' );

#-----------------------------+
# BIOPERL                     |
#-----------------------------+
# Bioperl packages use throughout
use_ok( 'Bio::SearchIO' );
use_ok( 'Bio::SeqIO' );
use_ok( 'Bio::Tools::Genemark' );
use_ok( 'Bio::Tools::HMMER::Results' );
use_ok( 'Bio::Tools::Run::StandAloneBlast' );
use_ok( 'Bio::Location::Simple' );

#-----------------------------+
# CONNECTIVITY                |
#-----------------------------+
# used for TE-NEST fetching
use_ok('LWP::Simple');
use_ok('LWP::UserAgent');
# used for push_file
#use_ok('Net::SFTP');

#-----------------------------+
# PRETTY PRINTING             |
#-----------------------------+
use_ok( 'Text::Wrap' );

#-----------------------------+
# Dont Know why               |
#-----------------------------+
use_ok( 'utf8' );
