#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# extract_pod.pl - Extract POD to html.                     |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/07/2009                                       |
# UPDATED: 04/07/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Extract the POD documentation from the scripts to html   |
#                                                           |
#-----------------------------------------------------------+

package DAWGPAWS;

my $outdir = "/Users/jestill/code/dawgpaws/branches/release-1.0/docs/html/";
my $css_location = "../css/dawgpaws_man.css";

# boolean to add the google analytics code to the resulting html
$do_analytics = 0;

my $analytics_code = "<script ".
    "src=\"http://www.google-analytics.com/urchin.js\"\n".
    "type=\"text/javascript\">\n".
    "</script>\n".
    "<script type=\"text/javascript\">\n".
    "_uacct = \"UA-2207628-5\"\;\n".
    "urchinTracker()\;\n".
    "</script>\n".
    "</body>\n";

my @programs = 
    ("cnv_gff2game",
     "batch_eugene",
     "batch_geneid",
     "batch_genemark",
     "batch_genscan",
     "batch_findltr",
     "batch_ltrfinder",
     "batch_ltrseq",
     "batch_findmite",
     "batch_trf",
     "batch_hmmer",
     "batch_blast",
     "batch_repmask",
     "batch_tenest",
     "cnv_fgenesh2gff",
     "cnv_genemark2gff",
     "cnv_findltr2gff",
     "cnv_ltrfinder2gff",
     "cnv_ltrseq2gff",
     "cnv_ltrstruc2gff",
     "cnv_repseek2gff",
     "cnv_blast2gff",
     "cnv_repmask2gff",
     "cnv_tenest2gff",
     "cnv_gff2game",
     "cnv_game2gff3",
     "batch_hardmask",
     "dir_merge",
     "vennseq",
     "batch_findgaps",
     "clust_write_shell",
     "cnv_seq2dir",
     "fasta_merge",
     "fasta_shorten",
     "fetch_tenest",
     "gff_seg",
     "ltrstruc_prep",
     "seq_oligocount",
    );

#my @programs = ("vennseq");

for my $program (@programs) {

    # RUN POD2HTML
    my $perl_program = $program.".pl";
    my $html_file = $outdir.$program.".html";
    print STDERR "Processing $perl_program\n";
    my $cmd = "pod2html --infile $perl_program".
	" --css $css_location".
	" --title $perl_program".
	" --header".
	" --outfile $html_file";
    system($cmd) || 
	warn "Could not run system command $cmd\n";

    print STDERR "\n";

    # ADD GOOGLE ANALYTICS CODE TO RESULTING HTML
    if ($do_analytics) {

	if (-e $html_file) {
	    open (INFILE, $html_file) ||
		warn "Count not open HTML file for editing\n";
	    my @html_content = <INFILE>;
	    close INFILE;
	    
	    # ITERATE THROUGH THE RESULTING HTML TO ADD THE ANALYTICS CODE
	    open (FILEOUT, ">$html_file") ||
		die "Can not open html file:\n$html_file\n";
	    foreach $line (@html_content) {
		$line =~ s/\<\/body\>/$analytics_code/;
		print FILEOUT $line;
	    }
	    close FILEOUT;

	}

    } # End of do_analytics

} # End of for each program in programs array

1;
__END__

# 04/07/2009
# -wrote program to do pod2thml for batch of perl programs
# -added code to add google analytics code if desired
