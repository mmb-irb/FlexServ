#!/usr/bin/perl
#
#       Manuel Rueda /11/06 
#
#       Program for generation of cartesian coordinates via mc-eigen.pl's output
#
#       Last modified 24/11/06
#

use strict;

sub usage {
	printf STDERR "usage:  $0 [-flags] \n";
	printf STDERR "Flags: -pdb reference (avg) \n";
	printf STDERR "       -evec eigenvectors file\n";
	printf STDERR "       -pout projected ASCII 'amber-crd-like' trajectory\n";
	printf STDERR "       -i input file (output from Montecarlo)\n";
        printf STDERR "       -n #vectors\n";  #needed for avoiding read all evecfile
	exit 1;
}
#
	my ($pdb,$evec,$traj,$mcfile,$nevec);
	if  ( scalar @ARGV != 10  ) {&usage();}
	else{
       	while (@ARGV) {
  		if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    		&usage();
  		 } elsif ($ARGV[0] eq "-pdb") {
    		shift @ARGV;
    		$pdb= shift @ARGV;
  		 } elsif ($ARGV[0] eq "-evec") {
    		shift @ARGV;
    		$evec=shift @ARGV;
  		 } elsif ($ARGV[0] eq "-pout") {
    		shift @ARGV;
    		$traj=shift @ARGV;
   		 } elsif ($ARGV[0] eq "-i") {
    		shift @ARGV;
    		$mcfile=shift @ARGV;
                  } elsif ($ARGV[0] eq "-n") {
                shift @ARGV;
                $nevec=shift @ARGV;
  		 } elsif ($ARGV[0] !~ /-pdb|-evec|-pout|-i|-n/) {
    		printf STDERR "invalid option\n";
    		&usage();
  		 } 
   	}
       }

#	Open pdb reference file
	my (@ref,$acount); # Take care since @ref starts at 0, not 1

	open(PDB,"$pdb") or die "failed to open $pdb";
	while(<PDB>){
  	next unless /^ATOM/;
  	$acount++;
  	my ($atom) = readPDBLine($_);
        push @ref,$atom->{x},$atom->{y},$atom->{z} ;
    	}
	close PDB;
        my $N=$acount;   # Number of atoms
	my $N3=3*$N;     # Number of elements

#	Read ptraj eigenvector file
	my $vec=read_evec($evec,$nevec);

#	Parse the input file (output from mc-eigen.pl
        my ($proj,$inumber,$start);
	open (MC,"$mcfile") or die "Cannot open $mcfile";
	while (<MC>){
        chomp;
        if (/\*\*\*\*/){ $inumber++; $start=0};
        $start++ if /^\-\-\-/;
        push @{$proj->[$inumber]} ,grep { /[0-9]/ } (split / +/,$_) if ($inumber && $start) ;
         }
      	close MC;

#	Printing trajectory
	print "Outputting trajectory to file $traj ...\n" ;
  	open(PROJECTION,">$traj") ||die "cannot open $traj";
     	print PROJECTION " MC generated trajectory \n";
       	projection();
   	close PROJECTION;
	print "-----------------------\n";
	print "Successfully ending\n";

#
sub projection {
	my $crd=0;
        my $j;
        for my $snap (1..$inumber) {
             #print "Generation of snapshot:$snap, using $nevec evecs\n";
                my $warning=0;
                my $format=0;
               # 10f8.3 is the trajectory format
                for my $i (0..($N3-1)){
                       #print "ATOM $i\n";
                      $j=$i+2; #since evec has index and eval
                           $format++;
                           $crd=$ref[$i];
                           for my $ndx (1..$nevec){
                           	$crd=$crd+($vec->[$ndx][$j]*$proj->[$snap][$ndx-1]);
                                #print"($vec->[$ndx][$j]#*PROJ#$proj->[$snap][$ndx-1]#)\n";
                           }
                           printf PROJECTION "%8.3f", $crd;
                              if ( $format == 10) { print PROJECTION "\n"; $format=0; $warning=1;}
                              else {$warning=0;}
                }
        #take care because if 10 columns at last line,we have a blank line that ptraj
        #will understand as a box!!!!!
         if ( $warning != 1 ) { print PROJECTION "\n";}
         }
}
#
sub readPDBLine {
        my ($line) = shift;
        my $newAt = {};
        my $newRes = {};
        $newAt->{name}=substr($line, 12,4);
        $newAt->{altloc}=substr($line, 16,1);
        $newAt->{residueId}=substr($line,17,9);
        $newAt->{x}=substr($line,30,8);
        $newAt->{y}=substr($line,38,8);
        $newAt->{z}=substr($line,46,8);
        $newAt->{occ}=substr($line, 54,6);
        $newAt->{Bfact}=substr($line, 60,6);
        $newAt->{charge}=$newAt->{occ};
        $newAt->{type}=substr($newAt->{name},1,1);
        return $newAt;
}
sub read_evec {
        my ($file,$nevec)=@_;
        my $header;
        my $vec=();
        open (EVEC,"$file") || die "Cannot open $file\n" if $file !~/\.gz$/ ;
        open (EVEC,"zcat $file |") || die "Cannot open $file\n" if $file =~/\.gz$/ ;
        my @maxevec= grep {/[0-9]/} ( map {split / +/,$_} (map scalar <EVEC>, 1..2) );
        my $max=$maxevec[0]-6;
        if ( $nevec > $max ) { print "$nevec is greater than $max\n"; exit }
        while(<EVEC>){
                chomp;
                $header++ if /\*\*\*\*/ ;
                last if ($header > $nevec);
                push @{$vec->[$header]} ,grep { /[0-9]/ } (split / +/,$_) if ($header) ; #$VEC start from 1, elements from 0
        }
        close EVEC;
        return $vec;
}

