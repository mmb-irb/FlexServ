#!/usr/bin/perl 
#
#
#      Manuel Rueda /09/06
#
#      Tirion's Hessian generation based on Harmonic Potential
#
#      last 27/06/07

use strict;
use FindBin '$Bin';

my $debug  = 0;
my $amberk = 0;
my $lorell = 1;
my $pdb    = $ARGV[0];     #input pdb file
my $out    = $ARGV[1];     #output hessian file
#my $kfile  = $ARGV[2];     #input k amber file only if ($amberk=1 or $lorell=1)
my $method = $ARGV[2];
my $cutoff = $ARGV[3];            #cutoff for bonded atoms
my $cutoff2 = $cutoff**2;  #avoids unncessary calculations
my $kforce = $ARGV[4];           # In Kcal/mol.A2 units
my $pablo  = 0;            # C=(3.8/Rab)^6
my $cabond = 0;            # In case we want to bond consecutive atoms
my $kbond  = 50;          # Bond Constant for consecutive atoms
my $kfile;

if ($#ARGV != 4) {
	print "Usage: $0 pdb hessian.dat method:0,1,2 cutoff fcte\n";
	exit;
}
die "ERROR cannot find the pdb file $pdb\n"if !-e $pdb;

$cutoff=999999 if ($amberk || $lorell);#########no hi ha cutoff si treballem amb kte d'amber o lorell
$cutoff2=$cutoff*$cutoff if ($amberk|| $lorell) ;#########no hi ha cutoff si treballem amb kte d'amber o lorell

#die "ERROR cannot find the amber cte file $kfile\n"if ($amberk||$lorell) && !-e $kfile;

#
# Parse pdb
# 
my (@x,@y,@z,$acount); # Take care since arrays start at 1, not 0
my $matk=undef if ($amberk || $lorell);

open(PDB,"$pdb") or die "failed to open $pdb";
while(<PDB>){
  next unless /^ATOM/;
  $acount++;
  my ($atom) = readPDBLine($_);
        $x[$acount]=$atom->{x};
        $y[$acount]=$atom->{y};
        $z[$acount]=$atom->{z};
    }
close PDB;
#
my  $N=$acount;  # Number of atoms
print "ATOMS: $N\n" if $debug;
my $N3=3*$N; # Number of elements

# Choose the method to use and setup the inputs
if ($method==0) {
	# Standard Linear NMA
	$pablo=0;
	$lorell=0;
} elsif ($method==1) {
	#Kovacs method
	$pablo=1;
	$lorell=0;
} elsif ($method==2) {
	# LOrell method
	$pablo=0;
	$lorell=1;

	# We need to compute intermediate force files
	$kfile = "fcte.dat";
	# passem coadena i num res a coords.dat
	#system("grep -G 'ATOM  ' $pdb | cut -c13-16,30-54 | grep CA | cut -c6- > coords.dat");
	system("grep -G 'ATOM  ' $pdb | cut -c13-27,30-54 | grep CA | cut -c6- > coords.dat");
	system("$Bin/lorellnma $N coords.dat $kfile");
	#unlink("coords.dat");
}

#
$matk=get_kmatrix($kfile) if $amberk;
$matk=get_lkmatrix($kfile) if $lorell;

hessian();

###################################################
sub get_kmatrix{###llegeix un fitxer de kte derivades de amber traj
#### ex /hosts/raid5/MoDEL/1OPC/1OPC_AMBER8_P99-T3P_0/ANALYSIS.ALL.ALL/CA-CA/1opc.ca.tab
###header of file
#CA-CA DISTANCES (ca) IN 1opc
#id   [    resdon    resacc ]   maskdon   maskacc    occ     kte     kstd [ error% ]   eqdst     avg [ error% ] (estimat   sigma ) (    std )
#00000 [    91-GLY    26-THR ]    :91@CA    :26@CA 100.00    2.70    2.52 [   6.73 ]  28.645  28.766 [   0.42 ] (  4.769   0.470 ) (  0.487 )
#00001 [    33-LEU    19-GLU ]    :33@CA    :19@CA 100.00    1.23    1.25 [   1.47 ]  18.777  18.756 [   0.11 ] (  6.965   0.695 ) (  0.690 )
###and so on for every ca-ca interaction
####
  my ($file)=shift;
  my $matrix=undef;
  my $max=0;
  open IN, "<$file" || die "cannot open kte file $file\n";
  while(<IN>){
    next if /^#/;
    my @tmp=split;
    my $tmp1=$tmp[5];
    my $tmp2=$tmp[6];
    my $kte=$tmp[8];
    $tmp1=~s/://;
    $tmp1=~s/\@CA\s*//;
    $tmp2=~s/://;
    $tmp2=~s/\@CA\s*//;
    if ($tmp1>$max){
      $max=$tmp1;
    }
 print "$tmp1 $tmp2 $kte\n";
    ${${$matrix}[$tmp1]}[$tmp2]=$kte;
    ${${$matrix}[$tmp2]}[$tmp1]=$kte;
  }
  $max++;
  for (my $i=1; $i<$max; $i++){ 
     for (my $j=1; $j<$max; $j++){ 
       printf "%5.2f ", ${${$matrix}[$i]}[$j]; 
     }
     print "\n";
  }
  return $matrix;
}
###################################################
sub get_lkmatrix{ #reads lorell kte files
# 1 2  9.75444221
# 1 3  8.89429379
	my $file=shift;
	my $matrix=undef;
	open (IN, "$file") || die "cannot open kte file $file\n";
	while(<IN>){
	next unless /\d/;
		my @tmp=split;
    		my $tmp1=$tmp[0];
    		my $tmp2=$tmp[1];
    		my $kte=$tmp[2];
    	${${$matrix}[$tmp1]}[$tmp2]=$kte;
    	#${${$matrix}[$tmp2]}[$tmp1]=$kte;
	}
	return $matrix;
}
###################################################
sub hessian {
##################################################

# See also Physical Review Letters, 77, 9, 1905 (1996) Tirion.

my ($count,$n,$m,$dx,$dy,$dz,$drdx,$drdy,$drdz,$r,$rr,@H,$k,$trace);

$count=0;

for $n (1..$N){
  for $m ($n+1..$N){
   # Upper right matrix

    $dx = $x[$n] - $x[$m];
    next if( abs($dx) > $cutoff);
    $dy = $y[$n] - $y[$m];
    next if( abs($dy) > $cutoff );
    $dz = $z[$n] - $z[$m];
    next if( abs($dz) > $cutoff );
    $rr = $dx*$dx + $dy*$dy + $dz*$dz;
    next if( $rr > $cutoff2);

    $r = sqrt($rr);
    $k=$kforce                     if (!$pablo);
    $k=((3.8/$r)**6) * $kforce     if $pablo;
    $k=$kbond                      if ( $cabond && $m == $n+1 );
    $k=${${$matk}[$n]}[$m]         if ($amberk || $lorell) ;###loading k from amber or lorell
    print "N:$n M:$m KTE:$k\n"	   if $debug;

    $drdx = $dx / $r;
    $drdy = $dy / $r;
    $drdz = $dz / $r;
    
    $H[3*$n-2][3*$n-2] += $k*$drdx*$drdx;
    $H[3*$n-2][3*$n-1] += $k*$drdx*$drdy;
    $H[3*$n-2][3*$n  ] += $k*$drdx*$drdz;
    $H[3*$n-1][3*$n-1] += $k*$drdy*$drdy;
    $H[3*$n-1][3*$n  ] += $k*$drdy*$drdz;
    $H[3*$n  ][3*$n  ] += $k*$drdz*$drdz;
    
    $H[3*$m-2][3*$m-2] += $k*$drdx*$drdx;
    $H[3*$m-2][3*$m-1] += $k*$drdx*$drdy;
    $H[3*$m-2][3*$m  ] += $k*$drdx*$drdz;
    $H[3*$m-1][3*$m-1] += $k*$drdy*$drdy;
    $H[3*$m-1][3*$m  ] += $k*$drdy*$drdz;
    $H[3*$m  ][3*$m  ] += $k*$drdz*$drdz;
    
    $H[3*$n-2][3*$m-2] += -$k*$drdx*$drdx;
    $H[3*$n-2][3*$m-1] += -$k*$drdx*$drdy;
    $H[3*$n-2][3*$m  ] += -$k*$drdx*$drdz;
    $H[3*$n-1][3*$m-2] += -$k*$drdy*$drdx;
    $H[3*$n-1][3*$m-1] += -$k*$drdy*$drdy;
    $H[3*$n-1][3*$m  ] += -$k*$drdy*$drdz;
    $H[3*$n  ][3*$m-2] += -$k*$drdz*$drdx;
    $H[3*$n  ][3*$m-1] += -$k*$drdz*$drdy;
    $H[3*$n  ][3*$m  ] += -$k*$drdz*$drdz;
    
    $count++;
  }
}

print "Number of interactions $count\n";

$count=0;
$trace=0;

#  Requirements for external diagonalization program
#  a matrix in the "i j non-zero-ij-matrix-element" format
#  left-bottom corner will be completed by external code

open(OUT,">$out") or die "Cannot open $out\n";

for $n ( 1 ..  $N3  ){
  $trace+=$H[$n][$n] if( defined $H[$n][$n] && $H[$n][$n] != 0 );
  for $m ( $n  ..  $N3 ){
      if( defined $H[$n][$m] && $H[$n][$m] != 0 ){
      printf(OUT "%6d %5d %15.12f\n",$n,$m,$H[$n][$m]);
      $count++;
    }
  }
}
close OUT;
print "Number of non-zero elements $count\n";
print "Matrix trace $trace\n";
}
###################################################
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

