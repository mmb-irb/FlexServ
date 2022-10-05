#!/usr/bin/perl
#
#       Manuel Rueda /09/06 (Adapted from M.Orozco code)
#       
#       Program for MC simulation in eigenvectors space
#
#       Last modified 27/02/07       
#	
#
use strict;

# Input Parameters.
my $long=@ARGV;
if ($long > 2 or $long < 1){
		    print "Usage: perl $0 <eigenvec.dat> [numSnapshots]\n";
		    exit(0);
}
my ($fforce,$nsnaps) =@ARGV;

$nsnaps = 1000 if(!$nsnaps);

#if  ($#ARGV != 0) {
#	print "$0 eval-force\n";
#	exit;
#}

#        my $fforce=$ARGV[0];
#
#       $nmc is the number of MC movements
#        my $nmc=125000; # $nmc=($buffer*<#snapshots_desired>)/$aratio

	my $nmc = 125 * $nsnaps;

#       Initial seed
	srand (12379987497); #if we need to reproduce numbers 

#       Number of modes used
        my $itot=50;

#       Buffer of lines
        my $buffer=50;

#       fundamental constants
        my $BOLTZ=0.00198717;
#       Temp
        my $TEMP=300;
#
        my $KbT=$BOLTZ*$TEMP;
	my $BETA = 1/$KbT;

#       read force constants (Kcal/mol.A^2) ->@stiff
        my @stiff=read_force($fforce);  #start at 0

#       xconf0 has the equilibrium parameters for each mode, set to 0.0

        my @xconf0=(); for my $i (0 .. ($itot-1)){ $xconf0[$i]=0; } # (although not necessary since empty=0, forced to be 0)

        my @xconft=@xconf0; # Initial configuration of the try (at equilibrium)
 
#       agr stands for the agressivity, the value
#       of the maximum allowed displacement along the first mode. Adjust it to have
#       around 40% acceptance. The others displacements are adjusted automatically
        my $agr=2.75*sqrt($KbT/$stiff[0]);
    
#       scale the movements

        my @scal=scaling();

#       Energy for all modes at equilibrium = 0
 
        my @ener=();
       	for my $i (0 .. ($itot-1)) {
		$ener[$i]=0;
	} # (although not necessary since empty=0, forced to be 0)
 
#       ener0 starting energy
        
        my $ener0=0;

#       Total energy at the begining

        my $enertot=$ener0;

###################################################
#       Now starts MC Metropolis run
#       Select Random 1 mode
###################################################
#       
my ($dx,$dd,$tdx,$itt,$ittp,$delta,$enerbase,$enex,$etry,$newline,$iflag,$counter,$accepted); #Logical value=0
for my $iter ( 1 .. $nmc) { #Start LOOP
        # $itt from 0 to $itot-1 
        # Problems if $itt=$itot (almost imposible, no cases although $nmc=1.000.000.000)
#       itt is the randomly selected node to be changed
        $itt=int(rand($itot));
#       note we sample both + and - displacements
#       $dx number between -1 and 1 
        $dx=rand(2);
        $dx=$dx-1;

#       we multiply it by the scaling factor
        $dx=$dx*$scal[$itt];

#	Energy of all the modes except the changed
        $enerbase=$enertot-$ener[$itt]; # 0 in the first iteration
#       Displacement of the try, taking into account the previous for the same mode
        $tdx=$xconft[$itt]+$dx;
#       enex has the new energy of the step after changing itt variable
	$dd=$tdx-$xconf0[$itt]; # x-x0 -> $xconf0[$it] shoukd be 0 in this case.
        $enex=$stiff[$itt]*$dd*$dd;

#       This is the total energy of the try
        $etry=$enerbase+$enex;
#	Do metropolis test
#	Stochastic processes use a random sampling procedure to search conformational space
#	The score (energy) is calculated in each step and compared to the previous. If the new energy is lower
# 	the step is accepted, otherwise the result is treated probabilistically by a Boltzmann mechanism
# 	The higher the temperature the higher the likehood the step is accepted
#       	PS: Another way to generate the DX foreach mode could be from 'normal random values'
#       	(see Box-Muller transform and Ziggurat algorithm) without any further Metropolis test.
#	(better if avoiding subroutines/sending variables-references->slower code)
        $delta=$etry-$enertot;
        if ( $delta < 0 || exp(-$BETA*$delta) > rand()  ) {
#       accepted configuration
		$iflag = 1;
        	$accepted++; 
        	$counter++;
                #New total energy
                $enertot=$etry;
                #New energy for the mode
                $ener[$itt]=$enex;
                #New DX for the mode
                $xconft[$itt]=$tdx;
                # Print
                if ($counter == $buffer ){
		        $ittp=$itt+1; #since arrays start at 0
		        print "\*\*\*\*\n";
		        print "ITERATION:$iter,MODE:$ittp,ETOT:$enertot,DX:$dx,CUM-DX($ittp):$tdx\n";
		        print "DX for $itot Modes\n";
		        print "-------------------------------------------------------------------\n";
		        $newline=0;
		        #DX Foreach mode
		        for (@xconft) {
				$newline++;
				if ( $newline != 10 ) {
					printf "%10.4f",$_;
				} else {
					printf "%10.4f\n",$_; $newline=0
				}
			}
		        $counter=0;
		       	print "\n\n" if ($newline != 0);
		        print "\n" if ($newline == 0);
               	}
		#
      } else {
            	$iflag = 0;
      }
}	#END LOOP
#
my $aratio=($accepted/$nmc)*100;
print "-------------------\n";
print "RATIO(%):$aratio\n";
print "-------------------\n";
############################################
# END MC
############################################

sub scaling {
	my @scal;
	$scal[0]=$agr;
	for my $i ( 1 .. ($itot-1) ) {
		my $rat=$stiff[0]/$stiff[$i];
		$scal[$i]=sqrt($agr*$agr*$rat);
	}
	return @scal;
}

sub read_force {
	my @force;
	my $file=shift;
	my @line;
	my $counter;
	open (EVAL, "$file") ||die "cannot open $file\n";
	while(<EVAL>){
		chomp;
		$counter =1 if (/\*\*\*\*/);
		if ( $counter eq 1 && /[0-9]/ ){
			#just in case, keep only numbers no blank spaces
			@line= grep { /[0-9]/ } (split (/ +/,$_));
			#Force constants divided by 2
			#E=0.5*K*x^2
			#print "$line[-1]\n";
			push @force,$line[-1]/2; 
			$counter = 0;
		}
	}
	close EVAL;
	return @force; #starts at 0
}
#
