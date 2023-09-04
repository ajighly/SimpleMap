#!/usr/bin/perl -w


######################################################################################################################

############################## SimpleMap 
############################## Script: Before_Mapping.pl 
############################## Written by: Abdulqader Jighly
############################## 		   abdulqader.jighly@depi.vic.gov.au
############################## 		   a.jighly@cgiar.org
############################## Usage: perl Before_Mapping.pl GenotypeFileName NoMarkers NoGenotypes RepulsionThreshold

######################################################################################################################


$filename1    = $ARGV[0];
$NoMarkers    = $ARGV[1];
$NoGenotypes  = $ARGV[2];
$MaxRepuls    = $ARGV[3];

unless(open(DATAFILE,"< $filename1")){						# Print error massege if the input file is not exist
	print "Unable to open $filename1: $!\n";
	exit 1; 
}
$x='';
open(CONVERT, ">Converted_Haplotype");
while(<DATAFILE>) {
	chomp($_);
	if (/(\S+)\t(.+)/){ 							# To diffrenciate the marker names from their haplotypes
		$Marker = $1;
		$Haplo  = $2;
	}
	$Haplo =~ s/\t//g; 							# Remove the Tab char
	$Haplo =~ tr/abcdhABCDH/1234512345/; 					# Unify the allele chars for the comparison
	print CONVERT "$Marker	$Haplo\n";
}
close CONVERT;
close DATAFILE;

open(CONVERTEDFILE1,"< Converted_Haplotype");

while(<CONVERTEDFILE1>) {							# The first loop

	$lineNoA++;
	if ($lineNoA == $NoMarkers) {goto PRINT;}				# Print (line 119) when you reach the last marker because this one will not be compared with any other marker
	chomp($_);

	if (/(\S+)\t(.+)/){ 							# To diffrenciate the marker names from their haplotypes
		$Marker1 = $1;
		$Haplo1  = $2;
	}

	if (scalar grep $Marker1 eq $_, @Compare) { next; } 			# If the marker is already existed in the compare list (which contain the doublicated markers) then go to the next marker

	@chars1 = split("", $Haplo1);						# Save the haplotype in a list

	open(CONVERTEDFILE2,"<Converted_Haplotype");

	while(<CONVERTEDFILE2>) { 						# The second loop; To compare with the folowing markers in the original data file

		$lineNoB++;
		if ($lineNoB < $lineNoA + 1) { next; } 				# To skip the previous markers
		chomp($_);

		if (/(\S+)\t(.+)/){ 						# To diffrenciate the marker names from their haplotypes
			$Marker2 = $1;
			$Haplo2  = $2;
		}

		@chars2 = split("", $Haplo2);					# Save the haplotype in a list
		
		$Repulsion = 0;

		for ($i=0; $i < $NoGenotypes; $i++) {				# Start the comparison for each pair of markers

			if ($chars1[$i] eq '-' or $chars2[$i] eq '-') { next; }	# Skip the genotype if it has a missing value in one of the genotypes
			if ($chars1[$i] == $chars2[$i])  { next; }
			if ($chars1[$i] == 3 && ($chars2[$i] == 2)) { next; }
			if ($chars1[$i] == 4 && ($chars2[$i] == 1)) { next; }
			if ($chars2[$i] == 3 && ($chars1[$i] == 2)) { next; }
			if ($chars2[$i] == 4 && ($chars1[$i] == 1)) { next; }	
			if ($chars1[$i] == 3 && ($chars2[$i] == 5)) { next; }
			if ($chars1[$i] == 4 && ($chars2[$i] == 5)) { next; }
			if ($chars2[$i] == 3 && ($chars1[$i] == 5)) { next; }
			if ($chars2[$i] == 4 && ($chars1[$i] == 5)) { next; }


			if ($chars1[$i] != $chars2[$i]) { $Repulsion++; }	# Increase the repulsion counter in the case of different haplotypes

			if ($Repulsion > $MaxRepuls) {				# End the comparison if the markers have more recombinations that the accepted threshold without saving them in the "ClustMarkers" hash
				goto ENDCOUNTING;				# Go to line 110
			}

		}

		if ($Marker2 ne $Marker1) {					# Save Marker2 as a key for the hash ClustMarkers to be saved in the Compare list and excluded from the first loop
			if (scalar grep $Marker2 eq $_, @Compare) {
				($x, $Repul) = split(/\t/, $ClustMarkers{$Marker2});
				if ($Repulsion >= $Repul) {
					goto ENDCOUNTING;
				} else {
					$ClustMarkers{$Marker2} = $Marker1 . "\t" . $Repulsion;
					goto ENDCOUNTING;
				}

			} else {
			$ClustMarkers{$Marker2} = $Marker1 . "\t" . $Repulsion;	
			}
		}
		ENDCOUNTING:
	}

	close CONVERTEDFILE2;
	@Compare = keys %ClustMarkers;						# Save the second marker of the highly linked markers in the compare list to pass it from the first loop if existed

}
close CONVERTEDFILE1;

PRINT:

open(OUTFILE,">Repulsions.txt");
foreach (keys %ClustMarkers) {							# The printing loop for the highly linked markers and the number of recombination(s) between them
    print OUTFILE "$_	$ClustMarkers{$_}\n";
}

open(DATAFILE,"< $filename1");
open(FORMAPPING,">For_Mapping.txt");

while(<DATAFILE>) {
	chomp($_);
	if (/(\S+)\t(.+)/){ 							# To diffrenciate the marker names from their haplotypes
		$Marker = $1;
		$Haplo  = $2;
	}

	if (scalar grep $Marker eq $_, @Compare) { 				# If the marker is not existed in the compare list (which means that this marker is unique) then print it in the file "For_Mapping.txt"
		next;
	} else {
		print FORMAPPING "$Marker	$Haplo\n";
	}
}
