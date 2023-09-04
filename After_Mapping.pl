#!/usr/bin/perl -w



#######################################################################################

############################## SimpleMap 
############################## Script: After_Mapping.pl 
############################## Written by: Abdulqader Jighly
############################## 		   abdulqader.jighly@depi.vic.gov.au
############################## 		   a.jighly@cgiar.org
############################## Usage: perl After_Mapping.pl InputMapFile OutputMapFile

#######################################################################################



$InputFile  = $ARGV[0];
$OutputFile = $ARGV[1];

unless(open(DATAFILE,"< $InputFile")){												# Print error massege if the input file is not exist
	print "Unable to open $InputFile: $!\n";
	exit 1; 
}
$i = 0;
while(<DATAFILE>) {														# Save the map file into an array and filter inappropriate lines and chars

	chomp($_);

	if ($_ =~ /^;/ or $_ =~ /^$/) {

		next;

	} elsif ($_ =~ /group\s+(\S+)/) {

		$map[$i] = 'group';
		$pos[$i] = $1;	
		$i++;
		next;

	} elsif ($_ =~ /(\S+)\s+(\d+\.*\d*)/) {

		$map[$i] = $1;
		$pos[$i] = $2;
		$i++;

	}
}
close DATAFILE;

@NewMap = @map;
@NewPos = @pos;
open(DATAFILE,"<Repulsions.txt");

while(<DATAFILE>) {														# Save the linked markers and the number of recombinations between them into a hash

	chomp($_);
	($Marker1, $Marker2, $Repulsion) = split(/\t/, $_);
	$LinkedMarkers{$Marker1} = $Marker2;
	$Repulsions{$Marker1}    = $Repulsion;
}
close DATAFILE;

open (DATAFILE,"<Converted_Haplotype");

while(<DATAFILE>) {														# Save the the marker haplotypes into a hash
	chomp($_);
	($MarkForHaplo, $haplotype) = split(/\t/, $_);
	$haplotypeHash{$MarkForHaplo} = $haplotype;
}
close DATAFILE;

foreach $KeyMarker (keys %LinkedMarkers) {						

	if (scalar grep $LinkedMarkers{$KeyMarker} eq $_, @map) {								# If the marker that we will compare with (which located in the "For_Mapping.txt" file) is mapped, then do the follow

		@i  = grep { $map[$_] eq $LinkedMarkers{$KeyMarker} } 0 .. $#map;						# Get the ID of the marker from the array @map
		$id = $i[0];
		@j  = grep { $NewMap[$_] eq $LinkedMarkers{$KeyMarker} } 0 .. $#NewMap;						# Get the ID of the marker from the array @map
		$NewId = $j[0];

		if ($Repulsions{$KeyMarker} == 0) {										# If there were no line that crossed over between both markers, then the position of the unmapped marker will equal to the position of the mapped marker. Continue next line...

			splice @NewMap, $NewId, 0, $KeyMarker;									# Add the new marker to the array @map
			splice @NewPos, $NewId, 0, $pos[$id];									# Add the position of the new marker to the array @pos	

		} else {													# If there were line(s) that crossed over between both markers, then do the follow
			$HaploMap = $haplotypeHash{$LinkedMarkers{$KeyMarker}};							# Get the haplotype of the mapped marker linked to $KeyMarker
			@charsMap = split("", $HaploMap);									# Save the haplotype in a list
			$HaploKey = $haplotypeHash{$KeyMarker};									# Get the haplotype of the $KeyMarker
			@charsKey = split("", $HaploKey);									# Save the haplotype in a list
			$Repulsion = $Repulsions{$KeyMarker};
			$NoGenotypes = @charsMap;
			$x = 0;
			$y = 0;

			if ($map[$id-1] eq 'group') {										# If the mapped marker was the first marker in the group

				$Haplo2 = $haplotypeHash{$map[$id+1]};								# Get the haplotype of the following marker
				@chars2 = split("", $Haplo2);									# Split the haplotype into the array @chars2

				for ($i=0; $i < $NoGenotypes; $i++) {								# Start the comparison for each pair of markers

					if ($charsKey[$i] eq '-' or $chars2[$i] eq '-' or $charsMap[$i] eq '-') { next; }	# Skip the genotype if it has a missing value in one of the genotypes
					if ($charsKey[$i] == $chars2[$i] && $charsKey[$i] == $charsMap[$i]) { next; }		# Next in the case of similar haplotypes
					if ($charsKey[$i] != $chars2[$i]) { $x++; }						# Increase the repulsion counter in the case of different haplotypes
					if ($charsMap[$i] != $chars2[$i]) { $y++; }
				}

				if ($x > $y) {											# If the distance between the mapped marker and the following one is less than this following marker and the $keymarker 

					$FinalPos = - ((100 * $Repulsion) / $NoGenotypes);					# Calculate the position of the marker
					splice @NewMap, $NewId, 0, $KeyMarker;							# Add the new marker to the array @map
					splice @NewPos, $NewId, 0, $FinalPos;							# Add the position of the new marker to the array @pos

				} elsif ($x < $y) {										# If the distance between the mapped marker and the following one is bigger than this following marker and the $keymarker 

					$FinalPos = $pos[$id+1] - ($x * ($pos[$id+1] / $y));					# Calculate the position of the marker
					splice @NewMap, $NewId+1, 0, $KeyMarker;						# Add the new marker to the array @map
					splice @NewPos, $NewId+1, 0, $FinalPos;							# Add the position of the new marker to the array @pos

				} elsif ($x == $y) {										# Could happens in the presence of missing data

					$FinalPos = - ((100 * $Repulsion) / $NoGenotypes);					# Calculate the position of the marker
					splice @NewMap, $NewId, 0, $KeyMarker;							# Add the new marker to the array @map
					splice @NewPos, $NewId, 0, $FinalPos;							# Add the position of the new marker to the array @pos
				}

			} elsif ($map[$id+1] eq 'group') {									# If the mapped marker was the last marker in the group

				$Haplo2 = $haplotypeHash{$map[$id-1]};								# Get the haplotype of the following marker
				@chars2 = split("", $Haplo2);									# Split the haplotype into the array @chars2

				for ($i=0; $i < $NoGenotypes; $i++) {								# Start the comparison for each pair of markers

					if ($charsKey[$i] eq '-' or $chars2[$i] eq '-' or $charsMap[$i] eq '-') { next; }	# Skip the genotype if it has a missing value in one of the genotypes
					if ($charsKey[$i] == $chars2[$i] && $charsKey[$i] == $charsMap[$i]) { next; }		# Next in the case of similar haplotypes
					if ($charsKey[$i] != $chars2[$i]) { $x++; }						# Increase the repulsion counter in the case of different haplotypes
					if ($charsMap[$i] != $chars2[$i]) { $y++; }
				}

				if ($x > $y) {

					$FinalPos = $pos[$id] + ((100 * $Repulsion) / $NoGenotypes);				# Calculate the position of the marker
					splice @NewMap, $NewId+1, 0, $KeyMarker;						# Add the new marker to the array @map
					splice @NewPos, $NewId+1, 0, $FinalPos;							# Add the position of the new marker to the array @pos

				} elsif ($x < $y) {

					$FinalPos = $pos[$id-1] + (($x * ($pos[$id] - $pos[$id-1])) / $y);			# Calculate the position of the marker
					splice @NewMap, $NewId, 0, $KeyMarker;							# Add the new marker to the array @map
					splice @NewPos, $NewId, 0, $FinalPos;							# Add the position of the new marker to the array @pos

				} elsif ($x == $y) {										# Could happens in the presence of missing data

					$FinalPos = $pos[$id] + ((100 * $Repulsion) / $NoGenotypes);				# Calculate the position of the marker
					splice @NewMap, $NewId+1, 0, $KeyMarker;						# Add the new marker to the array @map
					splice @NewPos, $NewId+1, 0, $FinalPos;							# Add the position of the new marker to the array @pos
				}

			} elsif ($map[$id+1] ne 'group' and $map[$id-1] ne 'group') {						# If the marker linked to the $Keymarker is in the middle of the linkage group

				$Haplo1 = $haplotypeHash{$map[$id-1]};								# Get the haplotype of the previous marker
				@chars1 = split("", $Haplo1);									# Split the haplotype into the array @chars1
				$Haplo2 = $haplotypeHash{$map[$id+1]};								# Get the haplotype of the following marker
				@chars2 = split("", $Haplo2);									# Split the haplotype into the array @chars2

				for ($i=0; $i < $NoGenotypes; $i++) {

					if ($charsKey[$i] eq '-' or $chars1[$i] eq '-') { next; }				# Skip the genotype if it has a missing value in one of the genotypes
					if ($charsKey[$i] == $chars1[$i]) { next; }						# Next in the case of similar haplotypes
					if ($charsKey[$i] != $chars1[$i]) { $x++; }						# Increase the repulsion counter in the case of different haplotypes

				}

				for ($i=0; $i < $NoGenotypes; $i++) {

					if ($charsKey[$i] eq '-' or $chars2[$i] eq '-') { next; }				# Skip the genotype if it has a missing value in one of the genotypes
					if ($charsKey[$i] == $chars2[$i]) { next; }						# Next in the case of similar haplotypes
					if ($charsKey[$i] != $chars2[$i]) { $y++; }						# Increase the repulsion counter in the case of different haplotypes

				}

				if ($x > $y) {

					$FinalPos = $pos[$id] + ((($pos[$id+1]-$pos[$id]) * $Repulsion) / ($Repulsion + $y));	# Calculate the position of the marker
					splice @NewMap, $NewId+1, 0, $KeyMarker;						# Add the new marker to the array @map
					splice @NewPos, $NewId+1, 0, $FinalPos;							# Add the position of the new marker to the array @pos

				} elsif ($x < $y) {

					$FinalPos = $pos[$id-1] + ((($pos[$id]-$pos[$id-1]) * $x) / ($Repulsion + $x));		# Calculate the position of the marker
					splice @NewMap, $NewId, 0, $KeyMarker;							# Add the new marker to the array @map
					splice @NewPos, $NewId, 0, $FinalPos;							# Add the position of the new marker to the array @pos

				} elsif ($x == $y) {										# Could happens in the presence of missing data

					$FinalPos = $pos[$id-1] + (($pos[$id+1]-$pos[$id-1]) / 2);				# Calculate the position of the marker
					splice @NewMap, $NewId, 0, $KeyMarker;							# Add the new marker to the array @map
					splice @NewPos, $NewId, 0, $FinalPos;							# Add the position of the new marker to the array @pos
				}
			}
		}
	} else { next; }													# If the marker that we will compare with (which located in the "For_Mapping.txt" file) is unmapped, then go to the next key in the hash %LinkedMarkers
}

$MappedMarkers = @NewMap;
SORTAGAIN:

for ($i = 2; $i < $MappedMarkers - 1; $i++) {											# To resort marker positions within each group

	if ($NewMap[$i] eq 'group') { next; }
	elsif ($NewMap[$i+1] eq 'group') { next;}
	else {
		if ($NewPos[$i] > $NewPos[$i+1]) {
			$Change       = $NewPos[$i];
			$NewPos[$i]   = $NewPos[$i+1];
			$NewPos[$i+1] = $Change;
			$Change       = $NewMap[$i];
			$NewMap[$i]   = $NewMap[$i+1];
			$NewMap[$i+1] = $Change;
		}
	}
}

for ($i = 2; $i < $MappedMarkers - 1; $i++) {											# To check the sort the markers within linkage groups

	if ($NewMap[$i] eq 'group') { next; }
	elsif ($NewMap[$i+1] eq 'group') { next;}
	else {
		if ($NewPos[$i] > $NewPos[$i+1]) {
			goto SORTAGAIN;												# Goto line 208 to repeat the sorting step
		}
	}
}

$rescale = 0;
for ($i = 1; $i < $MappedMarkers; $i++) {											# To let all linkage group positions start from 0

	if ($NewMap[$i] eq 'group') { 
		$rescale = 0; 
		next;
	}
	if ($NewPos[$i] < 0) {
		$rescale = abs ($NewPos[$i]);
	}
	$NewPos[$i] = $NewPos[$i] + $rescale;
}

$j = 0;
open(OUTMAP,">$OutputFile");													# Print the final map
foreach (@NewMap) {
	print OUTMAP "$NewMap[$j]	$NewPos[$j]\n";
	$j++;
}
close OUTMAP;
