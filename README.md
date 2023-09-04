# SimpleMap V2
A tool to streamline high density genetic linkage group constrcution

SimpleMap was previously available on https://sourceforge.net/projects/simplemap-aj/ but this is a new version with improved speed for the first step of the analysis

SimpleMap is consisted of two scripts. To use the first one “Before_Mapping.r” for the attached example, write on the command line:

Rscript Before_Mapping.r geno.txt 1463 163 1

Where geno.txt is the genotype file as JoinMap format; 1463 is the number of markers; 163 is the number of lines and 1 is the repulsion threshold.

This will end with three files: Repulsion.txt; Converted_Haplotype; and For_Mapping.txt. Use the file For_Mapping.txt with any mapping software, Generate the map and export is in a text file. (In our example the exported map file is “InputMap.Map”
PS. Don’t delete the Repulsion.txt and Converted_Haplotype files until you apply the second script “After_Mapping.pl.”
To use the second script “After_Mapping.pl”, write on the command line:

perl After_Mapping.pl InputMap.Map OutputMap.map


Citation:
Jighly A, Joukhadar R, Alagu M (2015) SimpleMap: A Pipeline to Streamline High Density Linkage Map Construction. doi: 10.3835/plantgenome2014.09.0056
