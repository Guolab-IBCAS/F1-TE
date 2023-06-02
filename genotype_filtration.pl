#Step 1: Filter out heterozygous loci in natural accessions as C. rubella is a highly selfing species
open OUT,">./all.filtered_snps.geno.update.content.filterParent1.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.vcf"); #Input file is the combined SNP matrix. 
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1; #Split the genotype into two alleles, below are the same
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P86_11,$P86_12)=split /\//,$P86_1;
($P86_21,$P86_22)=split /\//,$P86_2;
($P86_31,$P86_32)=split /\//,$P86_3;

($P9_11,$P9_12)=split /\//,$P9_1;
($P9_21,$P9_22)=split /\//,$P9_2;
($P9_31,$P9_32)=split /\//,$P9_3;

($PM_11,$PM_12)=split /\//,$PM_1;
($PM_21,$PM_22)=split /\//,$PM_2;
($PM_31,$PM_32)=split /\//,$PM_3;

#Skip this line if the two alleles of a parent were different
next if ($P8_11 != $P8_12 || $P8_21 != $P8_22 || $P8_31 != $P8_32 ||$P86_11 != $P86_12|| $P86_21 != $P86_22|| $P86_31 != $P86_32 || $P9_11 != $P9_12 || $P9_21 != $P9_22 ||$P9_31 != $P9_32||$PM_11 != $PM_12 || $PM_21 != $PM_22|| $PM_31 != $PM_32);
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 2: Filter out missing genotype loci in natural accessions as these allelic expression would not be counted, and keep the genotypes at a locus across three replicates of 879 the same
open OUT,">./all.filtered_snps.geno.update.content.filterParent2.1.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterParent1.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P86_11,$P86_12)=split /\//,$P86_1;
($P86_21,$P86_22)=split /\//,$P86_2;
($P86_31,$P86_32)=split /\//,$P86_3;

($P9_11,$P9_12)=split /\//,$P9_1;
($P9_21,$P9_22)=split /\//,$P9_2;
($P9_31,$P9_32)=split /\//,$P9_3;

($PM_11,$PM_12)=split /\//,$PM_1;
($PM_21,$PM_22)=split /\//,$PM_2;
($PM_31,$PM_32)=split /\//,$PM_3;

#Skip this line if the genotype of a parent contain missing allele
next if ($P8_1=~/\.\/\./ || $P8_2=~/\.\/\./ || $P8_3=~/\.\/\./ || $P86_1=~/\.\/\./ || $P86_2=~/\.\/\./ || $P86_3=~/\.\/\./ || $P9_1=~/\.\/\./ || $P9_2=~/\.\/\./ || $P9_3=~/\.\/\./ || $PM_1=~/\.\/\./ || $PM_2=~/\.\/\./ || $PM_3=~/\.\/\./);

#Print this line where the genotypes of three replicates of 879 were the same
if ($P8_1 eq $P8_2 && $P8_1 eq $P8_3){
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 3: Keep the genotypes at a locus across three replicates of 86IT1 the same
open OUT,">./all.filtered_snps.geno.update.content.filterParent2.2.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterParent2.1.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P86_11,$P86_12)=split /\//,$P86_1;
($P86_21,$P86_22)=split /\//,$P86_2;
($P86_31,$P86_32)=split /\//,$P86_3;

($P9_11,$P9_12)=split /\//,$P9_1;
($P9_21,$P9_22)=split /\//,$P9_2;
($P9_31,$P9_32)=split /\//,$P9_3;

($PM_11,$PM_12)=split /\//,$PM_1;
($PM_21,$PM_22)=split /\//,$PM_2;
($PM_31,$PM_32)=split /\//,$PM_3;

#Print this line where the genotypes of three replicates of 86IT1 were the same
if ($P86_1 eq $P86_2 && $P86_1 eq $P86_3){
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 4: Keep the genotypes at a locus across three replicates of 928 the same
open OUT,">./all.filtered_snps.geno.update.content.filterParent2.3.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterParent2.2.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P86_11,$P86_12)=split /\//,$P86_1;
($P86_21,$P86_22)=split /\//,$P86_2;
($P86_31,$P86_32)=split /\//,$P86_3;

($P9_11,$P9_12)=split /\//,$P9_1;
($P9_21,$P9_22)=split /\//,$P9_2;
($P9_31,$P9_32)=split /\//,$P9_3;

($PM_11,$PM_12)=split /\//,$PM_1;
($PM_21,$PM_22)=split /\//,$PM_2;
($PM_31,$PM_32)=split /\//,$PM_3;

#Print this line where the genotypes of three replicates of 928 were the same
if ($P9_1 eq $P9_2 && $P9_1 eq $P9_3){
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 5: Keep the genotypes at a locus across three replicates of MTE the same
open OUT,">./all.filtered_snps.geno.update.content.filterParent2.4.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterParent2.3.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P86_11,$P86_12)=split /\//,$P86_1;
($P86_21,$P86_22)=split /\//,$P86_2;
($P86_31,$P86_32)=split /\//,$P86_3;

($P9_11,$P9_12)=split /\//,$P9_1;
($P9_21,$P9_22)=split /\//,$P9_2;
($P9_31,$P9_32)=split /\//,$P9_3;

($PM_11,$PM_12)=split /\//,$PM_1;
($PM_21,$PM_22)=split /\//,$PM_2;
($PM_31,$PM_32)=split /\//,$PM_3;

#Print this line where the genotypes of three replicates of MTE were the same
if ($PM_1 eq $PM_2 && $PM_1 eq $PM_3){
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 6: Filter out the loci that the three replicates of 879 × 86IT1 contain genotypes both absent in 86IT1 and 879.
#e.g. if genotype at a locus in 879 × 86IT1 was "1/2", while in 879 was "0/0" and in 86IT1 was "1/1", the genotype of 879 × 86IT1 at this locus should be "0/1" rather than "1/2", and this locus would thus be filtered out.
open OUT,">./all.filtered_snps.geno.update.content.filterFill3.1.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterParent2.4.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P8_86_11,$P8_86_12)=split /\//,$P8_86_1;
($P8_86_21,$P8_86_22)=split /\//,$P8_86_2;
($P8_86_31,$P8_86_32)=split /\//,$P8_86_3;

($P86_11,$P86_12)=split /\//,$P86_1;
($P86_21,$P86_22)=split /\//,$P86_2;
($P86_31,$P86_32)=split /\//,$P86_3;

#Skip this line if the two alleles of any replicate of 879 × 86IT1 contain genotype both absent in 86IT1 and 879
next if (($P8_86_11 != $P8_11 and $P8_86_11 != $P86_11) || ($P8_86_12 != $P8_11 and $P8_86_12 != $P86_11) || ($P8_86_21 != $P8_11 and $P8_86_21 != $P86_11) || ($P8_86_22 != $P8_11 and $P8_86_22 != $P86_11) || ($P8_86_31 != $P8_11 and $P8_86_31 != $P86_11) || ($P8_86_32 != $P8_11 and $P8_86_32 != $P86_11)); 
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 7: Filter out the loci that the three replicates of 86IT1 × 879 contain genotypes both absent in 86IT1 and 879
open OUT,">./all.filtered_snps.geno.update.content.filterFill3.2.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterFill3.1.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P86_8_11,$P86_8_12)=split /\//,$P86_8_1;
($P86_8_21,$P86_8_22)=split /\//,$P86_8_2;
($P86_8_31,$P86_8_32)=split /\//,$P86_8_3;

($P86_11,$P86_12)=split /\//,$P86_1;
($P86_21,$P86_22)=split /\//,$P86_2;
($P86_31,$P86_32)=split /\//,$P86_3;

#Skip this line if the two alleles of any replicate of 86IT1 × 879 contain genotype both absent in 86IT1 and 879
next if (($P86_8_11 != $P8_11 and $P86_8_11 != $P86_11) || ($P86_8_12 != $P8_11 and $P86_8_12 != $P86_11) || ($P86_8_21 != $P8_11 and $P86_8_21 != $P86_11) || ($P86_8_22 != $P8_11 and $P86_8_22 != $P86_11) || ($P86_8_31 != $P8_11 and $P86_8_31 != $P86_11) || ($P86_8_32 != $P8_11 and $P86_8_32 != $P86_11));
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 8: Filter out the loci that the three replicates of 879 × 928 contain genotypes both absent in 928 and 879
open OUT,">./all.filtered_snps.geno.update.content.filterFill3.3.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterFill3.2.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P8_9_11,$P8_9_12)=split /\//,$P8_9_1;
($P8_9_21,$P8_9_22)=split /\//,$P8_9_2;
($P8_9_31,$P8_9_32)=split /\//,$P8_9_3;

($P9_11,$P9_12)=split /\//,$P9_1;
($P9_21,$P9_22)=split /\//,$P9_2;
($P9_31,$P9_32)=split /\//,$P9_3;

#Skip this line if the two alleles of any replicate of 879 × 928 contain genotype both absent in 928 and 879
next if (($P8_9_11 != $P8_11 and $P8_9_11 != $P9_11) || ($P8_9_12 != $P8_11 and $P8_9_12 != $P9_11) || ($P8_9_21 != $P8_11 and $P8_9_21 != $P9_11) || ($P8_9_22 != $P8_11 and $P8_9_22 != $P9_11) || ($P8_9_31 != $P8_11 and $P8_9_31 != $P9_11) || ($P8_9_32 != $P8_11 and $P8_9_32 != $P9_11));
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 9: Filter out the loci that the three replicates of 928 × 879 contain genotypes both absent in 928 and 879
open OUT,">./all.filtered_snps.geno.update.content.filterFill3.4.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterFill3.3.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P9_8_11,$P9_8_12)=split /\//,$P9_8_1;
($P9_8_21,$P9_8_22)=split /\//,$P9_8_2;
($P9_8_31,$P9_8_32)=split /\//,$P9_8_3;

($P9_11,$P9_12)=split /\//,$P9_1;
($P9_21,$P9_22)=split /\//,$P9_2;
($P9_31,$P9_32)=split /\//,$P9_3;

#Skip this line if the two alleles of any replicate of 928 × 879 contain genotype both absent in 928 and 879
next if (($P9_8_11 != $P8_11 and $P9_8_11 != $P9_11) || ($P9_8_12 != $P8_11 and $P9_8_12 != $P9_11) || ($P9_8_21 != $P8_11 and $P9_8_21 != $P9_11) || ($P9_8_22 != $P8_11 and $P9_8_22 != $P9_11) || ($P9_8_31 != $P8_11 and $P9_8_31 != $P9_11) || ($P9_8_32 != $P8_11 and $P9_8_32 != $P9_11));
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 10: Filter out the loci that the three replicates of 879 × MTE contain genotypes both absent in MTE and 879
open OUT,">./all.filtered_snps.geno.update.content.filterFill3.5.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterFill3.4.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P8_M_11,$P8_M_12)=split /\//,$P8_M_1;
($P8_M_21,$P8_M_22)=split /\//,$P8_M_2;
($P8_M_31,$P8_M_32)=split /\//,$P8_M_3;

($PM_11,$PM_12)=split /\//,$PM_1;
($PM_21,$PM_22)=split /\//,$PM_2;
($PM_31,$PM_32)=split /\//,$PM_3;

#Skip this line if the two alleles of any replicate of 879 × MTE contain genotype both absent in MTE and 879
next if (($P8_M_11 != $P8_11 and $P8_M_11 != $PM_11) || ($P8_M_12 != $P8_11 and $P8_M_12 != $PM_11) || ($P8_M_21 != $P8_11 and $P8_M_21 != $PM_11) || ($P8_M_22 != $P8_11 and $P8_M_22 != $PM_11) || ($P8_M_31 != $P8_11 and $P8_M_31 != $PM_11) || ($P8_M_32 != $P8_11 and $P8_M_32 != $PM_11));
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 11: Filter out the loci that the three replicates of MTE × 879 contain genotypes both absent in MTE and 879
open OUT,">./all.filtered_snps.geno.update.content.filterFill3.6.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterFill3.5.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($PM_8_11,$PM_8_12)=split /\//,$PM_8_1;
($PM_8_21,$PM_8_22)=split /\//,$PM_8_2;
($PM_8_31,$PM_8_32)=split /\//,$PM_8_3;

($PM_11,$PM_12)=split /\//,$PM_1;
($PM_21,$PM_22)=split /\//,$PM_2;
($PM_31,$PM_32)=split /\//,$PM_3;

#Skip this line if the two alleles of any replicate of MTE × 879 contain genotype both absent in MTE and 879
next if (($PM_8_11 != $P8_11 and $PM_8_11 != $PM_11) || ($PM_8_12 != $P8_11 and $PM_8_12 != $PM_11) || ($PM_8_21 != $P8_11 and $PM_8_21 != $PM_11) || ($PM_8_22 != $P8_11 and $PM_8_22 != $PM_11) || ($PM_8_31 != $P8_11 and $PM_8_31 != $PM_11) || ($PM_8_32 != $P8_11 and $PM_8_32 != $PM_11));
print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_1\t$P8_86_2\t$P8_86_3\t$P8_9_1\t$P8_9_2\t$P8_9_3\t$P8_M_1\t$P8_M_2\t$P8_M_3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_1\t$P86_8_2\t$P86_8_3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_1\t$P9_8_2\t$P9_8_3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_1\t$PM_8_2\t$PM_8_3\n";
}
close(CLUSTER); 
close(FASTA);
close(OUT);

#Step 12: Distinguish maternal and paternal alleles in each F1 hybrid according to the genotypes of its parents, i.e. this step is for genotype phasing. #e.g. if the genotypes at a locus in 879 and 86IT1 were "1/1" and "0/0", respectively, the genotypes in 86IT1 × 879 should be "0|1", and in 879 × 86IT1 should be "1|0"
open OUT,">./all.filtered_snps.geno.update.content.geno.vcf";
open(FASTA,"./all.filtered_snps.geno.update.content.filterFill3.6.vcf"); #Input file is the output file of the last step
while(<FASTA>){
chomp;
my @fields = split /\t/;
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$P8_1,$P8_2,$P8_3,$P8_86_1,$P8_86_2,$P8_86_3,$P8_9_1,$P8_9_2,$P8_9_3,$P8_M_1,$P8_M_2,$P8_M_3,$P86_1,$P86_2,$P86_3,$P86_8_1,$P86_8_2,$P86_8_3,$P9_1,$P9_2,$P9_3,$P9_8_1,$P9_8_2,$P9_8_3,$PM_1,$PM_2,$PM_3,$PM_8_1,$PM_8_2,$PM_8_3) = @fields;

($P8_11,$P8_12)=split /\//,$P8_1;
($P8_21,$P8_22)=split /\//,$P8_2;
($P8_31,$P8_32)=split /\//,$P8_3;

($P8_86_11,$P8_86_12)=split /\//,$P8_86_1;
($P8_86_21,$P8_86_22)=split /\//,$P8_86_2;
($P8_86_31,$P8_86_32)=split /\//,$P8_86_3;

($P86_8_11,$P86_8_12)=split /\//,$P86_8_1;
($P86_8_21,$P86_8_22)=split /\//,$P86_8_2;
($P86_8_31,$P86_8_32)=split /\//,$P86_8_3;

($P86_11,$P86_12)=split /\//,$P86_1;
($P86_21,$P86_22)=split /\//,$P86_2;
($P86_31,$P86_32)=split /\//,$P86_3;

($P8_9_11,$P8_9_12)=split /\//,$P8_9_1;
($P8_9_21,$P8_9_22)=split /\//,$P8_9_2;
($P8_9_31,$P8_9_32)=split /\//,$P8_9_3;

($P9_8_11,$P9_8_12)=split /\//,$P9_8_1;
($P9_8_21,$P9_8_22)=split /\//,$P9_8_2;
($P9_8_31,$P9_8_32)=split /\//,$P9_8_3;

($P9_11,$P9_12)=split /\//,$P9_1;
($P9_21,$P9_22)=split /\//,$P9_2;
($P9_31,$P9_32)=split /\//,$P9_3;

($P8_M_11,$P8_M_12)=split /\//,$P8_M_1;
($P8_M_21,$P8_M_22)=split /\//,$P8_M_2;
($P8_M_31,$P8_M_32)=split /\//,$P8_M_3;

($PM_8_11,$PM_8_12)=split /\//,$PM_8_1;
($PM_8_21,$PM_8_22)=split /\//,$PM_8_2;
($PM_8_31,$PM_8_32)=split /\//,$PM_8_3;

($PM_11,$PM_12)=split /\//,$PM_1;
($PM_21,$PM_22)=split /\//,$PM_2;
($PM_31,$PM_32)=split /\//,$PM_3;

$P8_86_F1 = (($P8_86_11 != $P8_86_12) and ($P8_86_11 eq $P86_11))? "$P8_86_12\/$P8_86_11" : "$P8_86_1";
$P8_86_F2 = (($P8_86_21 != $P8_86_22) and ($P8_86_21 eq $P86_11))? "$P8_86_22\/$P8_86_21" : "$P8_86_2";
$P8_86_F3 = (($P8_86_31 != $P8_86_32) and ($P8_86_31 eq $P86_11))? "$P8_86_32\/$P8_86_31" : "$P8_86_3";

$P86_8_F1 = (($P86_8_11 != $P86_8_12) and ($P86_8_11 eq $P8_11))? "$P86_8_12\/$P86_8_11" : "$P86_8_1";
$P86_8_F2 = (($P86_8_21 != $P86_8_22) and ($P86_8_21 eq $P8_11))? "$P86_8_22\/$P86_8_21" : "$P86_8_2";
$P86_8_F3 = (($P86_8_31 != $P86_8_32) and ($P86_8_31 eq $P8_11))? "$P86_8_32\/$P86_8_31" : "$P86_8_3";

$P8_9_F1 = (($P8_9_11 != $P8_9_12) and ($P8_9_11 eq $P9_11))? "$P8_9_12\/$P8_9_11" : "$P8_9_1";
$P8_9_F2 = (($P8_9_21 != $P8_9_22) and ($P8_9_21 eq $P9_11))? "$P8_9_22\/$P8_9_21" : "$P8_9_2";
$P8_9_F3 = (($P8_9_31 != $P8_9_32) and ($P8_9_31 eq $P9_11))? "$P8_9_32\/$P8_9_31" : "$P8_9_3";

$P9_8_F1 = (($P9_8_11 != $P9_8_12) and ($P9_8_11 eq $P8_11))? "$P9_8_12\/$P9_8_11" : "$P9_8_1";
$P9_8_F2 = (($P9_8_21 != $P9_8_22) and ($P9_8_21 eq $P8_11))? "$P9_8_22\/$P9_8_21" : "$P9_8_2";
$P9_8_F3 = (($P9_8_31 != $P9_8_32) and ($P9_8_31 eq $P8_11))? "$P9_8_32\/$P9_8_31" : "$P9_8_3";

$P8_M_F1 = (($P8_M_11 != $P8_M_12) and ($P8_M_11 eq $PM_11))? "$P8_M_12\/$P8_M_11" : "$P8_M_1";
$P8_M_F2 = (($P8_M_21 != $P8_M_22) and ($P8_M_21 eq $PM_11))? "$P8_M_22\/$P8_M_21" : "$P8_M_2";
$P8_M_F3 = (($P8_M_31 != $P8_M_32) and ($P8_M_31 eq $PM_11))? "$P8_M_32\/$P8_M_31" : "$P8_M_3";

$PM_8_F1 = (($PM_8_11 != $PM_8_12) and ($PM_8_11 eq $P8_11))? "$PM_8_12\/$PM_8_11" : "$PM_8_1";
$PM_8_F2 = (($PM_8_21 != $PM_8_22) and ($PM_8_21 eq $P8_11))? "$PM_8_22\/$PM_8_21" : "$PM_8_2";
$PM_8_F3 = (($PM_8_31 != $PM_8_32) and ($PM_8_31 eq $P8_11))? "$PM_8_32\/$PM_8_31" : "$PM_8_3";

print  OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$P8_1\t$P8_2\t$P8_3\t$P8_86_F1\t$P8_86_F2\t$P8_86_F3\t$P8_9_F1\t$P8_9_F2\t$P8_9_F3\t$P8_M_F1\t$P8_M_F2\t$P8_M_F3\t$P86_1\t$P86_2\t$P86_3\t$P86_8_F1\t$P86_8_F2\t$P86_8_F3\t$P9_1\t$P9_2\t$P9_3\t$P9_8_F1\t$P9_8_F2\t$P9_8_F3\t$PM_1\t$PM_2\t$PM_3\t$PM_8_F1\t$PM_8_F2\t$PM_8_F3\n";
}
close(CLUSTER); 
close(FASTA);
close(OUT);