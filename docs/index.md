---
# You don't need to edit this file, it's empty on purpose.
# Edit theme's home layout instead if you wanna make some changes
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults

title: Overview
layout: "default"
weight: 1
---

<h1><p style="text-align: center">Overview</p></h1>

-----
<br>
USDA APHIS Veterinary Services (VS) Mycobacterium tuberculosis complex, mainly M. bovis, and Brucella sp. genotyping from whole genome sequence (WGS) outputting BAM, VCF, SNP tables and phylogentic trees. 


Using BWA alignments VCFs are created for each isolate using GATK's haplotype caller.  The overview below discribes a two step pipeline, briefly: step 1 - outputing VCFs and step 2 - gathering those VCFs to produce SNP alignments, tables and phylogenetic trees of grouped isolates.  This two step pipeline has been developed to provide routine diagnostics results rapidly, with quality assurance and easy to intrepret reporting.

# Script 1
Step 1 is fairly straight forward.  Our main workflow includes Mycobacterium tuberculosis complex and Brucella sp. and therefore the script has been optimized for such.  The script begins by selecting the best reference and determining the spoliogtype for TB complex isolates and MLST for Brucella spieces.  The script then goes on to using what is now a fairly standard bacterial SNP calling pipeline where BWA aligns and GATK calls SNPs.  However, in addition to what is output by the HaplotypeCaller, zero coverage positions are added to VCFs.  This is not a standard GATK option and is done using python.  Including zero coverage positions allows a more accurate SNP summary to be represented in step 2.

# Script 2
Step 2 is called on VCFs output from step 1.  References chosen in script 1 have been selected because they have been found to be relatively close to the isolate.  The closer the reference is to the isolate the less overall SNP calling error is seen.  VCFs analyzed in step 2 must all be output from the same reference.  Obviously VCFs analyzed using different references can not be used in the same comparison.

In addition to choosing a closely related reference, which minimizes SNP calling error, there are three additional external files, or dependencies, used to create high quality, informative SNP alignments.  As shown in bovis_dependency_view.jpg and suis1_dependency_view.jpg the reference used to build VCFs must be reflected in the three dependencies.  The three dependent files are: filter file, defining SNPs, and gbk file.
