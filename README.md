# RTDmaker

Pipeline for the creation of High-Quality Reference Transcriptome Datasets.

----------------------------
# Table of Contents
----------------------------

   * [Overview](#overview)
   * [Installation](#installation)
   * [Modules](#modules)
   * [ShortReads](#shortreads)
      * [Input files](#input-files)
      * [Command and options](#command-and-options)
      * [Output files](#output-files)
   * [FAQ](#faq)
   * [Future work](#future-work)
   * [License](#license)
   * [Contact](#contact)


----------------------------
# Overview
----------------------------

RTDmaker is a computational pipeline to generate High-Quality transcriptome annotations, known as Reference-Transcript-Datasets (RTDs). Currently, RTDmaker consist of one module (ShortReads) to process transcript assemblies from RNA-seq data through a stringent filtering process with the option to incorporate available reference annotations (Fig. 1).

This module identifies and removes: 

1. Redundant transcripts 
   - Transcripts with identical combination of intron combinations of another transcript. The start coordinate of the first exon or the end coordinate of the last exon could be different. (Fig. 2)

2. Transcripts with ambiguous information
   - Unstranded transcripts for stranded RNA-seq data
   - Ttranscripts with annotations of wrong scaffold/chromosome names  comparing to the provided genome reference. 
   - Mono-exonic fragments that are antisense of a longer gene, which could be artifacts of strand labelling process (Fig. 3)

3. Transcripts poorly supported by reads
   - Transcripts containing Splice-Junctions with an insufficient number of uniquely mapped reads.
   - Transcripts with low TPM (transcript per million reads) abundance in a certain number of samples.

4. Fragmentary and chimeric transcripts
   - Fragments: transcripts whose length is much shorter than the encompassing locus.
   - Chimeric: transcripts which overlap two or more distinct loci (Fig. 4)


![Fig1_RTDmaker_Readme.PNG](https://drive.google.com/uc?export=view&id=1WArvkKiOamBOauXrjqH4IzFDoj4J1DNl)
**Figure 1.** **RTDmaker ShortReads workflow.** (**A**) Assembly QC checks the quality and read support of newly assembled annotations. (**B**) Merge QC integrates information from curated new assemblies and available references by applying a series of QC checks to remove redundant, fragmentary and chimeric transcripts. Finally, **Merge QC** also performs the necessary post-processing to generate the final RTD (i.e. gene re-annotation).



![Fig2_RTDmaker_Readme.PNG](https://drive.google.com/uc?export=view&id=1XbFl_KoeLfNg2Fbiz5utr0TTwFDBmjfE)
**Figure 2.**  **Filtering of redundant transcripts and transcript fragments.** (**A**) Assembled transcripts with identical splice junctions – shorter transcripts removed; (**B**) shorter transcripts of mono-exonic genes removed; (**C**) assembled fragments of transcripts with identical or new splice junctions removed.



![Fig3_RTDmaker_Readme_PGSC0003DMG401012750.png](https://drive.google.com/uc?export=view&id=1uMriO64596MvObyhXvqPzeHbQu1PeiLF)
**Figure 3.**  **Artefactual antisense mono-exonic genes.** The read distribution plot from a potato CHX treatment experiment (top, black) supports the reference PGSC0003DMG401012750 gene (conserved gene of unknown function) (middle, green) and its newly assembled isoforms (bottom, blue). However, the transcriptome assemblers also generate “novel” mono-exonic genes on the opposite strand – in this case in the region of the first and last exon. The PGSC reference also reported the mono-exonic transcript PGSC0003DMT400033204 (right, green) which is also annotated to the same gene PGSC0003DMG401012750, even though it is located in the opposite strand. The read distribution on this locus and the match of the mono-exonic “genes” to exons of the transcripts in the opposite strand, strongly suggest that these mono-exonic models are most likely to be mis-assembled models (both newly assembled and already annotated).



![Fig4_RTDmaker_Readme.PNG](https://drive.google.com/uc?export=view&id=1G_B51paA-6LnG5DyqMr0aDJvxHnKERTY)
**Figure 4.** **Filtering of fragmentary and chimeric transcripts.** The pipeline groups transcripts into loci with similar start-end exons (**white**). Then, for each locus it removes (as fragments) encompassed transcripts that are much shorter than the locus (**yellow**). Next, it identifies the group of transcripts overlapping two distinct loci (**grey**), and if the total number of overlapping transcripts are less than the total number of transcripts of all distinct locis, they are removed as chimeric (**grey**). Conversely, if the overlapping transcripts are more than those of the distinct groups (**white**), then these are removed as fragments and the longer, overlapping transcripts are kept.


----------------------------
# Installation
----------------------------

RTDmaker has been developed in Python 3.6 

RTDmaker requires the following packages:
- [Salmon v1.4.0](https://anaconda.org/bioconda/salmon)
- [Gffread v0.12.1](https://anaconda.org/bioconda/gffread)

Please beware that the installation of the specific version must be explicit to avoid errors, for example:
```
conda install salmon=1.4.0
```

RTDmaker dependencies require Linux OS.


RTDmaker is ready to use. The compressed file can be directly downloaded from the [GitHub repository](https://github.com/anonconda). Once decompressed, RTDmaker can be used directly from the command line by specifying the path to the main executable `RTDmaker.py`


----------------------------
# Modules
----------------------------

RTDmaker works with a module / module-options structure:

```
RTDmaker.py module options
```

Each module can be executed as follow:

```
python /path/to/RTDmaker.py module options
```

For example, to observed the help documentation of ShortReads module:

```
python /path/to/RTDmaker.py ShortReads --help
```



----------------------------
**ShortReads**
==============

----------------------------

RTDmaker 'ShortReads' module process newly assembled transcriptome annotations from RNA-seq data by using assembly tools, such as StringTie (Pertea et al., 2015) and Scallop (Shao et al., 2017), to generate a High-Quality RTD annotation (Fig. 1). 

This module identifies and removes: 
1. Redundant transcripts 
2. Transcripts with ambiguous location information
3. Transcripts poorly supported by reads 
4. Fragmentary and chimeric transcripts


## Input files
ShortReads takes the following inputs:

1. A folder containing available transcript annotations to integrate into the RTD, in GTF format (optional).
2. A folder containing newly assemble transcript annotations to merge and analyze, in GTF format.
3. A folder containing paired-end RNA-seq reads FASTQ files from which the new transcript assemblies were generated. 
4. A folder containing the splice junction SJ.out.tab files generate by the read aligner STAR. If the two-pass mode is used in STAR, the files from are automatically ignored.
5. The reference genome sequence FASTA file. 

- **Notes**: 

-  The program recursively parse the folders to get the GTF files. Thus, it is possible to populate these input folders with symbolic links of the original folders containing the GTF annotations to analyze.
-  Because most of the analysis are based on transcripts co-ordinates, it is important for the accuracy of the analyses that all the transcript annotations (new assemblies and available transcript references) must refer to the same Genome.
-  The FASTQ files will be autmatically paired to quantify the transcripts. A table reporting how the files are paired will be located in the report subfolder (See Output files below). This implies 2 requirements: (**1**) the folder MUST contain ONLY paired FASTQ files, and (**2**) the name of 2 PAIRED files must be alphabetically consecutive (Ex: EtOH4h-1_S1_R1_001.fq, EtOH4h-1_S1_R2_001.fq, etc). This is almost always already the case, just check the table ./report/fastq_table.csv in the report subfolder to confirm how the files has been paired.


## Command and options
Command to run ShortReads analysis:
```
python RTDmaker.py ShortReads [options]
```
```
python RTDmaker.py ShortReads --assemblies </path/to/assembly-folder> --references </path/to/references-folder> --SJ-data </path/to/SJ-folder> --SJ-reads <2,1> --genome <genome.fasta> --fastq </path/to/fastq-folder> --tpm <1, 1> --fragment-len <0.7> --antisense-len <0.5> --add <unstranded, intronic> --keep <intermediary, removed> --ram <8> --outpath </path/for/output-folder> --outname <outname> --prefix <prefix>
```

List of options available:
-**h**, **--help**:   show this help message and exit  

--**assemblies**: Path of the folder containing the new transcriptome assemblies (GTF format) to be analyzed. 

--**references**: Path of the folder containing reference annotations (GTF format) to be integrated into the RTD.

--**SJ-data**:  Path of the folder containing the splice-junction data ('SJ.out.tab' files) generated by STAR.

--**SJ-reads**:  Minimum number of uniquely-map reads (X) observed in a minimum number of samples (Y) to consider a splice junction (SJ) supported. Default: 2 1 (2 uniquely-map reads in at least 1 sample).

--**genome**: Path of the Genome FASTA file to extract the transcripts sequence to perform transcript quantification 

--**fastq**:  Path of the folder containing the FASTQ files to perform transcript quantification. Supported extensions: 'fq.gz', 'fq', 'fastq.gz', 'fastq'.  

--**tpm**:  Minimum TPM abundance (N), observed in a minimum number of samples (M), to consider a transcript as supported by quantification.Default: 1 1 (1 TPM in at least 1 sample).

--**add**:  Potentially non-coding transcripts to include into the annotation. Options: 'unstranded', 'intronic'. Default: exclude them.

--**fragment-len**: Fragment transcript in the gene model will be filtered if the length is less than this percentage of the gene length. Value must be within 0 and 1. Default: 0.7 

--**antisense-len**: Monoexonic-antisense fragment transcript in the gene model will be filtered if the transcript is in the opposite strand and the length is less than this percentage of the gene length. Value must be within 0 and 1. Default: 0.5 

--**ram**:  Maximum size (in Gb) of temporary files. Change this value only if your system cannot handle uploading into memory files with the current file size. Default: 8 Gb.

--**keep**:  Intermediary files to keep after the analysis for debug. Options: 'intermediary', 'removed'. Default: remove them.

--**prefix**:  Prefix to assign to the genes and transcripts IDs in the RTD. Default: outname.

--**outpath**:  Path of the output folder to save the results. 

--**outname**:  Prefix of the output file names.


- **NOTE**: The pipeline avoids the re-analysis of a QC step (Fig. 1) if it detects that the intermediary file resulting from that QC analysis already exist. Thus, if the intermediary files are kept (--**keep** intermediary), re-analysis of large datasets can be performed quickly by avoiding re-doing time-consuming steps (such as the **Redunancy QC** or the **transcript quantification**). This is useful to quickly re-start an analysis in the case of an abrupt interruption or to try a different value for a downstream QC step (i.e. no need to re-do transcripts quantification if you only want to try a different fragmentary length threshold). **HOWEVER**, this requires that the user **must manually delete** all the intermediary files that needs to be re-analyzed (which include **all the files downstream** of the file the user wish to re-analyze). 


Example:

```
python RTDmaker.py ShortReads --assemblies ./test_dataset/test_assemblies --references ./test_dataset/test_references --SJ-data ./test_dataset/test_SJdata --SJ-reads 2 1 --genome ./test_dataset/potato_dm_v404_all_pm_un.renamed.fasta --fastq ./test_dataset/test_fastq --tpm 0.1 1 --fragment-len 0.7 --antisense-len 0.5 --add intronic --keep intermediary --ram 8 --outpath ./test_dataset/outfolder --outname test_run --prefix MyRTD 
```

## Output files

RTDmaker ShortReads automatically generates an output folder (and subfolder tree) to store the files of the analysis:
/**&lt;outpath&gt;**/**&lt;outname&gt;**/


RTDmaker ShortReads generates the following output files:
1. A High-Quality RTD in *GTF* format and corresponding transcriptome FASTA
2. A padded version of the RTD to use for transcript quantification (*GTF* format) and corresponding transcriptome FASTA
3. Folder containing multiple files reporting information on: **1) the pipeline analysis** (i.e.: Nº of models per QC-step, logfiles of commands executed, etc), **2) input files** (i.e.: summary of references annodations, SJ dataset, etc), and **3) the final RTD**
4. Folder containing the intermediary files of the analysis (*GTF* files of accepted models at each QC-step, transcripts quantification files, etc)
5. Folder containing the rejected models at each quality control step (*GTF* format)


The name of the output files is generated as follows:
1. /**&lt;outpath&gt;**/**&lt;outname&gt;**\.gtf
2. /**&lt;outpath&gt;**/**&lt;outname&gt;**\_padded.gtf
3. /**&lt;outpath&gt;**/**&lt;outname&gt;**/report/
4. /**&lt;outpath&gt;**/**&lt;outname&gt;**/intermediary/
5. /**&lt;outpath&gt;**/**&lt;outname&gt;**/removed/


----------------------------
# FAQ
----------------------------

#### Can I run the pipeline only with newly assembled annotations / no references?

Yes, the --**references** argument is optional.


#### Can I run the pipeline only for one or more available transcript references / no newly assembled annotations?

Yes, if the --**assemblies** argument is not provided, the pipeline will proceed directly into the Merge QC subpipeline (Fig. 1B). 
As no assemblies are provided, it is not necessary to specify none of the "read support" arguments, i.e.:

```
python RTDmaker.py ShortReads --references </path/to/references-folder> --genome <genome.fasta> --fragment-len <0.7> --antisense-len <0.5> --add <unstranded, intronic> --keep <intermediary, removed> --ram <8> --outpath </path/for/output-folder> --outname <outname> --prefix <prefix>
```

In this case, the --**genome** becomes optional as it is only necessary to generate the transcriptome FASTA of the resulting RTD.


#### I don't have any 'SJ.out.tab' files / I use a different read aligner

RTDmaker ShortReads **requires** 'SJ.out.tab' files for the Splice-Junction support analysis.


#### I only have unpaired RNA-seq data

At the present time, RTDmaker ShortReads only support quantification with paired FASTQ data.


#### How is the transcript quantification run?

RTDmaker ShortReads perform an automated transcript quantification using Salmon (v1.4.0) as follow:

```
salmon quant -i <transcriptome_index> -l A -1 <fastq_P1.fq> -2 <fastq_P2.fq> -o <outfile> --useVBOpt --seqBias --threads 4
``` 

A logfile containing the list of executed commands is present in the 'report' subfolder.


#### Can I use different arguments for the Salmon transcript quantification?

Sure, if you are willing to modify the RTDmaker source code. The Salmon quantification command is located inside the '**generate_salmon_commands**' function on the file '*./RTDmaker/lib/tools/quantify_transcripts.py*'



----------------------------
# License
----------------------------

RTDmaker is released under the [MIT license](https://opensource.org/licenses/MIT)


----------------------------
# Contact
----------------------------

For any further enquiries please contact the main developer at <e.entizne@dundee.ac.uk>, <Juan.Carlos.Entizne@hutton.ac.uk>

