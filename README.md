<!-- ABOUT THE PROJECT -->
<a name="readme-top"></a>
## About PathoGD

PathoGD is a bioinformatic pipeline for the rapid and high-throughput design of RPA primers and guide RNAs (gRNA) for CRISPR-Cas12a-based nucleic acid detection. It incorporates two complementary modules- pangenome and *k*-mer- targeting the protein-coding and whole genome, respectively, for selection of a target region and subsequent guide RNA and primer design. PathoGD was initially developed for and tested on bacterial pathogens, but should work on any organism whose genome sequences are available.


## Introduction

The increasing availability of bacterial draft and whole-genome sequences provide an opportunity to identify alternative diagnostic markers from previously unexplored genomic regions. PathoGD was developed to facilitate the design of highly specific RPA primers and gRNAs for CRISPR-Cas12a-based pathogen detection. It provides a streamlined workflow to accelerate the *in silico* aspect of a CRISPR-Cas12a assay design, a task which can be very manual and labor-intensive. 

PathoGD incorporates multiple bioinformatics tools for performing the individual steps in the pipeline. It requires two databases- target and non-target- containing the genomes of the pathogen to be detected and closely-related organisms potentially resulting in cross-reactivity, respectively, for the selection of a highly-specific target region to be amplified.


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<a name="installation-top"></a>
## Installation

1. Clone the repo into your working directory
   ```sh
   git clone https://github.com/sjlow23/pathogd.git
   ```
2. Create a conda environment using the `pathogd.yml` file provided.

   _Installing to default conda environment path:_
   ```sh
   conda env create -n pathogd -f pathogd.yaml
   ```

   _Installing to user-specified directory:_
   ```sh
   conda env create --prefix /dir/to/pathogd -f pathogd.yaml
   ```
3. Move scripts into bin directory of conda environment
   ```sh
   mv pathogd $CONDA_PREFIX/bin/
   mv scripts/* $CONDA_PREFIX/bin/
   ```

4. Make scripts executable
   ```sh
   chmod u+x $CONDA_PREFIX/bin/pathogd $CONDA_PREFIX/bin/*.R
   ```

5. You will also need Prokka and Roary installed on your system and in your $PATH. See [Prokka](https://github.com/tseemann/prokka "Prokka") and [Roary](https://github.com/sanger-pathogens/Roary "Roary") for installation instructions.


<p align="right">(<a href="#installation-top">back to top</a>)</p>


<a name="modules-top"></a>
## Modules and workflows

Modules refer to the approach that is used for primer and gRNA design. The `-m` parameter is used for specifying which module is to be used. This parameter is mandatory only for workflows (see below) that require primer and gRNA design. Two modules are available in PathoGD:

1. `pangenome` - primers and gRNAs are designed from protein-coding regions only
2. `k-mer` - primers and gRNAs are designed from entire genome including non-coding regions

Workflows refer to the set of functions to be performed in the PathoGD pipeline. The `-w` parameter is used for specifying which workflow is to be used. This parameter is mandatory in the `pathogd` command. 


The following workflows perform specific functions only; no primer and gRNA design is executed:

1. `check` - check number of NCBI assemblies available for target and non-target taxa
2. `download_target` - only download NCBI target genome assemblies
3. `download_nontarget` - only download NCBI non-target genome assemblies


The following workflows are available for primer and gRNA design:

1. `ncbi_all_subsample` - download target and non-target genomes from NCBI; run subsampling of target genomes
2. `ncbi_all_nosubsample` - download target and non-target genomes from NCBI; use all target genomes
3. `user_target_subsample` - user-provided target genomes; download non-target genomes from NCBI; run subsampling of target genomes
4. `user_target_nosubsample` - user-provided target genomes; download non-target genomes from NCBI; use all target genomes
5. `user_all_subsample` - user-provided target and non-target genomes; run subsampling of target genomes
6. `user_all_nosubsample` - user-provided target and non-target genomes; use all target genomes


#### Integrated tools

| **Software/Tool**                                                  | **Function**         | **Module**            |
|--------------------------------------------------------------------|----------------------|-----------------------|
| ncbi-genome-download                                               | Genome download      | pangenome and *k*-mer |
| Prodigal, Prokka                                                   | Gene annotation      | pangenome             |
| Roary                                                              | Pangenome estimation | pangenome             |
| genometester4                                                      | *k*-mer enumeration  | pangenome and *k*-mer |
| MAFFT, BBMap                                                       | Sequence alignment   | pangenome and *k*-mer |
| cd-hit                                                             | Sequence clustering  | pangenome and *k*-mer |
| Primer3                                                            | Primer design        | pangenome and *k*-mer |
| isPcr                                                              | *in silico* PCR      | pangenome and *k*-mer |
| bedtools, cdbtools, csvtk,   taxonkit, emboss, samtools, seqkit, R | Data processing      | pangenome and *k*-mer |


<p align="right">(<a href="#modules-top">back to top</a>)</p>


<a name="config-top"></a>
## Configuration file

The configuration file is specified using the `-c` parameter in the `pathogd` command, and is mandatory for running PathoGD. This file is used for specifying most of the parameters and options for running the pipeline, including target and non-target organisms. A template of the file is provided (`config.txt`) above. 

Do not modify the file other than to include your specifications. Input specifications after the `|` sign. The following table defines which fields need to be completed for running the PathoGD pipeline. Asterisks indicate fields that are optional, but would be good to include for logging purposes. These do not have to be accurate as they are not explicitly used for running the pipeline.


| **Parameter/Workflow** | **check** | **download_target** | **download_nontarget** | **ncbi_all_subsample**     | **ncbi_all_nosubsample**   | **user_target_subsample**  | **user_target_nosubsample** | **user_all_subsample**     | **user_all_nosubsample**   |
|------------------------|-----------|---------------------|------------------------|----------------------------|----------------------------|----------------------------|-----------------------------|----------------------------|----------------------------|
| db                     | no        | yes                 | yes                    | yes                        | yes                        | yes                        | yes                         | no                         | no                         |
| target                 | no        | yes                 | no*                    | yes                        | yes                        | no                         | no                          | no                         | no                         |
| target_taxid           | yes       | yes                 | no*                    | yes                        | yes                        | no                         | no                          | no                         | no                         |
| offtarget              | no        | no*                 | yes                    | yes                        | yes                        | yes                        | yes                         | no                         | no                         |
| offtarget_taxid        | yes       | no*                 | yes                    | yes                        | yes                        | yes                        | yes                         | no                         | no                         |
| domain                 | yes       | yes                 | yes                    | yes                        | yes                        | yes                        | yes                         | no                         | no                         |
| assembly_level         | no        | yes                 | yes                    | yes                        | yes                        | yes                        | yes                         | no                         | no                         |
| kmer                   | no        | no                  | no                     | yes                        | yes                        | yes                        | yes                         | yes                        | yes                        |
| threshold              | no        | no                  | no                     | yes                        | yes                        | yes                        | yes                         | yes                        | yes                        |
| mismatch               | no        | no                  | no                     | yes                        | yes                        | yes                        | yes                         | yes                        | yes                        |
| primer_design          | no        | no                  | no                     | yes                        | yes                        | yes                        | yes                         | yes                        | yes                        |
| ref_annot              | no        | no                  | no                     | yes                        | yes                        | yes                        | yes                         | yes                        | yes                        |
| subsample              | no        | no                  | no                     | yes                        | no                         | yes                        | no                          | yes                        | no                         |
| prokka_env             | no        | no                  | no                     | yes for 'pangenome' module | yes for 'pangenome' module | yes for 'pangenome' module | yes for 'pangenome' module  | yes for 'pangenome' module | yes for 'pangenome' module |
| roary_env              | no        | no                  | no                     | yes for 'pangenome' module | yes for 'pangenome' module | yes for 'pangenome' module | yes for 'pangenome' module  | yes for 'pangenome' module | yes for 'pangenome' module || no                       |


<p align="right">(<a href="#config-top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
<a name="usage-top"></a>
## Getting started

The command for running PathoGD is:

   ```sh
   pathogd [-c <config file>] [-m <method (optional)>] [-w <workflow>] [-o <output directory>] [-r <reference genbank file (optional)>] [-t <number of cpus (default 8)>]
   ```

To see all options, run `pathogd -h`. Some examples of usage are described below:

### A. Check number of available NCBI GenBank and RefSeq genome assemblies for taxon of interest using `check` workflow

Use this for determining the number of genomes available for your target and non-target taxa before running pipeline. This workflow only works if your target and non-target taxa consist of distinct species.

To use this workflow, you will need to identify the NCBI taxonomy identifiers for your target and non-target taxa. Follow the steps below to obtain the NCBI taxids.

If your target or non-target taxa is a single species:

   **a) Identify species taxid of target taxon (single species)**

   Go to [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi "Taxonomy"), and search for your species of interest by name. Clicking on the species name will take you to a new page with the Taxonomy ID.
   For example, the species taxid for *Mycoplasmoides genitalium* is **2097**.

If your target or non-target taxa includes multiple species, *eg.* all species belonging to a genus:

   **b) Identify all species taxids for a parent taxon comprising multiple species**
   
  Go to [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi "Taxonomy"), and search for your genus (or other taxonomic rank) of interest by name. Clicking on the genus name will take you to a new page with the Taxonomy ID.
  For example, the taxid for *Mycoplasmoides* is **2995234**. 

  Specifying just the parent taxid will be sufficient as the pipeline will search for all children taxa.

   **c) Input taxids appropriate fields in config file**

   **d) Run the following command:**

   ```sh
   pathogd -c config.txt -w check -o pathogd_output
   ```

   **Important notes:**

   1. If target is a single species, ensure that the taxid provided is for the 'species' rank. Some organisms may have been reported at a sub-species or strain level, in which case the taxid will be different from the species taxid.

   2. To obtain genomes at taxonomic ranks below species level, eg. sub-species or strain level, please download genomes manually. See [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download "ncbi-genome-download") or [ncbi-datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ "ncbi_datasets") for options to download genomes using NCBI taxids and/or taxon names. 

   The outputs of the `check` workflow are:

   * Number of target and non-target genomes available (specified in log file)
   * Scripts for downloading target and non-target genomes for GenBank and/or RefSeq assemblies. Random subsampling is automatically performed to retain a maximum of 1000 and 100 genomes for each target and non-target species, respectively.


<p align="right">(<a href="#usage-top">back to top</a>)</p>


### B. Download genomes for target and/or non-target taxa of interest using the `download_target` or `download_nontarget` workflow

Use this to automatically download genome assemblies for target and/or non-target taxa by specifying their species names/taxids in the config file. Similar to the `check` workflow, this workflow only works if your target and non-target taxa consist of distinct species.

   _Download target genomes_

   ```sh
   pathogd -c config.txt -w download_target -o pathogd_output
   ```

   _Download non-target genomes_

   ```sh
   pathogd -c config.txt -w download_nontarget -o pathogd_output
   ```

   **Important notes:**

   1. No subsampling is performed for both target and non-target taxa using this workflow. If specifying species with a large number of genome assemblies, *eg.* *Escherichia coli*, all genomes will be downloaded which may not be ideal if there is insufficient storage.

   2. To prevent this, use the `check` workflow.

   3. If both species name and species taxids are provided in the config file, only the **species taxids** will be used for downloading genome assemblies.


<p align="right">(<a href="#usage-top">back to top</a>)</p>


### C. Design primers and guide RNAs using any of the 6 `ncbi_*` or `user_*` workflows for primer and gRNA design

Examples using *Mycoplasmoides genitalium* for either pangenome or *k*-mer module:

   **a) Download target and non-target genome assemblies from NCBI, no subsampling of target genomes**

   Pangenome module:

   ```sh
   pathogd -c config.txt -m pangenome -w ncbi_all_nosubsample -o pathogd_output
   ```

   *k*-mer module:

   ```sh
   pathogd -c config.txt -m kmer -w ncbi_all_nosubsample -o pathogd_output
   ```

   The following examples are for the pangenome module. To use the *k*-mer module, replace `-m pangenome` with `-m kmer` in the commands below: 

   **b) Download genome assemblies from NCBI, subsample target genomes to total specified in config file**

   ```sh
   pathogd -c config.txt -m pangenome -w ncbi_all_subsample -o pathogd_output
   ```


   **c) Provide own target and non-target genomes, no subsampling of target genomes**
   
   ```sh
   pathogd -c config.txt -m pangenome -w user_all_nosubsample -o pathogd_output
   ```


   **d) Provide own target genomes, download non-target genomes from NCBI, no subsampling of target genomes**
   
   ```sh
   pathogd -c config.txt -m pangenome -w user_target_nosubsample -o pathogd_output
   ```


   **Important notes:**

   1. If using the `ncbi_all_nosubsample` workflow, all target and non-target genomes will be downloaded which may not be ideal if the species is overrepresented, eg. *E. coli*. 

   2. To prevent this, use the `check` workflow, or provide own genomes.


<p align="right">(<a href="#usage-top">back to top</a>)</p>


## Output

For the `ncbi_*` and `user_*` workflows with primer and gRNA design, the most relevant output is a tab-delimited file `pathogd_primers_stats.tsv` containing information on all the primers and gRNAs that were designed. The table below explains the contents of each field:

| **Field**                         | **Remarks**                                                                                                                                   |
|-----------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|
| guide_set                         | Guide RNA id                                                                                                                                  |
| guide_seq                         | Guide RNA sequence                                                                                                                            |
| pam                               | Protospacer adjacent motif (sequence upstream of gRNA)                                                                                        |
| gene                              | Gene id                                                                                                                                       |
| gene_annotation                   | Gene name                                                                                                                                     |
| primer_set                        | Primer set id for a given guide RNA id (maximum of 5 sets per gRNA; name not unique)                                                          |
| fwd_primer                        | Forward primer sequence                                                                                                                       |
| rev_primer                        | Reverse primer sequence                                                                                                                       |
| primer_cas13a_compatibility       | Compatibility of primers with Cas13a (if 'yes', can use as is; if 'swap primers', swap forward and reverse primers). A T7 RNA polymerase promoter sequence is added on the 5' end of the forward primer for T7 transcription.                         |
| fwd_primer_length                 | Length of forward primer                                                                                                                      |
| rev_primer_length                 | Length of reverse primer                                                                                                                      |
| guide_length                      | Length of guide RNA                                                                                                                           |
| fwd_primer_gc                     | Forward primer GC content                                                                                                                     |
| rev_primer_gc                     | Reverse primer GC content                                                                                                                     |
| guide_gc                          | Guide RNA GC content                                                                                                                          |
| product_size                      | Amplicon size (based on consensus sequence)                                                                                                   |
| primerset_new                     | Unique primer set id (based on combination of primer and gRNA sequence; should occur only once)                                             |
| primerset_uniq                    | Primer set id (based on forward and reverse primer sequence; can occur more than once if multiple gRNAs are compatible with the primer set) |
| primer_target_prevalence_0mm      | Percentage of target genomes with perfect matches to primer sequence                                                                          |
| primer_target_prevalence_2mm      | Percentage of target genomes with up to 2 mismatches to primer sequence                                                                       |
| avg_amplicon_num_target           | Average number of amplicons in target genomes                                                                                                 |
| avg_amplicon_length_target        | Average amplicon size in target genomes                                                                                                       |
| primer_nontarget_prevalence_0mm   | Percentage of non-target genomes with perfect matches to primer sequence                                                                      |
| primer_nontarget_prevalence_2mm   | Percentage of non-target genomes with up to 2 mismatches to primer sequence                                                                   |
| primer_nontarget_prevalence_8mm   | Percentage of non-target genomes with up to 8 mismatches to primer sequence                                                                   |
| avg_amplicon_num_nontarget        | Average number of amplicons in non-target genomes                                                                                             |
| avg_amplicon_length_nontarget     | Average amplicon size in non-target genomes                                                                                                   |
| guide_average_copy_number_target  | Average number of guide RNAs in target genomes                                                                                                |
| guide_target_prevalence           | Percentage of target genomes with perfect matches to guide RNA sequence                                                                       |
| guide_offtarget_prevalence_max5mm | Percentage of non-target genomes with 3≤n≤5 mismatches to guide RNA sequence                                                                  |


<p align="right">(<a href="#readme-top">back to top</a>)</p>



## Citation

If you use this tool in your research, please cite https://github.com/sjlow23/pathogd.