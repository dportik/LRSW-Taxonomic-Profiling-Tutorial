# Taxonomic profiling for long-read shotgun metagenomics datasets <a name="TOP"></a>

## **Table of Contents** 

**Overview**
+ [**Introduction**](#INTRO)
+ [**Example Dataset**](#EXDA)

**Analyses**
+ [**Kraken2 & Bracken**](#KRBR)
+ [**Centrifuge**](#CENT)
+ [**MMSeqs2**](#MMSEQ)
+ [**BugSeq**](#BS) <- top performing method
+ [**Diamond & MEGAN-LR**](#MEG) <- top performing method

**Summary**
+ [**Comparative Analysis**](#COMP)


## Introduction <a name="INTRO"></a>

This tutorial was created for the [**Long Read Sequencing Workshop**](https://www.jax.org/education-and-learning/education-calendar/2022/may/long-read-sequencing-workshop) (2022) at The Jackson Laboratory for Genomic Medicine. It is an introduction to several methods that can be used to identify the number of species and associated relative abundances from a long-read metagenomic dataset. The analyses covered here are a subset of those performed in the following benchmarking study:

[**Evaluation of taxonomic profiling methods for long-read shotgun metagenomic sequencing datasets**](https://www.biorxiv.org/content/10.1101/2022.01.31.478527v1)

The above study evaluated a total of 9 methods using four mock community datasets. To keep the tutorial simple, I will demonstrate how to use 6 of these methods to analyze one mock community dataset:

+ **Kraken2**
+ **Bracken**
+ **Centrifuge**
+ **MMseqs2**
+ **BugSeq**  <- top performing method
+ **Diamond + MEGAN-LR**  <- top performing method

The top-performing methods were **BugSeq** and **Diamond + MEGAN-LR**, both in terms of precision/recall and relative abundance estimates. **MMSeqs2** performed moderately well, whereas **Kraken2**, **Bracken**, and **Centrifuge** performed less well. There is plenty of room for the development of newer and better tools for profiling long-read datasets. However, these existing tools serve an important function and provide a baseline for future improvements.

A small tutorial is provided for each of the 6 methods above. Each tutorial includes instructions for installing the program and associated reference database, and the commands to run the program.

The simplest way to install programs and their dependencies into a single environment is to use [Anaconda](https://docs.anaconda.com/anaconda/)/[Conda](https://docs.conda.io/projects/conda/en/latest/index.html). Instructions for conda installation are therefore provided for each program. Manual installation is also possible; links to the program sources are provided. Databases will also need to be obtained for each program. 

The main objective is to obtain a kraken-style report (kreport) for each method. This format is highly useful because it contains cumulative counts and level counts across the complete hierarchical taxonomy for each taxon assigned. The level count is the number of reads specifically assigned to a taxon, whereas the cumulative count is the sum of the level counts for a taxon plus its descendants. For example, the cumulative count of a genus is the level count for that genus plus the level counts of all species and strains contained in that genus. The following columns are present in the kreport format:

+ `Percent Reads Contained in Taxon`: The cumulative percentage of reads for this taxon and all descendants.
+ `Number of Reads Contained in Taxon`: The cumulative number of reads for this taxon and all descendants.
+ `Number of Reads Assigned to Taxon`: The number of reads assigned directly to this taxon (not a cumulative count of all descendants).
+ `Rank Code`: (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. 
+ `NCBI Taxon ID`: Numerical ID from the NCBI taxonomy database.
+ `Scientific Name`: The scientific name of the taxon.

Example contents of a kreport:

```
88.41	2138742	193618	K	2	Bacteria
0.16	3852	818	P	201174	  Actinobacteria
0.13	3034	0	C	1760	    Actinomycetia
0.13	3034	45	O	85009	      Propionibacteriales
0.12	2989	1847	F	31957	        Propionibacteriaceae
0.05	1142	352	G	1912216	          Cutibacterium
0.03	790	790	S	1747	            Cutibacterium acnes
10.07	243594	9919	P	976	  Bacteroidetes
9.66	233675	996	C	200643	    Bacteroidia
9.62	232679	65248	O	171549	      Bacteroidales
6.92	167431	1059	F	171551	        Porphyromonadaceae
6.88	166372	36434	G	836	          Porphyromonas
5.37	129938	129938	S	837	            Porphyromonas gingivalis
```

You are welcome to run some or all of the analyses. Outputs from each analysis are provided in the [**Example_Outputs**](https://github.com/dportik/LRSW-Taxonomic-Profiling-Tutorial/tree/main/Example_Outputs), and can be used to run the comparative analysis.

[**Back to top**](#TOP)

## Example Dataset <a name="EXDA"></a>

For this section we will use [**seqtk**](https://github.com/lh3/seqtk) for file format conversion. It can be installed using the following options:

conda
```
$ conda create --name seqtk-env -c bioconda seqtk 

# activate the environment
$ conda activate seqtk-env
```

github
```
$ git clone https://github.com/lh3/seqtk.git
$ cd seqtk; make
```

The example dataset used for the tutorial is [**ATCC MSA-1003**](https://www.atcc.org/products/msa-1003), sequencing with PacBio HiFi. The ATCC MSA-1003 mock community contains 20 bacteria species with the following staggered design:

|Number of species| Relative abundance|
|----:|----:|
|5|18.00%|
|5|1.80%|
|5|0.18%|
|5|0.02%|

The dataset was generated using the Sequel II System and contains 2.4 million HiFi reads with a median length of 8.3 kb, for a total of 20.5 Gbp of data.

This dataset is available on NCBI ([**SRR9328980**](https://www.ncbi.nlm.nih.gov/sra/SRR9328980)). We can download it as a `fastq.gz` file that is ~15 Gb in size. We will need to convert this to `fasta` format for several programs. The following commands can be used to download the data, rename the file, and convert to fasta format:

```
# make directory for data
$ mkdir data
$ cd data

# download from NCBI
$ wget https://sra-pub-src-1.s3.amazonaws.com/SRR9328980/m64015_190521_115634.Q20.fastq.gz.1

# rename file
$ mv m64015_190521_115634.Q20.fastq.gz.1 ATCC.fastq.gz

# convert to fasta using seqtk
$ seqtk seq -a ATCC.fastq.gz > ATCC.fasta

$ cd ..
```

The directory should now contain the following files:

```

├── data/
│	├── ATCC.fasta
│	└── ATCC.fastq.gz
```

[**Back to top**](#TOP)


## Kraken2 & Bracken <a name="KRBR"></a>

[**Kraken2**](https://ccb.jhu.edu/software/kraken2/) is a kmer-based approach, and [**Bracken**](https://github.com/jenniferlu717/Bracken) is a companion program that provides Bayesian refinement of abundance estimates. 

### Program & Database Installation

We will first install **Kraken2** and **Bracken**, and then obtain a suitable database. 

Conda installation:

```
$ conda create --name krakenbracken-env -c bioconda kraken2 bracken 

# activate the environment
$ conda activate krakenbracken-env
```

For manual installation, follow instructions for **Kraken2** ([**here**](https://github.com/DerrickWood/kraken2)) and for **Bracken** ([**here**](https://github.com/jenniferlu717/Bracken)).

**Kraken2** and **Bracken** rely on the same style of database. Several pre-compiled databases are available [**here**](https://benlangmead.github.io/aws-indexes/k2). We will use the `PlusPF` database from 2021, which contains archaea, bacteria, viral, plasmid, human, UniVec_Core, protozoa & fungi sequences. It has all files required for both **Kraken2** and **Bracken**. The unzipped size is roughly ~58 Gb.

```
# create a directory for the analysis
$ mkdir krakenbracken_analysis
$ cd krakenbracken_analysis

# make a directory to hold the database files
$ mkdir kraken_db
$ cd kraken_db

# download and unpack the database
$ wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20210517.tar.gz
$ tar -xvzf k2_pluspf_20210517.tar.gz
$ rm k2_pluspf_20210517.tar.gz
$ cd ..
``` 

The directory should now contain the following:

```
├── data/
│	├── ATCC.fasta
│	└── ATCC.fastq.gz
│
├── krakenbracken_analysis/
│	│
│	└── kraken_db/
│	    ├── hash.k2d
│	    ├── opts.k2d
│	    ├── taxo.k2d
│	    ├── seqid2taxid.map
│	    ├── inspect.txt
│	    ├── database100mers.kmer_distrib
│	    ├── database150mers.kmer_distrib
│	    ├── database200mers.kmer_distrib
│	    ├── database250mers.kmer_distrib
│	    ├── database300mers.kmer_distrib
│	    ├── database50mers.kmer_distrib
│	    └── database75mers.kmer_distrib
```


### Analysis Instructions

From the `krakenbracken_analysis/` directory, run **Kraken2** the following general command to generate the kreport file. Change the number of threads and path to fasta file as necessary. This analysis should take <30 min.

```
$ kraken2 --db kraken_db --threads 12 --output ATCC.kraken --report ATCC.kraken.kreport.txt /data/ATCC.fasta
```

We can then run **Bracken** to generate the refined kreport. Change the number of threads as necessary. This analysis will be nearly instantaneous.

```
$ bracken -d kraken_db -t 12 -i ATCC.kraken.kreport.txt -o ATCC.bracken -l S 

# rename output
$ mv ATCC.kreport_bracken_species.txt ATCC.bracken.kreport.txt
```

The directory should now contain two kreport files with the taxonomic profiles for the sample, one from Kraken2 and one from Bracken:

```
├── data/
│	├── ATCC.fasta
│	└── ATCC.fastq.gz
│
├── krakenbracken_analysis/
│	├── ATCC.bracken
│	├── ATCC.bracken.kreport.txt   <- kreport from Bracken
│	├── ATCC.kraken
│	├── ATCC.kraken.kreport.txt    <- kreport from Kraken2
│	└── kraken_db/
```

[**Back to top**](#TOP)

## Centrifuge <a name="CENT"></a>

[**Centrifuge**](http://www.ccb.jhu.edu/software/centrifuge/index.shtml) is a kmer-based profiling approach. 

### Program & Database Installation

We will first install **Centrifuge**, and then obtain a suitable database. 

Conda installation:

```
$ conda create --name centrifuge-env -c bioconda centrifuge

# activate the environment
$ conda activate centrifuge-env
```

For manual installation, follow instructions for **Centrifuge** on github [**here**](https://github.com/DaehwanKimLab/centrifuge) or on the website [**here**](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#obtaining-centrifuge).

For **Centrifuge**, a few older pre-compiled databases are available [**here**](https://benlangmead.github.io/aws-indexes/k2). We will use a pre-compiled database from 2018 (`Refseq: bacteria, archaea (compressed)`) which contains archaea and bacteria sequences. Note that this database is quite outdated. Creating a new database is essential for real analysis - instructions for doing so are [**here**](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building). However, the process is problematic as many scripts are outdated, which is why we use a pre-compiled database for the tutorial. The unzipped database size is roughly ~8.5 Gb.

```
# create a directory for the analysis
$ mkdir centrifuge_analysis
$ cd centrifuge_analysis

# download and unpack the database
$ wget https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed_2018_4_15.tar.gz
$ tar -xvzf p_compressed_2018_4_15.tar.gz
$ rm p_compressed_2018_4_15.tar.gz
```

The directory should contain the following:

```
├── data/
│	├── ATCC.fasta
│	└── ATCC.fastq.gz
│
├── centrifuge_analysis/
│	├── p_compressed.1.cf
│	├── p_compressed.2.cf
│	├── p_compressed.3.cf
│	└── p_compressed.4.cf
```


### Analysis Instructions

From the `centrifuge_analysis/` directory, run the following general command for **Centrifuge**. Change the number of threads and path to fasta file as necessary. This analysis should take <30 min.

```
$ centrifuge --threads 12 -f -k 20 -t -x p_compressed -U /data/ATCC.fasta -S ATCC.txt --report-file ATCC.centrifuge.tsv
```

We can then generate the kreport format:
```
$ centrifuge-kreport -x p_compressed ATCC.txt > ATCC.centrifuge.kreport.txt
```

In the benchmarking study we explored "long read" settings in **Centrifuge**. Note that this did not produce particularly good results. 
To run this, simply include the parameter `--min-hitlen 500`:
```
$ centrifuge --threads 12 --min-hitlen 500 -f -k 20 -t -x p_compressed -U /data/ATCC.fasta -S ATCC.500.txt --report-file ATCC.500.centrifuge.tsv
$ centrifuge-kreport -x p_compressed ATCC.500.txt > ATCC.centrifuge500.kreport.txt
```


The directory should now contain a kreport file with the taxonomic profile for the sample, for each parameter set used:

```
├── data/
│	├── ATCC.fasta
│	└── ATCC.fastq.gz
│
├── centrifuge_analysis/
│	├── ATCC.500.centrifuge.tsv
│	├── ATCC.500.txt
│	├── ATCC.centrifuge.kreport.txt     <- kreport; default settings
│	├── ATCC.centrifuge.tsv
│	├── ATCC.centrifuge500.kreport.txt  <- kreport; "long read" settings
│	├── ATCC.txt
│	├── p_compressed.1.cf
│	├── p_compressed.2.cf
│	├── p_compressed.3.cf
│	└── p_compressed.4.cf
```

[**Back to top**](#TOP)


## MMSeqs2 <a name="MMSEQ"></a>

[**MMSeqs2**](https://github.com/soedinglab/mmseqs2) is a protein alignment approach. 

### Program & Database Installation

We will first install **MMSeqs2**, and then obtain a suitable database. 

Conda installation:

```
$ conda create --name mmseqs-env -c bioconda mmseqs2 

# activate the environment
$ conda activate mmseqs-env
```

For manual installation, follow instructions for **MMSeqs2** on github [**here**](https://github.com/soedinglab/MMseqs2/wiki).

**MMSeqs2** has a module that can be used to download a number of different databases, which is explained [**here**](https://github.com/soedinglab/MMseqs2/wiki#downloading-databases). For the sake of time and space we will use the UniRef50 protein database (requires ~25 Gb space), but in the benchmarking study the NR database was used (requires ~215 Gb space).

```
# create a directory for the analysis
$ mkdir mmseqs2_analysis
$ cd mmseqs2_analysis

# download & prepare database (adjust threads as necessary)
$ mmseqs databases UniRef50 db_uniref50 tmp --threads 24

# remove the temporary directory
$ rm tmp

# to download NCBI NR instead, you could use the following commands:
$ mkdir ncbi_nr
$ mmseqs databases NR ncbi_nr tmp --threads 24
$ rm tmp
```

The directory should contain the following:

```
├── data/
│	├── ATCC.fasta
│	└── ATCC.fastq.gz
│
├── mmseqs2_analysis/
│	├── db_uniref50
│	├── db_uniref50.dbtype
│	├── db_uniref50_h
│	├── db_uniref50_h.dbtype
│	├── db_uniref50_h.index
│	├── db_uniref50.index
│	├── db_uniref50.lookup
│	├── db_uniref50_mapping
│	├── db_uniref50.source
│	├── db_uniref50_taxonomy
│	└── db_uniref50.version
```

[**Back to top**](#TOP)

### Analysis Instructions

From the `mmseqs2_analysis/` directory, run the following general command for **MMSeqs2**. Change the number of threads and path to fasta file as necessary. Note that the `tmp` directory specifies the temporary directory to write to, and this can be changed to another location if needed. This analysis will take ~5-10 hrs depending on resources available.

```
$ mmseqs easy-taxonomy /data/ATCC.fasta db_uniref50 ATCC tmp --threads 48 
```

After completion, we will rename the relevant output file:

```
$ mv ATCC_report ATCC.mmseqs.kreport.txt
```

The directory will contain a kreport file with the taxonomic profile for the sample:

```
├── data/
│	├── ATCC.fasta
│	└── ATCC.fastq.gz
│
├── mmseqs2_analysis/
│	├── ATCC_lca.tsv
│	├── ATCC.mmseqs.kreport.txt  <- kreport from mmseqs2
│	├── ATCC_tophit_aln
│	├── ATCC_tophit_report
│	├── db_uniref50
│	├── db_uniref50.dbtype
│	├── db_uniref50_h
│	├── db_uniref50_h.dbtype
│	├── db_uniref50_h.index
│	├── db_uniref50.index
│	├── db_uniref50.lookup
│	├── db_uniref50_mapping
│	├── db_uniref50.source
│	├── db_uniref50_taxonomy
│	├── db_uniref50.version
│	└── tmp/
```

[**Back to top**](#TOP)
     
## BugSeq <a name="BS"></a>

[**BugSeq**](https://bugseq.com) is a cloud-based analysis platform that relies on nucleotide alignments for matching.

### Program & Database Installation

No installation required! **BugSeq** is a cloud-based web service.

### Analysis Instructions

Create a free account with **BugSeq**: https://bugseq.com/app/register

Start a `New Analysis`. 

Upload the `ATCC.fastq.gz` file. 

Select the following options:

+ `Platform`: PacBio
+ `Device & Chemistry`: HiFi
+ `Metagenomic Database`: NCBI nt
+ `Sample Type`: Generic
+ `Sequenced Material`: DNA

Submit the analysis and wait for the results by email (typically <24 hrs).
   
[**Back to top**](#TOP)
    
## Diamond + MEGAN-LR <a name="MEG"></a>

[**Diamond**](https://github.com/bbuchfink/diamond) is used to perform translation alignments to a protein database and [**MEGAN-LR**](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/) is used to interpret the alignments to make taxonomic assignments.

Unlike the other programs listed above, there is a PacBio workflow available for this method. It is called `Taxonomic-Functional-Profiling-Protein` and is available on github in the [**pb-metagenomics-tools**](https://github.com/PacificBiosciences/pb-metagenomics-tools) repo. 

You will need access to an HPC to successfully run the analysis, because it requires substantial resources to run. This analysis can take multiple days to finish running. For complete instructions on how to run the workflow, please see the tutorial [**here**](https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/master/docs/Tutorial-Taxonomic-Functional-Profiling-Protein.md). The NCBI nr database must be downloaded and indexed by **Diamond**, and you must download **MEGAN-LR** and the associated mapping file before beginning. Instructions for doing this are included in the tutorial. 

Following completion of the workflow, the `MEGAN-RMA-Summary` workflow can be used to generate a taxonomic report. Complete instructions for running this workflow are available [**here**](https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/master/docs/Tutorial-MEGAN-RMA-summary.md). This workflow could be run locally, as it is not as resource intensive.

**NOTE**: The PacBio profiling pipeline is under development and is expected to change into a single workflow over the next couple months. 

[**Back to top**](#TOP)


## Comparative Analysis <a name="COMP"></a>

Outputs from each analysis are provided in the [**Example_Outputs**](https://github.com/dportik/LRSW-Taxonomic-Profiling-Tutorial/tree/main/Example_Outputs).
The results from the different methods can be compared using the **Jupyter** notebook available from: https://osf.io/uzk64/

This notebook will produce:
+ read utilization metrics (how many reads are assigned, and to which taxonomic rank)
+ precision, recall, F-measures at the species and genus level
+ visualizations for dataset characteristics

To use the notebook, you will need to have **Jupyter** installed. Installation instructions can be found [**here**](https://jupyter.org/install). The notebook also requires having the following Python libraries installed: `pandas`, `seaborn`, and `matplotlib`. These could be installed in a conda environment or using `pip` (as in the jupyter installation).

Then, run the following to start a session:
```
$ jupyter notebook
```
Navigate to the downloaded notebook location and select it to begin the session. 

Within the notebook, you will only need to change the locations for the kreport files and the contents of the file list. This notebook can be used to replicate the analyses from the benchmarking study. However, keep in mind that the results are expected to be different because we used reference databases that differ from those in the paper.


[**Back to top**](#TOP)