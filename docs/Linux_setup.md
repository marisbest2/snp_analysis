---
layout: "default"
title: Linux setup
weight : 2
#permalink: /linux_setup/
---

<h1><p style="text-align: center">Linux setup</p></h1>

-----
<br>

JEEVES INSTALLATION
=================

## Python environment setup

Instructions shown are to install without root privileges.  Use root privileges if available or with another user PATH setup if desired.

Script is written in Python and must be ran using Python 3.  

Anaconda is a highly trusted Python package distrubution platform.  If running Python 2 a new environment can be set without disrupting your current Python environment.  See option below for installing an additional Anaconda environment.  

Install Anaconda if not already installed.  Tested using Anaconda3-5.0.0.

    ~$ wget https://repo.continuum.io/archive/Anaconda3-5.0.0-Linux-x86_64.sh        
    ~$ bash Anaconda3-5.0.0-Linux-ppc64le.sh
    
When installing from the command line use Anaconda's default installation except when asked to, "prepend to PATH", choose yes.
    
Once Anaconda is installed close and reopen terminal.

---

<strong>Optional:</strong>
If you are currently using Python 2, and wish to keep it as your default, a virtual environment can be built.  `jeeves.py` can run in this virtual environment.  See notes at the bottom of this page.  It is recommend to not run in a virtual environment if possible.

---

<br>
Setup Bioconda channels.  Add them in the order shown below.  Order is important.

    ~$ conda config --add channels conda-forge
    ~$ conda config --add channels defaults # If warning, it can be ignored
    ~$ conda config --add channels r
    ~$ conda config --add channels bioconda
    
Install programs.
    
    ~$ conda install pyvcf biopython bwa samtools picard abyss gatk raxml newick_utils xlrd xlsxwriter gitpython regex cairosvg pandas

When gatk is downloaded using Anacoda it still needs to be registered.  GATK has a way to do this.  Go to GATK's website, download the GATK package, unzip it, and run:

    ~$ gatk-register /path/to/Downloads/GenomeAnalysisTK*/GenomeAnalysisTK.jar.  
    
After `gatk-register` is ran, GATK just downloaded from the GATK website, can be deleted.  The download was only needed to register the Anaconda GATK package.

## Script and dependents
Clone script: 

    ~$ git clone https://github.com/USDA-VS/snp_analysis.git
    
If git is unavailable, `~$ conda install git` will make it available.

Put `snp_analysis` in your $PYTHONPATH, or easier run lines below to put script in your anaconda PATH.

    ~$ ln -s ~/snp_analysis/jeeves.py ~/anaconda*/bin/

## Test step 1

Test files can be downloaded at:

    ~$ git clone https://github.com/USDA-VS/fastq_data_set-tb_complex
    ~$ git clone https://github.com/USDA-VS/fastq_data_set-brucella
    
Files have been cut to 200,000 reads, which will give around 20X coverage.  This file size is convenent for downloading and testing.  They should not be added to any currated database or used in reporting.  The complete sequence files are available in SRA.

Test by making directory containing FASTQs the working directory.

    ~$ cd ~/fastq_data_set-tb_complex

To aid in testing make a backup of the files

    $ mkdir original test; cp *gz original; mv *gz test; cd test; ls

`jeeves.py` must only see `*fastq.gz` files in the working directory.  It will exit if any other file type is found.  `jeeves` will batch FASTQs based on available computer resources.

    $ jeeves.py

If an error occurs it may have to do with running multiple samples and downloading dependencies.  Reset `test` directory and restart.

    $ rm -r ../test/*; cp ../original/* ./; ls
    $ jeeves.py

## Test step 2

In step 1 jeeves saw there were FASTQ files in the working directory and ran appropriately.  In step 2 it will look for VCF files.  `jeeves.py` must only see `*vcf` files in the working directory.  It will exit if any other file type is found.  

Test using VCF test files, or better yet use the VCFs you just produced from step 1 above.  Use VCF files in the `alignment` directory ending in *-zc.vcf.  Make a working directory containing only those VCFs and call `jeeves.py`.  
    
If using VCF test files

    $ cd ~
    
    ~$ git clone https://github.com/USDA-VS/vcf_test_files.git
    
    ~$ cd vcf_test_files/bovis

    $ jeeves.py
    
For list of options:
    
    $ jeeves.py -h
    
<br>

---

<br>

<strong>Optional:</strong>
If you are currently using Python 2, and wish to keep it as your default, a virtual environment can be built.  `jeeves.py` can run in this virtual environment.  Skip this step, however, if you are not concerned with the Python version you are running.

To setup virtual environment:

    ~$ conda create -n snp_analysis python=3.5 # Optional:  skip if wanting to use your current environment
    ~$ source activate snp_analysis # Optional:  skip if wanting to use your current environment

Proceed to setting up channels...

When optional environment is used:  To deactivate virtual environment, use:
    
    > source deactivate <env name> # skip if new environment was not activated
    
