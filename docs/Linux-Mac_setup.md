---
layout: "default"
title: Linux/Mac setup
weight : 2
#permalink: /linux-mac_setup/
---

<h1><p style="text-align: center">Linux/Mac setup</p></h1>

-----
<br>

INSTALLATION
=================

## Python environment setup

Instructions shown are to install with user, not root, privileges.  Use root privileges if available or with another user PATH setup if desired.

Script 2 is written in Python and must be ran using Python 3.  

Anaconda is a highly trusted Python package distrubution platform.  If running Python 2 a new environment can be set without disrupting your current Python environment.  See note below for installing an additional Anaconda environment.  

Install Anaconda if not already installed.  Tested using Anaconda3-4.3.1.

On Mac OS X
        
    ~$ wget https://repo.continuum.io/archive/Anaconda3-4.3.1-MacOSX-x86_64.sh
    ~$ bash Anaconda3-4.3.1-MacOSX-x86_64.sh

If `wget` is unavailable Anaconda can be installed by using the graphical installer for Python 3.  See the Anaconda download page.

    https://www.continuum.io/downloads

On Linux

    ~$ wget https://repo.continuum.io/archive/Anaconda3-4.3.1-Linux-x86_64.sh        
    ~$ bash Anaconda3-4.3.1-Linux-x86_64.sh
    
When installing from the command line use Anaconda's default installation except when asked if to, "prepend to PATH", choose yes.
    
Once Anaconda is installed close and reopen terminal.

Optional: Create a virtual environment to install the pipeline and it's requirement without breaking your other software requirements, or skip these two steps to setup environment as default:

    ~$ conda create -n snp_analysis python=3.5 # skip if wanting to use your current environment
    ~$ source activate snp_analysis # skip if wanting to use your current environment

Setup Bioconda channels.  Add them in the order shown below.  Order is important.

    ~$ conda config --add channels conda-forge
    ~$ conda config --add channels defaults # If warning, it can be ignored
    ~$ conda config --add channels r
    ~$ conda config --add channels bioconda
    
Install specific version of Python modules.

    ~$ conda install python=3.5
    
    ~$ conda install pyvcf biopython bwa samtools picard abyss gatk raxml xlrd xlsxwriter gitpython regex pandas=0.18.1
    ~$ conda update pyvcf biopython bwa samtools picard abyss gatk raxml xlrd xlsxwriter gitpython

When gatk is downloaded using Anacoda it still needs to be registered.  GATK has a way to do this.  Go to GATK's website, download GATK package, unzip it, and run:

    ~$ gatk-register /path/to/Downloads/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar.  
    
After this is ran GATK just downloaded from the GATK website can be deleted.  The download was only needed to register the Anaconda GATK package.

## Scripts and dependents
Clone scripts: 

    ~$ git clone https://github.com/USDA-VS/snp_analysis.git
    
If git is unavailable, `~$ conda install git` will make it available.

Put `snp_analysis` in your $PYTHONPATH, or easier run lines below to put script in your anaconda PATH.

    ~$ ln -s ~/snp_analysis/loopwrapper.py ~/anaconda*/bin/
    ~$ ln -s ~/snp_analysis/script1.py ~/anaconda*/bin/
    ~$ ln -s ~/snp_analysis/script2.py ~/anaconda*/bin/

The script will look in your home directory for dependencies.  These will be installed automatically when script 1 or 2 are ran if they are not already available.  

    $ ~/dependencies

If dependencies are not available they will be cloned from Github when script 1 and script 2 are ran.  Check this repo periodically, using `git pull`, or simply delete the `dependencies` directory in the home directory and rerun the script.  When the script does not find this directory it will download it anew with the latest updates.

    https://github.com/stuber/dependencies.git
    

## Test script 1

Test files can be downloaded at:

    ~$ git clone https://github.com/USDA-VS/fastq_data_set-tb_complex
    ~$ git clone https://github.com/USDA-VS/fastq_data_set-brucella
    
Files have been cut to 200,000 reads, which will give around 20X coverage.  This file size is convenent for downloading and testing.  They should not be added to any currated database or used in reporting.  The complete sequence files are available in SRA.

Test by making directory containing FASTQs the working directory.

    ~$ cd ~/fastq_data_set-tb_complex

To aid in testing make a backup of the files

    $ mkdir original test; cp *gz original; mv *gz test; cd test; ls

Calling loopwrapper.py will batch FASTQs based on available computer resources.

    $ loopwrapper.py

If an error occurs it may have to do with running multiple samples and downloading dependencies.  Reset `test` directory and restart.

    $ rm -r ../test/*; cp ../original/* ./; ls
    $ loopwrapper.py

## Test script 2

Test using VCF files bundled with dependencies, or better yet use the VCFs you just produced from script 1 above.  Use VCFs in `alignment` directory ending with *-zc.vcf.  Make working directory containing only those VCFs and call script.  In `test` directory:

    $ mkdir vcfs; find . -name "*zc.vcf" -exec cp {} vcfs \;; cd vcfs; ls
   
If running bovis VCFs, run the following:

    $ script2.py -s bovis
    
For list of options:
    
    $ script2.py -h
    
To deactivate virtual environment, use:
    
    > source deactivate anaconda400 # skip if new environment was not activated
    
