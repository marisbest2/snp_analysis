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

Install Anaconda if not already installed.  Tested using Anaconda3-4.3.1, or try the latest at: https://www.continuum.io/downloads

On Mac OS X
        
    ~$ wget https://repo.continuum.io/archive/Anaconda3-4.3.1-MacOSX-x86_64.sh
    ~$ bash Anaconda3-4.3.1-MacOSX-x86_64.sh

If `wget` is unavailable Anaconda can be installed using the graphical installer for Python 3.

    https://www.continuum.io/downloads

On Linux

    ~$ wget https://repo.continuum.io/archive/Anaconda3-4.3.1-Linux-x86_64.sh        
    ~$ bash Anaconda3-4.3.1-Linux-x86_64.sh
    
Use Anaconda's default installation except when asked if to, "prepend to PATH", choose yes.
    
Once Anaconda is installed close and reopen terminal.

Optional: Create a virtual environment to install the pipeline and it's requirement without breaking your other software requirements, or skip these two steps to setup environment as default:

    ~$ conda create -n snp_analysis python=3.5 # skip if wanting to use your current environment
    ~$ source activate snp_analysis # skip if wanting to use your current environment

Setup Bioconda channels.  Add them in the order shown below.  Order is important.

    ~$ conda config --add channels conda-forge
    ~$ conda config --add channels defaults # If warning, it can be ignored
    ~$ conda config --add channels r
    ~$ conda config --add channels bioconda
    
As of Anaconda3-4.3.1 ete3 requires python version < 3.6

    ~$ conda install python=3.5
    
    ~$ conda install ete3 pyvcf biopython bwa samtools picard abyss gatk raxml xlrd xlsxwriter gitpython regex pandas=0.18.1 scipy=0.17.1
    ~$ conda update ete3 pyvcf biopython bwa samtools picard abyss gatk raxml xlrd xlsxwriter gitpython

When gatk is downloaded using Anacoda it still needs to be registered.  GATK has a way to do this.  Go to GATK's website, register if not already, download GATK package, unzip it, and run:

    ~$ gatk-register /path/to/Downloads/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar.  
    
After this is ran GATK just downloaded from the GATK website can be deleted.  The download was only needed to register the Anaconda GATK package.

Xvfb (short for X virtual framebuffer) is a display server implementing the X11 display server protocol and must be in your environment to generate PDF and SVG tree files.  If X11 is not in your environment, xvfbwrapper will install, but not work.  Xvfb is likely already install but beware.  Root privileges will be needed if it is not yet in your environment.  It's available on Mac OS X via XQuartz and Linux via your package manager.

Note, Xvfb is finicky.  Using Xvfb is an unfortuante necessity.  Any feedback to improve the portablity generating PDF and SVG files will be appreciated.   On Mac OS X install using pip and on Linux use conda.  On Mac after following instructions below if error “Xvfb did not start”, then restart XQuartz application before running the script again.

On Mac OS X

Install XQuartz https://www.xquartz.org

On Linux
    
    ~$ sudo apt-get install xvfb  # Debian (Ubuntu)   
    ~$ sudo yum install xorg-x11-server-Xvfb  # CentOS
    
Install on Mac and Linux

    ~$ conda install xvfbwrapper

RAxML must be in your PATH as: raxmlHPC-SSE3, raxmlHPC-PTHREADS-AVX2, or raxml.  It seems raxmlHPC-SSE3 is the most system universal but if you have the correct computer architecture running raxmlHPC-PTHREADS-AVX2 is faster.  The script will first look for raxmlHPC-PTHREADS-AVX2.  If it is not found it will look for raxmlHPC-SSE3, then raxml.  If none are found in your PATH the script will fail.

## Scripts and dependents
Clone scripts: 

    ~$ git clone https://github.com/stuber/VCFs_to_SNP_alignment.git
    
If git is unavailable, `~$ conda install git` will make it available.

Put `VCFs_to_SNP_alignment` in your $PYTHONPATH, or easier run lines below to put script in your anaconda PATH.

    ~$ ln -s ~/VCFs_to_SNP_alignment/loopwrapper.py ~/anaconda3/bin/
    ~$ ln -s ~/VCFs_to_SNP_alignment/script1.py ~/anaconda3/bin/
    ~$ ln -s ~/VCFs_to_SNP_alignment/script2.py ~/anaconda3/bin/

The script will look in your home directory for dependencies.  These will be installed automatically when script 1 or 2 are ran if they are not already available.  

    $ ~/dependencies

If dependencies are not available they will be cloned from Github when script 1 and script 2 are ran.  Check this repo periodically, using `git pull`, or simply delete the `dependencies` directory in the home directory and rerun the script.  When the script does not find this directory it will download it anew with the latest updates.

    https://github.com/stuber/dependencies.git
    

## Test script 1

Test files can be downloaded at:

    ~$ git clone https://github.com/USDA-VS/fastq_data_set-tb_complex
    ~$ git clone https://github.com/USDA-VS/fastq_data_set-brucella
    
Files have been cut to 200,000 reads, which will give around 20X coverage.  This file size is convenent for downloading and testing.  They should not be added to any currated database or used in diagnostic reporting.  The complete sequence files are available in SRA.

Test by making directory containing FASTQs the working directory.

    ~$ cd ~/VCFs_to_SNP_alignment

To add in testing make a backup of the files

    $ mkdir original test; cp *gz original; mv *gz test; cd test; ls

Calling loopwrapper.py will batch FASTQs based on available computer resources.

    ~$ loopwrapper.py

If an error occurs it have to do with running multiple samples and downloading dependencies.  If so dependencies have been downloaded but the script needs to be restarted.  Reset `test` directory and start again.

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
    
