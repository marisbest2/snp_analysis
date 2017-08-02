---
layout: "default"
title: FAQ
weight: 7
permalink: /faq/
---

<h1><p style="text-align: center">Frequently Asked Questions</p></h1>

-----
<br>

## Step 1

-----

## GATK fails

Did `gatk-register` run successfully?  See setup instructions.  After running command a message similar to "jar file specified matches expected version" should be output.

## Samtools unable to access libraries

`samtools: error while loading shared libraries: libbz2.so.1.0: cannot open shared object file: No such file or directory`

The Samtools Anaconda install is unable to access libraries.  Fix: uninstall conda installation.  Use installation available via package manager such as apt-get or yum, or via a systems module setup.

`~$ conda uninstall samtools`

restart terminal

`~$ which samtools`
`~$ samtools`

The above should now print samtools list of commands.

### Step 2

-----
