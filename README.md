# hicmap
This program can construct raw Hi-C contact map accurately from reads and reference depend on unique kmer.
## Install
```
git clone https://github.com/ruoyu1123/hicmap.git
cd hicmap
make
```
## Usage
This program request **jellyfish** to obtained the unique *k-mer* and  shell command as follow: 

```
jellyfish count -m 21 -s 1G -t 16 -C reference.fasta 
jellyfish dump -c -U 1 mer_counts.jf > unique.kmer 
``` 
>NOTE: The parameter -c is necessary for rafilter because of input format. 

Reference format of fasta must be two lines, thr first line is start with ">", and the second line are seq line, such as follow:

 
    >S1_1
    CTAGCTCCAGTCCCACCCCGGCCTGCAGAGTGGCTGGGCTGCAGGCATACCCCA
    >S1_2
    TGCCACAGCGGAGCTTGGATGAGCAAAAGAGGAAGTGGAGCATCTGAACTCTT
    >S1_3
    GATCTCGAACCTACTCATCTTGTGTAACAAAACTTTATACCCTTTGAACAGTCACC
    >S1_4
    AATGTTTTTAAAGTGGCCATACTGCCCAAAGCAGTTTATAGATTCAATGCTATTCCT
    ...
 

When the unique.kmer and the two-line reference are prepared, You can use the following command to run the program.
```
Usage:
        HiCmap [-t <threads>] [-k <k-size>] [-o <outpath>] <kmer_dump> <reference> <read1.fq> <read2.fq>
Parameter:
        -t, --threads         Number of threads [8]
        -o                    Output path [./]
        <kmer_dump>           Dump file of jellyfish result
        <reference>           fasta file of reference
        <read1.fq>            Hi-C reads1
        <read2.fq>            Hi-C reads2
``` 
Example:

    hicmap -t 32 -o ./ reference.fasta read1.fastq read2.fastq
    
# Statement  

Copyright (C) [2023/7/10] [yangjinbao]

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses/>.
In accordance with the GNU General Public License, you should include the above copyright and license notice in the appropriate place throughout the software and its related documentation.

Note: This statement is based on an example declaration provided for the GNU General Public License (GPL) version 3 and may need modifications to fit your specific situation. Additionally, I am not a legal expert, so it is advisable to consult professional legal advice to ensure compliance with applicable laws and license requirements.




        