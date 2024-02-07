date
date
date
(base) ibg-4@ibg4-LIFEBOOK-A3510:~$ ls
AS        CL13           HS02      iliad.txt            project1
CL08      CL13.log       HS02.log  iliad.txt.1          Public
CL08.log  course         HS03      multiqc_data         R
CL10      Desktop        HS03.log  multiqc_report.html  rrr
CL10.log  Documents      HS04      Music                snap
CL11      Downloads      HS04.log  mydate.sh            Software
CL11.log  example.fastq  HS05      outf                 SRR8668759.fastq.gz
CL12      HS01           HS05.log  outf1                Templates
CL12.log  HS01.log       igv       Pictures             Videos
(base) ibg-4@ibg4-LIFEBOOK-A3510:~$ cd ~/AS/
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ ls
SRR2006791.fastq.gz  SRR2006795.fastq.gz
SRR2006792.fastq.gz  SRR2006796.fastq.gz
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ fastqc -o ~/AS/ -t 8 SRR2006791.fastq.gz 
application/gzip
Started analysis of SRR2006791.fastq.gz
Approx 5% complete for SRR2006791.fastq.gz
Approx 10% complete for SRR2006791.fastq.gz
Approx 15% complete for SRR2006791.fastq.gz
Approx 20% complete for SRR2006791.fastq.gz
Approx 25% complete for SRR2006791.fastq.gz
Approx 30% complete for SRR2006791.fastq.gz
Approx 35% complete for SRR2006791.fastq.gz
Approx 40% complete for SRR2006791.fastq.gz
Approx 45% complete for SRR2006791.fastq.gz
Approx 50% complete for SRR2006791.fastq.gz
Approx 55% complete for SRR2006791.fastq.gz
Approx 60% complete for SRR2006791.fastq.gz
Approx 65% complete for SRR2006791.fastq.gz
Approx 70% complete for SRR2006791.fastq.gz
Approx 75% complete for SRR2006791.fastq.gz
Approx 80% complete for SRR2006791.fastq.gz
Approx 85% complete for SRR2006791.fastq.gz
Approx 90% complete for SRR2006791.fastq.gz
Approx 95% complete for SRR2006791.fastq.gz
Analysis complete for SRR2006791.fastq.gz
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ fastqc -o ~/AS/ -t 8 SRR2006792.fastq.gz application/gzip
Started analysis of SRR2006792.fastq.gz
Approx 5% complete for SRR2006792.fastq.gz
Approx 10% complete for SRR2006792.fastq.gz
Approx 15% complete for SRR2006792.fastq.gz
Approx 20% complete for SRR2006792.fastq.gz
Approx 25% complete for SRR2006792.fastq.gz
Approx 30% complete for SRR2006792.fastq.gz
Approx 35% complete for SRR2006792.fastq.gz
Approx 40% complete for SRR2006792.fastq.gz
Approx 45% complete for SRR2006792.fastq.gz
Approx 50% complete for SRR2006792.fastq.gz
Approx 55% complete for SRR2006792.fastq.gz
Approx 60% complete for SRR2006792.fastq.gz
Approx 65% complete for SRR2006792.fastq.gz
Approx 70% complete for SRR2006792.fastq.gz
Approx 75% complete for SRR2006792.fastq.gz
Approx 80% complete for SRR2006792.fastq.gz
Approx 85% complete for SRR2006792.fastq.gz
Approx 90% complete for SRR2006792.fastq.gz
Approx 95% complete for SRR2006792.fastq.gz
Analysis complete for SRR2006792.fastq.gz
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ fastqc -o ~/AS/ -t 8 SRR2006795.fastq.gz application/gzip
Started analysis of SRR2006795.fastq.gz
Approx 5% complete for SRR2006795.fastq.gz
Approx 10% complete for SRR2006795.fastq.gz
Approx 15% complete for SRR2006795.fastq.gz
Approx 20% complete for SRR2006795.fastq.gz
Approx 25% complete for SRR2006795.fastq.gz
Approx 30% complete for SRR2006795.fastq.gz
Approx 35% complete for SRR2006795.fastq.gz
Approx 40% complete for SRR2006795.fastq.gz
Approx 45% complete for SRR2006795.fastq.gz
Approx 50% complete for SRR2006795.fastq.gz
Approx 55% complete for SRR2006795.fastq.gz
Approx 60% complete for SRR2006795.fastq.gz
Approx 65% complete for SRR2006795.fastq.gz
Approx 70% complete for SRR2006795.fastq.gz
Approx 75% complete for SRR2006795.fastq.gz
Approx 80% complete for SRR2006795.fastq.gz
Approx 85% complete for SRR2006795.fastq.gz
Approx 90% complete for SRR2006795.fastq.gz
Approx 95% complete for SRR2006795.fastq.gz
Analysis complete for SRR2006795.fastq.gz
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ fastqc -o ~/AS/ -t 8 SRR2006796.fastq.gz application/gzip

(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ multiqc -d .

  /// MultiQC 🔍 | v1.16

|           multiqc | Prepending directory to sample names
|           multiqc | Search path : /home/ibg-4/AS
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 12/12  
|            fastqc | Found 4 reports
|           multiqc | Report      : multiqc_report.html
|           multiqc | Data        : multiqc_data
|           multiqc | MultiQC complete


(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ kallisto index -i Solanum_lycopersicum.SL3.0.cdna.all.index Solanum_lycopersicum.SL3.0.cdna.all.fa.gz

(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ kallisto quant -i Solanum_lycopersicum.SL3.0.cdna.all.index -o k_SRR2006791 --single -l 180 -s 20 SRR2006791.fastq.gz

[quant] fragment length distribution is truncated gaussian with mean = 180, sd = 20
[index] k-mer length: 31
[index] number of targets: 34,658
[index] number of k-mers: 49,290,711
[index] number of equivalence classes: 68,193
[quant] running in single-end mode
[quant] will process file 1: SRR2006791.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 11,762,204 reads, 6,489,601 reads pseudoaligned
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 903 rounds

(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ kallisto quant -i Solanum_lycopersicum.SL3.0.cdna.all.index -o k_SRR2006792 --single -l 180 -s 20 SRR2006792.fastq.gz

[quant] fragment length distribution is truncated gaussian with mean = 180, sd = 20
[index] k-mer length: 31
[index] number of targets: 34,658
[index] number of k-mers: 49,290,711
[index] number of equivalence classes: 68,193
[quant] running in single-end mode
[quant] will process file 1: SRR2006792.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 14,138,812 reads, 7,682,756 reads pseudoaligned
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 922 rounds

(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ kallisto quant -i Solanum_lycopersicum.SL3.0.cdna.all.index -o k_SRR2006795 --single -l 180 -s 20 SRR2006795.fastq.gz

[quant] fragment length distribution is truncated gaussian with mean = 180, sd = 20
[index] k-mer length: 31
[index] number of targets: 34,658
[index] number of k-mers: 49,290,711
[index] number of equivalence classes: 68,193
[quant] running in single-end mode
[quant] will process file 1: SRR2006795.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 12,687,296 reads, 6,365,177 reads pseudoaligned
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 934 rounds

(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ kallisto quant -i Solanum_lycopersicum.SL3.0.cdna.all.index -o k_SRR2006796 --single -l 180 -s 20 SRR2006796.fastq.gz

[quant] fragment length distribution is truncated gaussian with mean = 180, sd = 20
[index] k-mer length: 31
[index] number of targets: 34,658
[index] number of k-mers: 49,290,711
[index] number of equivalence classes: 68,193
[quant] running in single-end mode
[quant] will process file 1: SRR2006796.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 14,057,885 reads, 7,139,820 reads pseudoaligned
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 938 rounds

(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ 
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$  multiqc -d .

  /// MultiQC 🔍 | v1.16

|           multiqc | Prepending directory to sample names
|           multiqc | Search path : /home/ibg-4/AS
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 44/44  
|            fastqc | Found 4 reports
|           multiqc | Existing reports found, adding suffix to filenames. Use '--force' to overwrite.
|           multiqc | Report      : multiqc_report_1.html
|           multiqc | Data        : multiqc_data_1
|           multiqc | MultiQC complete


#i did that before but forget to added 
base) ibg-4@ibg4-LIFEBOOK-A3510:~/LJ$ cd ~/AS/
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ ls
'bash code date.sh'   DiffGenes.txt   multiqc_data_1         'Mumen Uthman'  'R  steps'
 CG01                'linux file'     multiqc_report_1.html   myDGEList       TG01
 CG02                 multiqc_data    multiqc_report.html     problems        TG02
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ fastqc -o ~/AS/ -t 8 SRR2006796.fastq.gz application/gzip
Skipping 'SRR2006796.fastq.gz' which didn't exist, or couldn't be read
Skipping 'application/gzip' which didn't exist, or couldn't be read
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ fastqc -o ~/AS/ -t 8 SRR2006796.fastq.gz application/gzip
Skipping 'application/gzip' which didn't exist, or couldn't be read
application/gzip
Started analysis of SRR2006796.fastq.gz
Approx 5% complete for SRR2006796.fastq.gz
Approx 10% complete for SRR2006796.fastq.gz
Approx 15% complete for SRR2006796.fastq.gz
Approx 20% complete for SRR2006796.fastq.gz
Approx 25% complete for SRR2006796.fastq.gz
Approx 30% complete for SRR2006796.fastq.gz
Approx 35% complete for SRR2006796.fastq.gz
Approx 40% complete for SRR2006796.fastq.gz
Approx 45% complete for SRR2006796.fastq.gz
Approx 50% complete for SRR2006796.fastq.gz
Approx 55% complete for SRR2006796.fastq.gz
Approx 60% complete for SRR2006796.fastq.gz
Approx 65% complete for SRR2006796.fastq.gz
Approx 70% complete for SRR2006796.fastq.gz
Approx 75% complete for SRR2006796.fastq.gz
Approx 80% complete for SRR2006796.fastq.gz
Approx 85% complete for SRR2006796.fastq.gz
Approx 90% complete for SRR2006796.fastq.gz
Approx 95% complete for SRR2006796.fastq.gz
Analysis complete for SRR2006796.fastq.gz
(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ 








(base) ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ 

i use this for indexing#  ibg-4@ibg4-LIFEBOOK-A3510:~/AS$ kallisto quant -i Solanum_lycopersicum.SL3.0.cdna.all.index -o k_SRR2006791 --single -l 180 -s 20 SRR2006791.fastq.gz
