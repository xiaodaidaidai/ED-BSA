# ED
Using BSA mixed pool sequencing data for ED analysis and plotting

## Table of Contents
* Backgroud
* Install
* Usage
* Example
* Output
* License
## Backgroud

## Install
### Github source code:
https://github.com/xiaodaidaidai/ED.git
## Usage
usage: BSA_ED.py [-h] -vcf VCF -bulk1 BULK1 -bulk2 BULK2 [-power POWER] [-minDP MINDP] [-maxDP MAXDP] [-window WINDOW] [-step STEP] -outpre OUTPRE

Calculate ED power from VCF file

options:
  -h, --help      show this help message and exit
  -vcf VCF        Input VCF file
  -bulk1 BULK1    Sample name for bulk1
  -bulk2 BULK2    Sample name for bulk2
  -power POWER    ED power (default is 5)
  -minDP MINDP    Minimum depth of variation (per sample) (default is 4)
  -maxDP MAXDP    Maximum depth of variation (per sample) (default is 250)
  -window WINDOW  Window size (kilobases) (default is 2000)
  -step STEP      Step size (kilobases) (default is 100)
  -outpre OUTPRE  Output file prefix

usage: AtacAmp.py [-h] [--bam BAM] [--name NAME] [--isize_value ISIZE]
                  [--interval_size INTERVAL] [--mapq MAQP] [--mode {0,1,2}]
                  [--discbk DISCBK] [--type LIB] [--gtf GTF]
                  [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit  
  --bam BAM             input the bam file  
  --name NAME, -n NAME  prefix of output files  
  --isize_value ISIZE, -i ISIZE  
                        judge a pair of reads whether is discordant  
  --interval_size INTERVAL, -s INTERVAL  
                        size of interval when compute breakpoint nearby coverage  
  --mapq MAQP, -q MAQP  
  reads maqp threshold  
  --mode {0,1,2}, -m {0,1,2}  
                        choose the analysis mode,0/1/2  
  --discbk DISCBK, -d DISCBK  
                        if you choose mode 1,you need to input discbk file  
  --type LIB          
  choose library:sc/bulk  
  --gtf GTF             gtf file  
  --threads THREADS     threads  
  ## Example

  ## Output
  ### example output(.result)

  ## License
  
