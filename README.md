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


  ## Example

  ## Output
  ### example output(.result)

  ## License
  
