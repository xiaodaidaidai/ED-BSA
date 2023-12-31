# ED-BSA
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
## BSA_ED.py        
```
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
```

## BSA_ED_sliding.R
```
usage: BSA_ED_sliding.R [--] [--help] [--opts OPTS] ED window chr
       threshold outpre

Mapping the distribution of ED correlation values on chromosomes

positional arguments:                          
  ED          Input: ED.tsv                  
  window      Input: sliding_window.tsv                          
  chr         Input: Chromosome list file (chr.len)                        
  threshold   Input: Threshold of confidence interval (99 or 95)                        
  outpre      Output prefix                          

flags:
  -h, --help  show this help message and exit                  

optional arguments:                        
  -x, --opts  RDS file containing argument values  
```

## BSA_ED_loess.R    
```
usage: BSA_ED_loess.R [--] [--help] [--opts OPTS] ED chr threshold                                                                        
       outpre                                                

Generate ED figure                                                                        

positional arguments:                                                            
  ED          Input: ED.tsv                                                                        
  chr         Input: Chromosome length file (chr.len)                                        
  threshold   Input: Threshold of confidence interval (99 or 95)                                                    
  outpre      Output prefix                                            

flags:
  -h, --help  show this help message and exit                    

optional arguments:
  -x, --opts  RDS file containing argument values     
```

## vcf_processor.py
```
usage: vcf_processor.py [-h] input_file output_file                      

Process VCF file and generate output                                    

positional arguments:                            
  input_file   Input VCF file                                  
  output_file  Output file                            

options:
  -h, --help   show this help message and exit                                    
```

## Example 1：Calculate ED power from VCF file
```
python BSA_ED.py \
  -vcf all.clean.snp.qtlseq.vcf \
  -bulk1 S-pool -bulk2 T-pool \
  -power 5 \
  -minDP 4 -maxDP 100 \
  -outpre all.clean.snp.BSA-ED
```

## Example 2：Mapping the distribution of ED correlation values on chromosomes
```
Rscript BSA_ED_sliding.R \
    all.clean.snp.BSA-ED.snp_EDpower5.tsv \
    all.clean.snp.BSA-ED.sliding_EDpower5.tsv \
    chr.len \
    99 \
    all.clean.snp.BSA-ED5
```

## Example 3：Generate ED figure 
```
Rscript BSA_ED_loess.R \
    all.clean.snp.BSA-ED.snp_EDpower5.tsv \
    chr.len \
    99 \
    all.clean.snp.BSA-ED5
```

## Output
### example output(.result)
## License
  
