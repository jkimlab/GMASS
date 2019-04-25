GMASS
=================
The GMASS score is a novel measure for representing structural similarity between two assemblies. It represents the structural similarity of a pair of assemblies based on the length and number of similar genomic regions defined as consensus segment blocks (CSBs) in the assemblies.

![GMASS_diagram](https://user-images.githubusercontent.com/19523543/56632912-fa2ea800-6696-11e9-939d-7b9906c16aa0.jpg)

Quick start
----------------
``` 
git clone https://github.com/jkimlab/GMASS.git 
cd GMASS
perl setup.pl install
perl calculateGMASS.pl -p example/params.txt -o outdir
```
Install package
----------------
### System requirement
   - Linux x64 (Tested in CentOS 6.9, CentOS 6.10, Ubuntu 16.04 and Ubuntu 18.04)  
   - Perl >= 5.10
   - Perl modules  
      - Switch
      - Parallel::ForkManager
   - GCC >= 4.4.7
   - zlib >= 1.2.7
   - glib 2.14
   - glib 2.17
   
To install the package for calculating GMASS score,
```
git clone https://github.com/jkimlab/GMASS.git
cd GMASS
perl setup.pl install
``` 
User can test whether the package is installed properly as
```
perl calculateGMASSS.pl -p example/params.txt -o outdir
```
User can also uninstall this package as
```
perl setup.pl uninstall
```
Run package
---------------
### Usage
```
calculateGMASS.pl [options] -f1 <fasta1> -f2 <fasta2> -r <resolutions> -s <dist>
   
-Inputs:
  -f1/-f2		Uncompressed sequence files in fasta format
  -r|--resolution		Comma-separated resolution list
  -s|--strict		Alignment strictness [self|near|medium|far] (Default: near)

-Options:
  -c|--core		Core number  (Default: 1)
  -o|--outdir		Path of output directory  (Default: Current directory)
  -h|--help		Print help message
```
You can offer the input parameters to the package using a file as
```
  calculateGMASS.pl [options] -p <params>
```
### Inputs
  - -f1/-f2: \<fasta\>  
    Uncompressed sequence files of a pair of genome assemblies compared. The files must be written in fasta format. The details of fasta files can be confirmed from https://en.wikipedia.org/wiki/FASTA_format.
    
    **Providing assembly by -f1 does not always mean the assesmbly is used as reference assembly**. The assembly used as reference or target assembly will be automatically decided by N50 size. N50 size is calculated by the package, so user don't have to do. 

  - -r: \<resolutions\>  
    Comma-seperated resolution list. Each resolution will be used for setting minimum CSB size. The package requires resolution values at least two.

  - -s: \<dist\>  
    User can adjust the strictness for alignment. There are 4 options, 'self', 'near', 'medium', 'far', and the default is 'near'.
    > Example of usage
      - self: Comparng different versions of human genome assemblies
      - near: Comparing human genome assembly to chimpanzee genome assembly
      - medium: Comparing human genome assembly to mouse genome assembly
      - far: Comparing human genome assembly to chicken genome assembly

  - -p: \<param\>  
    A file containing input parameters.
    ```
      # Assemblies compared
      ### file format should be fasta format
      $ASSEMBLY1=example/Assembly1.fa
      $ASSEMBLY2=example/Assembly2.fa
      
      # Resolution lists (Comma-separated)
      $RESOLUTION=10000,20000,30000,40000,50000

      # Alignment strictness 
      ### Option: self, near(default), medium, far
      $STRICT=near
    ```

### Outputs
If the packages run successfully, The output directory looks like this.
```
  aln.log  assembly_CSB.stats.txt  chainNet/  CSB/  data/  scores.txt
```
  - aln.log  
    Log file of alignment run
    
  - assembly_CSB.stats.txt  
    The file containing the features of assemblies and CSB in each resolution
    >Description
    - AS_count: The number of scaffolds (chromosomes) in reference/target assembly  
    - AS_length: Total size of scaffolds (chromosomes) in reference/target assembly  
    - usedAS_count: The number of scaffolds being used for CSBs in reference/target assembly  
    - usedAS_length: Total size of scaffolds being used for CSBs in reference/target assembly  
    - CSB_count: The number of CSBs constructed between the assembly pair  
    - CSB_length: Total size of CSBs constructed between the assembly pair  
    
    >Example
    ```
      #stats  10000   20000   30000   40000   50000
      AS_count(ref)   19  19  19  19  19
      AS_count(tar)   17  17  17  17  17
      AS_length(ref)  2880676 2880676 2880676 2880676 2880676
      AS_length(tar)  2862930 2862930 2862930 2862930 2862930
      CSB_count(ref)  5   5   5   4   4
      CSB_count(tar)  5   5   5   4   4
      CSB_length(ref) 2837645 2837645 2837645 2799734 2799734
      CSB_length(tar) 2828992 2828992 2828992 2791108 2791108
      usedAS_count(ref)   5   5   5   4   4
      usedAS_count(tar)   5   5   5   4   4
      usedAS_length(ref)  2855326 2855326 2855326 2817019 2817019
      usedAS_length(tar)  2828992 2828992 2828992 2791108 2791108
    ```    
  - chainNet/  
    Directory containing alignment results
  
  - CSB/  
    Directory containing constructed CSB in each resolution
    
  - data/  
    Direcotory for package input data
    
  - scores.txt  
    The file containing GMASS score as well as Ci, Si, and Li scores in each resolution
    >Example
    ```
    #Resolution Li score    Ci score    Si score
    10000   0.996889512514958   1   0.996889512514958
    20000   0.996889512514958   1   0.996889512514958
    30000   0.996889512514958   1   0.996889512514958
    40000   0.996917865804394   1   0.996917865804394
    50000   0.996917865804394   1   0.996917865804394

    --
    GMASS   0.996900853830732
    ```

Third party tools
-------------------
- Bedtools (http://bedtools.readthedocs.io/en/latest/)  
- KentUtils (http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads)  
- LASTZ (http://www.bx.psu.edu/~rsharris/lastz/)  
- inferCars (http://www.bx.psu.edu/miller_lab/car/)   

Contact
-------------------
bioinfolabkr@gmail.com
