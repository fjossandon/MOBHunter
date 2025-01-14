# MOBHunter
Genomic islands and Mobile Genetic Elements predictor pipeline.

# Usage
MOBHunter is a complex pipeline that needs the installation of several tools, where each is executed using different languages like Perl, Phyton, and Javascript.
Because of this, using the specially prepared webserver to run it is highly recommended: https://informatica.utem.cl/mobhunter/,
which also has results-displaying capabilities. See the citation publication for details.

# Installation
MOBHunter core script
If you still want to use the MOBHunter core algorithm locally, you need to have available:
MOBHunter core script
If you still want to use the MOBHunter core algorithm locally, you need to have available:
* Perl 5.22
* BioPerl 1.6.924
* Alien Hunter 1.8 (executes in Javascript and Perl)
* Aragorn 1.2.36 (executes in C)
* IslandPath-DIMOB v0.3  (executes in Perl)
* PAI-DA (executes in C and Perl)
* tRNAscan-SE 1.3.2 (executes in C and Perl)
* Phage Finder v2.2 (executes in Perl)
* PhiSpy v2.3 (executes in Python)
* CONJscan HMM profiles
* TnPREd HMM profiles

Then you need to edit the script to update the hardcoded paths for each program and database near the beginning, so the script can find them.
After that, execute it in the command line specifying the working dir containing the input files:
```pl
perl MOBHunter.pl "workdir"
```

# How to Cite
If you use MOBHunter in your research, please cite our work as follows:
```
Rojas-Villalobos, C.*, Ossandon, F. J.*, Castillo-Vilcahuaman, C., Sepúlveda-Rebolledo, P., Castro, D., Zapata, A., Arisan, D., Issotta, F., & Quatrini, R.†, Moya-Beltrán, A.†
(2024). MOBHunter: A data integration platform for Mobile Genetic Elements identification in microbial genomes.
* Co-first authors
† Co-corresponding authors
```
