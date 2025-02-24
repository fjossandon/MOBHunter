# MOBHunter

Genomic islands and Mobile Genetic Elements predictor pipeline.

# Usage

MOBHunter is a complex pipeline that needs the installation of several tools,
where each is executed using different languages like Perl, Phyton, and
JavaScript. Because of this, using the specially prepared webserver to run it is
highly recommended: https://informatica.utem.cl/mobhunter/, which also has
results-displaying capabilities. See the citation publication for details.

# Installation

If you still want to use the MOBHunter core algorithm script locally, you need
to have the following available:

- Perl 5.22
- BioPerl 1.6.924
- Alien Hunter 1.8 (executes in Javascript and Perl)
- Aragorn 1.2.36 (executes in C)
- IslandPath-DIMOB v0.3 (executes in Perl)
- PAI-DA (executes in C and Perl)
- tRNAscan-SE 1.3.2 (executes in C and Perl)
- Phage Finder v2.2 (executes in Perl)
- PhiSpy v2.3 (executes in Python)
- CONJscan HMM profiles
- TnPREd HMM profiles

Then you need to edit the script to update the hardcoded paths for each program
and database near the beginning, so the script can find them. After that,
execute it in the command line specifying the working dir containing the input
files:

```pl
perl MOBHunter.pl "workdir"
```

# Extra

The `extra` folder contains auxiliar scripts used in the publication.

- mge_genome_comparison.pl: Used to align the Clostridium difficile 630 genome
  sequence (NC_009089) to other species from the same genus. All genomes must be
  placed in the same folder. The `blastn` version used was `2.12.0+`. The
  following accessions were used:
  - NC_009089.1: Clostridium difficile 630, complete genome.
  - NC_013315.1: Clostridium difficile CD196 chromosome, complete genome.
  - NC_013316.1: Clostridium difficile R20291 chromosome, complete genome.
  - NC_013974.1: Clostridium difficile BI9 chromosome.
  - NC_017173.1: Clostridium difficile CF5, complete genome.
  - NC_017174.1: Clostridium difficile M120, complete genome.
  - NC_017175.1: Clostridium difficile M68, complete genome.
  - NC_017178.1: Clostridium difficile complete genome, strain 2007855.
  - NC_017179.1: Clostridium difficile BI1, complete genome.
- mge_detectors_neg_control.pl: Uses the combined mapping file generated by the
  above script, (`Blastn_megablast.1_combined.ident_95`) along with the
  NC_009089 Genbank file, and the Genbank files of the tools' predictions that
  will be compared, and prints a performance comparison table.

All input and output files of both scripts mentioned above are available here:
https://doi.org/10.6084/m9.figshare.28464356

# How to Cite

If you use MOBHunter in your research, please cite our work as follows:

```
Rojas-Villalobos, C.*, Ossandon, F. J.*, Castillo-Vilcahuaman, C., Sepúlveda-Rebolledo, P., Castro, D., Zapata, A., Arisan, D., Issotta, F., & Quatrini, R.†, Moya-Beltrán, A.†
(2024). MOBHunter: A data integration platform for Mobile Genetic Elements identification in microbial genomes.
* Co-first authors
† Co-corresponding authors
```
