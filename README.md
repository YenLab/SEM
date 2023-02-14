# Size-based Expectation Maximum (SEM)<a name="sem"></a>

- [Introduction](#intro)
- [Installation](#install)
- [Quick Start](#quick-start)
- [Usage](#usage)
  * [Detecting nucleosome subtypes](#subtype)
  * [Running SEM](#run)
  * [Restrict nucleosome finding regions](#restrict)

## Introduction<a name="intro"></a>

The SEM algorithm is designed for nucleosome subtype finding, which uses a Bayesian probablistic model to model each individual nucleosome's contribution to each observed MNase-seq read pair, expectation maximum (EM) is used to find the maximum posteriori probability (MAP) of the nucleosome parameters, including:

- Dyad location
- Occupancy
- Fuzziness
- Nucleosome subtype mixture probability

The nucleosome subtype is defined as a Normal Distribution describing the probability of observing the particular protected fragment size from a nucleosome subtype, i.e., a canonical nucleosome should protect ~147bp DNA with some variation under extensive MNase digestion.


## Installation<a name="install"></a>

SEM is available on [bioconda](https://bioconda.github.io/). To install, simply run:

```bash
$ conda install -c bioconda sem
```

Also you can download the latest runnable JAR file from [SEM releases](https://github.com/YenLab/SEM/releases).

## Quick Start<a name="quick-start"></a>

Run SEM on a single MNase-seq experiment:

```bash
$ java -jar sem_latest.jar --threads 15 --geninfo mm10.fa.fai --out test_run/ --expt mES_MNase-seq.bam --format SAM 
```

The outputs will be inside the `test_run` directory in this example.

## Usage<a name="usage"></a>

```
Required:
	--out <output file prefix>
	--geninfo <genome info file>
	--expt <file name> AND --format <SAM/BED/SCIDX>
	OR
	--design <experiment design file name to use instead of --expt and --ctrl>
General:
	--threads <number of threads to use (default=1)>
	--verbose [flag to print intermediate files and extra output]
	--config <config file: all options here can be specified in a name<space>value text file, over-ridden by command-line args>
Genome:
	--providedPotenialRegions <bed file to restrict nucleosome detection regions> 
Loading Data:
	--nonunique [flag to use non-unique reads]
	--mappability <fraction of the genome that is mappable for these experiments (default=0.8)>
	--nocache [flag to turn off caching of the entire set of experiments (i.e. run slower with less memory)]
Detecting nucleosome type:
	--numClusters <number of nucleosome types> 
	--providedBindingSubtypes <custom binding subtypes (format: mean variance weight, sum of weights = 1)> 
Running SEM:
	--r <max. model update rounds, default=3>
	--alphascale <alpha scaling factor(default=1.0>
	--fixedalpha <minimum number of fragments a nucleosome should have (default=1, must >= 1)>
	--exclude <file of regions to ignore, usually blacklist regions>
	--consensusExclusion <consensus exclusion zone>
```

### Detecting nucleosome subtypes<a name="subtype"></a>

In SEM, Gaussian Mixture Model (GMM) is used on MNase-seq fragment size distribution to find the potential nucleosome subtypes, each nucleosome subtype is represented by a Normal Distribution with parameters mean and variance. When `--numClusters` is set as a positive integer, a finite GMM will be used to find out the parameters of each nucleosome subtype.

It's recommended to use `Picard CollectInsertSizeMetrics` first to check the distribution of fragment size distribution to decide the number of clusters. When there is no prior knowledge on the number of nucleosome subtypes, `--numClusters` can also be set as `-1` to let SEM decide it by a Dirichlet Process Mixture Model (DPMM).

Users can also provide their own nucleosome subtype information instead of using SEM built-in functions by setting `--providedBindingSubtypes`, file should be in the below format with tab delimited:

```
mean_1	var_1	weight_1
mean_2	var_2	weight_2
...
```

The sum of weights should be equal to 1

### Running SEM<a name="run"></a>

`--r` controls how many rounds of EM will be used, the default 3 rounds is able to return precise enough nucleosome predictions. This number can be increased to get more precise predictions but at the risk of overfitting.

`--fixedalpha` decides the threshold for nucleosome occupancy, during EM, all nucleosomes below this threshold will be terminated, which ensures all the remaining nucleosomes have occupancy >= `fixedalpha`

`--consensusExclusion` decides the exclusion zone between nucleosomes, the spacing between nucleosomes will be >= this threshold.

### Restrict nucleosome finding regions<a name="restrict"></a>

Since nucleosomes are everywhere on the genome, it's both computational intensive and time consuming to do EM on all the nucleosomes, besides, not all nucleosomes are of interest sometimes. `--providedPotenialRegions` can accept a bed file to restrict the regions where SEM will do nucleosome finding in, an example of potential regions could be candidate cis-regulatory regions (ccREs) from [ENCODE SCREEN project](https://screen.encodeproject.org/).

A potential issue when running SEM in this way is that large number of small pieces of regions could lead to imprecise nucleosome finding at the edge of region. So it's recommended to expand each region then merge the regions close to each other, below is the comparison between original ccRE regions and merged ccRE regions by `bedtools slop -b 500` and `bedtools merge -d 500`

![example region of user provided regions](https://raw.githubusercontent.com/YenLab/SEM/master/images/igv_comparison_expanded_merged_ccREs.png)

> Note: even provided with a potential region file, SEM will still use all MNase-seq fragments on the genome to decide nucleosome subtypes, if you want a different behavior, please refer to [Detecting nucleosome subtypes](#subtype) for how to provide custom nucleosome subtypes information.







