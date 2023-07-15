# parazitCUB

# Introduction

- Viruses exhibit discernible evolutionary patterns in their adaptation and reproductive strategies within host organisms. Theoretical considerations suggest that the codon usage of viruses reflects evolutionary changes that enable them to optimize their codon usage to their host organisms. While existing software tools can analyze codon usage in organisms, their performance can be further enhanced, as they do not specifically address the examination of co-adaptation in codon usage between viruses and their hosts.

- In 2019, we developed the vhcub (https://f1000research.com/articles/8-2137/v1) R package, a pivotal tool utilized for analyzing the co-adaptation of codon usage between a virus and its host, incorporating various indices and graphical representations.

- Building upon this foundation, we have recently introduced an updated version of the vhcub package, now named parazitCUB. This new version expands the scope of investigation from virus-host interactions to encompass the interplay between parasites and their hosts in codon usage bias (CUB).

# The update from vhcub to parazitCUB package:

- rectifying any bugs present in vhcub, ensuring improved functionality.
- The optimization of function codes enables faster processing speeds. 
- Increased capacity to handle more extensive datasets compared to its predecessor.
- Analyzing multiple viruses simultaneously.
- Providing a user-friendly interface that allows users with limited R programming skills to effectively utilize the package. However, 
- For users proficient in R programming, parazitCUB offers capabilities that extend beyond its primary purpose, enabling more advanced analyses and applications.

## The previous vhcub package incorporates the following measures: 

-	The effective number of codons (with two different implementations), ENc.
-	Synonymous codon Usage orderliness, SCUO.
-	Codon adaptation index, CAI.
-	Relative codon deoptimization index, RCDI.
-	Relative synonymous codon usage, RSCU.
-	Similarity index, SiD, 
-	Also, it provides a statistical dinucleotide over- and underrepresentation with three different models.
-	GC-overall, GC1, GC2, and GC3 content.
-	ENc-GC3 plot.
-	PR2 plot.

## The updated parazitCUB package implements the following measures:

-	All the old measures and features implemented in vhcub included in parazitCUB 
-	Maximum likelihood codon bias, MCB.
-	Measure Independent of Length and Composition, MILC.
-	MILC-based Expression Level Predictor, MELP.
-	A measure of codon bias is termed B.
-	A related measure of expression, E.
-	Frequency of optimal codons, Fop.
-	ENc-GC3 plot for genes per virus.
-	ENc-GC3 plot for all viruses.
-	PR2 plot for genes per virus.
-	PR2 plot for all viruses.
-	Neutrality plot.
-	CU Heatmap using dinucleotide over- and underrepresentation.
-	CU Heatmap using RSCU.
-	GC content Boxplot.


# Installation

### Bioconductor packages, parazitCUB dependencies:

```{r eval=FALSE}

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("coRdon")
BiocManager::install("Biostrings")

```

### Using devtools:

```{r eval=FALSE}

if (!requireNamespace("devtools", quietly=TRUE)){
        install.packages("devtools")}
devtools::install_github('AliYoussef96/parazitCUB')

```


# Example of Using parazitCUB

**Note; In this example, we will only examine the package functions without any biological interpretation.**

###

In this particular illustration, we will showcase the package's capacity and demonstrate its capabilities to conduct a CUB analysis specifically on Viruses that have the potential to induce respiratory illnesses in human hosts. Notably, this encompasses notable viral groups, such as the Influenza viruses from the Orthomyxoviridae family, and the coronaviruses from the Coronaviridae family. On one hand, Within the realm of influenza viruses, there exist four prominent types: influenza A, influenza B, influenza C, and influenza D. Among these, influenza A and B viruses prevail as the primary culprits behind seasonal flu outbreaks in the human population. Influenza C viruses, although capable of instigating respiratory ailments, typically induce milder symptoms in affected individuals. Conversely, influenza D viruses predominantly impact cattle and have yet to exhibit any established capacity to infect or elicit illness in human hosts. On the other hand, Within the scope of coronaviruses, we encounter MERS (Middle East Respiratory Syndrome) and SARS-CoV-2 (Severe Acute Respiratory Syndrome Coronavirus 2). MERS-CoV, initially identified in Saudi Arabia in 2012, primarily spreads from infected dromedary camels to humans, with the possibility of human-to-human transmission. This viral strain leads to severe respiratory symptoms, accompanied by a high fatality rate. SARS-CoV-2, responsible for the ongoing COVID-19 pandemic, emerged towards the end of 2019. It primarily transmits through respiratory droplets released during coughing, sneezing, talking, or breathing by infected individuals. The clinical presentation of COVID-19 can range from mild to severe, with most individuals experiencing moderate symptoms. However, advanced age and underlying health conditions increase the risk of severe respiratory complications.

**Note; It is important to clarify that in this specific instance, our focus will solely be on exploring the functions within the package, devoid of any biological interpretation.

The coding sequences of all the viruses used in this example were obtained from the NCBI virus portal. https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/

# 1- Loading sequences

The parazitCUB package has the capability to process multiple organism FASTA files, encompassing parasites (e.g., viruses, bacteriophages, fungi, leeches, protozoa, lice, tapeworms, etc.).
Additionally, it requires a single host organism in FASTA format for analysis.

## 1.1- Loading viruses FASTA files

In instances  where the header of the FASTA files is excessively long for figures, employing the ``sep`` function can facilitate the separation of the header, thereby enhancing its visualization.

```{r}
fasta.files <- list.files("flu fasta/", pattern = ".fasta", full.names = T) #parasite.fasta

list.virus <- read.virus(fasta.files, sep = "|")

head(list.virus,1)
```

The provided function will yield a list that encompasses the sequence information of all the viruses in the form of a DNAStringSet object. 

Each DNAStringSet object will be assigned the virus name derived from the corresponding file name. It can be invoked in the following manner:

```{r}
list.virus[[1]]@virus.name

list.virus[[9]]@virus.name
```

# Quality Control boxplot: outlier examination with specific cutoff. 

During codon usage bias analysis, coding sequences that are excessively long or short in length are excluded from the analysis to ensure unbiased estimation of the utilized indexes. Although there is no predefined cutoff, this package introduces a straightforward yet effective method. It involves plotting the lengths of all coding sequences (considering 3 nucleotides as 1 unit) using a boxplot and subsequently removing any outliers using the following function.

```{r}
QC.boxplot(list.virus)

```


In this particular case, coding sequences exceeding a length of 4000 and falling below a length of 100 will be excluded. This step is implemented to ensure that the codon usage bias analysis focuses on coding sequences within a specific length range, thereby maintaining the integrity of the analysis.

```{r}

list.virus <- QC.cutoff(list.virus, cut.off.up = 4000, cut.off.down = 100)
QC.boxplot(list.virus)
```

## 1.2- Loading host FASTA files

The same for the host (only one host allowed)

```{r}

theHost <- read.host("Human.fasta", sep = "|")
head(theHost,1)

theHost[[1]]@host.name

```

# 2- Calculate CU bias

Following the successful loading of both parasite and host sequences, we are now poised to commence the Codon Usage Bias (CUB) analysis.

Codon Usage Bias (CUB) analysis can be performed at both the gene level (Section 2.1 & 2.2) and the codon level (Section 2.3).

## 2.1 CUB indices that rely solely on the coding sequence of the parasite without the need for the inclusion of the host sequence.

### GC content:

As a demonstration, we can determine the GC content for each of the viruses examined in this study:

```{r}
GC.list <- GC.content(list.virus)
head(GC.list$`Influenza A H`)
```

### the effective number of codons (ENc):

This package incorporates two distinct methods for calculating the effective number of codons (ENc) in a scientifically sophisticated manner. The first method, known as Wright's equation, utilizes the formula ENc = (2 / ∑(f_i^2)), where ENc denotes the effective number of codons and f_i represents the frequency of the i-th synonymous codon. The summation (∑) is performed across all synonymous codons. This computation is implemented within the package as the ``ENc.values.old`` function.

```{r}

enc.list <- ENc.values.old(list.virus)
head(enc.list$`Influenza A NH`)
```

Furthermore, the package includes an alternative approach introduced by researcher Mathieu Joron Novembre to calculate ENc. This approach, known as Novembre's equation, employs the formula ENc = (2 / (Σ((f_i)^2 + 1/Nc))). Here, ENc represents the effective number of codons, f_i represents the frequency of the i-th synonymous codon, Σ denotes the summation across all synonymous codons, and Nc represents the total number of codons in the coding sequence. The implementation of this calculation is available as ``ENc.values.new`` within the package. This function has a parameter ``genetic.code`` to identify a specific genetic code for different organisms. 

By incorporating both Wright's equation and Novembre's equation, the package provides researchers with a comprehensive toolkit for evaluating and quantifying the effective number of codons, enabling detailed analyses of codon usage bias in various biological contexts.

**Moreover, by utilizing the following parameters, researchers have complete control over the analysis process.**

``df.fasta``   A list returned from read.virus function or read.host function.

``genetic.code``  A single string that uniquely identifies the genetic code to extract. Should be one of the values in the id or name2 columns of GENETIC_CODE_TABLE.

``threshold`` Optional numeric, specifying sequence length, in codons, used for filtering.

``stop.rm`` Logical, whether to remove stop codons. Default is FALSE.
alt.init Logical, whether to use alternative initiation codons. Default is TRUE.

``filtering`` Character vector, one of c("none", "soft", "hard"). Specifies whether sequences shorther than some threshold value of length (in codons), len.threshold, should be excluded from calculations. If "none" (default), length of sequences is not checked, if "soft", a warrning is printed if there are shorter sequences, and if "hard", these sequences are excluded from calculation.

**The same parameters can also be applied for other function such as;**

``MILC.values`` To calculate the Measure Independent of Length and Composition.

``B.values`` To calculate A measure of codon bias, termed B.

``MCB.values`` To calculate the Maximum likelihood codon bias.


```{r}

MILC.list.virus <- MILC.values(list.virus)
head(MILC.list.virus$`Influenza A NH`)

MILC.list.host <- MILC.values(theHost)
head(MILC.list.host$Human)

```

**NOTE;**

- All of these functions (``MILC.values`` , ``B.values``, and ``MCB.values``) can be utilized to calculate CUB of parasites by employing the host coding sequences as a reference set. However, they are also capable of operating solely on the coding sequences of parasites, independent of any host coding sequences serving as a reference set. Additionally, they can be applied to host coding sequences without the need for any external reference set.

## 2.2 CUB indices depend on the host coding sequence as a reference set.

These indices rely on the presence of the host coding sequence in order to analyze the CUB of parasites. Consequently, these indices cannot be employed independently without the inclusion of the host sequence. 

These indices include:

``CAI.values`` To calculate the Codon Adaptation Index.

``MELP.values`` To calculate the MILC-based Expression Level Predictor.

``FOP.values`` To calculate the Frequency of optimal codons.

``E.values`` To calculate the Related measure of expression.


```{r}

cai.list <- CAI.values(list.virus, host = theHost)
head(cai.list$`Influenza A H`)

melp.list <- MELP.values(list.virus, host = theHost)
head(melp.list$`Influenza D NH`)
```

**NOTE;**

- These functions have the same parameters as the previous functions in section 2.1. 

## 2.3 indices performing CUB analysis on the Codon Level:

``RSCU.values`` To calculate the Relative synonymous codon usage. This function can be used for the parasites and the host.

```{r}
rscu.virus <- RSCU.values(list.virus)

rscu.host <- RSCU.values(theHost)
```

``SiD.list`` To calculate similarity index between the RSCU of the parasites and the host. This function can be used only for parasites.

```{r}
SiD <-  SiD.list(RSCU.host = rscu.host, RSCU.virus = rscu.virus ) 
head(SiD)
```
```RCDI.calc``` To calculate Relative codon deoptimization index. This function requires the usage of ``ENc.values.old``.
```{r}
enc.host <- ENc.values.old(theHost)
rcdi <- RCDI.calc(list.virus , theHost, rscu.host, enc.host)
```

``dinuc.base`` A statistical dinucleotide over- and underrepresentation using base model.

```{r eval=FALSE}
base <- dinuc.base(list.virus, permutations = 100)

```

``dinuc.codon`` A statistical dinucleotide over- and underrepresentation using codon model.

```{r eval=FALSE}
codon <- dinuc.codon(list.virus, permutations = 100)
```

``dinuc.syncodon`` A statistical dinucleotide over- and underrepresentation using synonymous codons model.

```{r eval=FALSE}
syncodon <- dinuc.syncodon(list.virus, permutations = 100)
```


# 3- Visualisation of CUB

The parazitCUB package encompasses various visualization methods CUB analysis.

It utilizes ggplot2, enabling the generated plots to be treated as ggplot2 objects, thus facilitating further customization using the functions provided by ggplot2.


``GC.boxplot`` To generate Box plots for the GC content.

```{r}
GC.boxplot(GC.list) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```


``ENc.GC3plot`` To generate ENc-GC3 plot for specific virus.


```{r}
ENc.GC3plot(enc.list[["Influenza A H"]] , GC.list[["Influenza A H"]])
ENc.GC3plot(enc.list[["Influenza A NH"]] , GC.list[["Influenza A NH"]])
```

```{r}
ENc.GC3plot(enc.list[["Influenza B H"]] , GC.list[["Influenza B H"]])
ENc.GC3plot(enc.list[["Influenza B NH"]] , GC.list[["Influenza B NH"]])

```

```{r}
ENc.GC3plot(enc.list[["Influenza C H"]] , GC.list[["Influenza C H"]])
ENc.GC3plot(enc.list[["Influenza C NH"]] , GC.list[["Influenza C NH"]])

```

```{r}
ENc.GC3plot(enc.list[["Influenza D NH"]] , GC.list[["Influenza D NH"]])

```


``ENc.GC3plot.group`` To generate ENc-GC3 plot for all virus CUB average.


```{r}
ENc.GC3plot.group(enc.list, GC.list) + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "red" , "blue", "green" , "black", "grey", "#CCFF99",
                              "#660099", "#AD1457", "#004D40" , "#BBDEFB"))

```


``PR2.plot`` To generate PR2 plot for specific virus.


```{r}
PR2.plot(list.virus[[1]], fold4 = TRUE)
PR2.plot(list.virus[[2]], fold4 = TRUE) 
PR2.plot(list.virus[[3]], fold4 = TRUE) 
PR2.plot(list.virus[[4]], fold4 = TRUE) 
PR2.plot(list.virus[[5]], fold4 = TRUE) 
PR2.plot(list.virus[[6]], fold4 = TRUE) 
PR2.plot(list.virus[[7]], fold4 = TRUE) 
```


```{r}
list.virus[[1]]@virus.name
list.virus[[2]]@virus.name
list.virus[[3]]@virus.name
list.virus[[4]]@virus.name
list.virus[[5]]@virus.name
list.virus[[6]]@virus.name
list.virus[[7]]@virus.name
```


``Neutrality.plot`` To generate Neutrality plot for specific virus.


```{r}
Neutrality.plot(GC.list[["Influenza A H"]])
Neutrality.plot(GC.list[["Influenza A NH"]])
```


```{r}
Neutrality.plot(GC.list[["Influenza B H"]])
Neutrality.plot(GC.list[["Influenza B NH"]])
```


```{r}
Neutrality.plot(GC.list[["Influenza C H"]])
Neutrality.plot(GC.list[["Influenza C NH"]])
```


```{r}
Neutrality.plot(GC.list[["Influenza D NH"]])
```


## Advanced CUB analysis and visualization


```{r}
codon <- dinuc.codon(list.virus, permutations = 10)

cub.heatmap(codon, cluster_rows = T, aver = T)
```


``rscu.pca`` PCA using the RSCU values from the host and the viruses


```{r}
rscu.pca(rscu.virus , rscu.host, codons.exclude = c("ATG", "TAA", "TAG", "TGA", "TGG"))

```


``rscu.cluster`` Cluster analysis using the RSCU values from the host and the viruses


```{r}
rscu.cluster(rscu.virus, rscu.host, k = 4, rank = 2)

```


```{r}
rscu.cluster(rscu.virus, rscu.host, k = 2, rank = 2)

```

# Contribution Guidelines

For bugs and suggestions, the most effective way is by raising an issue on the github issue tracker. 
Github allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.

If you wish to contribute some changes to the code then you should submit a [pull request](https://github.com/AliYoussef96/parazitCUB/pulls)
How to create a Pull Request? [documentation on pull requests](https://help.github.com/en/articles/about-pull-requests)


# Citation

