---
title: 翻译 | A cross-package Bioconductor workflow for analysing methylation array
  data
author: zhengyanhua
date: '2023-08-08'
slug: []
categories:
  - workflow
tags:
  - bioconductor workflow
Description: ''
Tags: []
Categories: []
DisableComments: no
---

### Abstract

Methylation in the human genome is known to be associated with development and disease. The Illumina Infinium methylation arrays are by far the most common way to interrogate methylation across the human genome. This paper provides a Bioconductor workflow using multiple packages for the analysis of methylation array data. Specifically, we demonstrate the steps involved in a typical differential methylation analysis pipeline including: quality control, filtering, normalization, data exploration and statistical testing for probe-wise differential methylation. We further outline other analyses such as differential methylation of regions, differential variability analysis, estimating cell type composition and gene ontology testing. Finally, we provide some examples of how to visualise methylation array data.

已知人类基因组中的甲基化与发育和疾病有关。Illumina Infinium甲基化阵列是迄今为止在人类基因组中询问甲基化的最常见方法。本文提供了一个使用多个软件包分析甲基化阵列数据的Bioconductor工作流程。具体来说，我们演示了典型的差异甲基化分析管道所涉及的步骤，包括：质量控制、过滤、归一化、数据探索和探针差异甲基化的统计测试。我们进一步概述了其他分析，例如区域的差异甲基化，差异变异性分析，估计细胞类型组成和基因本体测试。最后，我们提供了一些如何可视化甲基化阵列数据的示例。

### Introduction

DNA methylation, the addition of a methyl group to a CG dinucleotide of the DNA, is the most extensively studied epigenetic mark due to its role in both development and disease (Bird [2002](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Bird2002); Laird [2003](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Laird2003)). Although DNA methylation can be measured in several ways, the epigenetics community has enthusiastically embraced the Illumina HumanMethylation450 (450k) array (Bibikova et al. [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Bibikova2011)) as a cost-effective way to assay methylation across the human genome. More recently, Illumina has increased the genomic coverage of the platform to \>850,000 sites with the release of their MethylationEPIC (850k) array. As methylation arrays are likely to remain popular for measuring methylation for the foreseeable future, it is necessary to provide robust workflows for methylation array analysis.

Measurement of DNA methylation by Infinium technology (Infinium I) was first employed by Illumina on the HumanMethylation27 (27k) array (Bibikova et al. [2009](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Bibikova2009)), which measured methylation at approximately 27,000 CpGs, primarily in gene promoters. Like bisulfite sequencing, the Infinium assay detects methylation status at single base resolution. However, due to its relatively limited coverage the array platform was not truly considered "genome-wide" until the arrival of the 450k array. The 450k array increased the genomic coverage of the platform to over 450,000 gene-centric sites by combining the original Infinium I assay with the novel Infinium II probes. Both assay types employ 50bp probes that query a [C/T] polymorphism created by bisulfite conversion of unmethylated cytosines in the genome, however, the Infinium I and II assays differ in the number of beads required to detect methylation at a single locus. Infinium I uses two bead types per CpG, one for each of the methylated and unmethylated states (Figure 1a). In contrast, the Infinium II design uses one bead type and the methylated state is determined at the single base extension step after hybridization (Figure 1b). The 850k array also uses a combination of the Infinium I and II assays but achieves additional coverage by increasing the size of each array; a 450k slide contains 12 arrays whilst the 850k has only 8.

DNA甲基化，将甲基添加到DNA的CG二核苷酸中，由于其在发育和疾病中的作用，是研究最广泛的表观遗传标记（Bird 2002;莱尔德2003）。虽然DNA甲基化可以通过多种方式测量，但表观遗传学界已经热情地接受了Illumina HumanMethylation450（450k）阵列（Bibikova等人，2011）作为一种经济有效的方法来分析整个人类基因组的甲基化。最近，Illumina通过发布其MethylationEPIC（850k）阵列，将该平台的基因组覆盖范围增加到\>850，000个位点。由于甲基化阵列在可预见的未来可能仍然很受欢迎，因此有必要为甲基化阵列分析提供强大的工作流程。

通过Infinium技术（Infinium I）测量DNA甲基化首先由Illumina用于HumanMethylation27（27k）阵列（Bibikova等人，2009），该阵列测量了大约27，000个CpG的甲基化，主要是在基因启动子中。与亚硫酸氢盐测序一样，Infinium测定法以单碱基分辨率检测甲基化状态。然而，由于其相对有限的覆盖范围，阵列平台在450k阵列到来之前并没有真正被认为是"全基因组的"。450k阵列通过将原始的Infinium I检测与新型Infinium II探针相结合，将平台的基因组覆盖率提高到超过450，000个以基因为中心的位点。两种检测类型均采用50bp探针，可查询由基因组中未甲基化胞嘧啶的亚硫酸氢盐转化产生的[C / T]多态性，但是，Infinium I和II检测在单个位点检测甲基化所需的磁珠数量不同。Infinium I每个CpG使用两种磁珠类型，每种用于甲基化和非甲基化状态（图1a）。相比之下，Infinium II设计使用一种磁珠类型，并且在杂交后的单碱基延伸步骤中测定甲基化状态（图1b）。850k阵列还使用Infinium I和II检测的组合，但通过增加每个阵列的尺寸来实现额外的覆盖范围;450K 幻灯片包含 12 个阵列，而 850K 只有 8 个阵列。

Regardless of the Illumina array version, for each CpG, there are two measurements: a methylated intensity (denoted by M) and an unmethylated intensity (denoted by U). These intensity values can be used to determine the proportion of methylation at each CpG locus. Methylation levels are commonly reported as either beta values (β=M/(M+U)) or M-values (Mvalue=log2(M/U)). For practical purposes, a small offset, α, can be added to the denominator of the β value equation to avoid dividing by small values, which is the default behaviour of the `getBeta` function in *minfi*. The default value for α

is 100. It may also be desirable to add a small offset to the numerator and denominator when calculating M-values to avoid dividing by zero in rare cases, however the default `getM` function in *minfi* does not do this. Beta values and M-values are related through a logit transformation. Beta values are generally preferable for describing the level of methylation at a locus or for graphical presentation because percentage methylation is easily interpretable. However, due to their distributional properties, M-values are more appropriate for statistical testing (Du et al. [2010](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Du2010)).

In this workflow, we will provide examples of the steps involved in analysing methylation array data using R (R Core Team [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-RCoreTeam2014)) and Bioconductor (Huber et al. [2015](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Huber2015)), including: quality control, filtering, normalization, data exploration and probe-wise differential methylation analysis. We will also cover other approaches such as differential methylation analysis of regions, differential variability analysis, gene ontology analysis and estimating cell type composition. Finally, we will provide some examples of useful ways to visualise methylation array data.

无论Illumina阵列版本如何，对于每个CpG，都有两种测量值：甲基化强度（用M表示）和未甲基化强度（用U表示）。这些强度值可用于确定每个CpG位点的甲基化比例。甲基化水平通常报告为β值（β=M/（M+U））或M值（M值=log2（M/U））。出于实际目的，可以在β值方程的分母中添加一个小偏移量α以避免除以小值，这是 minfi 中 getBeta 函数的默认行为。α的默认值为 100。在计算 M 值时，可能还需要向分子和分母添加一个小的偏移量，以避免在极少数情况下除以零，但是 minfi 中的默认 getM 函数不会这样做。Beta 值和 M 值通过 logit 变换相关联。β值通常更适合描述位点的甲基化水平或图形表示，因为甲基化百分比易于解释。然而，由于其分布特性，M值更适合统计测试（Du等人，2010）。

在此工作流程中，我们将提供使用R（R Core Team 2014）和Bioconductor（Huber等人，2015）分析甲基化阵列数据所涉及的步骤示例，包括：质量控制，过滤，归一化，数据探索和探针差异甲基化分析。我们还将介绍其他方法，例如区域的差异甲基化分析，差异变异性分析，基因本体分析和估计细胞类型组成。最后，我们将提供一些可视化甲基化阵列数据的有用方法示例。

### Differential methylation analysis

#### Obtaining the data

The data required for this workflow has been bundled with the R package that contains this workflow document. Alternatively, it can be obtained from [figshare](https://figshare.com/s/7a37f43c0ca2fec4669e). If you choose to download it seperately, once the data has been downloaded, it needs to be extracted from the archive. This will create a folder called `data`, which contains all the files necessary to execute the workflow.

Once the data has been downloaded and extracted, there should be a folder called `data` that contains all the files necessary to execute the workflow.

此工作流所需的数据已与包含此工作流文档的 R 包捆绑在一起。或者，可以从共享获得。如果您选择单独下载，则下载数据后，需要从存档中提取数据。这将创建一个名为 data 的文件夹，其中包含执行工作流所需的所有文件。下载并提取数据后，应该有一个名为 data 的文件夹，其中包含执行工作流所需的所有文件。

``` r
# set up a path to the data directory
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
# list the files
list.files(dataDirectory, recursive = TRUE)
```

To demonstrate the various aspects of analysing methylation data, we will be using a small, publicly available 450k methylation dataset ([GSE49667](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49667))(Zhang et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Zhang2013)). The dataset contains 10 samples in total: there are 4 different sorted T-cell types (naive, rTreg, act_naive, act_rTreg, collected from 3 different individuals (M28, M29, M30). For details describing sample collection and preparation, see Zhang et al. ([2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Zhang2013)). An additional `birth` sample (individual VICS-72098-18-B) is included from another study ([GSE51180](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51180))(Cruickshank et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Cruickshank2013)) to illustrate approaches for identifying and excluding poor quality samples.

There are several R Bioconductor packages available that have been developed for analysing methylation array data, including *minfi* (Aryee et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Aryee2014)), *missMethyl* (Phipson, Maksimovic, and Oshlack [2016](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Phipson2016)), *wateRmelon* (Pidsley et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Pidsley2013)), *methylumi* (Davis et al. [2015](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Davis2015)), *ChAMP* (Morris et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Morris2014)) and *charm* (Aryee et al. [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Aryee2011)). Some of the packages, such as *minfi* and *methylumi* include a framework for reading in the raw data from IDAT files and various specialised objects for storing and manipulating the data throughout the course of an analysis. Other packages provide specialised analysis methods for normalisation and statistical testing that rely on either *minfi* or *methylumi* objects. It is possible to convert between *minfi* and *methylumi* data types, however, this is not always trivial. Thus, it is advisable to consider the methods that you are interested in using and the data types that are most appropriate before you begin your analysis. Another popular method for analysing methylation array data is *limma* (Ritchie et al. [2015](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Ritchie2015)), which was originally developed for gene expression microarray analysis. As *limma* operates on a matrix of values, it is easily applied to any data that can be converted to a `matrix` in R. For a complete list of Bioconductor packages for analysing DNA methylation data, one can search for "DNAMethylation" in BiocViews (<https://www.bioconductor.org/packages/release/BiocViews.html#___DNAMethylation>) on the [Bioconductor website](https://www.bioconductor.org/).

We will begin with an example of a **probe-wise** differential methylation analysis using *minfi* and *limma*. By **probe-wise** analysis we mean each individual CpG probe will be tested for differential methylation for the comparisons of interest and p-values and moderated t-statistics (Smyth [2004](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-smyth2004ebayes)) will be generated for each CpG probe.

为了演示分析甲基化数据的各个方面，我们将使用一个小型的，公开可用的450k甲基化数据集（GSE49667）（Zhang等人，2013）。该数据集总共包含 10 个样本：有 4 种不同的分类 T 细胞类型（幼稚、rTreg、act_naive、act_rTreg，从 3 个不同的个体（M28、M29、M30）收集。有关描述样品收集和制备的详细信息，请参阅Zhang等人（2013）。另一项研究（GSE51180）（Cruickshank等人，2013年）包括了额外的出生样本（个体VICS-72098-18-B），以说明识别和排除劣质样本的方法。

已经开发了几种用于分析甲基化阵列数据的R Bioconductor软件包，包括minfi（Aryee等人，2014年），missMethyl（Phipson，Maksimovic和Oshlack 2016年），wateRmelon（Pidsley等人，2013年），methylumi（Davis等人，2015年），ChAMP（Morris等人，2014年）和charm（Aryee等人，2011年）。一些软件包，如minfi和methylumi，包括一个框架，用于从IDAT文件读取原始数据，以及用于在整个分析过程中存储和操作数据的各种专用对象。其他软件包为依赖于 minfi 或 methylumi 对象的归一化和统计测试提供了专门的分析方法。可以在 minfi 和 methylumi 数据类型之间进行转换，但是，这并不总是微不足道的。因此，建议在开始分析之前考虑您有兴趣使用的方法和最合适的数据类型。分析甲基化阵列数据的另一种流行方法是limma（Ritchie等人，2015），它最初是为基因表达微阵列分析而开发的。由于 limma 在值矩阵上运行，因此它很容易应用于任何可以转换为 R 矩阵的数据。有关用于分析DNA甲基化数据的Bioconductor软件包的完整列表，可以在Biooperd网站上的BiocViews（[https://www.bioconductor.org/packages/release/BiocViews.html#\_\_\_DNAMethylation）中搜索"DNAMethylation"。](https://www.bioconductor.org/packages/release/BiocViews.html#___DNAMethylation）中搜索“DNAMethylation”。)

我们将从一个使用minfi和limma进行探针差异甲基化分析的示例开始。通过探针分析，我们的意思是将测试每个单独的CpG探针的差异甲基化，以比较感兴趣的和p值，并且将为每个CpG探针生成调节的t统计量（Smyth 2004）。

#### Loading the data

It is useful to begin an analysis in R by loading all the packages that are likely to be required.

通过加载可能需要的所有包在 R 中开始分析非常有用。

``` r
# load packages required for analysis
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
```

The minfi, IlluminaHumanMethylation450kanno.ilmn12.hg19, IlluminaHumanMethylation450kmanifest, missMethyl, minfiData and DMRcate are methylation specific packages, while RColorBrewer and Gviz are visualisation packages. We use limma for testing differential methylation, and matrixStats and stringr have functions used in the workflow. The IlluminaHumanMethylation450kmanifest package provides the Illumina manifest as an R object which can easily be loaded into the environment. The manifest contains all of the annotation information for each of the CpG probes on the 450k array. This is useful for determining where any differentially methylated probes are located in a genomic context.

minfi，IlluminaHumanMethylation450kanno.ilmn12.hg19，IlluminaHumanMethylation450kmanifest，missMethyl，minfiData和DMRcate是甲基化特定的软件包，而RColorBrewer和Gviz是可视化软件包。我们使用limma来测试差异甲基化，matrixStats和stringr具有工作流程中使用的函数。IlluminaHumanMethylation450kmanifest 包将 Illumina manifest 作为 R 对象提供，可以轻松加载到环境中。清单包含 450k 阵列上每个 CpG 探测器的所有注释信息。这对于确定任何差异甲基化探针在基因组环境中的位置很有用。

``` r
# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
```

As for their many other BeadArray platforms, Illumina methylation data is usually obtained in the form of Intensity Data (IDAT) Files. This is a proprietary format that is output by the scanner and stores summary intensities for each probe on the array. However, there are Bioconductor packages available that facilitate the import of data from IDAT files into R (Smith et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Smith2013)). Typically, each IDAT file is approximately 8MB in size. The simplest way to import the raw methylation data into R is using the *minfi* function `read.metharray.sheet`, along with the path to the IDAT files and a sample sheet. The sample sheet is a CSV (comma-separated) file containing one line per sample, with a number of columns describing each sample. The format expected by the `read.metharray.sheet` function is based on the sample sheet file that usually accompanies Illumina methylation array data. It is also very similar to the targets file described by the *limma* package. Importing the sample sheet into R creates a `data.frame` with one row for each sample and several columns. The `read.metharray.sheet` function uses the specified path and other information from the sample sheet to create a column called `Basename` which specifies the location of each individual IDAT file in the experiment.

至于他们的许多其他BeadArray平台，Illumina甲基化数据通常以强度数据（IDAT）文件的形式获得。这是一种专有格式，由扫描仪输出，并存储阵列上每个探头的汇总强度。但是，有一些Bioconductor软件包可以帮助将数据从IDAT文件导入R（Smith等人，2013）。通常，每个 IDAT 文件的大小约为 8MB。将原始甲基化数据导入R的最简单方法是使用minfi函数read.metharray.sheet，以及IDAT文件和示例表的路径。示例表是一个 CSV（逗号分隔）文件，每个示例包含一行，其中有许多列描述每个示例。read.metharray.sheet 函数所需的格式基于通常伴随 Illumina 甲基化数组数据的示例工作表文件。它也与 limma 包描述的目标文件非常相似。将示例工作表导入 R 会创建一个 data.frame，其中每个示例有一行和几列。read.metharray.sheet 函数使用示例表中的指定路径和其他信息来创建一个名为 Basename 的列，该列指定实验中每个单独 IDAT 文件的位置。

``` r
# read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
targets
```

Now that we have imported the information about the samples and where the data is located, we can read the raw intensity signals into R from the IDAT files using the `read.metharray.exp` function. This creates an `RGChannelSet` object that contains all the raw intensity data, from both the red and green colour channels, for each of the samples. At this stage, it can be useful to rename the samples with more descriptive names.

现在我们已经导入了有关样本和数据所在位置的信息，我们可以使用 read.metharray.exp 函数将原始强度信号从 IDAT 文件读取到 R 中。这将创建一个 RGChannelSet 对象，其中包含每个样本的所有原始强度数据，包括来自红色和绿色通道的数据。在此阶段，使用更具描述性的名称重命名示例可能很有用。

``` r
# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)
rgSet
# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet
```

#### Quality control

Once the data has been imported into R, we can evaluate its quality. Firstly, we need to calculate detection p-values. We can generate a detection p-value for every CpG in every sample, which is indicative of the quality of the signal. The method used by *minfi* to calculate detection p-values compares the total signal (M+U)

for each probe to the background signal level, which is estimated from the negative control probes. Very small p-values are indicative of a reliable signal whilst large p-values, for example \>0.01, generally indicate a poor quality signal.

Plotting the mean detection p-value for each sample allows us to gauge the general quality of the samples in terms of the overall signal reliability (Figure 2). Samples that have many failed probes will have relatively large mean detection p-values.

将数据导入 R 后，我们可以评估其质量。首先，我们需要计算检测 p 值。我们可以为每个样本中的每个 CpG 生成一个检测 p 值，该值指示信号的质量。minfi 用于计算检测 p 值的方法比较总信号 （M+U）

对于每个探头的背景信号电平，这是从阴性对照探头估计的。非常小的 p 值表示信号可靠，而较大的 p 值（例如 \>0.01）通常表示信号质量较差。

绘制每个样本的平均检测p值，使我们能够根据整体信号可靠性来衡量样本的总体质量（图2）。具有许多失败探测器的样本将具有相对较大的均值检测 p 值。

``` r
# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")
```

The *minfi* `qcReport` function generates many other useful quality control plots. The *minfi* [vignette](http://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.pdf) describes the various plots and how they should be interpreted in detail. Generally, samples that look poor based on mean detection p-value will also look poor using other metrics and it is usually advisable to exclude them from further analysis.

minfi qcReport 函数可生成许多其他有用的质量控制图。小插曲描述了各种情节以及如何详细解释它们。通常，基于平均检测 p 值看起来较差的样本在使用其他指标时看起来也很差，通常建议将其从进一步分析中排除。

``` r
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")
```

Poor quality samples can be easily excluded from the analysis using a detection p-value cutoff, for example \>0.05. For this particular dataset, the `birth` sample shows a very high mean detection p-value, and hence it is excluded from subsequent analysis (Figure 2).

使用检测 p 值截止值（例如 \>0.05）可以轻松地从分析中排除质量差的样本。对于这个特定的数据集，出生样本显示出非常高的平均检测p值，因此将其排除在后续分析之外（图2）。

``` r
# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet
# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]
# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)
```

#### Normalisation

To minimise the unwanted variation within and between samples, various data normalisations can be applied. Many different types of normalisation have been developed for methylation arrays and it is beyond the scope of this workflow to compare and contrast all of them (Fortin et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Fortin2014); Wu et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Wu2014); Sun et al. [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Sun2011); Wang et al. [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Wang2012); Maksimovic, Gordon, and Oshlack [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Maksimovic2012); Mancuso et al. [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Mancuso2011); Touleimat and Tost [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Touleimat2012); Teschendorff et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Teschendorff2013); Pidsley et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Pidsley2013); Triche et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Triche2013)). Several methods have been built into *minfi* and can be directly applied within its framework (Fortin et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Fortin2014); Triche et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Triche2013); Maksimovic, Gordon, and Oshlack [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Maksimovic2012); Touleimat and Tost [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Touleimat2012)), whilst others are *methylumi*-specific or require custom data types (Wu et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Wu2014); Sun et al. [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Sun2011); Wang et al. [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Wang2012); Mancuso et al. [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Mancuso2011); Teschendorff et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Teschendorff2013); Pidsley et al. [2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Pidsley2013)). Although there is no single normalisation method that is universally considered best, a recent study by Fortin et al. ([2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Fortin2014)) has suggested that a good rule of thumb within the *minfi* framework is that the `preprocessFunnorm` (Fortin et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Fortin2014)) function is most appropriate for datasets with global methylation differences such as cancer/normal or vastly different tissue types, whilst the `preprocessQuantile` function (Touleimat and Tost [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Touleimat2012)) is more suited for datasets where you do not expect global differences between your samples, for example a single tissue. Further discussion on appropriate choice of normalisation can be found in (Hicks and Irizarry [2015](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-hicks2015quantro)), and the accompanying *quantro* package includes data-driven tests for the assumptions of quantile normalisation. As we are comparing different blood cell types, which are globally relatively similar, we will apply the `preprocessQuantile` method to our data (Figure 3). This function implements a stratified quantile normalisation procedure which is applied to the methylated and unmethylated signal intensities separately, and takes into account the different probe types. Note that after normalisation, the data is housed in a `GenomicRatioSet` object. This is a much more compact representation of the data as the colour channel information has been discarded and the M and U intensity information has been converted to M-values and beta values, together with associated genomic coordinates. Note, running the `preprocessQuantile` function on this dataset produces the warning: **'An inconsistency was encountered while determining sex'**; this can be ignored as it is due to all the samples being from male donors.

为了最大限度地减少样本内和样本之间不必要的变化，可以应用各种数据归一化。已经为甲基化阵列开发了许多不同类型的归一化，比较和对比所有这些方法超出了此工作流程的范围（Fortin等人，2014;吴等， 2014;孙等. 2011;王等. 2012;马克西莫维奇、戈登和奥什拉克 2012;曼库索等人，2011;图莱马特和托斯特2012;特申多夫等人，2013年;皮兹利等人，2013年;特里奇等人，2013 年）。minfi中内置了几种方法，可以直接在其框架内应用（Fortin等人，2014年;特里切等人，2013;马克西莫维奇、戈登和奥什拉克 2012;Touleimat和Tost 2012），而其他的则是甲基umi特异性的或需要自定义数据类型（Wu等人，2014;孙等. 2011;王等. 2012;曼库索等人，2011;特申多夫等人，2013年;皮兹利等人，2013 年）。虽然没有一种单一的归一化方法被普遍认为是最好的，但Fortin等人（2014）最近的一项研究表明，minfi框架内的一个好的经验法则是preprocessFunnorm（Fortin等人，2014）功能最适合具有全局甲基化差异的数据集，例如癌症/正常或截然不同的组织类型， 而预处理分位数函数（Touleimat和Tost 2012）更适合您不希望样本之间存在全局差异的数据集，例如单个组织。关于归一化适当选择的进一步讨论可以在（Hicks and Irizarry 2015）中找到，随附的quantro软件包包括分位数归一化假设的数据驱动测试。当我们比较全球相对相似的不同血细胞类型时，我们将对数据应用预处理分位数方法（图3）。此功能实现了分层分位数归一化程序，该程序分别应用于甲基化和非甲基化信号强度，并考虑了不同的探针类型。请注意，归一化后，数据存储在 GenomicRatioSet 对象中。这是更紧凑的数据表示，因为颜色通道信息已被丢弃，M和U强度信息已转换为M值和β值以及相关的基因组坐标。请注意，在此数据集上运行预处理分位数函数会产生警告："确定性别时遇到不一致";这可以忽略，因为所有样本都来自男性捐赠者。

``` r
# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) 
# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)
# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
```

#### Data exploration

Multi-dimensional scaling (MDS) plots are excellent for visualising data, and are usually some of the first plots that should be made when exploring the data. MDS plots are based on principal components analysis and are an unsupervised method for looking at the similarities and differences between the various samples. Samples that are more similar to each other should cluster together, and samples that are very different should be further apart on the plot. Dimension one (or principal component one) captures the greatest source of variation in the data, dimension two captures the second greatest source of variation in the data and so on. Colouring the data points or labels by known factors of interest can often highlight exactly what the greatest sources of variation are in the data. It is also possible to use MDS plots to decipher sample mix-ups.

多维缩放 （MDS） 图非常适合可视化数据，通常是浏览数据时应绘制的一些第一个图。MDS图基于主成分分析，是一种无监督方法，用于查看各种样本之间的异同。彼此更相似的样本应聚类在一起，而差异很大的样本应在图上相距更远。维度 1（或主成分 1）捕获数据中最大的变异源，维度 2 捕获数据中的第二大变异源，依此类推。按已知的感兴趣因素对数据点或标签进行着色通常可以准确突出显示数据中最大的变异来源。也可以使用 MDS 图来破译样本混淆。

``` r
# MDS plots to look at largest sources of variation
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
```

Examining the MDS plots for this dataset demonstrates that the largest source of variation is the difference between individuals (Figure 4). The higher dimensions reveal that the differences between cell types are largely captured by the third and fourth principal components (Figure 5). This type of information is useful in that it can inform downstream analysis. If obvious sources of unwanted variation are revealed by the MDS plots, we can include them in our statistical model to account for them. In the case of this particular dataset, we will include individual to individual variation in our statistical model.

检查该数据集的MDS图表明，最大的变异来源是个体之间的差异（图4）。较高的维度表明，细胞类型之间的差异主要由第三和第四主成分捕获（图5）。这种类型的信息很有用，因为它可以为下游分析提供信息。如果MDS图揭示了明显的不需要的变异来源，我们可以将它们包含在我们的统计模型中以解释它们。对于这个特定的数据集，我们将在我们的统计模型中包括个体到个体的变异。

``` r
# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
```

#### Filtering

Poor performing probes are generally filtered out prior to differential methylation analysis. As the signal from these probes is unreliable, by removing them we perform fewer statistical tests and thus incur a reduced multiple testing penalty. We filter out probes that have failed in one or more samples based on detection p-value.

性能不佳的探针通常在差异甲基化分析之前滤除。由于来自这些探针的信号不可靠，通过移除它们，我们执行的统计测试更少，从而减少了多次测试的损失。我们根据检测 p 值过滤掉一个或多个样本中失败的探针。

``` r
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
```

Depending on the nature of your samples and your biological question you may also choose to filter out the probes from the X and Y chromosomes or probes that are known to have common SNPs at the CpG site. As the samples in this dataset were all derived from male donors, we will not be removing the sex chromosome probes as part of this analysis, however example code is provided below. A different dataset, which contains both male and female samples, is used to demonstrate a [Differential Variability](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#Differential%20variability) analysis and provides an example of when sex chromosome removal is necessary (Figure 13).

根据样品的性质和生物学问题，您还可以选择从X和Y染色体或已知在CpG位点具有共同SNP的探针中过滤掉探针。由于此数据集中的样本均来自男性供体，因此作为此分析的一部分，我们不会删除性染色体探针，但是下面提供了示例代码。包含男性和女性样本的不同数据集用于演示差异变异性分析，并提供何时需要去除性染色体的示例（图 13）。

``` r
# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                        c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
```

There is a function in minfi that provides a simple interface for the removal of probes where common SNPs may affect the CpG. You can either remove all probes affected by SNPs (default), or only those with minor allele frequencies greater than a specified value.

minfi 中有一个功能，它提供了一个简单的接口，用于去除常见 SNP 可能影响 CpG 的探针。您可以移除所有受SNP影响的探针（默认），也可以仅删除那些具有次要等位基因频率大于指定值的探针。

``` r
# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt
```

We will also filter out probes that have shown to be cross-reactive, that is, probes that have been demonstrated to map to multiple places in the genome. This list was originally published by Chen et al. ([2013](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Chen2013)) and can be obtained from the authors' [website](http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48639-non-specific-probes-Illumina450k.xlsx).

我们还将过滤掉已被证明具有交叉反应性的探针，即已被证明可以映射到基因组中多个位置的探针。该列表最初由Chen等人（2013）发布，可以从作者的网站获得。

``` r
# exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)
```

Once the data has been filtered and normalised, it is often useful to re-examine the MDS plots to see if the relationship between the samples has changed. It is apparent from the new MDS plots that much of the inter-individual variation has been removed as this is no longer the first principal component (Figure 6), likely due to the removal of the SNP-affected CpG probes. However, the samples do still cluster by individual in the second dimension (Figure 6 and Figure 7) and thus a factor for individual should still be included in the model.

过滤和归一化数据后，重新检查MDS图以查看样本之间的关系是否发生变化通常很有用。从新的MDS图中可以明显看出，大部分个体间变异已被消除，因为这不再是第一个主成分（图6），可能是由于去除了受SNP影响的CpG探针。但是，样本仍然在二维中按个体聚类（图 6 和图 7），因此模型中仍应包括个体因子。

``` r
par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)])
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

par(mfrow=c(1,3))
# Examine higher dimensions to look at other sources of variation
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(1,3))
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(2,3))
legend("topright", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(3,4))
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")
```

The next step is to calculate M-values and beta values (Figure 8). As previously mentioned, M-values have nicer statistical properties and are thus better for use in statistical analysis of methylation data whilst beta values are easy to interpret and are thus better for displaying data. A detailed comparison of M-values and beta values was published by Du et al. ([2010](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Du2010)).

下一步是计算 M 值和 beta 值（图 8）。如前所述，M值具有更好的统计特性，因此更适合用于甲基化数据的统计分析，而β值易于解释，因此更适合显示数据。M值和β值的详细比较由Du等人（2010）发布。

``` r
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
```

#### Probe-wise differential methylation analysis

The biological question of interest for this particular dataset is to discover differentially methylated probes between the different cell types. However, as was apparent in the MDS plots, there is another factor that we need to take into account when we perform the statistical analysis. In the `targets` file, there is a column called `Sample_Source`, which refers to the individuals that the samples were collected from. In this dataset, each of the individuals contributes more than one cell type. For example, individual M28 contributes `naive`, `rTreg` and `act_naive` samples. Hence, when we specify our design matrix, we need to include two factors: individual and cell type. This style of analysis is called a paired analysis; differences between cell types are calculated *within* each individual, and then these differences are averaged *across* individuals to determine whether there is an overall significant difference in the mean methylation level for each CpG site. The *limma* [User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf) extensively covers the different types of designs that are commonly used for microarray experiments and how to analyse them in R.

We are interested in pairwise comparisons between the four cell types, taking into account individual to individual variation. We perform this analysis on the matrix of M-values in *limma*, obtaining moderated t-statistics and associated p-values for each CpG site. A convenient way to set up the model when the user has many comparisons of interest that they would like to test is to use a contrasts matrix in conjunction with the design matrix. A contrasts matrix will take linear combinations of the columns of the design matrix corresponding to the comparisons of interest.

Since we are performing hundreds of thousands of hypothesis tests, we need to adjust the p-values for multiple testing. A common procedure for assessing how statistically significant a change in mean levels is between two groups when a very large number of tests is being performed is to assign a cut-off on the false discovery rate (Benjamini and Hochberg [1995](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-benjamini1995fdr)), rather than on the unadjusted p-value. Typically 5% FDR is used, and this is interpreted as the researcher willing to accept that from the list of significant differentially methylated CpG sites, 5% will be false discoveries. If the p-values are not adjusted for multiple testing, the number of false discoveries will be unacceptably high. For this dataset, assuming a Type I error rate of 5%, we would expect to see 0.05\*439918=21996 statistical significant results for a given comparison, even if there were truly no differentially methylated CpG sites.

Based on a false discovery rate of 5%, there are 3021 significantly differentially methylated CpGs in the `naïve` vs `rTreg` comparison, while `rTreg` vs `act_rTreg` doesn't show any significant differential methylation.

这个特定数据集感兴趣的生物学问题是发现不同细胞类型之间的差异甲基化探针。但是，正如MDS图中所示，在执行统计分析时，我们还需要考虑另一个因素。在目标文件中，有一个名为 Sample_Source 的列，它指的是从中收集样本的个人。在这个数据集中，每个个体贡献了不止一种细胞类型。例如，单个 M28 提供朴素、rTreg 和act_naive样本。因此，当我们指定设计矩阵时，我们需要包括两个因素：个体和细胞类型。这种分析风格称为配对分析;计算每个个体内细胞类型之间的差异，然后将这些差异平均到个体之间，以确定每个CpG位点的平均甲基化水平是否存在总体显着差异。limma用户指南广泛涵盖了通常用于微阵列实验的不同类型的设计以及如何在R中分析它们。

我们对四种细胞类型之间的成对比较感兴趣，同时考虑到个体之间的差异。我们对limma中的M值矩阵进行此分析，获得每个CpG位点的调节t统计量和相关p值。当用户有许多感兴趣的比较需要测试时，设置模型的一种便捷方法是将对比矩阵与设计矩阵结合使用。对比矩阵将采用与感兴趣的比较相对应的设计矩阵列的线性组合。

由于我们正在执行数十万个假设检验，因此我们需要调整多个检验的 p 值。当进行大量测试时，评估两组之间平均水平变化的统计显着性的常用方法是根据错误发现率分配临界值（Benjamini and Hochberg 1995），而不是未调整的p值。通常使用5%的FDR，这被解释为研究人员愿意接受从显着差异甲基化CpG位点列表中，5%将是错误的发现。如果不针对多次检验调整 p 值，则错误发现的数量将高得令人无法接受。对于此数据集，假设I型错误率为5%，即使确实没有差异甲基化的CpG位点，我们也期望在给定的比较中看到0.05 \* 439918 = 21996的统计显着结果。

基于 5% 的错误发现率，在幼稚与 rTreg 比较中，有 3021 个显着差异化的甲基化 CpG，而 rTreg 与 act_rTreg 没有显示出任何显着的差异甲基化。

``` r
# this is the factor of interest
cellType <- factor(targets$Sample_Group)
# this is the individual effect that we need to account for
individual <- factor(targets$Sample_Source) 

# use the above to create a design matrix
design <- model.matrix(~0+cellType+individual, data=targets)
colnames(design) <- c(levels(cellType),levels(individual)[-1])
 
# fit the linear model 
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(naive-rTreg,
                           naive-act_naive,
                           rTreg-act_rTreg,
                           act_naive-act_rTreg,
                           levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

```

We can extract the tables of differentially expressed CpGs for each comparison, ordered by B-statistic by default, using the `topTable` function in *limma*. The B-statistic is the log-odds of differential methylation, first published by Lonnstedt and Speed (Lonnstedt and Speed [2002](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-lonnstedt2002replicated)). To order by p-value, the user can specify `sort.by="p"`; and in most cases, the ordering based on the p-value and ordering based on the B-statistic will be identical.The results of the analysis for the first comparison, `naive` vs. `rTreg`, can be saved as a `data.frame` by setting `coef=1`. The `coef` parameter explicitly refers to the column in the contrasts matrix which corresponds to the comparison of interest.

我们可以为每个比较提取差异表达的 CpG 表，默认情况下按 B 统计量排序，使用 limma 中的 topTable 函数。B统计量是差异甲基化的对数几率，最初由Lonnstedt和Speed（Lonnstedt和Speed 2002）发表。要按 p 值排序，用户可以指定 sort.by="p";在大多数情况下，基于 p 值的排序和基于 B 统计量的排序将是相同的。第一次比较的分析结果，朴素与rTreg，可以通过设置coef=1保存为data.frame。coef 参数显式引用对比矩阵中对应于感兴趣比较的列。

``` r
# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)
```

The resulting `data.frame` can easily be written to a CSV file, which can be opened in Excel.

``` r
write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)
```

It is always useful to plot sample-wise methylation levels for the top differentially methylated CpG sites to quickly ensure the results make sense (Figure 9). If the plots do not look as expected, it is usually an indication of an error in the code, or in setting up the design matrix. It is easier to interpret methylation levels on the beta value scale, so although the analysis is performed on the M-value scale, we visualise data on the beta value scale. The `plotCpg` function in *minfi* is a convenient way to plot the sample-wise beta values stratified by the grouping variable.

绘制最高差异甲基化 CpG 位点的样本甲基化水平图总是很有用的，可以快速确保结果有意义（图 9）。如果图看起来与预期不符，通常表明代码或设计矩阵设置有误。在β值范围内解释甲基化水平更容易，因此尽管分析是在M值范围内进行的，但我们在β值范围内可视化数据。minfi 中的 plotCpg 函数是绘制按分组变量分层的样本贝塔值的便捷方法。

``` r
# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})
```

#### Differential methylation analysis of regions

Although performing a *probe-wise* analysis is useful and informative, sometimes we are interested in knowing whether several proximal CpGs are concordantly differentially methylated, that is, we want to identify differentially methylated *regions*. There are several Bioconductor packages that have functions for identifying differentially methylated regions from 450k data. Some of the most popular are the `dmrFind` function in the [charm](http://www.bioconductor.org/packages/release/bioc/html/charm.html) package, which has been somewhat superseded for 450k arrays by the `bumphunter` function in [minfi](http://bioconductor.org/packages/release/bioc/html/minfi.html)(Jaffe et al. [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Jaffe2012); Aryee et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Aryee2014)), and, the recently published `dmrcate` in the [DMRcate](https://www.bioconductor.org/packages/release/bioc/html/DMRcate.html) package (Peters et al. [2015](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Peters2015)). They are each based on different statistical methods. In our experience, the `bumphunter` and `dmrFind` functions can be somewhat slow to run unless you have the computer infrastructure to parallelise them, as they use permutations to assign significance. In this workflow, we will perform an analysis using the `dmrcate`. As it is based on *limma*, we can directly use the `design` and `contMatrix` we previously defined.

Firstly, our matrix of M-values is annotated with the relevant information about the probes such as their genomic position, gene annotation, etc. By default, this is done using the `ilmn12.hg19` annotation, but this can be substituted for any argument compatible with the interface provided by the *minfi* package. The *limma* pipeline is then used for differential methylation analysis to calculate moderated t-statistics.

尽管进行探针分析是有用且信息丰富的，但有时我们有兴趣了解几个近端CpG是否一致地差异甲基化，也就是说，我们想要识别差异甲基化区域。有几种Bioconductor软件包具有从450k数据中识别差异甲基化区域的功能。一些最受欢迎的是魅力包中的 dmrFind 函数，它已被 minfi 中的 bumphunter 函数在某种程度上取代了 450k 数组（Jaffe 等人，2012 年;Aryee et al. 2014），以及最近在DMRcate包中发表的dmrcate（Peters et al. 2015）。它们各自基于不同的统计方法。根据我们的经验，碰撞猎手和 dmrFind 函数的运行速度可能有些慢，除非您有计算机基础设施来并行化它们，因为它们使用排列来分配重要性。在此工作流程中，我们将使用 dmrcate 执行分析。由于它基于 limma，我们可以直接使用我们之前定义的设计和 contMatrix。

首先，我们的M值矩阵用探针的相关信息进行注释，例如它们的基因组位置，基因注释等。默认情况下，这是使用 ilmn12.hg19 注释完成的，但这可以替换为与 minfi 包提供的接口兼容的任何参数。然后使用limma管道进行差异甲基化分析，以计算调节的t统计量。

``` r
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "naive - rTreg", arraytype = "450K")
str(myAnnotation)
```

Once we have the relevant statistics for the individual CpGs, we can then use the `dmrcate` function to combine them to identify differentially methylated regions. The main output table `DMRs$results` contains all of the regions found, along with their genomic annotations and p-values.

一旦我们获得了单个CpG的相关统计数据，我们就可以使用dmrcate函数将它们组合起来，以识别差异甲基化区域。主输出表 DMRs\$results 包含找到的所有区域，以及它们的基因组注释和 p 值。

``` r
#endif /* NEWSTUFF */
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
```

As for the probe-wise analysis, it is advisable to visualise the results to ensure that they make sense. The regions can easily be viewed using the DMR.plot function provided in the DMRcate package (Figure 10).

至于探针分析，建议将结果可视化以确保它们有意义。可以使用 DMRcate 包中提供的 DMR.plot 函数轻松查看这些区域（图 10）。

``` r
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]
# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")
```

#### Customising visualisations of methylation data

The Gviz package offers powerful functionality for plotting methylation data in its genomic context. The package vignette is very extensive and covers the various types of plots that can be produced using the Gviz framework. We will plot one of the differentially methylated regions from the DMRcate analysis to demonstrate the type of visualisations that can be created (Figure 11).

We will first set up the genomic region we would like to plot by extracting the genomic coordinates of one of the differentially methylated regions.

Gviz 软件包提供了强大的功能，用于在其基因组环境中绘制甲基化数据。该软件包的小插图非常广泛，涵盖了可以使用 Gviz 框架生成的各种类型的绘图。我们将绘制DMRcate分析中差异甲基化区域之一，以演示可以创建的可视化类型（图11）。

我们将首先通过提取其中一个差异甲基化区域的基因组坐标来设置我们想要绘制的基因组区域。

``` r
# indicate which genome is being used
gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 1
# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))
```

接下来，我们将添加一些感兴趣的基因组注释，例如 CpG 岛和 DNAseI 超敏位点的位置;这可以是您具有可用数据的任何感兴趣的特征或基因组注释。CpG群岛数据是使用Wu等人（2010）发表的方法生成的;DNaseI超敏位点数据是从UCSC基因组浏览器获得的。

``` r
# CpG islands
islandHMM <- read.csv(paste0(dataDirectory,
                             "/model-based-cpg-islands-hg19-chr17.txt"),
                      sep="\t", stringsAsFactors=FALSE, header=FALSE)
head(islandHMM)

islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
islandData
# DNAseI hypersensitive sites
dnase <- read.csv(paste0(dataDirectory,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                  sep="\t",stringsAsFactors=FALSE,header=FALSE)
head(dnase)
dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))),
                     data=dnase[,5])
dnaseData
```

Now, set up the ideogram, genome and RefSeq tracks that will provide context for our methylation data.

现在，设置表意文字，基因组和RefSeq轨道，为我们的甲基化数据提供背景。

``` r
iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)
```

Ensure that the methylation data is ordered by chromosome and base position.

``` r
ann450kOrd <- ann450kSub[order(ann450kSub$chr,ann450kSub$pos),]
head(ann450kOrd)
bValsOrd <- bVals[match(ann450kOrd$Name,rownames(bVals)),]
head(bValsOrd)
```

Create the data tracks using the appropriate track type for each data type.

``` r
# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(ann450kOrd$chr),
                   ranges=IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos),
                   strand=Rle(rep("*",nrow(ann450kOrd))),
                   betas=bValsOrd)
# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])

# methylation data track
methTrack <- DataTrack(range=cpgData, groups=targets$Sample_Group,genome = gen,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)
# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")

# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")

Set up the track list and indicate the relative sizes of the different tracks. Finally, draw the plot using the plotTracks function (Figure 11).

tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack,
               rTrack)
sizes <- c(2,2,5,2,2,2,3) # set up the relative sizes of the tracks
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
```

### Additional analyses

#### Gene ontology testing

Once you have performed a differential methylation analysis, there may be a very long list of significant CpG sites to interpret. One question a researcher may have is, "which gene pathways are over-represented for differentially methylated CpGs?" In some cases it is relatively straightforward to link the top differentially methylated CpGs to genes that make biological sense in terms of the cell types or samples being studied, but there may be many thousands of CpGs significantly differentially methylated. In order to gain an understanding of the biological processes that the differentially methylated CpGs may be involved in, we can perform gene ontology or KEGG pathway analysis using the `gometh` function in the *missMethyl* package (Phipson, Maksimovic, and Oshlack [2016](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Phipson2016)).

Let us consider the first comparison, naive vs rTreg, with the results of the analysis in the `DMPs` table. The `gometh` function takes as input a character vector of the names (e.g. cg20832020) of the significant CpG sites, and optionally, a character vector of all CpGs tested. This is recommended particularly if extensive filtering of the CpGs has been performed prior to analysis. For gene ontology testing (default), the user can specify `collection="GO”`. For testing KEGG pathways, specify `collection="KEGG”`. In the `DMPs` table, the `Name` column corresponds to the CpG name. We will select all CpG sites that have adjusted p-value of less than 0.05.

一旦您进行了差异甲基化分析，可能会有很长的重要 CpG 位点需要解释。研究人员可能遇到的一个问题是，"对于差异甲基化CpG，哪些基因途径被过度代表？在某些情况下，将顶部差异甲基化的CpG与在所研究的细胞类型或样品方面具有生物学意义的基因联系起来相对简单，但可能有数千个CpG显着差异甲基化。为了了解差异甲基化CpG可能参与的生物过程，我们可以使用missMethyl包中的gometh功能进行基因本体或KEGG通路分析（Phipson，Maksimovic和Oshlack 2016）。

让我们考虑第一个比较，朴素与rTreg，分析结果在DMP表中。gometh函数将重要CpG位点的名称（例如cg20832020）的字符向量作为输入，并可选地将所有测试的CpG的字符向量作为输入。如果在分析之前对CpG进行了广泛的过滤，则尤其建议这样做。对于基因本体测试（默认），用户可以指定collection="GO"。要测试 KEGG 途径，请指定 collection="KEGG"。在 DMP 表中，"名称"列对应于 CpG 名称。我们将选择调整后 p 值小于 0.05 的所有 CpG 站点。

``` r
# Get the 
significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
# First 10 significant CpGs
sigCpGs[1:10]
# Total number of significant CpGs at 5% FDR
length(sigCpGs)
# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
# Total number of CpG sites tested
length(all)
```

The `gometh` function takes into account the varying numbers of CpGs associated with each gene on the Illumina methylation arrays. For the 450k array, the numbers of CpGs mapping to genes can vary from as few as 1 to as many as 1200. The genes that have more CpGs associated with them will have a higher probability of being identified as differentially methylated compared to genes with fewer CpGs. We can look at this bias in the data by specifying `plot=TRUE` in the call to `gometh` (Figure 12).

gometh函数考虑了与Illumina甲基化阵列上每个基因相关的不同数量的CpG。对于 450k 阵列，映射到基因的 CpG 数量可以从少至 1 到多至 1200 不等。与具有较少CpG的基因相比，具有更多CpG的基因将具有更高的差异甲基化概率。我们可以通过在调用 gometh 时指定 plot=TRUE 来查看数据中的这种偏差（图 12）。

``` r
par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)
```

The `gst` object is a `data.frame` with each row corresponding to the GO category being tested. Note that the warning regarding multiple symbols will always be displayed as there are genes that have more than one alias, however it is not a cause for concern.

The top 20 gene ontology categories can be displayed using the `topGSA` function. For KEGG pathway analysis, the `topGSA` function will also display the top 20 enriched pathways.

gst 对象是一个 data.frame（数据帧），其中每一行都与正在测试的 GO 类别相对应。请注意，由于有些基因有多个别名，因此始终会显示有关多个符号的警告，但这并不令人担忧。

使用topGSA功能可以显示前20个基因本体类别。对于KEGG通路分析，topGSA功能还将显示前20个富集通路。

``` r
# Top 10 GO categories
topGSA(gst, number=10)
```

From the output we can see many of the top GO categories correspond to immune system and T cell processes, which is unsurprising as the cell types being studied form part of the immune system. Typically, we consider GO categories that have associated false discovery rates of less than 5% to be statistically significant. If there aren't any categories that achieve this significance it may be useful to scan the top 5 or 10 highly ranked GO categories to gain some insight into the biological system.

The `gometh` function only tests GO and KEGG pathways. For a more generalised version of gene set testing for methylation data where the user can specify the gene set to be tested, the `gsameth` function can be used. To display the top 20 pathways, `topGSA` can be called. `gsameth` accepts a single gene set, or a list of gene sets. The gene identifiers in the gene set must be Entrez Gene IDs. To demonstrate `gsameth`, we are using the curated genesets (C2) from the Broad Institute Molecular signatures [database](http://software.broadinstitute.org/gsea/msigdb). These can be downloaded as an `RData` object from the WEHI Bioinformatics [website](http://bioinf.wehi.edu.au/software/MSigDB/).

从输出中，我们可以看到许多顶级GO类别对应于免疫系统和T细胞过程，这并不奇怪，因为正在研究的细胞类型是免疫系统的一部分。通常，我们认为错误发现率低于 5% 的 GO 类别具有统计显著性。如果没有任何类别可以达到这种意义，那么扫描排名前5或10位的GO类别可能会很有用，以深入了解生物系统。

gometh函数仅测试GO和KEGG途径。对于甲基化数据的基因集测试的更通用版本，用户可以指定要测试的基因集，可以使用gsameth函数。要显示前 20 条通路，可以调用 topGSA。Gsameth接受单个基因集或基因集列表。基因集中的基因标识符必须是 Entrez 基因 ID。为了证明gsameth，我们使用了来自Broad Institute Molecular Signatures数据库中的精选基因集（C2）。这些可以从WEHI生物信息学网站下载RData对象。

``` r
# load Broad human curated (C2) gene sets
load(paste(dataDirectory,"human_c2_v5.rdata",sep="/"))
# perform the gene set test(s)
gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=Hs.c2)

# top 10 gene sets
topGSA(gsa, number=10)
```

While gene set testing is useful for providing some biological insight in terms of what pathways might be affected by abberant methylation, care should be taken not to over-interpret the results. Gene set testing should be used for the purpose of providing some biological insight that ideally would be tested and validated in further laboratory experiments. It is important to keep in mind that we are not observing gene level activity such as in RNA-Seq experiments, and that we have had to take an extra step to associate CpGs with genes.

虽然基因集测试有助于提供一些生物学见解，了解哪些途径可能受到异常甲基化的影响，但应注意不要过度解释结果。基因集测试应该用于提供一些生物学见解，理想情况下将在进一步的实验室实验中进行测试和验证。重要的是要记住，我们没有观察到像RNA-Seq实验那样的基因水平活性，并且我们不得不采取额外的步骤将CpG与基因联系起来。

#### Differential variability

Rather than testing for differences in mean methylation, we may be interested in testing for differences between group variances. For example, it has been hypothesised that highly variable CpGs in cancer may contribute to tumour heterogeneity (Hansen et al. [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Hansen2011)). Hence we may be interested in CpG sites that are consistently methylated in one group, but variably methylated in another group.

Sample size is an important consideration when testing for differentially variable CpG sites. In order to get an accurate estimate of the group variances, larger sample sizes are required than for estimating group means. A good rule of thumb is to have at least ten samples in each group (Phipson and Oshlack [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Phipson2014)). To demonstrate testing for differentially variable CpG sites, we will use a publicly available dataset on ageing [GSE30870](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30870), where whole blood samples were collected from 18 centenarians and 18 newborns and profiled for methylation on the 450k array (Heyn et al. [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Heyn2012)). The data (`age.rgSet`) and sample information (`age.targets`) have been included as an R data object in both the workflow package or the data archive you downloaded from [figshare](https://figshare.com/s/7a37f43c0ca2fec4669e). We can load the data using the `load` command, after which it needs to be normalised and filtered as previously described.

与其检验平均甲基化的差异，我们可能感兴趣的是检验组方差之间的差异。例如，有人假设癌症中高度可变的CpGs可能导致肿瘤异质性（Hansen等人，2011）。因此，我们可能对在一个组中始终甲基化但在另一组中可变甲基化的CpG位点感兴趣。

在测试差异可变的CpG位点时，样本量是一个重要的考虑因素。为了获得组方差的准确估计值，需要比估计组均值更大的样本数量。一个好的经验法则是每组中至少有十个样本（Phipson and Oshlack 2014）。为了演示对差异可变CpG位点的测试，我们将使用有关老化GSE30870的公开数据集，其中从18名百岁老人和18名新生儿收集全血样本，并在450k阵列上分析甲基化（Heyn等人，2012）。数据 （age.rgSet） 和示例信息 （age.targets） 已作为 R 数据对象包含在工作流包或从 figshare 下载的数据存档中。我们可以使用 load 命令加载数据，之后需要按照前面所述对其进行规范化和过滤。

``` r
load(file.path(dataDirectory,"ageData.RData"))

# calculate detection p-values
age.detP <- detectionP(age.rgSet)

# pre-process the data after excluding poor quality samples
age.mSetSq <- preprocessQuantile(age.rgSet)

# add sex information to targets information
age.targets$Sex <- getSex(age.mSetSq)$predictedSex

# ensure probes are in the same order in the mSetSq and detP objects
age.detP <- age.detP[match(featureNames(age.mSetSq),rownames(age.detP)),]
# remove poor quality probes
keep <- rowSums(age.detP < 0.01) == ncol(age.detP) 
age.mSetSqFlt <- age.mSetSq[keep,]

# remove probes with SNPs at CpG or single base extension (SBE) site
age.mSetSqFlt <- dropLociWithSnps(age.mSetSqFlt, snps = c("CpG", "SBE"))

# remove cross-reactive probes
keep <- !(featureNames(age.mSetSqFlt) %in% xReactiveProbes$TargetID)
age.mSetSqFlt <- age.mSetSqFlt[keep,] 
```

As this dataset contains samples from both males and females, we can use it to demonstrate the effect of removing sex chromosome probes on the data. The MDS plots below show the relationship between the samples in the ageing dataset before and after sex chromosome probe removal (Figure 13). It is apparent that before the removal of sex chromosome probes, the sample cluster based on sex in the second principal component. When the sex chromosome probes are removed, age is the largest source of variation present and the male and female samples no longer form separate clusters.

由于该数据集包含来自男性和女性的样本，我们可以使用它来证明去除性染色体探针对数据的影响。下面的MDS图显示了性染色体探针去除前后老化数据集中样本之间的关系（图13）。很明显，在去除性染色体探针之前，样本簇基于性别的第二主成分。当性染色体探针被移除时，年龄是存在的最大变异来源，男性和女性样本不再形成单独的簇。

``` r
# tag sex chromosome probes for removal
keep <- !(featureNames(age.mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                            c("chrX","chrY")])

age.pal <- brewer.pal(8,"Set1")
par(mfrow=c(1,2))
plotMDS(getM(age.mSetSqFlt), top=1000, gene.selection="common", 
        col=age.pal[factor(age.targets$Sample_Group)], labels=age.targets$Sex, 
        main="With Sex CHR Probes")
legend("topleft", legend=levels(factor(age.targets$Sample_Group)), 
       text.col=age.pal)

plotMDS(getM(age.mSetSqFlt[keep,]), top=1000, gene.selection="common", 
        col=age.pal[factor(age.targets$Sample_Group)], labels=age.targets$Sex, 
        main="Without Sex CHR Probes")
legend("top", legend=levels(factor(age.targets$Sample_Group)),
       text.col=age.pal)

# remove sex chromosome probes from data
age.mSetSqFlt <- age.mSetSqFlt[keep,]
```

We can test for differentially variable CpGs using the `varFit` function in the *missMethyl* package. The syntax for specifying which groups we are interested in testing is slightly different to the standard way a model is specified in `limma`, particularly for designs where an intercept is fitted (see *missMethyl* [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/missMethyl/inst/doc/missMethyl.pdf) for further details). For the ageing data, the design matrix includes an intercept term, and a term for age. The `coef` argument in the `varFit` function indicates which columns of the design matrix correspond to the intercept and grouping factor. Thus, for the ageing dataset we set `coef=c(1,2)`. Note that design matrices without intercept terms are permitted, with specific contrasts tested using the `contrasts.varFit` function.

我们可以使用 missMethyl 封装中的 varFit 函数来测试差分可变的 CpG。用于指定我们有兴趣测试的组的语法与在 limma 中指定模型的标准方式略有不同，特别是对于拟合截距的设计（有关更多详细信息，请参阅 missMethyl 小插图）。对于老化数据，设计矩阵包括截距项和年龄项。varFit 函数中的 coef 参数指示设计矩阵的哪些列对应于截距和分组因子。因此，对于老化数据集，我们设置 coef=c（1，2）。请注意，允许使用没有截距项的设计矩阵，并使用 contrasts.varFit 函数测试特定的对比度。

``` r
# get M-values for analysis
age.mVals <- getM(age.mSetSqFlt)

design <- model.matrix(~factor(age.targets$Sample_Group)) 
# Fit the model for differential variability
# specifying the intercept and age as the grouping factor
fitvar <- varFit(age.mVals, design = design, coef = c(1,2))

# Summary of differential variability
summary(decideTests(fitvar))

topDV <- topVar(fitvar, coef=2)
# Top 10 differentially variable CpGs between old vs. newborns
topDV
```

Similarly to the differential methylation analysis, is it useful to plot sample-wise beta values for the differentially variable CpGs to ensure the significant results are not driven by artifacts or outliers (Figure 14).

与差异甲基化分析类似，绘制差异可变 CpG 的样品 β 值是否有用，以确保显著结果不受伪影或异常值的影响（图 14）。

``` r
# get beta values for ageing data
age.bVals <- getBeta(age.mSetSqFlt)
```

``` r
par(mfrow=c(2,2))
sapply(rownames(topDV)[1:4], function(cpg){
  plotCpg(age.bVals, cpg=cpg, pheno=age.targets$Sample_Group, 
          ylab = "Beta values")
})
```

An example of testing for differential variability when the design matrix does not have an intercept term is detailed in the *missMethyl* [vignette](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/%22http://www.bioconductor.org/packages/release/bioc/vignettes/missMethyl/inst/doc/missMethyl.pdf%22).

#### Cell type composition

As methylation is cell type specific and methylation arrays provide CpG methylation values for a population of cells, biological findings from samples that are comprised of a mixture of cell types, such as blood, can be confounded with cell type composition (Jaffe and Irizarry [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Jaffe2014)). The *minfi* function `estimateCellCounts` facilitates the estimation of the level of confounding between phenotype and cell type composition in a set of samples. The function uses a modified version of the method published by Houseman et al. ([2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Houseman2012)) and the package `FlowSorted.Blood.450k`, which contains 450k methylation data from sorted blood cells, to estimate the cell type composition of blood samples.

由于甲基化是细胞类型特异性的，甲基化阵列为细胞群提供CpG甲基化值，因此由细胞类型混合物（例如血液）组成的样品的生物学发现可能与细胞类型组成混淆（Jaffe和Irizarry 2014）。minfi 函数 estimateCellCounts 有助于估计一组样本中表型和细胞类型组成之间的混淆水平。该函数使用Houseman等人（2012）发表的方法的修改版本和包含来自分选血细胞的450k甲基化数据的FlowSorted.Blood.450k包来估计血液样本的细胞类型组成。

``` r
# load sorted blood cell data package
library(FlowSorted.Blood.450k)
# ensure that the "Slide" column of the rgSet pheno data is numeric
# to avoid "estimateCellCounts" error
pData(age.rgSet)$Slide <- as.numeric(pData(age.rgSet)$Slide)
# estimate cell counts
cellCounts <- estimateCellCounts(age.rgSet)
```

``` r
# plot cell type proportions by age
par(mfrow=c(1,1))
a = cellCounts[age.targets$Sample_Group == "NewBorns",]
b = cellCounts[age.targets$Sample_Group == "OLD",]
boxplot(a, at=0:5*3 + 1, xlim=c(0, 18), ylim=range(a, b), xaxt="n", 
        col=age.pal[1], main="", ylab="Cell type proportion")
boxplot(b, at=0:5*3 + 2, xaxt="n", add=TRUE, col=age.pal[2])
axis(1, at=0:5*3 + 1.5, labels=colnames(a), tick=TRUE)
legend("topleft", legend=c("NewBorns","OLD"), fill=age.pal)
```

As reported by Jaffe and Irizarry ([2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Jaffe2014)), the preceding plot demonstrates that differences in blood cell type proportions are strongly confounded with age in this dataset (Figure 15). Performing cell composition estimation can alert you to potential issues with confounding when analysing a mixed cell type dataset. Based on the results, some type of adjustment for cell type composition may be considered, although a naive cell type adjustment is not recommended. Jaffe and Irizarry ([2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Jaffe2014)) outline several strategies for dealing with cell type composition issues.

正如Jaffe和Irizarry（2014）所报告的，上图显示，在该数据集中，血细胞类型比例的差异与年龄密切相关（图15）。在分析混合细胞类型数据集时，进行细胞组成估计可提醒您注意潜在的混杂问题。根据结果，可以考虑对细胞类型组成进行某种类型的调整，但不建议进行天真的细胞类型调整。Jaffe和Irizarry（2014）概述了处理细胞类型组成问题的几种策略。

### Discussion

Here we present a commonly used workflow for methylation array analysis based on a series of Bioconductor packages. While we have not included all the possible functions or analysis options that are available for detecting differential methylation, we have demonstrated a common and well used workflow that we regularly use in our own analysis. Specifically, we have not demonstrated more complex types of analyses such as removing unwanted variation in a differential methylation study (Maksimovic et al. [2015](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Maksimovic2015); Leek et al. [2012](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Leek2012); Teschendorff, Zhuang, and Widschwendter [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Teschendorff2011)), block finding (Hansen et al. [2011](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Hansen2011); Aryee et al. [2014](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Aryee2014)) or A/B compartment prediction (Fortin and Hansen [2015](http://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#ref-Fortin2015)). Our differential methylation workflow presented here demonstrates how to read in data, perform quality control and filtering, normalisation and differential methylation testing. In addition we demonstrate analysis for differential variability, gene set testing and estimating cell type composition. One important aspect of exploring results of an analysis is visualisation and we also provide an example of generating region-level views of the data.

在这里，我们介绍了一种基于一系列Bioconductor包的甲基化阵列分析的常用工作流程。虽然我们没有包括可用于检测差异甲基化的所有可能的功能或分析选项，但我们已经展示了一种常见且使用良好的工作流程，我们经常在自己的分析中使用。具体来说，我们还没有证明更复杂的分析类型，例如在差异甲基化研究中去除不需要的变异（Maksimovic等人，2015;Leek等人，2012;Teschendorff， Zhuang， and Widschwendter 2011）， 块查找 （Hansen et al. 2011;Aryee等人，2014）或A / B隔室预测（Fortin和Hansen 2015）。我们这里介绍的差异甲基化工作流程演示了如何读取数据、执行质量控制和过滤、归一化和差异甲基化测试。此外，我们还演示了差异变异性分析、基因集测试和细胞类型组成的估计。探索分析结果的一个重要方面是可视化，我们还提供了一个生成数据的区域级视图的示例。
