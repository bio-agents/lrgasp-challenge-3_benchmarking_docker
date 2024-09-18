# LRGASP challenge 3 mouse benchmarking Docker
The Long-read RNA-seq Genome Annotation Assessment Project \([LRGASP](https://www.gencodegenes.org/pages/LRGASP/)\) consortium is organizing a systematic evaluation of different transcript computational identification and quantification methods using long-read sequencing technologies such as [PacBio](https://www.pacb.com/) and [Oxford Nanopore](https://nanoporetech.com/). We are interested in characterizing the strengths and potential remaining challenges in using these technologies to annotate and quantify the transcriptomes of both model and non-model organisms.

The consortium will generate cDNA and direct RNA datasets using different platforms and protocols in [human](https://www.gencodegenes.org/human/), [mouse](https://www.gencodegenes.org/mouse/), and manatee samples. Participants will be provided with the data to generate annotations of expressed genes and transcripts and measure their expression levels. Evaluators from different institutions will determine which pipelines have the highest accuracy for different aspects, including transcription detection, quantification, and differential expression.

Challenge 3 of LRGASP was to evaluate the performance of different [mouse](https://www.gencodegenes.org/mouse/) and manatee samples pipelines to reconstruct a transcriptome **without** using previous annotation or a reference genome. However, for its evaluation, we will use a *de novo* ONT-based genome available here to obtain information about the submitted transcriptomes.

You can find the LRGASP data and submission guidelines at https://lrgasp.github.io/lrgasp-submissions.

This repository branch contains the docker image for running the LRGASP benchmarking https://openebench.bsc.es/benchmarking/OEBC010 on the OpenEBench platform. 

This repository branch hosts the Docker image for executing the LRGASP benchmarking,specifically tailored for **mouse** transcriptome analysis. For manatee benchmarking, please visit the dedicated GitHub repository at: https://github.com/TianYuan-Liu/lrgasp-challenge-3_manatee_benchmarking_docker 

## Usage
1. Install [OpenEBench VRE executor](https://github.com/inab/vre-process_nextflow-executor/blob/master/INSTALL.md) and go to the **tests** folder
```
cd vre-process_nextflow-executor
source .py3Env/bin/activate
cd tests
```

2. Clone the [LRGASP mouse workflow repository](https://github.com/TianYuan-Liu/lrgasp-challenge-3_benchmarking_workflow) and rename the folder to LRGASP
```
git clone https://github.com/TianYuan-Liu/lrgasp-challenge-3_benchmarking_workflow.git
mv lrgasp-challenge-3_benchmarking_workflow LRGASP
```
3. Materialize both the containers and datasets needed by the LRGASP test:
```
cd LRGASP
bash ./materialize-data.sh
```
4. Run the tests from LRGASP example
```
cd ../..
./test_VRE_NF_RUNNER.sh LRGASP
```