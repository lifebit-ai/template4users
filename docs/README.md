# bcftools-s3-nf

## Pipeline description

A nextflow pipeline with bcftools subcomamnds to run quick actions with VCF files directly from s3.

## Input

### `--input`

Input file is a VCF file with s3 path meantioned in each new line. Example - `--input file.vcf` (contents shown below)

```
s3://lifebit-featured-datasets/pipelines/vcf-aggregation/test_data/large_tests_vcfs/test_sample_1.vcf.gz
```

## Output

Output is a uncompressed VCF file format.

## Usage

```bash
nextflow run main.nf --input file.vcf --input_index file.vcf.csi -with-docker
```

## Options

### `--genomic_region`

The genomic range which need to be merged. Example - `--genomic_region "chr1:146421218-165998654"`

### `--bcftools_view_options`

Any addional `bcftools merge` options. Example - `--bcftools_view_options "--samples id1 --force-samples"`


<!-- For Sphinx doc, This option will be auto rendered help() section from Nextflow main.nf in the doc build -->


<!------------------
Build of this doc in github handle by - .github/workflows/build-deploy-doc.yml

To build this doc locally follow these steps.

Needs to have installed - 
1. sphinx
2. sphinx-rtd-theme
3. nextflow

Supposing your currently in base directory of the pipeline -
```
cd docs && bash src/pre-build.sh
cp README.md src
cd src && make html 
```
index.html will be generated in `docs/src/build/html` folder
-->
