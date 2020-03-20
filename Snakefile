
## Configuration file
# for exome
import os
if len(config) == 0:
  if os.path.isfile("./config.yaml"):
    configfile: "./config.yaml"
  else:
    sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Make sure that all expected variables from the config file are in the config dictionary
configvars = ['annotation', 'organism', 'build', 'release', 'txome', 'genome', 'gtf', 'readlength', 'fldMean', 'fldSD', 'metatxt', 'design', 'contrast', 'genesets', 'ncores', 'FASTQ', 'fqext1', 'fqext2', 'fqsuffix', 'output', 'useCondaR', 'Rbin', 'run_trimming', 'run_bwa']
for k in configvars:
	if k not in config:
		config[k] = None

## If any of the file paths is missing, replace it with ""
def sanitizefile(str):
	if str is None:
		str = ''
	return str
	
config['txome'] = sanitizefile(config['txome'])
config['gtf'] = sanitizefile(config['gtf'])
config['genome'] = sanitizefile(config['genome'])
config['bwaindex'] = sanitizefile(config['bwaindex'])
config['metatxt'] = sanitizefile(config['metatxt'])
config['metacsv'] = os.path.splitext(os.path.basename(config['metatxt']))[0]+'.csv'

## Read metadata
if not os.path.isfile(config["metatxt"]):
  sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

import pandas as pd
samples = pd.read_csv(config["metatxt"], sep='\t')

if not set(['names','type']).issubset(samples.columns):
  sys.exit("Make sure 'names' and 'type' are columns in " + config["metatxt"])


## Sanitize provided input and output directories
import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

proj_dir = os.path.abspath(os.path.normpath(getpath(config["proj_dir"])))
outputdir = os.path.abspath(getpath(config["output"])) + "/"
FASTQdir = getpath(config["FASTQ"])

## Define the conda environment for all rules using R
if config["useCondaR"] == True:
	Renv = "envs/environment_R.yaml"
else:
	Renv = "envs/environment.yaml"

## Define the R binary
Rbin = config["Rbin"]

def dbtss_output(wildcards):
	input = []
	input.extend(expand(outputdir + "dbtss_coverage/{sample}_dbtss_coverage_over_10.txt", sample = samples.names[samples.type == 'PE'].values.tolist()))
	return input
	
def jbrowse_output(wildcards):
  input = []
  input.extend(expand("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/samples/{sample}.bw", sample = samples.names[samples.type == 'PE'].values.tolist()))
  input.append("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/trackList.json")
  return input
  
	
## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
rule all:
	input:
		outputdir + "MultiQC/multiqc_report.html",
		# dbtss_output,
		jbrowse_output

rule setup:
	input:
		outputdir + "Rout/pkginstall_state.txt",
		outputdir + "Rout/softwareversions.done"

## Install R packages
rule pkginstall:
	input:
		script = "scripts/install_pkgs.R"
	output:
	  outputdir + "Rout/pkginstall_state.txt"
	params:
		flag = config["annotation"],
		ncores = config["ncores"],
		organism = config["organism"]
	priority:
		50
	conda:
		Renv
	log:
		outputdir + "Rout/install_pkgs.Rout"
	benchmark:
	  outputdir + "benchmarks/install_pkgs.txt"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args outtxt='{output}' ncores='{params.ncores}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## FastQC on original (untrimmed) files
rule runfastqc:
	input:
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_val_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_val_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## bwa alignment
rule runbwa:
	input:
		expand(outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.names.values.tolist()),
		expand(outputdir + "bwabigwig/{sample}_Aligned.sortedByCoord.out.bw", sample = samples.names.values.tolist())
		
## DBTSS coverage
rule rundbtss:
	input:
		expand(outputdir + "dbtss_coverage/{sample}/{sample}_dbtss_coverage.txt", sample = samples.names.values.tolist())

## List all the packages that were used by the R analyses
rule listpackages:
	log:
		outputdir + "Rout/list_packages.Rout"
	params:
		Routdir = outputdir + "Rout",
		outtxt = outputdir + "R_package_versions.txt",
		script = "scripts/list_packages.R"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args Routdir='{params.Routdir}' outtxt='{params.outtxt}'" {params.script} {log}'''

## Print the versions of all software packages
rule softwareversions:
	output:
		touch(outputdir + "Rout/softwareversions.done")
	conda:
		"envs/environment.yaml"
	shell:
		"echo -n 'ARMOR version ' && cat version; "
		"salmon --version; trim_galore --version; "
		"echo -n 'cutadapt ' && cutadapt --version; "
		"fastqc --version; bwa --version; samtools --version; multiqc --version; "
		"bedtools --version"

## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
	input:
		fastq = FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = outputdir + "FastQC"
	log:
		outputdir + "logs/fastqc_{sample}.log"
	benchmark:
		outputdir + "benchmarks/fastqc_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## FastQC, trimmed reads
rule fastqctrimmed:
	input:
		fastq = outputdir + "FASTQtrimmed/{sample}.fq.gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = outputdir + "FastQC"
	log:
		outputdir + "logs/fastqc_trimmed_{sample}.log"
	benchmark:
		outputdir + "benchmarks/fastqc_trimmed_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"
		
## Rseqc gene body coverage plot
rule genebodycoverage:
	input:
		bigwig = outputdir + "bwabigwig/{sample}_Aligned.sortedByCoord.out.bw"
	output:
		coverage_txt = outputdir + "rseqc/{sample}.geneBodyCoverage.txt"
	params:
		sample = outputdir + "rseqc/{sample}",
		bed = config["bed"]
	log:
		outputdir + "logs/genebodycoverage_{sample}.log"
	benchmark:
		outputdir + "benchmarks/genebodycoverage_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'geneBody_coverage.py version:\n' > {log}; geneBody_coverage.py --version >> {log}; "
		"geneBody_coverage2.py -r {params.bed} -i {input.bigwig}  -o {params.sample}"

# The config.yaml files determines which steps should be performed
def multiqc_input(wildcards):
	input = []
	input.extend(expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist()))
	input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()))
	input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()))
	if config["run_genebodycoverage"]:
	  input.extend(expand(outputdir + "rseqc/{sample}." + "geneBodyCoverage.txt", sample = samples.names[samples.type == 'PE'].values.tolist()))
	if config["run_trimming"]:
		input.extend(expand(outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.names[samples.type == 'SE'].values.tolist()))
		input.extend(expand(outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz", sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz", sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_val_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_val_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()))
	if config["run_bwa"]:
		input.extend(expand(outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.names.values.tolist()))
	return input

## Determine the input directories for MultiQC depending on the config file
def multiqc_params(wildcards):
	param = [outputdir + "FastQC"]
	if config["run_trimming"]:
		param.append(outputdir + "FASTQtrimmed")
	if config["run_bwa"]:
		param.append(outputdir + "bwa")
	return param

## MultiQC
rule multiqc:
	input:
		multiqc_input
	output:
		outputdir + "MultiQC/multiqc_report.html"
	params:
		inputdirs = multiqc_params,
		MultiQCdir = outputdir + "MultiQC"
	log:
		outputdir + "logs/multiqc.log"
	benchmark:
		outputdir + "benchmarks/multiqc.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdirs} -f -o {params.MultiQCdir}"


## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgaloreSE:
	input:
		fastq = FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	params:
		FASTQtrimmeddir = outputdir + "FASTQtrimmed"
	log:
		outputdir + "logs/trimgalore_{sample}.log"
	benchmark:
		outputdir + "benchmarks/trimgalore_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"

rule trimgalorePE:
	input:
		fastq1 = FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz",
		fastq2 = FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz",
		outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz"
	params:
		FASTQtrimmeddir = outputdir + "FASTQtrimmed"
	log:
		outputdir + "logs/trimgalore_{sample}.log"
	benchmark:
		outputdir + "benchmarks/trimgalore_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2}"

## ------------------------------------------------------------------------------------ ##
## bwa mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with bwa
rule bwaPE:
	input:
		fastq1 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz",
		fastq2 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz"
	output:
		bam = outputdir + "bwa/{sample}/{sample}_Aligned.out.bam"
	threads:
		config["ncores"]
	log:
		outputdir + "logs/bwa_{sample}.log"
	benchmark:
		outputdir + "benchmarks/bwa_{sample}.txt"
	params:
		bwaindex = config["bwaindex"],
		bwadir = outputdir + "bwa"
	threads:
	  config["ncores"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bwa' > {log}; "
		"bwa mem -M -t {threads} {params.bwaindex} {input.fastq1} {input.fastq2} | samtools view -Sbo {output.bam}"

# convert and sort sam files
rule bamsort:
	input:
		bam = outputdir + "bwa/{sample}/{sample}_Aligned.out.bam"
	output:
		sorted_bam = outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	log:
		outputdir + "logs/samtools_sort_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_sort_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools sort -O bam -o {output.sorted_bam} {input.bam}"

## Index bam files
rule bamindexbwa:
	input:
		bam = outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	log:
		outputdir + "logs/samtools_index_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_index_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

# section        

rule mark_duplicates:
  input:
    outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
  output:
    bam=temp(outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.deduped.bam"),
    metrics="qc/dedup/{sample}.metrics.txt"
  log:
    "logs/picard/dedup/{sample}.log"
  params:
    config["params"]["picard"]["MarkDuplicates"]
  wrapper:
    "0.26.1/bio/picard/markduplicates"

## Index bam files
rule bamindexdeduped:
	input:
	  bam = outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.deduped.bam"
	output:
		outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.deduped.bam.bai"
	log:
		outputdir + "logs/samtools_index_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_index_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"


rule recalibrate_base_qualities:
  input:
    bam = outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.deduped.bam",
    bai = outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.deduped.bam.bai",
    ref = config["genome"],
    known=config["known-variants"]
  output:
    bam=protected("recal/{sample}.bam")
  params:
  log:
    "logs/gatk/bqsr/{sample}.log"
  wrapper:
    "0.27.1/bio/gatk/baserecalibrator"

# section        

## Convert BAM files to bigWig
rule bigwigbwa:
	input:
		bam = "recal/{sample}.bam"
	output:
		outputdir + "bwabigwig/{sample}_Aligned.sortedByCoord.out.bw"
	params:
		bwabigwigdir = outputdir + "bwabigwig",
		chrl = config["chrom_sizes"]
	log:
		outputdir + "logs/bigwig_{sample}.log"
	benchmark:
		outputdir + "benchmarks/bigwig_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | LC_COLLATE=C sort -k1,1 -k2,2n > "
		"{params.bwabigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.bwabigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{params.chrl} {output}; rm -f {params.bwabigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"

## ------------------------------------------------------------------------------------ ##
## Input variable check
## ------------------------------------------------------------------------------------ ##
def geneset_param(wildcards):
	if config["run_camera"]:
                gs = config["genesets"].replace(" ", "") if config["genesets"] is not None else "NOTDEFINED"
		return "genesets='" + gs + "'"
	else:
		return ""


## check design matrix and contrasts
rule checkinputs:
    input:
        "config.yaml",
        script = "scripts/check_input.R"
    output:
        outputdir + "Rout/check_input.txt"
    log:
        outputdir + "Rout/check_input.Rout"
    benchmark:
    	outputdir + "benchmarks/check_input.txt"
    params:
        gtf = config["gtf"],
        genome = config["genome"],
        txome = config["txome"],
        fastqdir = config["FASTQ"],
        metatxt = config["metatxt"],
        design = config["design"].replace(" ", "") if config["design"] is not None else "NOTDEFINED",
        contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "NOTDEFINED",
        annotation = config["annotation"].replace(" ", "") if config["annotation"] is not None else "NOTDEFINED",
        genesets = geneset_param,
        fqsuffix = str(config["fqsuffix"]),
        fqext1 = str(config["fqext1"]),
        fqext2 = str(config["fqext2"]),
        organism = config["organism"]    
    conda:
	    Renv
    shell:
        '''{Rbin} CMD BATCH --no-restore --no-save "--args metafile='{params.metatxt}' design='{params.design}' contrast='{params.contrast}' outFile='{output}' gtf='{params.gtf}' genome='{params.genome}' fastqdir='{params.fastqdir}' fqsuffix='{params.fqsuffix}' fqext1='{params.fqext1}' fqext2='{params.fqext2}' txome='{params.txome}' organism='{params.organism}' {params.genesets} annotation='{params.annotation}'" {input.script} {log};
        cat {output}
        '''
       

## ------------------------------------------------------------------------------------ ##
## Differential expression
## ------------------------------------------------------------------------------------ ##
rule edgeR:
	input:
		outputdir + "Rout/pkginstall_state.txt",
		rds = outputdir + "seurat/unfiltered_seu.rds",
		script = "scripts/run_render.R",
		template = "scripts/edgeR_dge.Rmd"
	output:
		html = outputdir + "outputR/edgeR_dge.html",
		rds = outputdir + "outputR/edgeR_dge.rds"
	params:
		directory = outputdir + "outputR",
		organism = config["organism"],        
                design = config["design"].replace(" ", "") if config["design"] is not None else "",
                contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "",
		genesets = geneset_param
	log:
		outputdir + "Rout/run_dge_edgeR.Rout"
	benchmark:
		outputdir + "benchmarks/run_dge_edgeR.txt"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' organism='{params.organism}' design='{params.design}' contrast='{params.contrast}' {params.genesets} rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='edgeR_dge.html'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Differential transcript usage
## ------------------------------------------------------------------------------------ ##
## DRIMSeq
rule DRIMSeq:
	input:
	    outputdir + "Rout/pkginstall_state.txt",
		rds = outputdir + "outputR/edgeR_dge.rds",
		script = "scripts/run_render.R",
		template = "scripts/DRIMSeq_dtu.Rmd"
	output:
		html = outputdir + "outputR/DRIMSeq_dtu.html",
		rds = outputdir + "outputR/DRIMSeq_dtu.rds"
	params:
		directory = outputdir + "outputR",
		organism = config["organism"],
		ncores = config["ncores"],
                design = config["design"].replace(" ", "") if config["design"] is not None else "",
                contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else ""
	log:
		outputdir + "Rout/run_dtu_drimseq.Rout"
	benchmark:
		outputdir + "benchmarks/run_dtu_drimseq.txt"
	conda:
		Renv
	threads:
		config["ncores"]
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' design='{params.design}' contrast='{params.contrast}' ncores='{params.ncores}' rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='DRIMSeq_dtu.html'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## shiny app
## ------------------------------------------------------------------------------------ ##
def shiny_input(wildcards):
	input = [outputdir + "Rout/pkginstall_state.txt"]
	if config["run_bwa"]:
		input.extend(expand(outputdir + "bwabigwig/{sample}_Aligned.sortedByCoord.out.bw", sample = samples.names.values.tolist()))
	return input

def shiny_params(wildcards):
	param = ["outputdir='" + outputdir + "outputR'"]
	if config["run_bwa"]:
		param.append("bigwigdir='" + outputdir + "bwabigwig'")
	return param

## shiny
rule shiny:
	input:
		shiny_input,
		rds = "outputR/edgeR_dge.rds",
		script = "scripts/run_render.R",
		gtf = config["gtf"],
		template = "scripts/prepare_shiny.Rmd"
	output:
		html = outputdir + "outputR/prepare_shiny.html",
		rds = outputdir + "outputR/shiny_sce.rds"
	params:
		p = shiny_params
	log:
		outputdir + "Rout/prepare_shiny.Rout"
	benchmark:
		outputdir + "benchmarks/prepare_shiny.txt"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' gtffile='{input.gtf}' rmdtemplate='{input.template}' outputfile='prepare_shiny.html' {params.p}" {input.script} {log}'''
		
## ------------------------------------------------------------------------------------ ##
## dbtss coverage mapping
## ------------------------------------------------------------------------------------ ##
## compute coverage from dbtss
rule dbtss:
	input:
		sorted_bam = outputdir + "bwa/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		coverage_txt = outputdir + "dbtss_coverage/{sample}_dbtss_coverage_over_10.txt"
	threads:
		config["ncores"]
	log:
		outputdir + "logs/dbtss_{sample}.log"
	benchmark:
		outputdir + "benchmarks/dbtss_{sample}.txt"
	params:
		dbtss_bed = config["dbtss_bed"],
		genome_file = config["chrom_sizes"]
	conda:
		"envs/environment.yaml"
	shell:
	  "bedtools coverage -sorted -s -g {params.genome_file} -a {params.dbtss_bed} -b {input.sorted_bam} | awk '$4 > 10 {{print}}' > {output.coverage_txt}"

## ------------------------------------------------------------------------------------ ##
## configure jbrowse
## ------------------------------------------------------------------------------------ ##
## configure jbrowse

rule jbrowsemeta:
  input: 
    metatxt = config['metatxt']
  output:
    metacsv = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + config['metacsv']
  threads:
    config['ncores']
  conda: 
    "envs/environment.yaml"
  shell:
    "cp {input.metatxt} {output.metacsv}"

rule jbrowse:
	input:
		hisatbigwig = outputdir + "bwabigwig/{sample}_Aligned.sortedByCoord.out.bw"
	output:
	  bigwig_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/samples/{sample}.bw"
	threads:
		config["ncores"]
	params:
	  proj_name = os.path.basename(proj_dir)
	conda:
		"envs/environment.yaml"
	shell:
	  "ln -s {input.hisatbigwig} {output.bigwig_symlink}"
	  
rule jbrowsetracklist:
	input:
		script = "scripts/format_tracklist_json.py",
		refdir = config["refdir"],
		gff = config["gff"],
		gff_tbi = config["gff_tbi"],
		refseq = config["refseq"]
	output:
	  tracklist_json = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/trackList.json",
	  refdir_symlink = directory("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/reference/"),
	  gff_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/Homo_sapiens.GRCh38.87.sorted.gff3.gz",
	  gff_tbi_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/Homo_sapiens.GRCh38.87.sorted.gff3.gz.tbi",
	  refseq_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/seq/refSeqs.json"
	threads:
		config["ncores"]
	params:
	  proj_name = os.path.basename(proj_dir),
	  metacsv = os.path.basename(proj_dir) + "/" + config['metacsv']
	conda:
		"envs/environment.yaml"
	shell:
	  "{input.script} '{params.proj_name}' '{params.metacsv}';"
	  "ln -sr {input.refdir} {output.refdir_symlink};"
	  "ln -s {input.gff} {output.gff_symlink};"
	  "ln -s {input.gff_tbi} {output.gff_tbi_symlink};"
	  "ln -s {input.refseq} {output.refseq_symlink};"

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")
