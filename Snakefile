configfile: "example_config.yaml"
DATASETS = [d for d in config for s in config[d]]
TEMPLATES = [s for d in config for s in config[d]]

rule all:
    input: 
        expand("consensus/{dataset}/{template}_dUMI.fasta", zip, dataset = DATASETS, template = TEMPLATES) 

rule filter_fastq:
    input:
        "fastq/{dataset}.fastq"
    output:
        "fastq/{dataset}_filt.fastq"
    params:
        min_length = 1440,
        max_length = 3500,
        error_rate = 0.01 
    script:
        "src/filter.jl"

rule demux:
    input:
        "fastq/{dataset}_filt.fastq"
    output:
        directory("demux/{dataset}")
    params:
        config = lambda wc: config[wc.dataset]
    script:
        "src/demux.jl"

rule porpid:
    input:
        "demux/{dataset}"
    output:
        "tagged/{dataset}/{template}/{template}_UMI1_family_tags.csv",
        directory("tagged/{dataset}/{template}/UMI1/"),
        directory("tagged/{dataset}/{template}/UMI1_keeping/"),
        directory("tagged/{dataset}/{template}/UMI2/"),
        directory("tagged/{dataset}/{template}/dUMI/"),
        "tagged/{dataset}/{template}/{template}_dUMI_ranked.csv"   
    params:
        dir = "tagged/{dataset}",
        sUMI_primer = lambda wc: config[wc.dataset][wc.template]["sUMI_Primer_Sequence"],
        dUMI_primer = lambda wc: config[wc.dataset][wc.template]["dUMI_Primer_Sequence"]  
    script:
        "src/porpid.jl"
        
rule consensus:
    input:
        "tagged/{dataset}/{template}/dUMI/",
        "tagged/{dataset}/{template}/UMI1_keeping/"
    output:
        "consensus/{dataset}/{template}_dUMI.fasta",
        "consensus/{dataset}/{template}_sUMI.fasta"
    params:
        sUMI_primer = lambda wc: config[wc.dataset][wc.template]["sUMI_Primer_Sequence"],
        dUMI_primer = lambda wc: config[wc.dataset][wc.template]["dUMI_Primer_Sequence"]       
    script:
        "src/consensus.jl"
        
onsuccess:
    shell("tar -czf consensus.tar.gz consensus"),
    shell("tar -czf tagged.tar.gz tagged")
