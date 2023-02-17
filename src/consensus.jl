using Distributed, DataFrames, DataFramesMeta
@everywhere ENV["MPLBACKEND"] = "Agg"
include("../../src/functions.jl")

sUMI_primer = snakemake.params["sUMI_primer"] 
dUMI_primer = snakemake.params["dUMI_primer"]
template_name = snakemake.wildcards["template"]

### Calculate sUMI consensus sequences for each sUMI family ###
println("Generating sUMI consensus...")
t1 = time()
println("Processing $(template_name)")
base_dir = "$(snakemake.input[2])"
@time seq_collection, seqname_collection = sUMI_generateConsensusFromDir(base_dir, template_name)
write_fasta(snakemake.output[2],
    reverse_complement.(seq_collection),
    names = seqname_collection)
t2 = time()
println("Consensus generation took $(t2-t1) seconds.")

### Calculate dUMI consensus sequences for each dUMI family ###
println("Generating dUMI consensus...")
t1 = time()
println("Processing $(template_name)")
base_dir = "$(snakemake.input[1])"
@time seq_collection, seqname_collection = dUMI_generateConsensusFromDir(base_dir, template_name)
write_fasta(snakemake.output[1],
    reverse_complement.(seq_collection),
    names = seqname_collection)
t2 = time()
println("Consensus generation took $(t2-t1) seconds.")

