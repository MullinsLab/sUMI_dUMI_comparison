using NextGenSeqUtils

#filter .fastq to remove reads with low quality
filtered_path = snakemake.output[1]
println("Filtering .fastq file...")
@time fastq_filter(snakemake.input[1],
                   filtered_path, #path here
                   error_rate = snakemake.params["error_rate"],
                   min_length = snakemake.params["min_length"],
                   max_length = snakemake.params["max_length"])
 