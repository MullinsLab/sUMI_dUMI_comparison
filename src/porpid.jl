using Distributed, DataFrames, DataFramesMeta
@everywhere ENV["MPLBACKEND"] = "Agg"
include("../../src/functions.jl")

sUMI_primer = snakemake.params["sUMI_primer"] #fill in here.
dUMI_primer = snakemake.params["dUMI_primer"]
filtered_data_file = snakemake.input[1]*"/"*snakemake.wildcards["template"]*".fastq"
template_name = snakemake.wildcards["template"]

### Processing UMI 1 ###
println("Extracting UMI 1...")
t1 = time()
templates1 = Dict()
templates1["UMI1"] = replace(uppercase(sUMI_primer), "N" => "n")*"*"

cfg = Configuration()
cfg.files = [filtered_data_file]
cfg.filetype = fastq
cfg.start_inclusive = 0
cfg.end_inclusive = length(first(values(templates1))) + 3
cfg.try_reverse_complement = false #The sequences are already oriented
for (name, template) in templates1
    println("Using template $(template)")
    push!(cfg.templates, Template(name, template))
end

println("Writing UMI1 to $(snakemake.output[2])...")
dir_dict = Dict()
my_output_func(source_file_name, template, tag, output_sequence, score) = my_write_to_file_count_to_dict(dir_dict, source_file_name, template, tag, output_sequence, score, snakemake.params["dir"])
say_print_func = function(count)
    println("Processed $(count) sequences")
end
# This is the slow bit
extract_tags_from_file(cfg.files[1], cfg, my_output_func, print_every=5000, print_callback=say_print_func)
directories = collect(keys(dir_dict))


tag_dfs = []
	for (ix,template) in enumerate(templates1)
    local_template_name = basename(directories[ix])
    tag_dict = dir_dict[directories[ix]]
    delete!(tag_dict, "REJECTS")
    tag_counts = tag_dict
    tags = collect(keys(tag_dict))  
	

    #Builts matrix of conditional probabilities.
    tag_to_index, index_to_tag = tag_index_mapping(tags)
    pacbio_error_rate = 0.005
    recurse = 1
    probabilities_array = prob_observed_tags_given_reals(tag_to_index, PORPID.PacBioErrorModel(pacbio_error_rate), recurse)
    indexed_counts = index_counts(tag_counts, tag_to_index);
    #Runs the LDA inference
    most_likely_real_for_each_obs = LDA(probabilities_array, indexed_counts); 
     

    likely_real = []
    UMI = []
    bin_sizes = []
    tags = []
    probs = []
    #Filter and copy to "_keeping"
    tag_df = filterCCSFamilies(most_likely_real_for_each_obs, snakemake.output[2], index_to_tag, tag_counts, local_template_name, template[2])
    tag_df[!, :Sample] .= snakemake.wildcards["template"]
    push!(tag_dfs, tag_df)
end

sorted =  sort!(vcat(tag_dfs...), [:tags, :fs], rev = [false, true])
CSV.write(snakemake.output[1], sorted); #io step


t2 = time()
println("UMI 1 identification took $(t2-t1) seconds.")

### Processing UMI 2 ###
println("Extracting UMI 2...")
t1 = time()
templates2 = Dict()
templates2["UMI2"] = replace(dUMI_primer, "N" => "n")*"*" #need to add to config!

cfg = Configuration()
cfg.files = [filtered_data_file]
cfg.filetype = fastq
cfg.start_inclusive = 0
cfg.end_inclusive = length(first(values(templates2))) + 3
cfg.try_reverse_complement = true ## need to look at reverse complement for opposite UMI...
for (name, template) in templates2
    println("Using template $(template)")
    push!(cfg.templates, Template(name, template))
end

println("Writing UMI2 to $(snakemake.output[4])...")
dir_dict = Dict()
my_output_func(source_file_name, template, tag, output_sequence, score) = my_write_to_file_count_to_dict(dir_dict, source_file_name, template, tag, output_sequence, score, snakemake.params["dir"])
say_print_func = function(count)
    println("Processed $(count) sequences")
end

# This is the slow bit
extract_tags_from_file(cfg.files[1], cfg, my_output_func, print_every=5000, print_callback=say_print_func)

### Merging UMI calls ###

dir_UMI1 = "$(snakemake.output[3])";
dir_UMI2 = "$(snakemake.output[4])";

joined = join_UMIs_by_rname(dir_UMI1, dir_UMI2);
ranked = get_ranked_families(joined);
CSV.write("$(snakemake.output[6])",ranked);

write_duplex_output(dir_UMI1, ranked[ranked.fs .> 4, :], joined; outdir = snakemake.output[5])
