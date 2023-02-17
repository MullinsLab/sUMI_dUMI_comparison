using NextGenSeqUtils, StatsBase

#load required packages
function unique_not_substr(a)
    out = []
    for i in unique(a)
        res = true
        for j in unique(a)
            if occursin(i, j) & (i != j)
                res = false
            end
        end
        if res
            push!(out, i)
        end
    end
    return out
end

function iterative_primer_match(seqs,full_primers,window::Int,slide_by::Int;tol_one_error=true)
    if(slide_by + window - 1 > minimum(length.(full_primers)))
        @warn ("Matching window extends beyond shortest primer. This is ok, but check that you aren't matching something too short.")
    end
    primers = [p[1:min(window,minimum(length.(full_primers)))] for p in full_primers]
    filter = fast_primer_match(seqs,primers,tol_one_error=tol_one_error);
    for i in 2:slide_by
        unresolved = filter .== 0
        primers = [p[i:min(i + window - 1,minimum(length.(full_primers)))] for p in full_primers]
        filter[unresolved] = fast_primer_match(seqs[unresolved],primers,tol_one_error=tol_one_error);
    end
    return filter
end

function sliding_demux_dict(seqs,fwd_primers,window::Int,slide_by::Int; verbose = true, phreds = nothing, tol_one_error = true)
    fwd_matches = iterative_primer_match(seqs,fwd_primers,window,slide_by,tol_one_error=tol_one_error)
    rev_comp_bool = fwd_matches .< 0
    keepers = abs.(fwd_matches) .> 0
    fwd_matches = abs.(fwd_matches)
    pair_keeps = fwd_matches[keepers]
    pair_counts = countmap(pair_keeps)
    sorted_pairs = sort([(k,pair_counts[k]) for k in keys(pair_counts)])
    if verbose
        for s in sorted_pairs
            println(s[1], " => ", s[2])
        end
    end
    if phreds == nothing
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                d_key = fwd_matches[i]

                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],i))
                end
            end
        end
        return seq_dict
    else
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Vector{Int8},Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                d_key = fwd_matches[i]
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),reverse(phreds[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],phreds[i],i))
                end
            end
        end
        return seq_dict
    end
end

function longest_conserved_5p(seqs)
    for i in 1:length(seqs[1])
        if length(unique(getindex.(seqs,i))) != 1
            return seqs[1][1:i-1]
        end
    end
    return seqs[1]
end

SAMPLES = snakemake.params["config"]
mkdir(snakemake.output[1])

#load sequences
println("Reading filtered sequences...")
seqs,phreds,seq_names = read_fastq(snakemake.input[1])

t1 = time()

#forward
fwd_ends = [v["Forward_Primer_2ndRd_Sequence"] for (k,v) in SAMPLES]
unique_fwd_ends = unique_not_substr(fwd_ends)
fwd_end_group_arr = []
for e in fwd_ends
    matches = []
    for (j, group) in enumerate(unique_fwd_ends)
        if occursin(e, group)
            push!(matches, (group => j))
        end
    end
    push!(fwd_end_group_arr, sort(matches, by = x -> length(x[1]))[1][2])
end

println("Splitting by primers...")
rev_adapters = [v["Reverse_Primer_2ndRd_Sequence"] for (k,v) in SAMPLES]
rev_adapter = longest_conserved_5p(rev_adapters)
rev_adapter = String(split(rev_adapter,r"[a-z]")[1]) #if all contain sample ID keep sample ID

#sliding window demultiplex on forward primers
fwd_demux_dic = sliding_demux_dict(seqs,
                                   unique_fwd_ends,
                                   12,
                                   6, #must not extend into the sampleID
                                   verbose=false,
                                   phreds = phreds)
#iterate by each forward primer group
for j in unique(fwd_end_group_arr)
    #define templates
    template_names = [k for (k,v) in SAMPLES][fwd_end_group_arr .== j]
    templates = [SAMPLES[n]["sUMI_Primer_Sequence"] for n in template_names]
    sampleIDs = uppercase.([m.match for m in match.(r"[a-z]+", templates)])
    IDind2name = Dict(zip(collect(1:length(sampleIDs)),template_names));

    #retrieve from demux_dic
    seqs_fwd = [i[1] for i in fwd_demux_dic[j]];
    phreds_fwd = [i[2] for i in fwd_demux_dic[j]];
    seq_names_fwd = seq_names[[i[3] for i in fwd_demux_dic[j]]]

    println("$(length(seqs_fwd)) reads matching forward primer $(unique_fwd_ends[j])")

    #match to reverse adapter
    rev_matches = iterative_primer_match(seqs_fwd, [rev_adapter], 12, 6, tol_one_error=true);
    rev_keepers = rev_matches .< 0

    #filter to reverse adapter matches
    seqs_both = seqs_fwd[rev_keepers]
    phreds_both = phreds_fwd[rev_keepers]
    seq_names_both = seq_names_fwd[rev_keepers]
    println("$(length(seqs_both)) contain reverse adapter $(rev_adapter)")

    #trim and orient primers, phreds
    println("Trimming primers for read group $(j)...")
    trimmed = [double_primer_trim(seqs_both[i],
                                        phreds_both[i],
                                        uppercase(unique_fwd_ends[j]),
                                        uppercase(rev_adapter)) for i in 1:length(seqs_both)]

    #separate by donor (sample) IDs
    println("Splitting read group $(j) by donor ID...")
    split_donor_dic = demux_dict([s[1] for s in trimmed],
                                       sampleIDs,
                                       nothing,
                                       verbose=false,
                                       phreds=[s[2] for s in trimmed],
                                       tol_one_error=false)

    println("Writing individual donor files...")
    for i in 1:length(sampleIDs)
        if i in keys(split_donor_dic)
            println(IDind2name[i]," => ",length(split_donor_dic[i]))
            template = template_names[i]
            write_fastq(snakemake.output[1]*"/"*template*".fastq",
                        [i[1] for i in split_donor_dic[i]],
                        [i[2] for i in split_donor_dic[i]];
                        names = seq_names_both[[i[3] for i in split_donor_dic[i]]])
        else
            println(IDind2name[i]," => 0")
        end
    end
end

t2 = time()
println("Demultiplex took $(t2 - t1) seconds.")

#remove filtered fastq file after finished
dataset = snakemake.wildcards["dataset"]
run(`rm fastq/$(dataset)_filt.fastq`)
