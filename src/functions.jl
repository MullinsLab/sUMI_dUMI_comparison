using Distributed
@everywhere ENV["MPLBACKEND"] = "Agg"
using PORPID, PyPlot, StatsBase, HypothesisTests,
    IterTools, DataFrames, CSV, DataFramesMeta, BioSequences, FASTX
@everywhere using RobustAmpliconDenoising, NextGenSeqUtils


function read_fasta_with_everything(filename; seqtype=String)
    records = NextGenSeqUtils.read_fasta_records(filename)
    return seqtype[FASTA.sequence(seqtype, r) for r in records], [FASTA.identifier(r) for r in records], [FASTA.description(r) for r in records]
end

function read_fasta_with_descriptors_in_names(f)
   seqs,seqnames,descriptors = read_fasta_with_everything(f);
   seqnames = seqnames .* " " .* descriptors
   return seqnames,seqs
end

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

"""
Given results of the LDA, filter likely offspring UMIs and apply qc
"""
function filterCCSFamilies(most_likely_real_for_each_obs, path, index_to_tag, tag_counts, template_name, template)
    likely_real = []
    UMI = []
    bin_sizes = []
    tags = []
    probs = []
    
    if path[end] != '/'
        path = "$(path)/"
    end	

    for (observed_index_local, tuple_local) in enumerate(most_likely_real_for_each_obs)
        tagSeq = index_to_tag[observed_index_local]
        prob = tuple_local[2]
 
        push!(UMI, tagSeq)
        push!(bin_sizes, tag_counts[tagSeq])
        push!(probs, log(1-prob))

        #Testing for whether there is a drop in QC scores in the barcode region
        #This filters out sequencing errors, heteroduplexes	        
        read_seqs, read_phreds, read_names = read_fastq(path*tagSeq*".fastq");       
        averages = mean([Array{Float64}(read_phreds[i][1:50]) for i in 1:length(read_phreds)])
        umi_ix = findfirst(r"n+", template)
        p_val = 1
        if mean(averages[umi_ix])<mean(averages[umi_ix[end]+5:umi_ix[end]+25])
            p_val=pvalue(EqualVarianceTTest(averages[umi_ix],averages[umi_ix[end]+5:umi_ix[end]+25]))
        end
 
        if p_val < 1/(5*tag_counts[tagSeq]^2) && minimum(averages[umi_ix]) < mean(averages[umi_ix[end]+5:umi_ix[end]+25])/2
            push!(tags, "heteroduplex")
            continue
        end

        if prob < .995
            push!(tags, "LDA-rejects")
            continue
        end

        #Filter out if not length 8
        if length(tagSeq) != 8
            push!(tags, "UMI_len != 8")
            continue
        end

        thresh = 5
        if tag_counts[tagSeq] < thresh # set fs thresh here
            push!(tags, "fs<$(thresh)")
            continue
        end

        push!(tags, "likely_real")
        push!(likely_real, (prob, tagSeq))
    end

    #sort!(likely_real, by=x->tag_counts[x[2]])
    
    directory = "$(path[1:end-1])_keeping/"
    if isdir(directory) == false
    	mkdir(directory)
    end	

    for tag_tuple in likely_real
        run(`cp $(path)$(tag_tuple[2]).fastq $(directory)$(tag_tuple[2])_$(length(tag_tuple[2]))_$(round(tag_tuple[1],digits = 5))_$(tag_counts[tag_tuple[2]]).fastq`)
    end

    tag_df = DataFrame(Sample=template_name, UMI = UMI, fs = bin_sizes, tags = tags, probs = probs)
    return tag_df
end

function dUMI_generateConsensusFromDir(dir, template_name)#, fwd_ix, rev_ix)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq"]
    if length(files) > 0
        println("Generating consensus for $(length(files)) templates")
    else
        println("WARNING: no template families for $(template_name)")
        return [], []
        exit()
    end
    cons_collection = pmap(dUMI_ConsensusFromFastq, files)
    #println("$(cons_collection)")
    seq_collection = [i[1] for i in cons_collection]
    seqname_collection = [template_name*i[2] for i in cons_collection]
    #seqname_collection = [template_name*i[2]*" Index=$(fwd_ix)_$(rev_ix)" for i in cons_collection]
    return seq_collection, seqname_collection
end

function sUMI_generateConsensusFromDir(dir, template_name)#, fwd_ix, rev_ix)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq"]
    if length(files) > 0
        println("Generating consensus for $(length(files)) templates")
    else
        println("WARNING: no template families for $(template_name)")
        return [], []
        exit()
    end
    cons_collection = pmap(sUMI_ConsensusFromFastq, files)
    #println("$(cons_collection)")
    seq_collection = [i[1] for i in cons_collection]
    seqname_collection = [template_name*i[2] for i in cons_collection]
    #seqname_collection = [template_name*i[2]*" Index=$(fwd_ix)_$(rev_ix)" for i in cons_collection]
    return seq_collection, seqname_collection
end

function generateConsensusFromDirSingleThread(dir, template_name)#, fwd_ix, rev_ix)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq"]
    if length(files) > 0
        println("Generating consensus for $(length(files)) templates")
    else
        println("WARNING: no template families for $(template_name)")
        exit()
    end
    cons_collection = [ConsensusFromFastq(f) for f in files]
    seq_collection = [i[1] for i in cons_collection]
    seqname_collection = [template_name*i[2] for i in cons_collection]
    #seqname_collection = [template_name*i[2]*" Index=$(fwd_ix)_$(rev_ix)" for i in cons_collection]
    return seq_collection, seqname_collection
end

@everywhere function dUMI_ConsensusFromFastq(file)
    seqs,phreds,seq_names = read_fastq(file)
    draft = consensus_seq(seqs)
    draft2 = refine_ref(draft, seqs)
    final_cons = refine_ref(draft2,seqs)
    trimmed_collection = double_primer_trim(final_cons,uppercase(sUMI_primer), uppercase(dUMI_primer))
    alignments, maps, matches, matchContent = getReadMatches(trimmed_collection, seqs, 0)
    cons_name = split(basename(file),"_")[1]*"_"*split(basename(file),"_")[2]*" num_CCS=$(length(seqs)) min_agreement=$(round(minimum(matches); digits = 2))"
    return trimmed_collection, cons_name
end

@everywhere function sUMI_ConsensusFromFastq(file)
    seqs,phreds,seq_names = read_fastq(file)
    draft = consensus_seq(seqs)
    draft2 = refine_ref(draft, seqs)
    final_cons = refine_ref(draft2,seqs)
    trimmed_collection = double_primer_trim(final_cons,uppercase(sUMI_primer), uppercase(dUMI_primer))
    alignments, maps, matches, matchContent = getReadMatches(trimmed_collection, seqs, 0)
    cons_name = split(basename(file),"_")[1]*" num_CCS=$(length(seqs)) min_agreement=$(round(minimum(matches); digits = 2))"
    return trimmed_collection, cons_name 
end

@everywhere begin
"""
Returns an array of degapped coordinates, such that
coords(ref, read)[i] gives you the position the aligned read/ref
that matches the i'th ungapped position in ref.
"""
function coords(ref, read)
    if length(ref) != length(read)
        error("Aligned strings are meant to be the same length.")
    end
    degappedRef = degap(ref)
    coordMap = zeros(Int64, length(degappedRef))
    count = 1
    for i in 1:length(degappedRef)
        while ref[count] == '-'
            count += 1
        end
        coordMap[i] = count
        count += 1
    end
    return coordMap
end
end

@everywhere begin
"""
Return matches to a candidate reference from a set of reads.
"""
function getReadMatches(candidate_ref, reads, shift; degap_param = true, kmer_align = true)
    alignments = []
    if kmer_align
        alignments = map(i -> kmer_seeded_align(candidate_ref, i), reads)
    else
        alignments = map(i -> nw_align(candidate_ref, i), reads)
    end

    maps = [coords(i...) for i in alignments]

    if (degap_param)
        matchContent = [[degap(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
        matches = [freq(matchContent[k], degap(candidate_ref[k:k+shift])) for k in 1:length(matchContent)]
    else
       matchContent = [[(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
       matches = [freq(matchContent[k], candidate_ref[k:k+shift]) for k in 1:length(matchContent)]
    end
    return alignments, maps, matches, matchContent
end
end

"""
Returns UMIs, readnames and their corresponding files from a directory of .fastq files
"""

function collect_UMIs_names(dir)

    UMIs = []
    names = []
    files_out = []

    dir = strip(dir, '/')
    files = [f for f in readdir(dir) if f[end-5:end] == ".fastq"]

    for file in files
        if file == "REJECTS.fastq"
            continue
        end

        UMI = match(r"^.+?(?=_|.fastq)", file).match # very flexible -- be more stringent to only select 8

        open(dir*"/"*file) do f
            for (i,ln) in enumerate(eachline(f))
                if rem(i - 1,4) == 0
                    push!(UMIs, UMI)
                    push!(names, strip(ln,'@'))
                    push!(files_out, file)
                end
            end
        end
    end
    return UMIs, names, files_out
end

"""
Matches UMIs from the same read for dUMI
"""

function join_UMIs_by_rname(dir_f, dir_r)
	
    fwd_cols = collect_UMIs_names(dir_f)
    rev_cols = collect_UMIs_names(dir_r)

    fwd_df = DataFrame(UMI_1 = fwd_cols[1], name = fwd_cols[2], file = fwd_cols[3]);
    rev_df = DataFrame(UMI_2 = rev_cols[1], name = rev_cols[2]);

    joined = innerjoin(fwd_df, rev_df, on = :name)

    return joined
end

"""
Ranks UMI read families for dUMI processing
"""

function get_ranked_families(joined)
    if size(joined,1) == 0
        @warn "No dUMI families!"
        return DataFrame([String, String, Int, Int, Float64, String], [:UMI_1, :UMI_2, :fs, :rank, :prop, :names])
    else
        gdf = DataFrames.groupby(joined, [:UMI_1, :UMI_2]);
        duplex_fs = combine(x-> DataFrame(fs = size(x)[1], names = join(x[:,:name],";")), gdf);
        #duplex_fs = by(joined, [:UMI_1, :UMI_2], x-> DataFrame(fs = size(x)[1], names = join(x[:,:name],";")))
        sort!(duplex_fs, :fs, rev=true)

        # select largest family size for each tag_1 ID
        ranked = combine(DataFrames.groupby(duplex_fs, :UMI_1), x-> DataFrame( 
        								   UMI_2 = x[:,:UMI_2],
                                           fs=x[:,:fs],
                                           rank=range(1; stop=length(x[:,:UMI_2])),
                                           prop=x[:,:fs]/sum(x[:,:fs]),
                                           names = x[:,:names]))                                                              
        return ranked
    end
end

"""
Write duplex output to new files
"""

function write_duplex_output(dir, duplex_out, joined; outdir = ".")
	
    dir = strip(dir,'/')
    outdir = strip(outdir,'/')
     
    if isdir(outdir) == false
    	mkdir(outdir)
    end	

    for x in eachrow(duplex_out)

        out_seqs = []
        out_phreds = Array{Int8,1}[]
        out_names = []

        fh = "$(x[:UMI_1])_$(x[:UMI_2])_$(x[:fs])_$(round(x[:prop]; digits=2)).fastq"

        UMI_1 = x[:UMI_1]
        UMI_2 = x[:UMI_2]
		
        family = @linq joined |> where(:UMI_1 .== UMI_1, :UMI_2 .== UMI_2)
        
        seqs, phreds, names = read_fastq(dir*"/"*family[1,:file])

        in_family = [name in family[:,:name] for name in names]

        for (i,bool) in enumerate(in_family)
            if bool
                push!(out_seqs, seqs[i])
                push!(out_phreds, phreds[i])
                push!(out_names, names[i])
            end
        end
        write_fastq(outdir*"/"*fh, out_seqs, out_phreds; names = out_names)
    end
end

### writing PORPID output
function my_write_to_file(source_file_name, template, tag, output_sequence, score, outdir)
  source_file_name = basename(source_file_name)
  sample = split(basename(source_file_name),".fast")[1]
  output_file_name = "$(outdir)/$(sample)/$(template.name)/$(tag).fastq"
  mkpath("$(outdir)/$(sample)/$(template.name)")
  fo = open(output_file_name, "a")
  writer = FASTQ.Writer(fo)
  write(writer, output_sequence)
  close(writer)
end

function my_write_to_dictionary(dictionary, source_file_name, template, tag, output_sequence, score)
  directory = "$(source_file_name)/$(template.name)"
  if !haskey(dictionary, directory)
    dictionary[directory] = Dict()
  end
  directory_dict = dictionary[directory]
  if !haskey(directory_dict, tag)
    directory_dict[tag] = []
  end
  push!(directory_dict[tag], (score, output_sequence))
end

function my_write_to_file_count_to_dict(dictionary, source_file_name, template, tag, output_sequence, score, outdir)
  my_write_to_file(source_file_name, template, tag, output_sequence, score, outdir)
  directory = "$(source_file_name)/$(template.name)"
  if !haskey(dictionary, directory)
    dictionary[directory] = Dict()
  end
  directory_dict = dictionary[directory]
  directory_dict[tag] = get(directory_dict, tag, 0) + 1
end
