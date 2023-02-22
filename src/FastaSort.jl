module FastaSort

using BioSequences
using FASTX
using DataStructures

export fasta_sort, fasta_sort_file

function fasta_sort_file(infile, outfile)    
    @info "Read file $infile: $(stat(infile).size/1E9) GB"
    @time records = collect(FASTAReader(open(infile)))
    @info "Sort"
    @time sorted = fasta_sort(records)
    @info "Write output: $(length(sorted)) records."
    @time open(outfile,"w") do writer
        for l in sorted
            write(writer, ">"*string(l[1])*"\n"*string(l[2])*"\n")
        end
    end

end

function fasta_sort(records::Vector{FASTX.FASTA.Record})
#    t1 = (identifier(records[1]), sequence(records[1]))
#    data = typeof(t1)[]
    data = []
#    @info "Split records"
    time_split_records = @timed \
    for record in records
        push!(data, (identifier(record), sequence(record)))
    end
#    println(time_split_records)
#    @info "Processed $(length(data)) records in $(time_read_file.time) s using $(time_read_file.bytes/1E6) GB"
#    @info "Build Heap"
    time_build_heap = @timed h = BinaryMinHeap(data)
#    println(time_build_heap)
    
#    @info "Sort Heap"
    time_sort_heap = @timed s = extract_all!(h)
#    println(time_sort_heap)

#    @info "return $(length(s)) records"
    return(s)
end


end # module FastaSort
