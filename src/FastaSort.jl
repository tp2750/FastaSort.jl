module FastaSort

using BioSequences
using FASTX
using DataStructures

export fasta_sort, fasta_sort_file
export fasta_sort_file_0, fasta_sort_file_00, fasta_sort_file_000, fasta_sort_file_1, fasta_sort_file_2, fasta_sort_file_4
export records_split_loop_any, records_split_loop_type, records_split_map
export fasta_perm, fasta_perm_file


function split_record(record::FASTX.FASTA.Record)
    (identifier(record), sequence(record))
end

function fasta_sort_file_000(infile, outfile)
    @info "Process file $infile: $(stat(infile).size/1E9) GB"
    @time open(outfile,"w") do writer
        for l in extract_all!(BinaryMinHeap(map(x -> (identifier(x), sequence(x)), FASTAReader(open(infile)))))
            write(writer, ">"*string(l[1])*"\n"*string(l[2])*"\n")
        end
    end
end
function fasta_sort_file_00(infile, outfile)
    @info "Process file $infile: $(stat(infile).size/1E9) GB"
    @time open(outfile,"w") do writer
        for l in extract_all!(BinaryMinHeap(map(split_record, FASTAReader(open(infile)))))
            write(writer, ">"*string(l[1])*"\n"*string(l[2])*"\n")
        end
    end
end

function fasta_sort_file_0(infile, outfile)
    @info "Process file $infile: $(stat(infile).size/1E9) GB"
    @time open(outfile,"w") do writer
        for l in extract_all!(BinaryMinHeap(map(split_record, collect(FASTAReader(open(infile))))))
            write(writer, ">"*string(l[1])*"\n"*string(l[2])*"\n")
        end
    end
end


function fasta_sort_file_1(infile, outfile)
    @info "Process file $infile: $(stat(infile).size/1E9) GB"
    @time sorted = extract_all!(BinaryMinHeap(map(split_record, collect(FASTAReader(open(infile))))))
    @info "Write output: $(length(sorted)) records."
    @time open(outfile,"w") do writer
        for l in sorted
            write(writer, ">"*string(l[1])*"\n"*string(l[2])*"\n")
        end
    end
end


function fasta_sort_file_4(infile, outfile)
    @info "Read file $infile: $(stat(infile).size/1E9) GB"
    @time records = collect(FASTAReader(open(infile)))
    @info "Spilt"
    @time splitted = map(split_record, records)
    @info "Build Heap"
    @time h = BinaryMinHeap(splitted)
    @info "Sort"
    @time sorted = extract_all!(h)
    @info "Write output: $(length(sorted)) records."
    @time open(outfile,"w") do writer
        for l in sorted
            write(writer, ">"*string(l[1])*"\n"*string(l[2])*"\n")
        end
    end
end


function fasta_sort_file_2(infile, outfile)
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


function records_split_loop_any(records::Vector{FASTX.FASTA.Record})
    data = []
    for record in records
        push!(data, (identifier(record), sequence(record)))
    end
    data
end

function records_split_loop_type(records::Vector{FASTX.FASTA.Record})
    t1 = (identifier(records[1]), sequence(records[1]))
    data = typeof(t1)[]
    for record in records
        push!(data, (identifier(record), sequence(record)))
    end
    data
end


function records_split_map(records::Vector{FASTX.FASTA.Record})
    map(split_record, records)
end

function fasta_sort(records::Vector{FASTX.FASTA.Record})
    data = records_split_map(records)
    h = BinaryMinHeap(data)
    s = extract_all!(h)
    return(s)
end

function fasta_perm(records::Vector{FASTX.FASTA.Record})
    ids  = identifier.(records)
    sorter = sortperm(ids)
    records[sorter]
end

function fasta_perm_file(infile, outfile)
    @info "Process file $infile: $(stat(infile).size/1E9) GB"
    records = FASTAReader(open(infile))
    sorter = sortperm(identifier.(records))
    @time open(outfile,"w") do writer
        for l in records[sorter]
            write(writer, l)
        end
    end
end


end # module FastaSort
