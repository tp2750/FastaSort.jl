using Revise
using FastaSort
using BioSequences
using FASTX
using DataStructures
using BenchmarkTools

const medium_file = "data/hep1m.fasta"   ##   164MB. 500_000 records. copy 2.5 s
@time const medium  = collect(FASTAReader(open(medium_file)));

@time map = records_split_map(medium);             # 0.033937 seconds (2 allocations: 38.147 MiB)
@time perm = fasta_perm(medium); ## 0.595614 seconds (7 allocations: 26.703 MiB) So not better!

@time fasta_sort_file_0(medium_file, "/tmp/m1.fas") ## 2.5 sec
@time fasta_sort_file_00(medium_file, "/tmp/m2.fas") ## 2.2 sec
@time fasta_sort_file_000(medium_file, "/tmp/m3.fas") ## 2.2 sec

@btime fasta_sort_file_0(medium_file, "/tmp/m1.fas")   # 2.192 s (7523157 allocations: 1.16 GiB)
@btime fasta_sort_file_00(medium_file, "/tmp/m2.fas")  # 2.155 s (7523156 allocations: 1.16 GiB)
@btime fasta_sort_file_000(medium_file, "/tmp/m3.fas") # 2.117 s (7523156 allocations: 1.16 GiB)
