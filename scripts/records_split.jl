# Test splitting implementations
using Revise
using FastaSort
using BioSequences
using FASTX
using DataStructures


#small_file  = "data/hep.fasta"     ##      342 kB
medium_file = "data/hep1m.fasta"   ##   164MB. 500_000 records. copy 2.5 s
#large_file  = "data/hep100m.fasta" ## 17GB

# @time const small   = collect(FASTAReader(open(small_file)))
@time const medium  = collect(FASTAReader(open(medium_file)));
## 1.000142 seconds (1.51 M allocations: 221.956 MiB, 28.02% gc time, 29.64% compilation time)
# @time const large   = collect(FASTAReader(open(large_file)))
# 0.608767 seconds (1.50 M allocations: 221.479 MiB, 37.77% gc time)

## split records
@time loop_any = records_split_loop_any(medium); #   0.066569 seconds (500.01 k allocations: 51.004 MiB)
@time loop_type = records_split_loop_type(medium); # 0.072355 seconds (13 allocations: 45.523 MiB)
@time map = records_split_map(medium);             # 0.033937 seconds (2 allocations: 38.147 MiB)
