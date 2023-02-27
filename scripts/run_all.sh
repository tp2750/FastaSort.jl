# run from FastaSort.jl
# bash scripts/runall.sh > scripts/log_2023-02-27.out

/usr/bin/time -v  /opt/Julia/julia-1.8.2/bin/julia --project=. -O 3 scripts/fs_0.jl  data/hep25m.fasta /tmp/out25m0.fas
/usr/bin/time -v  /opt/Julia/julia-1.8.2/bin/julia --project=. -O 3 scripts/fs_1.jl  data/hep25m.fasta /tmp/out25m1.fas
/usr/bin/time -v  /opt/Julia/julia-1.8.2/bin/julia --project=. -O 3 scripts/fs_2.jl  data/hep25m.fasta /tmp/out25m2.fas
/usr/bin/time -v  /opt/Julia/julia-1.8.2/bin/julia --project=. -O 3 scripts/fs_4.jl  data/hep25m.fasta /tmp/out25m4.fas

/usr/bin/time -v /z/linux/bin/seqkit sort --by-name  --two-pass  data/hep25m.fasta > /tmp/s25_2.fas
/usr/bin/time -v /z/linux/bin/seqkit sort --by-name  data/hep25m.fasta > /tmp/s25_1.fas
