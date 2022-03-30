include("DDForDP.jl")

benchmarkSOPs(1024, "outputFiles/bnb_1024_2.txt",3600,numStart=22,numEnd=41,usePeel=false)
