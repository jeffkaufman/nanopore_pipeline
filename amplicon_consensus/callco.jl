#!/usr/bin/env julia
using BioAlignments
using FASTX
using BioSequences
#using Plots
using Random
#using StatsPlots
#using Statistics
using DelimitedFiles
using ArgParse
refseq=dna""


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--b1"
            arg_type = Int
            help = "size selection lower limit"
            default=100
        "--b2"
            help = "size selection upper limit"
            arg_type = Int
            default = 10000
	"--s1"
	    help="first sample #"
	    arg_type=Int
            default=1
	"--s2"
	    help="last sample #"
	    arg_type=Int
	    default=0
        "--out"
            help = "output file name (e.g. consensus_expid_requester.fasta)"
	    arg_type=String
	    default=""
	"--id"
            help="file containing user-provided sample names"
	    arg_type=String
	    default="./samples.csv"
        "path"
            help = "path of demultiplexed folders"
            required = true
    end

    return parse_args(s)
end

#%%

#%%
pargs = parse_commandline()
function getReads(path)
    files=readdir(path,join=true)
    allReads1=[]
    for file in files
        @show file
        @show reader=open(FASTQ.Reader,file)
        push!(allReads1,collect(reader))
        close(reader)
    end
    allReads1=reduce(vcat,allReads1)
    lenDist1=length.(sequence.(allReads1))
    frag1selidx1=pargs["b1"] .<lenDist1.<pargs["b2"]
    frag1reads1=allReads1[frag1selidx1]
    shuffle!(frag1reads1)
    if length(frag1reads1)>1000
	frag1reads1=frag1reads1[1:1000]
    end
    return  frag1reads1
end

function toFastaRec(fastqRec)
    return FASTA.Record(FASTQ.identifier(fastqRec),FASTQ.description(fastqRec),FASTQ.sequence(fastqRec))
end

function getSampleIds(sampleIdfile)
return readdlm(sampleIdfile,',')
end


function callconsensus(readsubset,ref)
    open(FASTA.Writer,"/dev/shm/refSeq.fa") do io
        write(io,FASTA.Record("ref","ref",ref))
    end
    open(FASTA.Writer,"/dev/shm/concat_seg_tmp.fa") do io
        for frag1read in readsubset
            write(io,toFastaRec(frag1read))
        end
    end
    blasrCmd=`blasr /dev/shm/concat_seg_tmp.fa /dev/shm/refSeq.fa -m 5 --affineOpen 10 --nCandidates 15 --nproc 16`
    pbdagconCmd=`pbdagcon /dev/stdin -m 100 -c 3 -j 16 -t 0`
    blasrCmd1=`blasr /dev/shm/concat_seg_tmp.fa /dev/shm/pbdagcontemp1.fa -m 5 --affineOpen 10 --nCandidates 15 --nproc 16`
    blasrCmd2=`blasr /dev/shm/concat_seg_tmp.fa /dev/shm/pbdagcontemp2.fa -m 5 --affineOpen 10 --nCandidates 15 --nproc 16`
    blasrCmd3=`blasr /dev/shm/concat_seg_tmp.fa /dev/shm/pbdagcontemp3.fa -m 5 --affineOpen 10 --nCandidates 15 --nproc 16`
    run(pipeline(blasrCmd,pbdagconCmd,"/dev/shm/pbdagcontemp1.fa"))
    run(pipeline(blasrCmd1,pbdagconCmd,"/dev/shm/pbdagcontemp2.fa"))
    #run(pipeline(blasrCmd2,pbdagconCmd,"/dev/shm/pbdagcontemp3.fa"))
    consesus=split(readchomp(pipeline(blasrCmd2,pbdagconCmd)),'\n')
    @show consesus=parse(LongSequence{DNAAlphabet{4}},consesus[2])
    return consesus
end

#%%
ref=dna""
conList=[]
lastSample=0
if pargs["s2"]==0
    lastSample=size(getSampleIds(pargs["id"]))[1]
    print(lastSample)
else
   lastSample=pargs["s2"]
end
for i in range(pargs["s1"],lastSample,step=1)
    fastqPath="bw_"*lpad(i-1,3,'0')
    fastqPath=(pargs["path"]*"/")*fastqPath
    reads=""
    try
    	reads=getReads(fastqPath)
    catch
    	#reads=getReads(fastqPath)
        push!(conList,"CCCC")
	continue
	end
    con=""
    try
    	con=callconsensus(reads,reads[10])
    catch
	try
    	con=callconsensus(reads,sequence(reads[15]))
	catch
	try
    	con=callconsensus(reads,sequence(reads[51]))
	catch
	con="AAAA"
	end
	end
    end
    push!(conList,con)
end
#%%
open(FASTA.Writer,pargs["out"],width=1000000) do io
    sampleIds=getSampleIds(pargs["id"])
    for (i,read) in enumerate(conList)
        write(io,FASTA.Record("$(sampleIds[pargs["s1"]+i-1])","$(pargs["s1"]+i-1)",read))
    end
end
#%%
