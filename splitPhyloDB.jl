#!/usr/bin/env julia

using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "fasta"
        help = "input fasta file to split"
        required = true
        arg_type = String
    end
    
    return parse_args(s)
end

function getFasta(fp)
    inFasta = false
    seq = ""
    for line in eachline(fp)
        if inFasta && line[1] == '>'
            produce(seq)
            seq = ""
        elseif line[1] == '>'
            inFasta = true
        end
        seq = seq * line
    end
    produce(seq)
end


function main()
    parsed_args = parse_commandline()
    fasta = parsed_args["fasta"]
    fp = open(fasta)
    ver = last(split(fasta, "_"))
    transcript = open("phylodb_transcript_" * ver, "w")
    kegg = open("phylodb_kegg_" * ver, "w")
    mito = open("phylodb_mito_" * ver, "w")
    chloro = open("phylodb_chloro_" * ver, "w")
    rest = open("phylodb_rest_" * ver, "w")
    fastaReader = Task(()->getFasta(fp))
    for seq in fastaReader
        header = first(split(seq, '\n', 2))
        flags = last(split(header, ' '))
        if ismatch(r"kegg", flags)
            print(kegg, seq)
        elseif ismatch(r"transcript", flags)
            print(transcript, seq)
        elseif ismatch(r"mito", flags)
            print(mito, seq)
        elseif ismatch(r"chloro", flags)
            print(chloro, seq)
        else
            print(rest, seq)
        end
    end
    close(fp)
    close(transcript)
    close(kegg)
    close(mito)
    close(chloro)
    close(rest)
end

main()
