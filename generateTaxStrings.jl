#!/usr/bin/env julia

using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "organisms.txt"
        help = "tab file from phylodb (taxon_id, taxonomy)"
        required = true
        arg_type = String
        "--aliases", "-a"
        help = "taxonomic aliases file (from, to)"
        arg_type = String
        "--silva", "-s"
        help = "silva taxonomy file (accession, start, end, taxonomy)"
        arg_type = String
        "--pr2", "-p"
        help = "pr2 qiime taxonomy file (entry, taxonomy)"
        arg_type = String
        "--viruses", "-v"
        help = "virus taxonomy file (taxon_id, taxonomy)"
        arg_type = String
        "--output", "-o"
        help = "output to file, not stdout"
        arg_type = String
        default = "taxstrings.txt"
    end
    
    return parse_args(s)
end

function loadSilva(file)
    info = Dict()
    parents = Dict()
    if file != nothing
        println(STDERR, "Loading silva...")
        fp = open(file)
        for line in eachline(fp)
            silva, start, stop, tx, sp = split(chomp(line), "\t")
            tx *= ";$(sp)"
            tx = split(replace(tx, "Candidatus ",""), ";")
            if sort(tx) != sort(unique(tx)) # duplicated rank
                newtx = String[]
                prev = ""
                for taxon in tx
                    suffix = "_"
                    if count(x->x==taxon, taxon) == 1
                        push!(newtx, taxon)
                    elseif !in(taxon, newtx)
                        push!(newtx, taxon)
                    else
                        suffix = suffix * "X"
                        taxon = taxon * suffix
                        push!(newtx, taxon)
                    end
                    if !haskey(parents, taxon)
                        parents[taxon] = String[]
                    end
                    if !in(prev, parents[taxon])
                        push!(parents[taxon], prev)
                    end
                    prev = taxon
                end
                tx = newtx
            end
            # these ranks make things too long
            badranks = ["Myxococcaceae", "Cystobacteraceae", "Polyangiaceae", 
                        "Haliangiaceae", "Sphaerobacteraceae", 
                        "Sphaerobacterineae", "Nannocystaceae"]
            tx = filter(x->!in(x, badranks), tx)
            sp = pop!(tx)
            if !in("Chloroplast", tx) && !in("Eukaryota", tx) && 
                !haskey(info, last(tx)) && !ismatch(r"phage|virus"i, sp)
                info[last(tx)] = tx
                if !haskey(info, sp) && sp != "uncultured"
                    info[sp] = tx 
                end
            end
        end
        close(fp)
    end
    return info
end

function loadPr2(file)
    info = Dict{String, Array{String}}()
    best = Dict{String, Integer}()
    count = Dict{Array{String}, Integer}()
    if file != nothing
        println(STDERR, "Loading pr2...")
        fp = open(file)
        for line in eachline(fp)
            entry, tx = split(chomp(line), '\t', 2)
            tx = split(tx, ';')
            pop!(tx)
            genus = last(tx)
            if genus == "marine" || in("Organelle", tx)
                continue
            end
            count[tx] = get(count, tx, 0) + 1
            if !haskey(best, genus) || best[genus] < count[tx]
                best[genus] = count[tx]
            end
        end
        for (tx, num) in count
            genus = last(tx)
            if num == best[genus]
              info[genus] = tx
            end
        end
        close(fp)
    end
    return info
end

function loadTaxonIdTaxonomy(file)
    info = Dict()
    if file != nothing
        println(STDERR, "Loading ", file, "...")
        fp = open(file)
        for line in eachline(fp)
            try
                taxon_id, tx = split(chomp(line), '\t')
                taxon_id = int(taxon_id)
                tx = split(tx, ";")
                info[taxon_id] = tx
            catch
            end
        end
        close(fp)
    end
    return info
end

function loadAliases(file)
    info = Dict()
    if file != nothing
        println(STDERR, "Loading ", file, "...")
        fp = open(file)
        for line in eachline(fp)
            try
                old, new = split(chomp(line), '\t')
                info[old] = new
            catch
            end
        end
        close(fp)
    end
    return info
end

function processContigs(file, silva, pr2, viruses, custom, aliases)
    fp = open(file)
    println("Processing Contigs...")
    tax = Dict{Int, Array{String}}()
    sps = Dict{Int, String}()
    for line in eachline(fp)
        tid, tx = split(chomp(line), '\t')
        sp = last(split(tx, ";"))
        if tid == "taxon_id"
            continue
        end
        tid = int(tid)
        tax[tid] = String[]
        if ismatch(r"Eukaryota", tx)
            txinfo = pr2
        elseif ismatch(r"Viruses", tx)
            txinfo = viruses
        else
            txinfo = silva
        end
        if length(split(sp, ' ')) == 1 # MMETSP no spaces
            sp = replace(sp, '-',' ')
        end
        sp = replace(sp, "Candidatus ","")
        genus = first(split(sp, ' '))
        genus = replace(genus, Set("[]()'"), "")
        if ismatch(r"phytoplasma"i, sp) || ismatch(r"phytoplasma"i, genus)
            genus = "Phytoplasma"
        end
        if ismatch(r"uncultured"i, genus)
            genus = split(sp, ' ')[2]
        end
        if ismatch(r"Candidate Division"i, genus)
            genus = split(sp, ' ')[3]
        end
        if ismatch(r"proteobacterium"i, sp)
            tax[tid] = ["Bacteria", "Proteobacteria"]
            pclass = match(r"alpha|beta|gamma|delta|epsilon|zeta"i, sp)
            if pclass != nothing
                push!(tax[tid], ucfirst(pclass.match)*"proteobacteria")
                push!(tax[tid], 
                      "Misc. "*ucfirst(pclass.match)*"proteobacteria")
                push!(tax[tid], 
                      "Environmental "*ucfirst(pclass.match)*"proteobacteria")
                push!(tax[tid], 
                      "Uncultured "*ucfirst(pclass.match)*"proteobacteria")
            end
        end
        if haskey(aliases, genus)
            genus = aliases[genus]
        end
        if haskey(aliases, sp)
            sp = aliases[sp]
        end
        if haskey(txinfo, sp)
            tax[tid] = txinfo[sp]
        else
            if haskey(txinfo, sp) && !haskey(txinfo, genus)
                genus = sp
            end
            if haskey(txinfo, genus)
                tax[tid] = txinfo[genus]
            end
        end
        if haskey(viruses, tid)
            tax[tid] = viruses[tid]
        end
        if haskey(custom, tid)
            tax[tid] = custom[tid]
        end
        if haskey(tax, tid)
            sps[tid] = sp
        end
    end
    close(fp)
    return tax, sps
end

function printTax(file, tax, sps)
    fp = open(file, "w")
    missing = 0
    println("Writing taxstrings...")
    for tid in sort(collect(keys(tax)))
        tx = tax[tid]
        if isempty(tx)
            println("Missing Taxonomy for ", tid)
            missing = missing + 1
            continue
        end
        sp = sps[tid]
        if (in("Eukaryota", tx) && length(tx) < 8) || 
            (!in("Eukaryota", tx) && length(tx) < 7)
            if sp == last(tx)
                sp = sp * " sp."
            end
            push!(tx, sp)
        else
            tx[length(tx)] = sp
        end
        taxstring = join(tx,";")
        println(fp, tid, '\t', taxstring)        
    end
    close(fp)
    println(STDERR, length(collect(keys(tax))), " taxa mapped, ",
            missing, " umapped.")
end

function main()
    parsed_args = parse_commandline()
    silva = loadSilva(parsed_args["silva"])
    pr2 = loadPr2(parsed_args["pr2"])
    viruses = loadTaxonIdTaxonomy(parsed_args["viruses"])
    custom = loadTaxonIdTaxonomy(parsed_args["organisms.txt"])
    aliases = loadAliases(parsed_args["aliases"])
    tax, sps = processContigs(parsed_args["organisms.txt"], 
                                       silva, pr2, viruses, custom, aliases)
    printTax(parsed_args["output"], tax, sps)
end

main()
