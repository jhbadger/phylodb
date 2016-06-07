#!/usr/bin/env julia

type Taxon
    taxon_id::Int
    name::ASCIIString
    parent_id::Int
    rank::ASCIIString
end

function loadNCBI(filename)
    println(STDERR, "Loading ncbi taxonomy...")
    ncbi = Dict{ASCIIString, Taxon}()
    fp = open(filename)
    for line in eachline(fp)
        tid, nm, parent, rk = split(chomp(line), '\t')
        tid = parse(Int, tid)
        parent = parse(Int, parent)
        ncbi[nm] = Taxon(tid, nm, parent, rk)
    end
    close(fp)
    return ncbi
end

function loadExisting(filename)
    println(STDERR, "Loading existing taxonomy...")
    existing = Dict{ASCIIString, Int}()
    fp = open(filename)
    for line in eachline(fp)
        tid, nm, parent, rk = split(chomp(line), '\t')
        tid = parse(Int, tid)
        existing[nm] = tid
    end
    close(fp)
    return existing
end

function findMinNum(existing, taxstrings)
    minNum = 1000000000000
    fp = open(taxstrings)
    for line in eachline(fp)
        num, rest = split(chomp(line), '\t', limit = 2)
        num = parse(Int, num)
        if num < minNum
            minNum = num
        end
    end
    close(fp)
    for key in keys(existing)
        if existing[key] < minNum
            minNum = existing[key]
        end
    end
    return minNum
end

function findTaxon(hierarchy, nameLookup, name, ncbi, existing, rname, minNum)
    if haskey(nameLookup, name)
        taxon = hierarchy[nameLookup[name]]
    elseif haskey(ncbi, name)
        taxon = Taxon(ncbi[name].taxon_id, name, 0, rname)
    elseif haskey(existing, name)
        taxon = Taxon(existing[name], name, 0, rname)
    else
        minNum = minNum - 1
        taxon = Taxon(minNum, name, 0, rname)
    end
    return taxon, minNum
end

function setParent(hierarchy, nameLookup, taxon, parent)
    if taxon.taxon_id == parent.taxon_id
        println(STDERR, "Error: Recursive taxonomy for ", taxon.taxon_id,
                " ", taxon.name)
    end
    taxon.parent_id = parent.taxon_id
    if haskey(nameLookup, taxon.name) && taxon.taxon_id != nameLookup[taxon.name]
        println(STDERR, "Error: Duplicate name for: ", taxon.taxon_id, 
                " ", taxon.name)
    else
        hierarchy[taxon.taxon_id] = taxon
        nameLookup[taxon.name] = taxon.taxon_id
    end
    return taxon
end

function buildHierarchy(filename, ncbi, existing, minNum)
    hierarchy = Dict{Int, Taxon}()
    nameLookup = Dict{ASCIIString, Int}()
    hierarchy[1] = Taxon(1, "root", 0, "")
    ranks = ["kingdom","phylum","class","order","family",
             "genus", "species", "subspecies"]
    fp = open(filename)
    for line in eachline(fp)
        tid, tx = split(chomp(line), '\t')
        tid = parse(Int, tid)
        tax = split(tx,';')
        parent = hierarchy[1]
        tlen = length(tax)
        if tlen > 8
            tlen = 8
        end
        for i in 1:tlen
            taxon, minNum = findTaxon(hierarchy, nameLookup, tax[i], 
                                      ncbi, existing, ranks[i], minNum)
            if i == tlen
                taxon.taxon_id = tid
            end
            taxon = setParent(hierarchy, nameLookup, taxon, parent)
            parent = taxon
        end
    end
    close(fp)
    return hierarchy
end

function printHierarchy(hierarchy)
    for key in sort(collect(keys(hierarchy)))
        println(hierarchy[key].taxon_id, '\t', hierarchy[key].name, '\t',
                hierarchy[key].parent_id, '\t', hierarchy[key].rank)
    end
end

if (length(ARGS) < 3)
    println(STDERR, "Usage: buildTaxonomyHierachy.jl taxstrings.txt ncbi_taxonomy.txt existing.txt")
else
    ncbi = loadNCBI(ARGS[2])
    existing = loadExisting(ARGS[3])
    minNum = findMinNum(existing, ARGS[1])
    hierarchy = buildHierarchy(ARGS[1], ncbi, existing, minNum)
    printHierarchy(hierarchy)
end
