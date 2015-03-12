# calculate seguid on string
class String
   def seguid
      require 'digest/sha1'
      final = self.upcase.scan(/[A-Z]/).join
      [Digest::SHA1.digest(final)].pack("m").chomp("=\n")
   end
end

# downloads PR2 taxonomy if needed
def getPr2
  if !File.exists?("pr2_genus.txt")
    `wget http://application.sb-roscoff.fr/redgenes//static/pr2_gb191_taxonomy.xls`
    # Problem -- need manual parsing of microsoft silliness as useful text file not provided
    STDERR << "Please save genus tab as text to pr2_genus.txt\n"
    exit(1)
  end
end

def loadSilva(file)
   STDERR << "Loading Silva taxonomy...\n"
   info = Hash.new
   parents = Hash.new
   File.new(file).each do |line|
      silva, tx = line.chomp.split(" ", 2)
      tx.gsub!("Candidatus ","")
      tx = tx.split(";")
      if tx.sort != tx.uniq.sort # duplicated rank that needs to be fixed
         newtx = []
         prev = nil
         tx.each do |taxon|
            suffix = "_"
            if tx.count(taxon) == 1
               newtx.push(taxon)
            elsif !newtx.include?(taxon)
               newtx.push(taxon)
            else
               suffix += "X"
               taxon = taxon + suffix
               newtx.push(taxon)
            end
            parents[taxon] = [] if !parents[taxon]
            parents[taxon].push(prev) if prev && !parents[taxon].include?(prev)
            prev = taxon
         end
         tx = newtx
      end
      # these ranks make things too long
      tx -= ["Myxococcaceae", "Cystobacteraceae", "Polyangiaceae", "Haliangiaceae", "Sphaerobacteraceae", "Sphaerobacterineae","Nannocystaceae"]
      sp = tx.pop
      if !tx.include?("Chloroplast") && !tx.include?("Eukaryota") && !info[tx[-1]] && sp !~/phage|virus/i
         info[tx[-1]] = tx
         info[sp] = tx if !info[sp] && sp != "uncultured"
      end
   end
   info
end

def loadPr2(file)
   info = Hash.new
   count = Hash.new
   File.new(file).each do |line|
      num, tx = line.chomp.split("\t", 2)
      num = num.to_i
      tx = tx.split("\t")
      genus = tx.last
      next if genus == "marine" || tx.include?("Organelle")
      count[genus] = 0 if !count[genus]
      if num > count[genus]
        info[genus] = tx
        count[genus] = num
      end
   end
   info
end

def loadViruses(file)
   info = Hash.new
   File.new(file).each do |line|
      sp, tx = line.chomp.split("\t")
      if !tx.nil?
         tx = tx.split(";")
         sp = tx.pop
         info[sp] = tx
      end
   end
   info
end

# make kegg taxonid file
def makeKeggTaxonId(file)
   out = File.new(file, "w")
   keggtaxurl = "http://www.genome.jp/kegg-bin/get_htext?query=&htext=Organisms&filedir=&option=-e&extend=F65F548C25-2E42-3&uploadfile=&format=&wrap=&length=&open=&close=&hier=18&oneclick=on"
   out = File.new("kegg_taxon_ids.txt", "w")
   File.new(open(keggtaxurl)).each do |line|
      if line =~/([^\>]*)<\/I>.*id=([0-9]*).*>([a-z]*)</
         sp, taxid, name = $1, $2, $3
         out.printf("%s\t%d\t%s\n", name, taxid, sp)
      end
   end
   out.close
end

# return Fasta Format
def to_fasta(seq, header, len = 60)
   ">#{header}\n#{seq.gsub("*","").gsub(Regexp.new(".{1,len}"), "\\0\n")}\n"
end

# load kegg taxonomy
def makeKeggTaxonomy(keggdir)
   tax = nil
   taxonomy = Hash.new
   lastlin = false
   File.new(keggdir + "/genes/genome/genome").each do |line|
      if line =~/TAX:([0-9]*)/
         tax = $1.to_i
      elsif line =~/LINEAGE ([^\n]*)\n/
         taxonomy[tax] = $1.lstrip.rstrip
         lastlin = true
      elsif line =~/;/ && lastlin
         taxonomy[tax] += " " + line.chomp.lstrip.rstrip
         lastlin = false
      end
   end
   taxonomy
end

# get and untar ncbi taxdump if needed
def getTaxDump
   if !File.exists?("taxdump.tar.gz")
      STDERR << "Getting taxdump...\n"
      require 'net/ftp'
      Net::FTP.open("ftp.ncbi.nih.gov") do |ftp|
         ftp.login
         files = ftp.chdir("/pub/taxonomy")
         ftp.getbinaryfile("taxdump.tar.gz", "taxdump.tar.gz", 1024)
      end
   end
   if !File.exists?("names.dmp")
      `tar xvfz taxdump.tar.gz`
   end
end

# get NCBI taxonomy if needed
def get_ncbi_taxonomy(exclude=nil)
   getTaxDump if !File.exists?("names.dmp") && !File.exists?("ncbi_taxonomy.txt")
   if !File.exists?("ncbi_taxonomy.txt")
      ranks = Hash.new
      parents = Hash.new
      seen = Hash.new
      tax = File.new("ncbi_taxonomy.txt", "w")
      STDERR << "Loading nodes...\n"
      File.new("nodes.dmp").each do |line|
         num, parent, rank = line.chomp.split("\t|\t")
         num = num.to_i
         parent = parent.to_i
         ranks[num] = rank
         parents[num] = parent
      end
      File.new("names.dmp").each do |line|
         num, name, foo, type = line.chomp.split("\t|\t")
         type = type.split("\t").first
         num = num.to_i
         if type == "scientific name"
            if seen[num]
               seen[num] = 1
            end
            name.gsub!("Candidatus ","")
            tax.print [num.to_s, name, parents[num].to_s, ranks[num]].join("\t") + "\n"
         end
      end
      File.unlink("citations.dmp", "delnodes.dmp", "division.dmp", "gencode.dmp", "merged.dmp",
      "names.dmp", "nodes.dmp", "gc.prt", "readme.txt", "taxdump.tar.gz")
   end
   excluded = Hash.new
   if exclude
      File.new(exclude).each do |line|
         num, rest = line.chomp.split("\t")
         excluded[num.to_i] = true
      end
   end
   ncbi = Hash.new
   ncbi_line = Hash.new
   STDERR << "Loading ncbi taxonomy...\n"
   File.new("ncbi_taxonomy.txt").each do |line|
      num, name, parent = line.chomp.split("\t")
      if !excluded[num.to_i]
         ncbi[name] = num.to_i
         ncbi_line[num.to_i] = [name, parent.to_i]
      end
   end
   [ncbi, ncbi_line]
end

# load existing taxonomy
def get_phylodb_taxonomy(existing, ncbi, ncbi_line)
  STDERR << "Loading existing taxonomy...\n"
  phylotax = ncbi.dup
  phylopar = Hash.new
  phyloused = Hash.new
  if existing
    File.new(existing).each do |line|
      num, name, parent = line.chomp.split("\t")
      if !phylotax[name]
        phylotax[name] = num.to_i
        phyloused[num.to_i] = name
      end
    end
  end
  phylotax["root"] = 1
  phylopar[1] = 0
  [phylotax, phylopar, phyloused]
end

# make seguid file for phylodb
def writePhyloSeguid(db)
   out = File.new("phylodb.seguid","w")
   db.query("SELECT proteins.name, seguid, contigs.taxon_id, contigs.species FROM proteins, contigs WHERE contigs.name=contig_name").each do |row|
      out.print row.join("\t")+"\n"
   end
   out.close
end

def loadPhyloDB(phylo)
   phyloseg = Hash.new
   spHash = Hash.new
   segHash = Hash.new
   STDERR << "Loading phylodb seguid...\n"
   File.new(phylo).each do |line|
      id, seguid, taxonid, sp = line.chomp.split("\t")
      taxonid = taxonid.to_i
      phyloseg[taxonid] = Hash.new if !phyloseg[taxonid]
      phyloseg[taxonid][id] = seguid
      spHash[taxonid] = sp if !spHash[taxonid]
      segHash[seguid] = true
   end
   [spHash, phyloseg, segHash]
end

# parse slash parsed tags to hash
def parseSlashDesc(desc)
   hash = Hash.new
   desc.split("/").each do |entry|
      key,value = entry.split("=",2)
      if !value.nil?
         value.gsub!("\"","")
         hash[key] = value.chomp(" ")
      end
   end
   hash
end

# fast(er) simple FASTA parser
class FastaParser
   def initialize(file)
      if file.index(".gz")
         @fp = IO.popen("gunzip -c #{file}")
      elsif file.index(".bz2")
         @fp = IO.popen("bunzip2 -c #{file}")
      else
         @fp = File.new(file)
      end
   end
   def each
      @fp.each("\n>") do |seq|
         seq.chomp!(">")
         if seq != ""
            header, seq = seq.split("\n", 2)
            name, desc = header.split(" ", 2)
            if !seq.nil?
               seq.tr!("\n","")
               name.gsub!(">","")
               yield [name, desc, seq]
            end
         end
      end
   end
end

# delete data in Phylodb based on description field
def deleteDescriptionPhyloDB(db, description)
  seen = Hash.new
  STDERR << "Processing " << description << "...\n"
  db.query("SELECT name FROM contigs WHERE description=\"#{description}\"").each do |contig|
    seen[contig.first] = true
  end
  infield = '("'+seen.keys.join('","') + '")'
  ["geneorders", "proteins", "rrnas", "transcripts"].each do |tbl|
    db.query("DELETE FROM #{tbl} WHERE contig_name IN #{infield}")
  end
  db.query("DELETE FROM contigs WHERE name IN #{infield}")
end

# delete data in Phylodb for given taxonid
def deleteTaxonFromPhyloDB(db, num)
   contigs = Hash.new
   db.query("SELECT name FROM contigs WHERE taxon_id = #{num}").each do |contig|
      contigs[contig.first] = true
   end
   contigs.keys.each do |key|
      ["contigs", "geneorders", "proteins", "rrnas", "transcripts"].each do |table|
         if (table == "contigs")
            name = "name"
         else
            name = "contig_name"
         end
         query = "DELETE FROM #{table} WHERE #{name} = '#{key}'"
         STDERR << query + "...\n"
         db.query(query)
      end
   end
end



# is the taxonomy too long?
def tax_too_long(tax)
   if tax[0] == "Eukaryota"
      goodLen = 8
   else
      goodLen = 7
   end
   tax.length > goodLen
end

# is the taxonomy  too long?
def tax_too_short(tax)
   if tax[0] == "Eukaryota"
      goodLen = 8
   else
      goodLen = 7
   end
   tax.length < goodLen
end

# check length of taxonomy array based on Euk or not
def check_taxonomic_length(tax)
   !tax_too_short(tax) && !tax_too_long(tax)
end

# make sure every taxon has one and only one parent
def check_taxonomic_consistency(file)
   parents = Hash.new
   levels = Hash.new
   prob = false
   short = File.new("short.txt", "w")
   long = File.new("long.txt", "w")
   mult = File.new("mult.txt", "w")
   scount = 0
   lcount = 0
   mcount = 0
   File.new(file).each do |line|
      tag, taxonomy = line.chomp.split("\t")
      taxonomy = taxonomy.to_s.split(";")
      if tax_too_short(taxonomy)
         scount += 1
         short.print line
         prob = true
      elsif tax_too_long(taxonomy)
         lcount += 1
         long.print line
         prob = true
      end
      rtax = taxonomy.reverse
      0.upto(rtax.size - 2) do |i|
         child, parent = rtax[i], rtax[i+1]
         if !parents[child]
            parents[child] = [parent]
            levels[child] = [parent]
         elsif !parents[child].include?(parent)
            parents[child].push(parent)
         end
      end
   end
   parents.keys.each do |child|
      if parents[child].size > 1
         mult << child << " has multiple parents " << parents[child].join("||") << "\n"
         mcount += 1
         prob = true
      end
   end
   if prob
      STDERR << scount << " too short. " << lcount << " too long. " << mcount << " with multiple parents\n"
   else
      STDERR << "Taxonomy looks good! Hooray!\n"
   end
   short.close
   long.close
   mult.close
end

# patch taxstring based on patch hash
def patch_taxonomy(taxstring, patch)
   patch.keys.each do |pat|
      if taxstring =~/#{pat}/
         taxstring.sub!(/#{pat}/, patch[pat])
      end
   end
   taxstring
end

def loadKegg(kegg, keggtaxon, spHash)
   keggtax = Hash.new
   keggseg = Hash.new
   STDERR << "Loading kegg seguid...\n"
   File.new(keggtaxon).each do |line|
      abbr, taxonid, sp = line.chomp.split("\t")
      taxonid = taxonid.to_i
      keggtax[abbr] = taxonid
      spHash[taxonid] = sp if !spHash[taxonid]
   end
   File.new(kegg).each do |line|
      id, seguid = line.chomp.split("\t")
      abbr = id.split(":").first
      taxonid = keggtax[abbr]
      if (taxonid.nil?)
         STDERR << "Problem: No taxon for " << abbr << "\n"
         next
      else
         keggseg[taxonid] = Hash.new if !keggseg[taxonid]
         keggseg[taxonid][id] = seguid
      end
   end
   keggseg
end

def compareGenome(segs1, segs2)
   count1 = segs1.values.uniq.size
   count2 = segs2.values.uniq.size
   if count1 > count2
      overlap = 1 - (segs1.values.uniq - segs2.values.uniq).size / count1.to_f
   else
      overlap = 1 - (segs2.values.uniq - segs1.values.uniq).size / count2.to_f
   end
   [count1, count2, overlap]
end

def comparePhyloDB(segs1, segs2)
   count1 = segs1.values.uniq.size
   count2 = 0
   segs1.values.uniq.each do |seguid|
      count2 += 1 if segs2[seguid]
   end
   [count2, count1, count2/count1.to_f]
end

def compareKeggPhyloDB(taxonomy, phylo, kegg, keggtaxon)
   spHash, phyloseg, segHash = loadPhyloDB(phylo)
   keggseg = loadKegg(kegg, keggtaxon, spHash)
   updates = File.new("kegg_updates.csv", "w")
   updates.print ["Taxon Id", "Species", "Genes in PhyloDB", "Genes in Kegg", "Percent Overlap", "Taxonomy"].to_csv
   newgenomes = File.new("kegg_new.csv", "w")
   newgenomes.print ["Taxon Id", "Species", "Genes in PhyloDB", "Genes in Kegg", "Overlap with PhyloDB", "Taxonomy"].to_csv
   STDERR << "Comparing Kegg and PhyloDB...\n"
   keggseg.keys.each do |taxonid|
      if phyloseg[taxonid]
         phylo, kegg, overlap = compareGenome(phyloseg[taxonid], keggseg[taxonid])
         if overlap < 0.99
            updates.print [taxonid, spHash[taxonid], phylo, kegg,(overlap*100).to_i/100.0, taxonomy[taxonid]].to_csv
         end
      else
         phylo, kegg, overlap = comparePhyloDB(keggseg[taxonid], segHash)
         if overlap < 0.99
            newgenomes.print [taxonid, spHash[taxonid], phylo, kegg,(overlap*100).to_i/100.0, taxonomy[taxonid]].to_csv
         end
      end
   end
   updates.close
   newgenomes.close
end

def recurseTaxonomy(db, current, count = 0)
   if db.class == Hash # not using db
      name, parent = db[current]
   else
      name, parent = db.query("SELECT name, parent_id FROM taxonomy WHERE taxon_id=#{current}").fetch_row
   end
   if (current == 1 || name.nil? || name == "Bacteria" || name == "Eukaryota" || name == "Viruses" || name == "Archaea" || count > 50)
      [name]
   else
      recurseTaxonomy(db, parent, count + 1).to_a + [name]
   end
end

