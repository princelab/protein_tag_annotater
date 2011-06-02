#!/usr/bin/env ruby

def upcaser(sequence, test_match, position_arr = nil)
	if position_arr.class == NilClass
	elsif position_arr.class == Array and not position_arr.empty?
		regex = test_match 
		position_arr.flatten.each do |i|
			a =  sequence.to_s.rpartition(/#{regex}/i)
			a[1][i.to_i-1] = a[1][(i.to_i-1)].upcase unless a[1][i.to_i-1].nil? || a[1][i.to_i-1].empty?
			sequence = a.shift + a.first + a.last 
		end
	else
		regex = test_match
		sequence = sequence.gsub(/#{regex}/i, "#{regex.upcase}")
	end
#	puts "Sequence: #{sequence}"
	sequence
end
def bolder(sequence)
	sequence.gsub(/([A-Z]+)/, '<strong>\1</strong>')
end
def markdown(sequence)
	sequence.gsub(/([A-Z]+)/, '*\1*')
end

def colorify(sequence)
	sequence.gsub(/([A-Z]+)/, '<span class="colorify">\1</span>')
end	
def color_misses(sequence)
	sequence.gsub(/(k)/, '<span class="miss">\1</span>')
end

if ARGV.size < 1
	puts "Usage: #{File.basename(__FILE__)} Database.fasta 1.pepxml 2.pepxml ... ?.pepxml "
	puts "Output: 1.coverage 2.coverage ... ?.coverage"
	exit
end
require 'bio'
require 'nokogiri'
#fasta_file = Bio::FlatFile.auto(ARGV.shift)
fasta_file = '/home/ryanmt/lab/DB/uni_bovin_var_100518_fwd.fasta'
fasta = Bio::FlatFile.auto(fasta_file)
prot_db = Hash.new { |h,k| h[k] = [] }
pep_db = Hash.new { |h,k| h[k] = [] }
fasta.each_entry do |entry|
	prot_db[entry.accession.to_sym] << entry.seq
	pep_db[entry.seq.to_sym] << entry.accession
end

ARGV.each do |file|
	doc = Nokogiri.XML(File.open(file))
	top_hits = doc.xpath('//xmlns:search_hit[@hit_rank="1"]')
	proteins = Hash.new { |h,k| h[k] = [] }
	pep_centric = Hash.new { |h,k| h[k] = [] }
	amidinated = []; non_amidinated = []
	ProteinEvidence = Struct.new(:accession, :peptide, :modification_arr, :evidence_arr )
	top_hits.each do |hit| 
		prot = hit.to_s[/protein="(\w*?)"/,1]
		peptide = hit.to_s[/peptide="([A-Z]*?)"/,1]
		next if peptide.nil? or prot.nil?
		amidination_site = hit.to_s.scan(/mod_aminoacid_mass\sposition="(\d*?)"\smass="169.121512"/)
		pep_centric[peptide.to_sym] << prot
		if amidination_site.empty?
			proteins[prot.to_sym] << ProteinEvidence.new(prot, peptide, nil, amidination_site)
		else
			proteins[prot.to_sym] << ProteinEvidence.new(prot, peptide, amidination_site, [])
		end
	end
######################################################################
#	outputing 
	# pep_db vs pep_centric
	# prot_db vs proteins
	output = []
	ProteinModMap = Struct.new(:accession, :modified_seq, :observed_norm_seq, :output)
	prot_mods = proteins.map do |key, evidence_structs_array|
#		p key
#		p evidence_structs_array
		output_arr = []
		sequence = prot_db[key.to_sym]
		if sequence.length > 1
			output_arr << "============"
			output_arr << "Protein: #{key} failed because of multiple prot_db sequences"
			output_arr << "============"
		else
			output_arr << "##ACCESSION: #{key}"
#			output_arr << "#{evidence_structs_array}"
			tmp = sequence.first.downcase
			modified_seq = tmp
			observed_norm_seq = tmp
			evidence_structs_array.each do |evidence_struct| # 	ProteinEvidence = Struct.new(:accession, :peptide, :modification_arr, :evidence_arr )
				modified_seq = upcaser(modified_seq, evidence_struct.peptide, evidence_struct.modification_arr)
				#puts "OBS: #{observed_norm_seq}"
				observed_norm_seq = upcaser(observed_norm_seq, evidence_struct.peptide, evidence_struct.evidence_arr)
			end
			output_arr << "Modified:\n"
			output_arr << "Observed:"
			output_arr << color_misses(bolder(modified_seq))
			output_arr << bolder(observed_norm_seq)
		end

		ProteinModMap.new(key, modified_seq, observed_norm_seq, output_arr)
	end
output << "#HEADER: \n##ACCESSION: example\n##Modified\n_Sequence goes here_\n##Unmodified\n_Sequence goes here_\n\n------------\n#Data\n"
#puts output_arr
prot_mods.each {|prot| output << prot.output.join("\n")}
outstring = output.join("\n")
#p prot_mods


file_outname = File.basename(file).gsub(File.extname(file), '') + ".markdown"
File.open(file_outname,'w'){|f|
	f.print outstring
}
end
	

# => [:definition, :definition=, :data, :data=, :entry_overrun, :entry, :query, :fasta, :blast, :seq, :comment, :length, :naseq, :nalen, :aaseq, :aalen, :to_biosequence, :to_seq, :identifiers, :entry_id, :gi, :accession, :accessions, :acc_version, :locus, :tags, :exists?, :get, :fetch]

