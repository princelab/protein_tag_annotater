#!/usr/bin/env ruby

require 'nokogiri'

def cleanup_modification_array(array)
	outs = []
	array.each_index do |i|
		if array[i][0] == array[i-1][0] 
			if array[i][1] == array[i-1][1]
				array[i][2] << array[i-1][2]
				array[i-1] = nil
			end
		end
		array[i][2] = array[i][2].flatten
		outs << array[i]
	end
	outs
end
def upcaser(sequence, test_match)
	regex = test_match
	sequence.downcase.gsub(/#{regex.downcase}/, "#{regex.upcase}")
end
def bolder(sequence)
	sequence.gsub(/([A-Z]+)/, '<span id="bold">\1</span>')
end

def colorify(sequence)
	sequence.gsub(/([A-Z]+)/, '<span id="colorify">\1</span>')
end	



if ARGV.size < 1
	puts "Usage: #{File.basename(__FILE__)} Database.fasta 1.pepxml 2.pepxml ... ?.pepxml "
	puts "Output: 1.coverage 2.coverage ... ?.coverage"
	exit
end
require 'bio'
#fasta = Bio::FlatFile.auto(ARGV.shift)
FASTA = '/home/ryanmt/lab/DB/uni_bovin_var_100518_fwd.fasta'
fasta = Bio::FlatFile.auto(FASTA)
prot_db = Hash.new { |h,k| h[k] = [] }
pep_db = Hash.new { |h,k| h[k] = [] }
fasta.each_entry do |entry|
	prot_db[entry.accession.to_sym] << entry.seq
	pep_db[entry.seq.to_sym] << entry.accession
end
pep_db.each {|k,v| puts v if v.length > 10}
ARGV.each do |file|
	doc = Nokogiri.XML(File.open(file))
	top_hits = doc.xpath('//xmlns:search_hit[@hit_rank="1"]')
	proteins = Hash.new { |h,k| h[k] = [] }
	pep_centric = Hash.new { |h,k| h[k] = [] }
	amidinated = []; non_amidinated = []
	top_hits.each do |hit| 
		prot = hit.to_s[/protein="(\w*?)"/,1]
		peptide = hit.to_s[/peptide="([A-Z]*?)"/,1]
		next if peptide.nil? or prot.nil?
		amidination_site = hit.to_s.scan(/mod_aminoacid_mass\sposition="(\d*?)"\smass="169.121512"/)
		pep_centric[peptide.to_sym] << prot
		proteins[prot.to_sym] << peptide
		if amidination_site.empty?
			non_amidinated << [prot, peptide,[amidination_site]]
		else
			amidinated << [prot, peptide,[amidination_site]]
		end
	end
	amidinated = amidinated.sort_by{|a| a[1]}
	amid_outs = []; non_amid_outs = []
	non_amid_outs, amid_outs = [non_amidinated, amidinated].map do |array|
		cleanup_modification_array(array)
	end
	p amid_outs.length
	p non_amid_outs.length

######################################################################
#	outputing 
	# pep_db vs pep_centric
	# prot_db vs proteins
	output = []
	output << "HEADER: \nACCESSION: example\nModified:Unmodified"
	proteins.each do |key, value|
		sequence = prot_db[key]
		if sequence.length > 1
			output << "Protein: #{key} failed because of multiple prot_db sequences" 
			next
		end
		


end


# => [:definition, :definition=, :data, :data=, :entry_overrun, :entry, :query, :fasta, :blast, :seq, :comment, :length, :naseq, :nalen, :aaseq, :aalen, :to_biosequence, :to_seq, :identifiers, :entry_id, :gi, :accession, :accessions, :acc_version, :locus, :tags, :exists?, :get, :fetch]

