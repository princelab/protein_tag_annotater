#!/usr/bin/env ruby

require 'nokogiri'

if ARGV.size < 1
	puts "Usage: #{File.basename(__FILE__)} 1.pepxml 2.pepxml ... ?.pepxml "
	puts "Output: 1.coverage 2.coverage ... ?.coverage"
	exit
end
FASTA = '/home/ryanmt/lab/DB/uni_bovin_var_100518_fwd.fasta'
fasta = File.open(FASTA, 'r')
ARGV.each do |file|
	doc = Nokogiri.XML(File.open(file))
	top_hits = doc.xpath('//xmlns:search_hit[@hit_rank="1"]')
	proteins = []
	pep_centric = []
	amidinated = top_hits.map do |hit| 
		prot = hit.to_s[/protein="(\w*?)"/,1]
		peptide = hit.to_s[/peptide="([A-Z]*?)"/,1]
		proteins << [prot, peptide]
		amidination_site = hit.to_s.slice(/mod_aminoacid_mass\sposition="(\d)"\smass="169.121512"/,1)
		[prot, peptide,[amidination_site]]
	end
	#amidinated = doc.xpath("//xmlns:search_hit[@hit_rank='1']/modification_info/mod_aminoacid_mass")#[@mass='169.121512']")
	p top_hits.size
	p amidinated.size
	amidinated = amidinated.sort_by{|a| a[1]}
	outs = []
	amidinated.each_index do |i|
		if amidinated[i][0] == amidinated[i-1][0] 
			if amidinated[i][1] == amidinated[i-1][1]
				amidinated[i][2] << amidinated[i-1][2]
				amidinated[i-1] = nil
			end
		end
		amidinated[i][2]= amidinated[i][2].flatten.compact.uniq
		outs << amidinated[i]
	end
p outs
		
end
