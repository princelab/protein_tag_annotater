require 'spec_helper'

describe "TaggingProject" do
	it 'upcases the string' do 
		seq = "YVSSKALQRQHSEGAAGKAPCILPIIENGK"
		test = "HSEGAAGK"
		upcaser(seq, test).should.equal "yvsskalqrqhseGAAGKAPCilpiiengk"
	end
	it 'bolds a sequence' do 
		seq = "yvsskalqrqhseGAAGKAPCilpiiengk"
		bolder(seq).should.equal "YVSSKALQRQ<span>HSEGAAGK</span>APCILPIIENGK"
	end
	it '"colorifies" a sequence' do 
		seq = "YVSSKALQRQHSEGAAGKAPCILPIIENGK"
		test2 = "GAAGKAPC"
		seq2 = upcaser(seq, test2)
		colorify(seq2).should.equal "yvsskalqrqhse<span id=\"colorify\">GAAGKAPC</span>ilpiiengk"
	end
end
