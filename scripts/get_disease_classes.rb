#! /usr/bin/env ruby

require 'net/http'
require 'uri'
require 'csv'

def load_disease_data(file)
	disease_classes = {}
	count = 0
	File.open(file).each do |line|
		count += 1
		next if count == 1
		line.chomp!
		fields = line.split("\t")
		dis_class = fields[2]
		genes = fields.last.scan(/\((\d+)\)/).flatten
		query = disease_classes[dis_class]
		if query.nil?
			disease_classes[dis_class] = genes
		else
			disease_classes[dis_class] = genes | query
		end
	end	
	return disease_classes
end



min_genes = ARGV[1].to_i
disease_classes = load_disease_data(ARGV[0])
disease_classes.select!{ | dis_class, genes| genes.length >= min_genes && dis_class != 'Unclassified'}
disease_classes.each do |dis_class, genes|
	puts "#{dis_class}\t#{genes.join(',')}"
end

=begin
all_genes = disease_classes.values.flatten.uniq

temp_file = 'temp_data'
if !File.exists?(temp_file)
	xml_template = '<?xml version="1.0" encoding="UTF-8"?>
	<!DOCTYPE Query>
	<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
		<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
			<Filter name = "entrezgene" value = "' + all_genes.join(',')  + '"/>
			<Attribute name = "ensembl_peptide_id" />
			<Attribute name = "entrezgene" />
		</Dataset>
	</Query>'
	url_string = 'http://www.ensembl.org/biomart/martservice' 
	xml_template.gsub(/[\t\n]/,'')
	uri = URI.parse(url_string)
	http = Net::HTTP.new(uri.host, uri.port)
	request = Net::HTTP::Post.new(uri.request_uri)
	request.body = 'query=' + xml_template.gsub(/[\t\n]/,'')
	response = http.request(request)
	query_data = response.body
	File.open(temp_file, 'w') {|f| f.print query_data}
else
	query_data = File.open(temp_file).read
end
puts CSV.parse(query_data, { :col_sep => "\t" }).select{|q| !q[0].nil?}.inspect
=end
