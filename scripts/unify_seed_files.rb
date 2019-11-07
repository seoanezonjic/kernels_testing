#! /usr/bin/env ruby

# @author Fernando Moreno Jabato (jabato<at>uma<dot>es) based on Pedro Seoane Zonjic code

# Load neceessary functionalities
require 'optparse'

# Define necessary functions
def load_input_list(file)
    return File.open(file).readlines.map!{|line| line.chomp}
end

def write_Hash(ofile,info,sep=',')
	File.open(ofile, 'w') do |f|
		info.each do |target,seeds|
			f.puts target + "\t" + seeds.join(sep)
		end
	end
end

def load_dictionary(id_relations_file)
    relations = {}
    File.open(id_relations_file).each do |line|
        line.chomp!
        to_id, from_ids = line.split("\t")
        from_ids.split("|").each do |from_id|
            relations[from_id] = to_id
        end
    end
    return relations
end

def convert_ids(seed_names, dictionary)
    return seed_names.map{|s| dictionary[s]}.compact
end


######################################################################################
##                                   OPTPARSE                                       ##
######################################################################################
options = {}

optparse = OptionParser.new do |opts|
    options[:seeds_path] = nil
    opts.on( '-s', '--seeds_path PATH', 'Path where seed files are stored. File name will be used as target and file content as seeds' ) do |opt|
        options[:seeds_path] = opt
    end

    options[:output_file] = nil
    opts.on( '-o', '--output_file PATH', 'Output unified file' ) do |opt|
        options[:output_file] = opt
    end

    options[:id_conversion] = nil
    opts.on( '-i', '--id_conversion PATH', '[OPTIONAL] Dictionary to translate IDs <To, From>' ) do |opt|
        options[:id_conversion] = opt
    end

    opts.banner = "Usage: #{File.basename(__FILE__)} [options] \n\n"

    opts.on( '-h', '--help', 'Display this screen' ) do
            puts opts
            exit
    end
end

optparse.parse!



######################################################################################
##                                 HANDLE DATA                                      ##
######################################################################################

# Obtain target files
seeds_files = Dir.glob(File.join(options[:seeds_path], '*')).select {|f| !File.directory? f}
# Load seed files
targets = []
seeds_files.each do |seed_file|
	# Load and store
	targets << [File.basename(seed_file), load_input_list( seed_file)]
end

# Check if translation are available
if !options[:id_conversion].nil?
    # Load dictionary
    dictionary = load_dictionary(options[:id_conversion])
    # Translate
    targets.map!{|sg| [dictionary[sg[0]], convert_ids(sg[1], dictionary).uniq] }
    # targets.map!{|sg| [sg[0], convert_ids(sg[1],dictionary)] }
end

# Write in te correct format
write_Hash(options[:output_file], targets)
