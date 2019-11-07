#! /usr/bin/env ruby

# @author Fernando Moreno Jabato (jabato<at>uma<dot>es) based on Pedro Seoane Zonjic code


# Load neceessary functionalities
require 'optparse'
require 'nmatrix'

# Define necessary functions
def load_input_list(file)
    return File.open(file).readlines.map!{|line| line.chomp}
end

def load_seed_groups(seed_group_file)
	seed_groups = []
	File.open(seed_group_file).each do |line|
		line.chomp!
		id, members = line.split("\t")
		members = members.split(',')
		seed_groups << [id, members]
	end
	return seed_groups
end

def load_matrix_file(input_file)
	dimension_elements = 0
	adjacency_vector = []
	File.open(input_file).each do |line|
 	   	line.chomp!
    		adjacency_vector.concat(line.split("\t").map{|c| c.to_f })
	    	dimension_elements += 1
	end
	matrix = NMatrix.new([dimension_elements, dimension_elements], adjacency_vector) # Create working matrix
	return matrix
end

def get_disjoint_matrix(nmatrix, cols)
	dj_matrix = nil
	cols.each do |col|
		nmatrix_col = nmatrix.column(col, :reference)
		if dj_matrix.nil?
			dj_matrix = nmatrix_col
		else
			dj_matrix = dj_matrix.concat(nmatrix_col)
		end
	end
	return dj_matrix
end

def get_ranked_list(score_matrix, seed_cols)
	integrated_scores = []
	count = 0
	seed_matrix =  get_disjoint_matrix(score_matrix, seed_cols)
	score_list = seed_matrix.sum(1) / seed_cols.length
	score_list.each do |score|
		integrated_scores << [count, score] if !seed_cols.include?(count)
		count += 1
	end
	integrated_scores.sort!{|is1, is2| is2[1] <=> is1[1]}
	return integrated_scores
end






def median(array, already_sorted=false)
	return nil if array.empty?
	array = array.sort unless already_sorted
	m_pos = array.size / 2
	return array.size % 2 == 1 ? array[m_pos] : array[m_pos-1..m_pos].mean 
end

module Enumerable

    def sum
      self.inject(0){|accum, i| accum + i }
    end

    def mean
      self.sum/self.length.to_f
    end

    def sample_variance
      m = self.mean
      sum = self.inject(0){|accum, i| accum +(i-m)**2 }
      sum/(self.length - 1).to_f
    end

    def standard_deviation
      Math.sqrt(self.sample_variance)
    end

end 


# True Positive Rate above a given threshold
def TPR(ranks,universe_length,threshold = 0.1)
	# Metric limit
	limit = threshold * universe_length
	# Obtain number of ranks which are found before LIMIT
	num_uppers = ranks.select{|rank| rank <= limit}.length
	# return rate
	return (num_uppers*100 / ranks.length)
end












######################################################################################
##                                   OPTPARSE                                       ##
######################################################################################
options = {}

optparse = OptionParser.new do |opts|
    options[:targets_file] = nil
    opts.on( '-s', '--seeds PATH', 'Targets and seeds file') do |opt|
        options[:targets_file] = opt
    end

    options[:matrix_file] = nil
    opts.on( '-m', '--matrix_file PATH', 'Input matrix file' ) do |opt|
        options[:matrix_file] = opt
    end

    options[:input_matrix_format] = 'bin'
    opts.on( '-f', '--input_matrix_format STRING', 'Input matrix file format. "bin" or "text"' ) do |opt|
        options[:input_matrix_format] = opt
    end

    options[:matrix_ids] = nil
    opts.on( '-n', '--node_ids PATH', 'Matrix element IDs file' ) do |opt|
	options[:matrix_ids] = opt
    end

    options[:output_file] = nil
    opts.on( '-o', '--output_file PATH', 'Output metrics file') do |opt|
        options[:output_file] = opt
    end

    opts.banner = "Usage: #{File.basename(__FILE__)} [options] \n\n"

    opts.on( '-h', '--help', 'Display this screen' ) do
            puts opts
            exit
    end
end



optparse.parse!



######################################################################################
##                               LOAD AND EVALUATE                                  ##
######################################################################################

# Load data
seed_groups = load_seed_groups(options[:targets_file])  # Targets
node_names = load_input_list(options[:matrix_ids])      # Kernel elements
if options[:input_matrix_format] == 'bin'               # Kernel
	score_matrix = NMatrix.read(options[:matrix_file])
elsif  options[:input_matrix_format] == 'text'
	score_matrix = load_matrix_file(options[:matrix_file])
end


# Tranform IDs to Matrix coordinates
seed_groups.map!{|gr| [gr[0], gr[1].map{|seed_name| node_names.index(seed_name)}.compact]}

# Obtain result info per each target
target_ranks = []
seed_groups.each do |target,seeds_coords|
	# Find target position
	t_indx = node_names.index(target)
	# Check
	next if t_indx.nil?
	# Rank for this target
	ranked_list = get_ranked_list(score_matrix,seeds_coords)
	# Find target into ranked list
	target_rank = nil
	ranked_list.each_with_index do |element, i|
		if t_indx == element[0]
			target_rank = i
			break
		end
	end
	# Add new target rank
	target_ranks << target_rank if !target_rank.nil?
end

raise 'Any TARGET (from seeds) have been found' if target_ranks.empty?

# Evalueate
metrics = {}

# Calculate ratio of each rank
pos_ratio = target_ranks.map{|rank| rank*100/node_names.length}


# Eval!
metrics['median'] = median(pos_ratio).round(2)
metrics['std']    = pos_ratio.standard_deviation.round(2)
metrics['lowest'] = pos_ratio.max.round(2)
metrics['tpr01']  = TPR(target_ranks, node_names.length, 0.01).round(1)
metrics['tpr05']  = TPR(target_ranks, node_names.length, 0.05).round(1)
metrics['tpr10']  = TPR(target_ranks, node_names.length, 0.10).round(1)
metrics['tpr15']  = TPR(target_ranks, node_names.length, 0.15).round(1)
metrics['tpr20']  = TPR(target_ranks, node_names.length, 0.20).round(1)
metrics['tpr30']  = TPR(target_ranks, node_names.length, 0.30).round(1)

# Write metrics
File.open(options[:output_file], 'w') do |f|
	metrics.each do |metric,value|
		f.puts metric + "\t" + value.to_s
	end
end
