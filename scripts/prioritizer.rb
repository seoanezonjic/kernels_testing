#! /usr/bin/env ruby
require 'optparse'
require 'numo/narray'
require 'numo/linalg'
require 'npy'

#require 'pp'

##############################################################################
## METHODS
##############################################################################
# Statistical methods for Unbiased way evaluation methods (FERNANDO)
# ----------------------------------------------------------------------------
module Enumerable
    def sum
      return self.inject(0){|accum, i| accum + i }
    end

    def mean
      return self.sum.fdiv(self.length)
    end

    def median
    	res = nil
    	if !self.empty?
	    	sorted_arr = self.sort 
	    	m_pos = self.size / 2 
	    	res = self.size % 2 == 1 ? sorted_arr[m_pos] : sorted_arr[m_pos-1..m_pos].mean
	    end
	    return res
    end

    def sample_variance
      m = self.mean
      sum = self.inject(0){|accum, i| accum +(i-m)**2 }
      return sum.fdiv(self.length - 1)
    end

    def standard_deviation
      return Math.sqrt(self.sample_variance)
    end
end



# I/O 
#------------------------------------------------------------------------------
def load_matrix_file(input_file, splitChar = "\t")
	matrix = nil
	counter = 0
	File.open(input_file).each do |line|
	    	line.chomp!
    		row = line.split(splitChar).map{|c| c.to_f }
    		if matrix.nil?
    			matrix = Numo::DFloat.zeros(row.length, row.length)
    		end
    		row.each_with_index do |val, i|
    			matrix[counter, i] = val 
    		end
    		counter += 1
	end
	return matrix
end

def load_input_list(file)
    return File.open(file).readlines.map!{|line| line.chomp}
end

def load_dictionary(id_relations_file)
	relations = {}
	File.open(id_relations_file).each do |line|
		line.chomp!
		from_ids, to_id = line.split("\t")
		from_ids.split("|").each do |from_id|
			relations[from_id] = to_id
		end
	end
	return relations
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

def write_list(array, output_file)
	File.open(output_file, 'w') do |f|
		array.each do |record|
			f.puts record.join("\t")
		end
	end
end

# General methods 
#------------------------------------------------------------------------------
def get_ranked_list(score_matrix, seed_cols)
	integrated_scores = []
	count = 0
	seed_matrix =  score_matrix[true, seed_cols]
	score_list = get_score_list(seed_matrix, seed_cols)
	score_list.each do |score|
		integrated_scores << [count, score] if !seed_cols.include?(count)
		count += 1
	end
	integrated_scores.sort!{|is1, is2| is2[1] <=> is1[1]}
	return integrated_scores
end

def convert_ids(seed_names, dictionary)
	return seed_names.map{|s| dictionary[s]}.compact
end

def get_score_list(seed_matrix, seed_cols)
	score_list = seed_matrix.sum(1) / seed_cols.length
	return score_list
end

# classic evaluation
#------------------------------------------------------------------------------
def get_rank_percentage_for_member(score_matrix, seed_cols, test_member)
	seed_matrix = score_matrix[true, seed_cols]
	score_list = get_score_list(seed_matrix, seed_cols)
	members_below_test = 0
	ref_score = score_list[test_member]
	score_list.each do |score|
		members_below_test += 1 if score > ref_score 
	end
	return members_below_test.fdiv(score_list.rows - 1)
end

def evaluate_priotirizer_classic_way(score_matrix, gr_seed_cols)
	full_rank = []
	gr_seed_cols.each do |gs_id, gs_members|
		puts gs_id
		gs_members.each do |seed_col|
			seed_subset = gs_members.dup
			test_member = seed_subset.delete(seed_col)
			full_rank << ['P', get_rank_percentage_for_member(score_matrix, seed_subset, test_member)]
		end
		non_related_members = get_non_related_members(gr_seed_cols, gs_id, gs_members)
		non_related_members.each do |member_col|
			full_rank << ['N', get_rank_percentage_for_member(score_matrix, gs_members, member_col)]
		end
		break
	end
	return full_rank
end

def get_non_related_members(gr_seed_cols, gs_id, gs_members)
	gs_id_number = gs_members.length
	non_related_members = []
	(gs_id_number/2).times do
		current_gs = nil
		current_gs_id = gs_id
		while current_gs_id == gs_id # pick random group that are not the current one
			current_gs = gr_seed_cols.sample
			current_gs_id = current_gs.first 
		end
		member = nil
		while member.nil? # extract random gene
			extract = current_gs.last.sample
			member = extract if !non_related_members.include?(extract) && !gs_members.include?(extract)
		end
		non_related_members << member
	end
	return non_related_members
end

# Unbiased way evaluation methods (FERNANDO)
# ----------------------------------------------------------------------------
def get_ranking_of_targets(seed_groups, score_matrix, node_names)	
	target_ranks = [] # Obtain result info per each target
	seed_groups.each do |target, seeds_coords|		
		t_indx = node_names.index(target) # Find target position
		next if t_indx.nil? # Check
		ranked_list = get_ranked_list(score_matrix, seeds_coords) # Rank for this target
		target_rank = get_position_element(ranked_list, t_indx)
		target_ranks << target_rank if !target_rank.nil? # Add new target rank
	end
	return target_ranks
end

def get_position_element(ranked_list, original_position)
	target_rank = nil # Find target into ranked list
	ranked_list.each_with_index do |record, i|
		previous_position, score = record
		if original_position == previous_position
			target_rank = i
			break
		end
	end
	return target_rank
end

# True Positive Rate above a given threshold
def get_TPR(ranks, universe_length, threshold = 0.1)
	# Metric limit
	limit = threshold * universe_length
	# Obtain number of ranks which are found before LIMIT
	num_uppers = ranks.select{|rank| rank <= limit}.length
	# return rate
	return (num_uppers*100 / ranks.length)
end

def compute_unbiased_way_metrics(target_ranks, item_number,targets_size)
	pos_ratio = target_ranks.map{|rank| rank*100/item_number} # Calculate ratio of each rank
	metrics = {}
	metrics['median'] = pos_ratio.median.round(2)
	metrics['mean'] = pos_ratio.mean.round(2)
	metrics['std']    = pos_ratio.standard_deviation.round(2)
	metrics['lowest'] = pos_ratio.max.round(2)
	metrics['auc'] = get_auc_from_ranked_targets(target_ranks, item_number)
	metrics['auc30'] = get_auc_from_ranked_targets(target_ranks, item_number, 0.3)
	metrics['tpr01']  = get_TPR(target_ranks, item_number, 0.01).round(1)
	metrics['tpr05']  = get_TPR(target_ranks, item_number, 0.05).round(1)
	metrics['tpr10']  = get_TPR(target_ranks, item_number, 0.10).round(1)
	metrics['tpr15']  = get_TPR(target_ranks, item_number, 0.15).round(1)
	metrics['tpr20']  = get_TPR(target_ranks, item_number, 0.20).round(1)
	metrics['tpr30']  = get_TPR(target_ranks, item_number, 0.30).round(1)
	metrics['Ranked_items'] = target_ranks.length
	metrics['Target_items'] = targets_size
	metrics['Total_items'] = item_number
	metrics['Ranked_ratio'] = target_ranks.length / targets_size * 100
	return metrics
end

def get_auc_from_ranked_targets(target_ranks, ranked_list_length, threshold = 1)
	aucs = target_ranks.map{ |position| single_target_AUROC(position, ranked_list_length, threshold) }
	return aucs.mean
end

def single_target_AUROC(position, ranked_list_length, threshold)
	auc = 0
	position_threshold = (ranked_list_length * threshold).to_i
	if position <= position_threshold
		auc = (position_threshold - position + 1).fdiv(position_threshold)
	end
	return auc
end

##############################################################################
## OPTPARSE
##############################################################################
options = {}

optparse = OptionParser.new do |opts|
    options[:matrix_file] = nil
    opts.on( '-m', '--matrix_file PATH', 'Input matrix file' ) do |opt|
        options[:matrix_file] = opt
    end

    options[:input_matrix_format] = 'bin'
    opts.on( '-f', '--input_matrix_format STRING', 'Input matrix file format. "bin" or "text"' ) do |opt|
        options[:input_matrix_format] = opt
    end


    options[:id_conversion] = nil
    opts.on( '-c', '--conversion_file PATH', 'Conversion file between seed ids and kernel ids' ) do |opt|
	options[:id_conversion] = opt
    end

    options[:output_file] = nil
    opts.on( '-o', '--output_file PATH', 'Output matrix file' ) do |opt|
        options[:output_file] = opt
    end

    options[:seed_file] = nil
    opts.on( '-s', '--seed_file PATH', 'Input seed file' ) do |opt|
        options[:seed_file] = opt
    end

    options[:axis_names_file] = nil
    opts.on( '-n', '--axis_names_file PATH', 'Input matrix axis names file' ) do |opt|
        options[:axis_names_file] = opt
    end

    options[:evaluation_mode] = nil
    opts.on( '-e', '--eval_mode STRING', 'None for normal mode, "cla" for classic evaluation and "unb" for unbiased way evaluation' ) do |opt| 
        options[:evaluation_mode] = opt
    end

    options[:byte_format] = :float64
  	opts.on( '-b', '--byte_format STRING', 'Format of the numeric values stored in matrix. Default: float64, warning set this to less precission can modify computation results using this matrix.' ) do |opt|
    	options[:byte_format] = opt.to_sym
  	end

    opts.banner = "Usage: #{File.basename(__FILE__)} [options] \n\n"

    opts.on( '-h', '--help', 'Display this screen' ) do
            puts opts
            exit
    end
end

optparse.parse!

##############################################################################
## MAIN
##############################################################################
dictionary_ids = load_dictionary(options[:id_conversion]) if !options[:id_conversion].nil?
node_names = load_input_list(options[:axis_names_file])
if options[:input_matrix_format] == 'bin'
	#score_matrix = Marshal.load(File.binread(options[:matrix_file]))
	score_matrix = Npy.load(options[:matrix_file])
elsif  options[:input_matrix_format] == 'text'
	score_matrix = load_matrix_file(options[:matrix_file])
end
if !options[:evaluation_mode].nil?
	seed_groups = load_seed_groups(options[:seed_file]) #[[group_id, [seed_ids]]]
	seeds_length = seed_groups.length
	seed_groups.map!{|sg| [sg[0], convert_ids(sg[1], dictionary_ids).uniq] } if !options[:id_conversion].nil?
	seed_groups.map!{|sg| [sg[0], sg[1].map{|seed_name| node_names.index(seed_name)}.compact]} # Tranform IDs to Matrix coordinates
	if options[:evaluation_mode] == 'cla'
		output_data = evaluate_priotirizer_classic_way(score_matrix, seed_groups)
		output_data.unshift(['tag', 'rank'])
	elsif options[:evaluation_mode] == 'unb'
		ranked_targets = get_ranking_of_targets(seed_groups, score_matrix, node_names)	
		if ranked_targets.empty?
			raise 'ERROR: Any TARGET (from seeds) have been found' 
		else		
			output_data = compute_unbiased_way_metrics(ranked_targets, node_names.length, seeds_length).to_a
		end
	end
else
	seed_names = load_input_list(options[:seed_file]).uniq
	seed_names = convert_ids(seed_names, dictionary_ids).uniq if !options[:id_conversion].nil?
	seed_cols = seed_names.map{|seed_name| node_names.index(seed_name)}.compact
	missed_ids = seed_names.length - seed_cols.length
	warn("#{missed_ids} have not been found in names file") if missed_ids > 0
	output_data = get_ranked_list(score_matrix, seed_cols)
end
write_list(output_data, options[:output_file])
