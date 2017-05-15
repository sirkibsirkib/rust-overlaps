

# a_match_len is the size (in chars) of the A string (from pattern)
# b_match_len is the size (in chars) of the B string (from text)

# index --block_id_lookup--> block_id 	(relative to beginning of pattern)
# filter_func --mode--> int				//how many errors allowed at this position
# ???? --candidate_condition--> bool	//returns whether can generate a candidate at this point

def forward(candidate_set, p, p_i_start, p_i_next, p_i_end, sp, ep, p_T_index, p_id,
			errors, filter_func, block_id_lookup, indel_balance, b_match_len, a_match_len):

	# dead end. No matches in index
	if sp > ep:
		return

	block_id : int        = block_id_lookup(p_i_next)
	permitted_error : int = filter_func(block_id-s_id, suff_len)
	spare_error : bool    = errors < permitted_error

	if self.arguments.indels and spare_error:
		# INSERTION
		if indel_balance >= 0:
			for a in self.sorted_alphabet:
				sp_ = self.C[a] + self.rank(a, sp)
				ep_ = self.C[a] + self.rank(a, ep + 1) - 1
				self.forward(candidate_set, p, p_i_start, p_i_next, p_i_end, sp_, ep_, p_T_index, p_id,
							 errors + 1, error_lookup, block_id_lookup, MATCHED + a.lower(), 1, b_match_len+1, a_match_len)


	# ADD SUFFIX CANDIDATES
	if self.condition_met_f(p_i_start, max(p_i_next, p_i_next-a_match_len+b_match_len), self.arguments.thresh, block_id_lookup[p_i_next-1-p_i_start], errors):
		sp_ = 0 + self.rank('$', sp)
		ep_ = 0 + self.rank('$', ep + 1) - 1
		a_ovr = a_match_len + p_i_start
		b_ovr = b_match_len + p_i_start
		debug_string = '' if not search_debug else 'MATCHED[{}] patt[{}]'.format(MATCHED, p[:p_i_end])
		for i in range(sp_, ep_ + 1):
			b_index = self.string_start_in_not_backwards_T(self.sSAT[i])
			x = self.new_candidate(a_index=p_T_index,
								   b_index=b_index,
								   a_ovr=a_ovr,
								   b_ovr=b_ovr,
								   b_tail=0,
								   debug_str=debug_string
								   )
			# if x[0] != x[1]:
			if x not in candidate_set:
				candidate_set.add(x)

	# no more characters in patt to match
	if p_i_next >= p_i_end:
		# INCLUSIONS!!!
		if self.arguments.inclusions and p_i_next > p_i_start:
			for i in range(sp, ep + 1):
				# print('\n\nINCLUSION')
				a_ovr = a_match_len + p_i_start
				b_ovr = b_match_len + p_i_start
				b_index, b_tail = self.index_inside_to_front_not_backwards(self.sSAT[i])
				x = self.new_candidate(a_index=p_T_index,
									   b_index=b_index,
									   a_ovr=a_ovr,
									   b_ovr=b_ovr,
									   b_tail=b_tail,
									   debug_str=debug_string
									   )
				if x[0] != x[1]:
					if x not in candidate_set:
						candidate_set.add(x)
					else:
						self.duplicate_candidate_count += 1
		return


	for a in self.sorted_alphabet:
		# SUBSTITUTION
		sp_ = self.C[a] + self.rank(a, sp)
		ep_ = self.C[a] + self.rank(a, ep + 1) - 1
		errors_ = errors if p[p_i_next] == a and a != 'N' else errors + 1
		if errors_ <= error_lookup[p_i_next-p_i_start]:
			self.forward(candidate_set, p, p_i_start, p_i_next+1, p_i_end, sp_, ep_, p_T_index, p_id,
					errors_, error_lookup, block_id_lookup, MATCHED+a, 0, b_match_len+1, a_match_len+1)

	if self.arguments.indels and errors < error_lookup[p_i_next-p_i_start]:
		# DELETION
		if indel_balance <= 0 and p_i_next < p_i_end-2:
			if indel_balance == -1 or p[p_i_next+1] != p[p_i_next]:
				self.forward(candidate_set, p, p_i_start, p_i_next + 1, p_i_end, sp, ep, p_T_index, p_id,
					 errors + 1, error_lookup, block_id_lookup, MATCHED + '_', -1, b_match_len, a_match_len+1)