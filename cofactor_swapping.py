import cobra
""" e. coli ijo1366 model """
swap_reactions = ["NADH16pp", "GAPD", "G6PDH2r", "GND", "PDH", "FADRx", "PGCD", "ICDHyr", "ASAD","MDH", "MTHFD", "KARA1", "HSDy", "IPMD", "SHK3Dr", "ACALD", "ALCD2x", "LCARR","LDH_D", "ME1", "ME2"]

def opt_yield(model, swap_reactions, target_reaction, biomass_reaction=None, min_biomass=0, max_swaps=1, copy_model=True):
    if copy_model:
        model = model.copy()

    if biomass_reaction is None:
        biomass_reaction = (reaction.id for reaction in model.reactions if reaction.objective_coefficient != 0).next()

    # initialize max swaps constraint
    max_swaps_constraint = cobra.core.Metabolite('__max_swaps_constraint')

    for r in swap_reactions:
        reaction = model.reactions.getbyId(r)
        
        #addConstr
        #add non native reactions
        non_native_reaction = cobra.core.Reaction(r + "__swap")

        new_mets = {}
        for met, stoich in reaction.metabolites.iteritems():
            if met.id == 'nad_c':
                new_mets[model.metabolites.getbyId('nadp_c')] = stoich
            elif met.id == 'nadh_c':
                new_mets[model.metabolites.getbyId('nadph_c')] = stoich
            else:
                new_mets[met] = stoich
        non_native_reaction.add_metabolites(new_mets)
        model.add_reaction(non_native_reaction)

        # s_d
        active_var = cobra.core.Reaction(r+'__active')
        active_var.lower_bound = 0
        active_var.upper_bound = 1
        active_var._variable_type = 'integer'
        # t_d
        swap_active_var = cobra.core.Reaction(r+'__swap_active')
        swap_active_var.lower_bound = 0
        swap_active_var.upper_bound = 1
        swap_active_var._variable_type = 'integer'

        # s_t + t_d = 1
        swap_constraint = cobra.core.Metabolite(r+'__swap_constraint')
        swap_constraint._constraint_sense = "E"
        swap_constraint._bound = 1
        swap_active_var.add_metabolites({swap_constraint: 1})
        active_var.add_metabolites({swap_constraint: 1})

        # s_d LB_x_d <= v_x_d 
        # LB_x_d * s_d - v_x_d <= 0
        active_lower_bound = cobra.core.Metabolite(r+'__active_lower_bound')
        active_lower_bound._constraint_sense = "LE" #look up L or LE
        active_lower_bound._bound = 0
        active_var.add_metabolites({active_lower_bound: reaction.lower_bound})
        reaction.add_metabolites({active_lower_bound: -1})
        # repeat for other 3 

        
        


    for reaction in model.reactions:
        reaction.objective_coefficient = (1 if reaction.id == target_reaction else 0)
    
    model.reactions.getbyId(biomass_reaction).lower_bound = min_biomass    
        
    solution =  model.optimize()
    # solution.swapped_reactions = ... # e.g. ['GAPD']
    # solution.not_swapped_reactions = ... # e.g. ["NADH16pp", "G6PDH2r", "GND", "PDH", "FADRx", "PGCD", "ICDHyr", "ASAD","MDH", "MTHFD", "KARA1", "HSDy", "IPMD", "SHK3Dr", "ACALD", "ALCD2x", "LCARR","LDH_D", "ME1", "ME2"]

    return solution
