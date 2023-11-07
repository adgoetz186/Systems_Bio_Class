import random
import numpy as np
import copy


class rxn_model:
	def __init__(self, rxn_dict=None,param_dict = None):
		if rxn_dict is not None:
			self.rxn_dict = rxn_dict
		else:
			self.rxn_dict = {}
		if param_dict is not None:
			self.param_dict = param_dict
		else:
			self.param_dict = {}
	
	def add_rxn(self,name,reactants,products,propensity):
		self.rxn_dict[name] = rxn(reactants,products,propensity)
		
	def add_param(self,name_or_dict,value=None):
		if isinstance(name_or_dict,dict):
			self.param_dict.update(name_or_dict)
		else:
			self.param_dict[name_or_dict] = value
	
	def Gillespie_simulate(self, t_vals,ivc = None,n = 1):
		species_values = {}
		propensity_values = {}
		for rxn in self.rxn_dict.keys():
			propensity_val = self.rxn_dict[rxn].propensity
			params = list(self.param_dict.keys())
			sorted_params = list(sorted(params, key = len,reverse=True))

			for param in sorted_params:
				propensity_val = propensity_val.replace(param,str(self.param_dict[param]))
			propensity_values[rxn] = propensity_val
			

		for rxn in self.rxn_dict.keys():
			species_values.update(self.rxn_dict[rxn].reactant_dict)
			species_values.update(self.rxn_dict[rxn].product_dict)
		
		# initializes all species, does this from reactions so species dont need to be explicitly declared
		for rxn in species_values.keys():
			species_values[rxn] = 0

		# sets initial value conditions
		if ivc is not None:
			species_values.update(ivc)

		
		initial_species = copy.copy(species_values)
		prop_names_ordered = list(propensity_values.keys())
		prop_values_ordered = list(propensity_values.values())

		# initializes trajectories for all species
		species_values_traj = copy.copy(species_values)
		for species in species_values_traj.keys():
			species_values_traj[species] = np.zeros((n,np.size(t_vals)))
			species_values_traj[species][:,0] = copy.copy(species_values[species])
		for repeat_ind in range(n):
			t = 0
			prop_values_ordered_eval = np.zeros(len(prop_values_ordered))
			species_values = copy.copy(initial_species)
			for t_val_ind in range(1,len(t_vals[1:])+1):
				while t_val_ind < len(t_vals) and t <= t_vals[t_val_ind]:

					# makes increasing list of propensities
					for ind_val in range(len(prop_values_ordered)):
						propen = copy.copy(prop_values_ordered[ind_val])
						
						species = list(species_values.keys())
						sorted_species = list(sorted(species, key=len, reverse=True))
						for species in sorted_species:
							propen = propen.replace(species,str(species_values[species]))
						if ind_val == 0:
							prop_values_ordered_eval[ind_val] = eval(propen)
						else:
							prop_values_ordered_eval[ind_val] = eval(propen) + prop_values_ordered_eval[ind_val-1]
					
					# randomly generates exponential timestep
					sum_val = np.max(prop_values_ordered_eval)
					rand_val = random.random()
					t += 1/sum_val*np.log(1/rand_val)
					
					# if the timestep goes over the time partition set the species values to be that of the previous iteration
					while t_val_ind < len(t_vals) and t >= t_vals[t_val_ind]:
						for species in species_values_traj.keys():
							species_values_traj[species][repeat_ind,t_val_ind] = species_values[species]
						t_val_ind+=1

					# selects a random reaction to fire (in weighted fashion) before running that reaction.
					prop_values_ordered_eval /= sum_val
					rand_val = random.random()
					rxn_to_use = prop_names_ordered[np.where((prop_values_ordered_eval > rand_val))[0][0]]
					for species in self.rxn_dict[rxn_to_use].reactant_dict.keys():
						species_values[species] -= self.rxn_dict[rxn_to_use].reactant_dict[species]
					for species in self.rxn_dict[rxn_to_use].product_dict.keys():
						species_values[species] += self.rxn_dict[rxn_to_use].product_dict[species]

		return species_values_traj


class rxn:
	def __init__(self, reactant_dict=None,product_dict = None,propensity = None):
		self.reactant_dict = reactant_dict
		self.product_dict = product_dict
		self.propensity = propensity
		

		
