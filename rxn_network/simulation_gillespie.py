import rxn_model
import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd
import seaborn as sns

# very simple example of the simulation at work

# obtained from https://stackoverflow.com/questions/12309269/how-do-i-write-json-data-to-a-file
# might be useful for later, as of now is not used
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

p_error_B_mean = np.zeros((10,10))
p_error_A_mean = np.zeros((10,10))

p_error_B_var = np.zeros((10,10))
p_error_A_var = np.zeros((10,10))

for i in np.logspace(-2,2,10):
    for j in np.logspace(-2,2,10):
        reaction_model = rxn_model.rxn_model()
        reaction_model.add_rxn("bind_A",{"DNA_A_e":1,"A":1},{"DNA_A_f":1},"A*DNA_A_e*kA_on")
        reaction_model.add_rxn("unbind_A",{"DNA_A_f":1},{"DNA_A_e":1,"A":1},"DNA_A_f*kA_off")
        
        reaction_model.add_rxn("bind_B",{"DNA_B_e":1,"B":1},{"DNA_B_f":1},"B*DNA_B_e*kB_on")
        reaction_model.add_rxn("unbind_B",{"DNA_B_f":1},{"DNA_B_e":1,"B":1},"DNA_B_f*kB_off")
        
        reaction_model.add_param({"kA_on":i,"kA_off":1,"kB_on":j,"kB_off":1})
        
        
        
        
        
        tvals = np.linspace(0,150,150)
        result_dict = reaction_model.Gillespie_simulate(tvals,ivc={"DNA_A_e": 1,"DNA_B_e": 1,"A":1,"B":1},n = 1000)

        p_error_B_mean[int(2+np.log10(i)),int(2+np.log10(j))] = np.abs(np.average(result_dict["DNA_A_f"][:,-1]) - j/(j+1) / (j/(j+1)))
        p_error_A_mean[int(2+np.log10(i)),int(2+np.log10(j))] = np.abs(np.average(result_dict["DNA_B_f"][:,-1]) - i/(i+1) / (i/(i+1)))

        p_error_B_var[int(2+np.log10(i)),int(2+np.log10(j))] = np.abs((j / (j + 1)**2 - np.var(result_dict["DNA_B_f"][:, -1]))/((j / (j + 1)**2)))
        p_error_A_var[int(2+np.log10(i)),int(2+np.log10(j))] = np.abs((i / (i + 1)**2 - np.var(result_dict["DNA_A_f"][:, -1]))/((i / (i + 1)**2)))
        print([i,j])
sns.heatmap(p_error_B_mean,vmin = 0,vmax = 1)
plt.title("Error in approximating B binding mean")
plt.ylabel("$q_A$")
plt.xlabel("$q_B$")
plt.show()
sns.heatmap(p_error_A_mean,vmin = 0,vmax = 1)
plt.title("Error in approximating A binding mean")
plt.ylabel("$q_A$")
plt.xlabel("$q_B$")
plt.show()
sns.heatmap(p_error_B_var,vmin = 0,vmax = 1)
plt.title("Error in approximating B binding variance")
plt.ylabel("$q_A$")
plt.xlabel("$q_B$")
plt.show()
sns.heatmap(p_error_A_var,vmin = 0,vmax = 1)
plt.title("Error in approximating A binding variance")
plt.ylabel("$q_A$")
plt.xlabel("$q_B$")
plt.show()
