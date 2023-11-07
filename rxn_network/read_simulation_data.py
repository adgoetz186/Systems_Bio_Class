import json
import numpy as np


with open('Data/result_dict.json', 'r') as f:
    json_load = json.load(f)

a_restored = np.asarray(json_load)

#print(np.average(np.asarray(json_load['DNA_A_e']),axis=0))
print(np.average(np.asarray(json_load['DNA_A_f']),axis=0))
print(np.var(np.asarray(json_load['DNA_A_f']),axis=0))

print(np.average(np.asarray(json_load['DNA_B_f']),axis=0))
#print(np.average(np.asarray(json_load['DNA_B_e']),axis=0))
print(np.var(np.asarray(json_load['DNA_B_f']),axis=0))

