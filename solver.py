import numpy as np
import dwave
import neal
from pyqubo import Binary, Array
import pyqubo
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite


"""
M = np.array([[0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],
              [1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,0],
              [1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0],
              [0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,1],
              [1,1,1,0,1,1,1,1,1,1,0,1,0,1,1,0,0],
              [0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0],
              [0,0,1,1,1,1,1,1,1,0,0,1,0,1,1,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
              [0,0,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0],
              [0,0,1,1,1,1,1,1,1,1,0,1,0,1,1,0,0],
              [0,0,1,1,1,0,1,1,0,1,0,0,0,0,0,0,0],
              [0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
              [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]])
"""

M = np.array([[1,0,1],[0,1,0],[1,0,1]])



array = Array.create('x', shape=np.shape(M), vartype='BINARY')


def penalty():
    f=0
    for i in range(len(M)):
        for y in range(len(M[0])):
            f+= (1 - 2 * np.sum(M, axis=1)[i]*array[i,y])
            for z in range(y):
                f+= 2*array[i,y]*array[i,z]

    for i in range(len(M[0])):
        for y in range(len(M)):
            f+= (1 - 2 * np.sum(M, axis=0)[i])*array[y,i]
            for z in range(y):
                f+= 2*array[y,i]*array[z,i]
    return f


sampler = EmbeddingComposite(DWaveSampler(solver ={'qpu' : True}))
#sampler = neal.SimulatedAnnealingSampler()
model = penalty().compile()
bqm = model.to_bqm()
sampleset = sampler.sample(bqm, num_reads =1000)
decoded_samples = model.decode_sampleset(sampleset)
best_sample = min(decoded_samples, key=lambda x: x.energy)

print(sampleset)

sol = np.zeros((len(M), len(M[0])))

res = best_sample.sample.items()

print('\n\n\n')

print('best sample :', best_sample,'\n')


for item in res:
    for i in range(len(sol)):
        for y in range(len(sol[0])):
            if int(item[0][2])==i and int(item[0][5]) == y :
                sol[i,y] = item[1]

print(sol)

print('\n\n\n')


print('worst sample :', max(decoded_samples, key=lambda x: x.energy),'·\n')

for item in max(decoded_samples, key=lambda x: x.energy).sample.items():
    for i in range(len(sol)):
        for y in range(len(sol[0])):
            if int(item[0][2])==i and int(item[0][5]) == y :
                sol[i,y] = item[1]

print(sol)

