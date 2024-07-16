#Libraries
import numpy as np
import dwave
import neal
from pyqubo import Binary, Array
import pyqubo
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from random import sample
import os

#random generator to create the matrices
def r():return sample([0,1],1)[0]

#Randomly creating the binary matrix with given size x*y
def Matrix(x,y):
    M = []
    for i in range(y):
        vect = []
        for z in range(x):
            vect.append(r())
        M.append(np.array(vect))
    return np.array(M)



#Penalty function for QUBO following the Ising model for ferromagnetic system to simulate
def penalty(M):
    array = Array.create('x', shape=np.shape(M), vartype='BINARY')

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


#Function for knowing the minimal energy that can be obtained and then get the accuracy
def optimal(M):
    f=0
    for i in range(len(M)):
        for y in range(len(M[0])):
            f+= (1 - 2 * np.sum(M, axis=1)[i]*M[i,y])
            for z in range(y):
                f+= 2*M[i,y]*M[i,z]

    for i in range(len(M[0])):
        for y in range(len(M)):
            f+= (1 - 2 * np.sum(M, axis=0)[i])*M[y,i]
            for z in range(y):
                f+= 2*M[y,i]*M[z,i]
    return f


#Solving function using the Dwave Leap Computer     (TOKEN NEEDED)
def solve(M,Qcomputer):

    if Qcomputer:
        sampler = EmbeddingComposite(DWaveSampler(solver ={'qpu' : True}, token="YOUR-TOKEN"))
    else:sampler = neal.SimulatedAnnealingSampler()
    model = penalty(M).compile()
    bqm = model.to_bqm()
    sampleset = sampler.sample(bqm, num_reads =1000)
    decoded_samples = model.decode_sampleset(sampleset)
    best_sample = min(decoded_samples, key=lambda x: x.energy)


    en,oc = np.array([]),np.array([])
    for i in sampleset.record:
        if i[1] not in en:
            en = np.append(en,i[1])
            oc = np.append(oc,i[2])
        else:
            for y in range(len(en)):
                if en[y]==i[1]:
                    oc[y]+=i[2]

    return optimal(M),en,oc


#Saving the results in the file chosen
def save(x,y,Qcomputer):

    path = ''
    filename = ''
    file = open('{}/{}.txt'.format(path,filename), 'a')

    emin,en,oc = solve(Matrix(x,y),Qcomputer)

    file.write(str('matrix of size {}*{}  \n'.format(x,y)))
    file.write(str('minimal energy : '+str(emin)+'\n'+'resuts:'+'\n'))

    for i in range(len(en)):
        file.write(str('energy:  '+str(en[i])+'   occurence:  '+str(oc[i])+'\n'))

    file.write(str('\n\n\n\n\n'))



#saving the results for several matrices in the choosen file
for x in range(2,11):
    for y in range(x,11):
        save(x,y,Qcomputer = False)
        
        #choose Qcomputer = True to run on the quantum computer of leap, for the simulator take Qcomputer = False
print(sol)

