import matplotlib.pyplot as plt
import numpy as np
from qutip import *

def fHamiltonian(N,w1,w2,w3,alpha1,alpha2,alpha3,g13,g12,g23):
    
    # N is the demension of each qubit,w123 denotes the frequencies,alpha123 denotes the anharmonicity,
    # g denotes XY coupling strength, dont need to multiply 2pi
    
    a = destroy(N)
    a_dagger = create(N)  
    
    H_w1= w1 * tensor(a_dagger * a, qeye(N), qeye(N)) 
    H_a1= alpha1 / 2 * tensor(a_dagger * a_dagger * a * a, qeye(N), qeye(N))
    
    H_w2=  w2 * tensor( qeye(N),a_dagger * a, qeye(N)) 
    H_a2= alpha2 / 2 * tensor( qeye(N),a_dagger * a_dagger * a * a, qeye(N))
    
    H_w3=  w3 * tensor( qeye(N), qeye(N),a_dagger * a) 
    H_a3= alpha3 / 2 * tensor( qeye(N), qeye(N),a_dagger * a_dagger * a * a)
    
    H_13 = g13 * tensor(a+a_dagger,qeye(N),a+a_dagger)
    H_12 = g12 * tensor(a + a_dagger, a + a_dagger, qeye(N))
    H_23 = g23 * tensor(qeye(N), a + a_dagger, a + a_dagger)

    H_com= (H_w1+H_a1+H_w2+H_a2+H_w3+H_a3+H_12+H_23+H_13)*np.pi*2
    
    return H_com

def feigen(H,n1,n2,n3,N): # Give eigenstate. 
    
    # N is the demension of each qubit,n1 n2 n3 denotes the photon number on qubit123,H is the Hamiltonian
    
    eigenlist=[np.abs(arr[n1*N*N+n2*N+n3][0][0].real) for arr in H.eigenstates()[1]]
    max_value=max(eigenlist)
    max_index = None
                  
    for index, value in enumerate(eigenlist):
        if value == max_value:
            max_index = index
            break
            
            ### output: state in qobj, energylevel magnitude, index number of the energylevel
            
    return [H.eigenstates()[1][max_index],H.eigenenergies()[max_index]/np.pi/2,max_index]
            

    
def HXYpulse(N,H_static,tg,t1,qn,ampl,freq):
    
    ampl=ampl*2*np.pi
    
    def mpulse(t, args): 
        if t<t1:
            return ampl * np.cos(freq*np.pi*2*t)* (1/2-1/2*np.cos(np.pi*t/t1))
        if t1<=t<=tg-t1:
            return ampl * np.cos(freq*np.pi*2*t)
        if tg-t1<t:
            return ampl * np.cos(freq*np.pi*2*t)* (1/2-1/2*np.cos(np.pi*(t-tg+2*t1)/t1))
        
    a = destroy(N)
    a_dagger = create(N)  
    
    # To which qubit 
    if qn==1:
        HXY=[H_static,[tensor(a+a_dagger,qeye(N), qeye(N)), mpulse]]
    if qn==2:
        HXY=[H_static,[tensor(qeye(N),a+a_dagger, qeye(N)), mpulse]]
    if qn==3:
        HXY=[H_static,[tensor(qeye(N), qeye(N),a+a_dagger), mpulse]]
        
    return HXY

def HZpulse(N,H_static,tg,qn,ampl):
    
    ### half sin shape

    def fpulse(t, args): 
        return ampl * np.sin(np.pi*t/tg)
    
    a = destroy(N)
    a_dagger = create(N)  
    
    if qn==1:
        HZ=[H_static,[tensor(a_dagger*a,qeye(N), qeye(N)), fpulse]]
    if qn==2:
        HZ=[H_static,[tensor(qeye(N),a_dagger*a, qeye(N)), fpulse]]
    if qn==3:
        HZ=[H_static,[tensor(qeye(N), qeye(N),a_dagger*a), fpulse]]

    return HZ


def coplist(N,G1up,G1down,G1Z,G2up,G2down,G2Z,G3up,G3down,G3Z):
    
    ao= destroy(N)
    ad= create(N)
    # ??
    cop1up = np.sqrt(G1up)*tensor(ad,qeye(N),qeye(N))
    cop1down = np.sqrt(G1down)*tensor(ao,qeye(N),qeye(N))
    cop1Z = np.sqrt(G1Z)*tensor(ad*ao,qeye(N),qeye(N))
    cop2up = np.sqrt(G2up)*tensor(qeye(N),ad,qeye(N))
    cop2down = np.sqrt(G2down)*tensor(qeye(N),ao,qeye(N))
    cop2Z = np.sqrt(G2Z)*tensor(qeye(N),qeye(N),ad*ao)
    cop3up = np.sqrt(G3up)*tensor(qeye(N),qeye(N),ad)
    cop3down = np.sqrt(G3down)*tensor(qeye(N),qeye(N),ao)
    cop3Z = np.sqrt(G3Z)*tensor(qeye(N),qeye(N),ad*ao)
    
    return [cop1up,cop1down,cop1Z,cop2up,cop2down,cop2Z,cop3up,cop3down,cop3Z]


    
    
def seoutput(H_d,tg,step,state):
    tlist=np.linspace(0, tg, step)
 
    option = Options(rtol=1e-8)
        
    output = mesolve(H_d, state, tlist, [], [], options=option)
    
    # result1 = np.angle(state.dag() * output.states[-1])
        
    return output

def meoutput(H_d,tg,step,state,coplist):
    tlist=np.linspace(0, tg, step)
 
    option = Options(rtol=1e-8)
        
    output = mesolve(H_d, state, tlist, coplist, [], options=option)

        
    return output



def dlist(output,step,state2):

    
    datalist=[]
    for ii in range(0,step):
        datalist.append(np.abs(((output.states[ii]).dag()*state2)[0][0][0])**2)
        
    return datalist


def dlist_density(output,step,state2):

    
    datalist=[]
    for ii in range(0,step):
        datalist.append(((output.states[ii] * output.states[ii].dag()) * (state2 * state2.dag())).tr())
        
    return datalist
    

