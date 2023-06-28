import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

n1 =int(input("Enter number of reactants: "))
n2 =int(input("Enter number of products: "))
li_stc=[]
li_cp=[]
li_c0=[] 
li_alpha = []
li_theta = []
for i in range(n1+n2):
    if i==0:
        li_stc.append(float(input("Enter coeff for limiting reagent with sign: ")))
        li_cp.append(float(input("Enter Cp value of limiting reagent: ")))
        li_c0.append(float(input("Enter initial concentration of limiting reagent: ")))
        if i<n1:
            li_alpha.append(float(input("Enter the order wrt this reactant: ")))
    else:
        li_stc.append(float(input("Enter stoichiometric coefficient with sign: ")))
        li_cp.append(float(input("Enter Cp value: ")))
        li_c0.append(float(input("Enter initial concentration: ")))

k0 =float(input("Enter the reference rate contsant wrt limiting reagent: "))
tr =float(input("Enter the reference temprature for k0: "))
V0 =float(input("Enter initial volumetric flow rate: "))
x = float(input("Enter the conversion: "))
phase =int(input("Enter 0 for liquid and 1 for gas: "))
tf = float(input("Enter feed temprature: "))
ea =float(input("Enter the activation energy of this reaction: "))
R =float(input("Enter the value of gas constant in corresponding unit: "))
delta_h =float(input("Enter the reference enthalpy of reaction: "))
trh = float(input("Enter the reference temprature for enthalpy: "))
n_li_stc = []
for i in range(0,n1+n2):
    n_li_stc.append(li_stc[i]/li_stc[0])
    li_theta.append(li_c0[i]/li_c0[0])


def inst_conc(n_li_stc, li_c0, x, li_theta,t):
    l= []
    if phase==0:
        for i in range(n1+n2):
            l.append(li_c0[i]*(li_theta[i]-x))
        return l
    else:
        delta = sum(li_theta)
        ct = sum(li_c0)
        eps = delta * li_c0[0]/ct
        z = tf/t
        for i in range(n1+n2):
            l.append(li_c0[0]*(li_theta[i]+(n_li_stc[i]*x))*z/(1-(eps*x)))
        return l  
    

def f1(v,x,T):
    z = inst_conc(n_li_stc,li_c0, x,li_theta,T)
    k = k0*math.exp((ea/R)*((1/tr)-(1/T)))
    fa0 = li_c0[0]*V0
    ans = k
    for i in range(0,n1):
        ans *= pow(z[i], li_alpha[i])
    #print("f1: ",ans,"  ",fa0)
    return ans/fa0

def f2(v,x,T):
    z = inst_conc(n_li_stc,li_c0, x,li_theta,T)
    k = k0*math.exp((ea/R)*((1/tr)-(1/T)))
    fa0 = li_c0[0]*V0
    ans = k
    ti_cpi =0 
    delta_cp =0
    for i in range(0,n1):
        ans *= pow(z[i], li_alpha[i])
    for i in range(0,n1+n2):
        ti_cpi += li_theta[i]*li_cp[i]
        delta_cp += li_stc[i]*li_cp[i]
        
    n_delta_h = delta_h + delta_cp*(T-tr)
    #print("f2: ",ans," ",n_delta_h," ",fa0," ",ti_cpi," ",x," ",delta_cp)
    return -ans*n_delta_h/(fa0*(ti_cpi+(x*delta_cp)))

def euler(f1,f2, v0, T0 ,x0, f_vol, h):
    v = np.arange(v0, f_vol+h, h)  # generate volume steps
    x = []  # allocate solution array
    T = []
    x.append(x0)  
    T.append(T0)
    for i in range(len(v)-1):
        k1 = f1(v[i], x[i],T[i])
        k1_ = f2(v[i], x[i],T[i])
        #k2 = f1(v[i] + 3/4*h, x[i] + 3/4*k1*h, T[i] + 3/4*k1_*h )
        #k2_ = f2(v[i] + 3/4*h, x[i] + 3/4*k1*h, T[i] + 3/4*k1_*h)
        #x[i+1] = x[i] + (1/3*k1 + 2/3*k2)1*h
        #T[i+1] = T[i] + (1/3*k1_ + 2/3*k2_)*h
        x.append(x[i] + k1*h) 
        T.append(T[i] +k1_*h) 
        if(x[i+1]>1 or x[i+1]<0):
            break   
    return (x,v,T)

(x,v,t) = euler(f1,f2, 0, tf ,0, 1000000, 0.01)
vf=[]
for i in range(len(x)):
    vf.append(v[i])
plt.plot(vf, t)
plt.xlabel('volume')
plt.ylabel('temperature')
plt.show()

plt.plot(vf,x)
plt.ylabel('conversion')
plt.xlabel('volume')
plt.show()


