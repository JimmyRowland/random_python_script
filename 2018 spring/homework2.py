from sympy import Symbol, solve
from math import *
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
e=constants.value(u'elementary charge')
g=constants.value(u'standard acceleration of gravity')
p0=constants.value(u'standard atmosphere')
RRR=constants.value(u'molar gas constant')
rhoW=1000
rhoMercury=13690
Avogadro=6.022*10**23
x = Symbol("x")
y = Symbol("y")
Ke = (4*pi*constants.epsilon_0)**(-1)
epsilon_0=constants.epsilon_0
eMass=9.10938356 / 10**31
protonMass=1.6726219/10**27
# Permeability=4*pi/10**7
Permeability_air=constants.mu_0
SPEED_OF_LIGHT=constants.c

def char1744():
    m0=0.2
    T1=-40
    T2=80
    T3=20
    c_water=4187
    c_ice=2108
    L=3.34*10**5
    print(solve((40*c_ice*m0+m0*L+20*c_water*m0-60*x*c_water),x))
char1744()
def lab2():
    f1=427
    df1=f1*0.005
    f2=384
    df2=f2*0.005
    v0=[331.4,259]
    const=[0.6,0.4]
    d1=[58,59,18,19]
    d2=[20,21,64,65]
    d3=[68,69,43,44]
    d4=[37,38,61,63]
    d3_cooked=[71,72,43,44]
    d4_cooked=[45,46,76,77]
    temp=[24, 24, 18, 18]
    temp_original = [24.1, 23.9, 21, 18.2]
    d=[d1,d2,d3_cooked,d4_cooked]
    l=[]
    v=[]
    dv=[]
    dv_theo=[]
    l_theo=[]
    v_theory=[]
    for x in d:
        l.append(abs((x[1]+x[0])/2-(x[3]+x[2])/2)/100)
    for y in range(0,2):
        v.append(l[y*2]*2*f1)
        print(l[y*2],f1,l[y * 2+1],f2)
        v.append(l[y * 2+1] * 2 * f2)
        dv.append(v[y * 2]*(df1/f1+(0.04)/l[y * 2] ))
        dv.append(v[y * 2+1] * (df2 / f2 + (0.04) / l[y * 2+1]))
    for z in range(0,2):
        v_theory.append(v0[z]+temp[z*2]*const[z])
        v_theory.append(v0[z] + temp[z*2+1]*const[z])
        l_theo.append(v_theory[z*2]/2/f1)
        l_theo.append(v_theory[z*2+1] / 2 / f2)
        dv_theo.append(const[z]*2)
        dv_theo.append(const[z]*2)
    error_percentage = []
    for z in range (0,4):
        error_percentage.append(abs((v[z]-v_theory[z])/v_theory[z]))
        print(v[z]-dv[z],v[z],v[z]+dv[z],v_theory[z])
    # # plt.errorbar(range(10),y=[v[0]],yerr=error_percentage[0])
    # # example data
    # # x = np.arange(0.1, 4, 0.5)
    # # y = np.exp(-x)
    # #
    # # # example variable error bar values
    # # yerr = 0.1 + 0.2 * np.sqrt(x)
    # # xerr = 0.1 + yerr
    #
    # # First illustrate basic pyplot interface, using defaults where possible.
    # u = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    #
    # plt.figure()
    # # plt.errorbar(x, y, xerr=0.2, yerr=0.4)
    # plt.errorbar([v[0]], [f1], xerr=error_percentage[0], yerr=0.005, fmt='r^')
    # plt.errorbar([v[1]], [f2], xerr=error_percentage[1], yerr=0.005,fmt='r^')
    # # plt.errorbar(xdata = [v[1]],xerr=error_percentage[1], fmt='r^')
    # plt.savefig("graph.svg")
    # plt.show()
    def v_the(gamma, temp, M):
        R=8.314
        v=sqrt(gamma*R*(temp+273.15)/M*1000)
        print(v)
    v_the(1.4,24,28.8)
    v_the(1.28, 18, 44)


    print(v,v_theory)
    print(l,l_theo)
    print(dv,dv_theo)
    print(df1,df2)
    print(error_percentage)
# lab2()
def char_16_10():
    T=20+273
    phia=4/v*w
    phib=3/v*w
    print((phia-phib)*180/pi)
    print(1/wl*2*180)
    Ab=4*pi*rb**2
    Ib=Pb/Ab
    Lb=10*log(Ib*10**12,10)
    Aa = 4 * pi * ra ** 2
    Ia = Pa / Aa
    La = 10 * log(Ia * 10 ** 12, 10)
    print(Ib,Lb)
    print((sqrt(Ia)-sqrt(Ib))**2)
    print(10 * log((sqrt(Ia)-sqrt(Ib))**2 * 10 ** 12, 10))
    print(Ia,La)

# char_16_10()
def char_16_9():
    v=344
    fs=80.8*1000
    f=83.7*1000
    vb=4.3
    vs=v-fs/f*(v+vb)

    print(solve(((v+vb)/(v+x)*y-f,(v-x)/(v-vb)*fs-y),x,y))
# char_16_9()
def char_16_8():
    T=77+273
    f=500
    l1=0.18
    l2=0.555
    l3=0.93
    wl=l3-l1
    v=wl*f
    R=8.314
    M=28.8/10**3
    r=v**2*M/T/R
    print(l2/l1,l3/l1)
    print(v)
    print(r)
    print(wl/4-l1)
# char_16_8()
def char_16_7():
    l=0.1075
    f=600
    v=344
    print(v/4/l)
# char_16_7()
def char_16_5():
    fs=520
    fl=490
    v=344
    vl=v-fl/fs*v
    print(vl)
# char_16_5()
def char_16_4():
    va=0
    vl=15
    vb=35
    v=344
    f=392
    fla=(v-vl)/v*f
    print(fla)
    flb=(v+vl)/(v+vb)*f
    print(flb)

# char_16_4()
def char_16_3():
    v=344
    l=0.14
    wl=4*l
    f=v/wl
    print(f)

# char_16_3()
def char_16_2():
    d1=0.35
    d2=1.6
    difference=10*log((d2**2/d1**2),10)
    print(difference)
# char_16_2()
def char_16_1():
    temp=20
    f=147
    A=5.5/10**6
    B=1.42*10**5
    v=343
    wl=v/f
    pmax=B*2*pi/wl*A
    print(pmax)
    rho=1.2
    intensity = sqrt(rho*B)*(2*pi*f)**2*A**2/2
    print(intensity)
    intensitylevel=10*log(intensity*10**12,10)
    print(intensitylevel)
# char_16_1()
def char_15_9():
    L=0.54
    r=1.18/2/1000
    f=311
    rho=7800
    A=r**2*pi
    V=A*L
    m=rho*V
    linearDensity=m/L
    v=2*f*L
    tension=4*L**2*f**2*pi*r**2*rho
    T=rho*pi*r**2*f**2*L**2
    print(tension,T)
# char_15_9()
def char_15_6():
    m=42/1000
    l=0.75
    f=64
    A=0.32/100
    linearDensity=m/l
    # F=m*g
    # v=sqrt(F/linearDensity)
    w=2*pi*f
    v=f*2*l
    tension=v**2*linearDensity
    vmax=w*A
    amax=w**2*A
    print(v,tension,vmax,amax)

# char_15_6()
def char_15_5():
    I1=9.2
    I2=1/10**6
    r1=30.2
    r2=sqrt(I1/I2*r1**2)
    print(r2)
    print(r2/4)
    P=4*pi*r2**2*I2
    print(P)
# char_15_5()
def char_15_4():
    f=40
    A=3/100
    linearDensity=5/100
    tension=5
    v=sqrt(tension/linearDensity)
    wl=v/f
    print(v,wl)
    w=2*pi*f
    amax=w**2*A
    print(amax)
# char_15_4()
def char_15_3():
    f=111
    m=1.55
    linearDensity=5.4/100
    T=m*g
    v=sqrt(T/linearDensity)
    print(v)
    wl=v/f
    print(wl)
# char_15_3()
def char_15_2():
    T=3.7/100
    f=1/T
    wl=0.29
    v=wl/T
    print(v)
    print(f)
    print()
# char_15_2()
def char_15_1():
    wl=800*1000
    T=1*3600
    v=wl/T
    V2=v*3.6
    print(v,V2)
# char_15_1()

def char_12_13():
    h1=8
    A2=0.048
    A3=0.016
    print(solve((g*h1-x**2/2),x))
    v3=12.5262284826679
    discharge = 12.5262284826679*A3
    print(discharge)
    print(solve(y * A2 - v3 * A3,),y)
    v2=4.17540949422262
    print(solve(x+v2**2/2*rhoW-p0-v3**2/2*rhoW),x)
    print(171061.177777778-p0)

# char_12_13()

def char_12_11():
    alDensity=2700
    w=48
    l=40
    buoyant=w-l
    v=buoyant/rhoW/g
    auDensity=19320
    print(solve((x+y-v,x*auDensity*g+y*alDensity*g-48),x,y))
    print(0.000161976582566891*auDensity*g)
# char_12_11()


def char_12_10():
    a=0.56
    b=0.21
    h=0.086
    density=700
    rhoLead=11340
    Vwood=a*b*h
    mwood=Vwood*density
    wwood=mwood*g
    print(solve(((x+Vwood)*rhoW-Vwood*density-rhoLead*x),x))
    vlead=0.000293431334622824
    m=vlead*rhoLead
    print(m)
# char_12_10()
def char_12_9():
    hw=0.15
    p1=rhoW*g*hw
    print(p1)
    hmer=p1/g/rhoMercury
    print(hmer)
    print(hw-hmer)
# char_12_9()
def char_12_7():
    b=5.9
    a=4.2
    h=2.8
    V=a*b*h
    W=V*g*rhoW
    A=a*h
    F2=A*rhoW*g*h/2
    print(W,F2)
# char_12_7()
def char_12_6():
    r1=0.22
    A1=r1**2*pi
    flow=1.7
    v1=flow/A1
    v2=2.55
    A2=flow/v2
    r2=sqrt(A2/pi)
    print(v1,r2)
# char_12_6()
def char_12_5():
    a=0.1
    rho_o=790
    A=a*a
    h1=0.025
    h2=h1+a
    P1=h1*rho_o*g
    P2=a*rho_o*g+h1*rhoW*g
    print(P1,P2)
    m=(P2-P1)*A/g
    density=m/A/a
    print(density)
    print(m)
# char_12_5()
def char_12_4():
    w=16.1
    tension=12.1
    buoyant=w-tension
    V=buoyant/rhoW/g
    density=w/g/V
    print(V,density)
# char_12_4()
def char_12_3():
    w=45
    r=0.15
    rho=850
    h=0.75
    A=r**2*pi
    p1=w/A
    w1=92+w
    p2 = w1 / A
    p_delta=p2-p1

    print(p1,p_delta)
# char_12_3()
def char_12_2():
    h=6.1
    p=g*h*rhoW
    print(p)
# char_12_2()
def chq4_12_1():
    h=1300
    g1=3.71
    rho=1000
    p=h*g1*rho

    h1=(p-p0)/g/rho
    print(p,h1,g)
    h2=500
    p=h2*g1*rhoW
    h3=p/rhoW/g
    print(p,h3)
# chq4_12_1()

