from sympy import Symbol, solve
from math import *
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
e=constants.value(u'elementary charge')
g=constants.value(u'standard acceleration of gravity')
rhoW=1000
rhoMercury=13690

x = Symbol("x")
y = Symbol("y")
z = Symbol("z")
Ke = (4*pi*constants.epsilon_0)**(-1)
epsilon_0=constants.epsilon_0
eMass=9.10938356 / 10**31
protonMass=1.6726219/10**27
# Permeability=4*pi/10**7
Permeability_air=constants.mu_0
SPEED_OF_LIGHT=constants.c
waterHeatOffusion=3.34*10**5
waterHeatCapacity=4.186*10**3
waterVapor=2.256*10**6
Kcopper=385
Ksteel=50.2
p0=constants.value(u'standard atmosphere')
gas_constant=constants.value('molar gas constant')
boltzmann = constants.value('Boltzmann constant')
Avogadro=6.022*10**23
def char191():
    n=2
    T0=30+273.15
    T1=110+273.15

    dT=110-30
    print(n*gas_constant*dT)
# char191()
def char186():
    T0=17+273.15
    m1=3.34/10**27
    m2=5.34/10**26
    print(solve((T0/m1-x/m2),x)[0]-273.15)
# char186()
def char185():
    h=3.5
    p=4.2*10**5
    h0=1
    p1=p/(4-2.9)*0.5
    h1=2.9
    print(solve((p+h*rhoW*g-p0-h0*rhoW*g-x**2/2*rhoW),x))
    print(solve((p1 + h1 * rhoW * g - p0 - h0 * rhoW * g - x ** 2 / 2 * rhoW), x))
    h3=2
    p3=p/2*0.5
    print(solve((p3 + h3 * rhoW * g - p0 - h0 * rhoW * g - x ** 2 / 2 * rhoW), x))
    print(solve(((4-h)*p-y*(4-x),y+x*rhoW*g-p0-h0*rhoW*g),x,y))
    print(solve(((4 - h) * p - y * (4 - h3), y + h3 * rhoW * g - p0 - h0 * rhoW * g-x**2/2*rhoW), x, y))
    print(sqrt(2/rhoW*(p3-p0)+2*g*(h3-h0)))
# char185()
def char184():
    p=1.5*p0
    T=94
    print(T-273.15)
    density=p/gas_constant/T*Avogadro
    print(density)
    density1=p0/gas_constant/(22+273.15)*Avogadro
    print(density1)

# char184()
def char183():
    M=18/1000
    cv=3*gas_constant/M
    print(cv)
# char183()
def char1810000():
    l=1.5
    r=45/100
    T=22+273.15
    p=21*p0
    M=32/1000
    V=r**2*l*pi
    n=p*V/gas_constant/T
    m=n*M
    print(n,m)
# char1810000()

def char181():
    l=1.6
    r=94/2/100
    T0=26+273.15
    p=21.5*p0
    m=32/10**3
    n=p*pi*r**2*l/gas_constant/T0
    print(n)
    print(n*m)
    print(gas_constant)
# char181()
def char1712():
    flowrate=0.6/60
    T=18
    V=120
    A=15
    print(solve((1800-waterHeatCapacity*flowrate*(x-T)),x))
# char1812()
def char1711():
    m0=0.44
    mice=0.095
    T0=0
    Ccopper=385
    msteam=0.035
    print(solve((-mice*waterHeatOffusion-mice*waterHeatCapacity*x-m0*Ccopper*x+msteam*waterVapor+waterHeatCapacity*(100-x)*msteam),x))
    print(mice+msteam)
# char1811()
def char1710():
    mice=29
    T=17
    T2=4
    k=0.0108
    a=0.55
    b=0.85
    h=0.55
    A=2 * (a * b + a * h + b * h)
    
    print(A)
    print(waterHeatCapacity*mice*6)
    print(waterHeatOffusion * mice)
    print(solve((waterHeatOffusion * mice + waterHeatCapacity * mice * T2 - 7 * 24 * 3600 * k * A * (T-T2) / x), x))

    print(solve((waterHeatOffusion * mice - y * k * 2 * (
                a * b + a * h + b * h) * T / x,waterHeatCapacity*mice*6-z*k*2*(a*b+a*h+b*h)*T/x,y+z-7*24*3600), x,y,z))


# char1810()
def char18100():
    mice=22
    k=0.0109
    a=0.55
    b=0.85
    h=0.5
    print(solve((waterHeatOffusion*mice+waterHeatCapacity*mice*6-7*24*3600*k*2*(a*b+a*h+b*h)*9/x),x))
    print(solve((waterHeatOffusion * mice - y * k * 2 * (
                a * b + a * h + b * h) * 15 / x,waterHeatCapacity*mice*6-z*k*2*(a*b+a*h+b*h)*9/x,y+z-7*24*3600), x,y,z))

def char179():
    n=3
    dT=200
    T0=27+273.15
    T1=227+273.15
    Q=(29.5*dT+8.2/10**3/2*(T1**2-T0**2))*n
    print(Q)
# char179()
# char18100()
def char178():
    m=250
    f=440
    dT=40
    # print(solve(f-1/2*sqrt(m*g/x),x))

    l=1
    copperLinearCo = 17 / 10 ** 6
    print(f*(1/sqrt(copperLinearCo*40+1)-1))
    print(sqrt(copperLinearCo*40+1))

    l1=l+l*copperLinearCo*dT

    mu1=m/l1
# char178()


def char177():
    stephanBoldzman=5.67/10**8
    T=2450
    e=0.35
    H=100
    A=H/e/T**4/stephanBoldzman
    print(A)
# char177()
def char176():
    l1=1
    A=4/10**4
    T=65
    print(solve((y-Kcopper*A*(100-65)/l1,y-Ksteel*A*(T)/x),x,y))
# char176()
def char175():
    dT=100
    l=0.45
    tgradiet=dT/l
    print(tgradiet)
    Kcopper=385
    A=1.25/10**4
    H=Kcopper*dT/l*A
    print(Kcopper*dT/l*A)
    dT1=0.12*H/Kcopper/A
    print(100-dT1)
# char175()

def char174():
    m=0.26
    dT=17
    Q=dT*m*waterHeatCapacity+m*waterHeatOffusion
    print(Q)
    print(Q/waterHeatCapacity)
    print(Q/1055)
# char174()

def char173():
    m=24000
    v=12.5
    a=65
    b=20
    h=12
    p=1.2
    CH=1020
    V=a*b*h
    E=v**2/2*m
    m_air=V*p
    dC=E/CH/m_air
    print(dC)
# char173()
def char1723():
    l=1.5
    print(solve((x*(420-20)*l-0.19),x))
    y=2*10**11
    a=3.2/10**5
    stress=y*a*400
    print(stress)
# char1723()
def char172():
    B=5.1/10**5
    print(solve((B*(x-20)-0.0015),x))
# char172()
def char171():
    d=4.5/10**3
    a=2.4/10**5
    print(d+a*d*(23+78))
# char171()
def lab1():
    d= 2.37/10**3
    dd=0.01/10**3
    l=[0.823,0.80,0.83,0.827,0.845,0.805,0.969,0.84,0.826]
    l1=[0.827,0.805,0.826]
    angle1=[55,66,35]
    delta_T1=[76,75,75]
    angle=[50,5,48,55,10,66,37.5,32,35]
    dl=0.01
    T0=[22,100]
    delta_T=[78,74,74,76,74,76,76,76,76,76]
    dT=2
    dangle=2
    a_theory=[19/10**6,24/10**6,13/10**6]
    percentage=[]
    a_ex=[]
    da=[]
    for x in range(len(l1)):
        a_ex.append(d/2*radians(angle1[x])/l1[x]/delta_T1[x])
        # a_theory.append(pi/360/l1[x]/delta_T1[x]*angle1[x]*d)
        da.append(a_ex[x]*(dd/d+dangle/angle1[x]+dl/l1[x]+dT*2/delta_T1[x]))
        percentage.append(abs(a_theory[x]-a_ex[x])/a_theory[x])

    for x in range(len(l1)):
        print(a_ex[x]-da[x],a_ex[x],a_ex[x]+da[x],a_theory[x])
    print(a_ex,a_theory)
    print(da,percentage)

# lab1()



def char1744():
    m0=0.2
    T1=-40
    T2=80
    T3=20
    c_water=4187
    c_ice=2108
    L=3.34*10**5
    print(solve((40*c_ice*m0+m0*L+20*c_water*m0-60*x*c_water),x))
# char1744()
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

