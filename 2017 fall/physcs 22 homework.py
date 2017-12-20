from sympy import Symbol, solve
from math import *
import numpy
from scipy import constants
e=constants.value(u'elementary charge')
g=9.8
Avogadro=6.022*10**23
x = Symbol("x")
Ke = (4*pi*constants.epsilon_0)**(-1)
epsilon_0=constants.epsilon_0
eMass=9.10938356 / 10**31
protonMass=1.6726219/10**27
# Permeability=4*pi/10**7
Permeability_air=constants.mu_0
SPEED_OF_LIGHT=constants.c
def Eforce(charge1,charge2,distance):
    return Ke*charge1*charge2/distance**2
def EforceMicro(charge1,charge2,distance):
    return Ke*charge1*charge2/10**12/distance**2
def EforceNano(charge1,charge2,distance):
    return Ke*charge1*charge2/10**18/distance**2
def potential(q,d):
    return Ke*q/d
def C_series(*args):
    C=0
    for x in args:
        C+=1/x
    return 1/C
def R_parellel(*args):
    C = 0
    for x in args:
        C += 1 / x
    return 1 / C
def Q_magnetic_field(q,v,r,angle=90):
    return Permeability_air/4/pi*q*v*sin(radians(angle))/r**2
def char_31_6():
    V=120
    P=20
    I=P/V
    print(I)
    R=V/I
    print(R)
char_31_6()
def char_31_5():
    f=120
    w=2*pi*f
    emf=87
    R=74
    L=115
    Z=115
    angle=atan(sqrt(Z**2-R**2)/R)
    I=emf/Z
    P=I*emf*cos(angle)
    print(P)
# char_31_5()
def char_31_4():
    R=175
    C=12.5/10**6
    L=8/10**3
    emf=25
    w_min=sqrt(1/L/C)
    print(w_min)
    Xl=w_min*L
    Xc=1/w_min/C
    Z=sqrt(R**2+(Xl-Xc)**2)
    Imax=emf/Z
    angle=atan((Xl-Xc)/R)
    V=emf*cos(pi/3+angle)
    i1=Imax/2
    V1=i1*Z
    VR=Imax*R/2
    Vl=Imax*Xl*cos(pi/2+pi/3)
    Vc=Imax*Xc*cos(-pi/2+pi/3)


    print(V1,V)
    print(V1,VR,Vl,Vc)
    print(Imax)
# char_31_4()

def char_31_3():
    R=160
    L=0.38
    emf=31
    w=270
    Xl=L*w
    Z=sqrt(R**2+Xl**2)
    print(Z)
    i=emf/Z
    print(i)
    VR=i*R
    print(VR)
    Vl=i*Xl
    print(Vl)
    angle=atan(Xl/R)
    print(degrees(angle))

# char_31_3()
def char_31_2():
    c=6/10**6
    L=0.4
    emf=30
    w=250
    Xl=w*L
    Xc=1/w/c
    Z=w*L-1/w/c
    print(Z)
    i=emf/Z
    print(i)
    Vl=Xl*i
    Vc=Xc*i
    print(Vl,Vc)

# char_31_2()
def char_30_13():
    L=65/1000
    C=250/10**6
    q0=6.5/10**6
    V0=q0/C
    print(V0)
    w=sqrt(1/L/C)
    Imax=w*q0
    print(Imax)
    U=V0**2*C/2
    print(U)
    q1=sqrt(3)/2*q0
    print(q1)
    U1=U-q1**2/C/2
    print(U1)
# char_30_13()
def char_30_12():
    C=15/10**9
    L=26/1000
    R=74
    w=sqrt(1/L/C-R**2/4/L**2)
    f=w/2/pi
    print(f)
    t01=log(0.1)/R*2*L
    print(t01)
    R_cri=sqrt(4*L/C)
    print(R_cri)
# char_30_12()
def char_30_11():
    emf=62
    R=220
    L=0.17
    I1=emf/R
    t=4.5/10000
    It=I1*2.7**(-t*R/L)
    print(I1)
    print(It)
    Vbc=-It*R
    print(Vbc)
    I2=I1/2
    t=log(I2/I1)*L/R
    print(t)
# char_30_11()
def char_30_10():
    L=2.3
    R=8
    V=6
    didt0=V/L
    I=0.55
    emf=V-I*R
    didt1=emf/L

    print(didt0)
    print(didt1)
    I1=V/R
    t=0.210
    It=I1*(1-2.7**(-R/L*t))
    print(I1)
    print(It)
# char_30_10()
def char_30_9():
    R1=100
    R4=50
    R3=75
    emf=50
    V1=50/3*2
    print(V1)
    V2=V1/2
    print(V2)
    I1=V1/R1
    print(I1)
    I3=V2/R4
    print(I3)
    R=R1+R_parellel(R3,R4)
    I11=emf/R
    V11=I11*R1
    print(V11,I11)
    V14=emf-V11
    I13=V14/R4
    print(V14,I13)
    I2=V14/R3
    print(I2)

# char_30_9()
def char_30_8():
    R=15
    L=12/1000
    emf=25
    I=emf/R
    print(I)
# char_30_8()
def char_30_7():
    emf=60
    R1=45
    R2=27
    L=0.308
    VR1=emf/R2*R1
    print(VR1)
    Vab=emf/R2*(R1+R2)
    print(Vab)
# char_30_7()
def char_30_6():
    P=205
    U=P*24*3600
    print(U)
    I=85
    L=2*U/I**2
    print(L)
# char_30_6()
def char_30_5():
    L=0.28
    didt=-1.5/100
    emf=-L*didt
    print(emf)
# char_30_5()
def char_30_4():
    r=15.5/100
    A=5.05/10000
    I=12.5
    E=0.395
    L=2*E/I**2
    n=sqrt(2*pi*r*L/Permeability_air/A)
    print(n)
# char_30_4()
def char_30_3():
    didt=6.45/100
    emf=1.64/100
    L=emf/didt
    print(L)
    n=400
    I=0.725
    flux=I*L/n
    print(flux)
# char_30_3()
def char_30_2():
    n1=690
    n2=420
    I1=6.5
    flux2=3.5/100
    M=flux2*n2/I1
    I2=2.55
    flux1=I2*M/n1

    print(M)
    print(flux1)
# char_30_2()
def char_30_1():
    M=3.3/10**4
    didt1=850
    emf2=didt1*M
    print(emf2)
# char_30_1()
def char_30_32():
    C=20/10**6
    emf=150
    L=0.28/10**3
    w0=sqrt(1/L/C)
    T=2*pi/w0
    f=1/T
    print(f)
    energy=emf**2*C/2
    print(energy)
    q0=C*emf
    I1=q0*w0*sin(w0*1.3/10**3)
    print(I1)
    energyL=I1**2*L/2
    print(energyL)
# char_30_32()
def char_30_39():
    L=0.45
    C=2.5/10**5
    R=0
    w0=sqrt(1/L/C)
    print(w0)
    w1=0.95*w0
    R1=sqrt(4*L**2*(1/L/C-w1**2))
    print(R1)
# char_30_39()
def char_29_8():
    r=6.3/2/100
    B=0.98
    a=pi*r**2
    t=0.27
    emf=a*B/t
    print(emf)
# char_29_8()
def char_29_7():
    C=20/10**6
    V=100
    n=25
    R0=10

    asmall=0.1
    bsmall=0.2
    c=0.05
    Asmall=asmall*bsmall
    n=25
    abig=2
    bbig=4
    t=200/10**6
    i=c*V/R0/c*2.7**(-t/R0/C)
    diOverDt=10/t*2.7**(-1)
    R1 = 1 * (2 * (asmall + bsmall)) * n
    emf = n*Permeability_air / 2 / pi * bsmall * diOverDt * log((asmall+c) / c)
    induce=emf/R1
    print(diOverDt)
    print(10/t*2.71828**(-1))
    print()
    print(emf)
    print(induce)
# char_29_7()


def char_29_6():
    B=0.8
    v0=7.5
    ab=0.5
    emf0=ab*v0*B
    print(emf0)
    R=1.5
    I=emf0/R
    F=I*ab*B
    print(F)
# char_29_6()
def char_29_3():
    a=0.3
    b=0.4
    B=1.2
    v=0.02
    emp1=b*v*B
    print(emp1)
# char_29_3()
def char_29_2():
    a=0.12
    b=0.36
    L=0.24
    dioverdt=9.6
    emf=Permeability_air/2/pi*L*dioverdt*log(b/a)
    print(emf)
# char_29_2()
def char_29_1():
    n=52
    a=0.25
    b=0.31
    B=1.4
    A=a*b
    phi=A*B*n
    print(phi)
    t=0.24
    E=phi/t
    print(E)
# char_29_1()
def char_28_14():


    l=4/100
    r=sin(radians(6))*l*2
    m = 1.2 / 100
    F = tan(radians(6)) * m * g
    I=sqrt(2*pi*r*F/Permeability_air)
    print(I)
# char_28_14()
def char_28_13():
    r=25/100
    I=2.3
    x=27/100/2

    B=2*Permeability_air * I  * r ** 2 / 2 / (x ** 2 + r ** 2)** (3 / 2)
    v=2800
    F=B*v*e
    print(B)
    print(F)
# char_28_13()

def char_28_12():
    B=3.6/10**2
    I=14.5
    L=42/100
    n=B/Permeability_air/I
    d=2.2/100*2
    length=pi*d*n*L
    print(n)
    print(length)
# char_28_12()
def char_28_11():
    line=3.01/10**4
    I=line/Permeability_air
    print(I)
# char_28_11()
def char_28_10():
    r=3.2/100/2
    n=620
    I=0.44
    print(Permeability_air*I*n/2/r)
    x=6.5/100
    print(Permeability_air * I * n*r**2 / 2 / (x**2+r**2)**(3/2))
# char_28_10()
def char_28_7():
    d=10/100
    I=4
    r1=20/100
    r2=30/100
    r=sqrt((5/100)**2+(20/100)**2)
    R1=5/100
    R2=20/100
    R3=sqrt(R1**2+R2**2)

    print(r)
    print(Permeability_air/2/pi*I*(1/r1+1/r2))
    print(Permeability_air/2/pi*I/(20/100)*2/r*(5/100))
    print(Permeability_air / 2 / pi * I / R3 *R1/R3*2)
# char_28_7()
def char_28_6():
    nvA=3.2*10**18
    r=4/100
    I=nvA*e
    print(Permeability_air/2/pi*I/r)
# char_28_6()
def practice_29_2():
    n=200
    a=12/10000
    t=0.04
    B=6/10**5
    flux0=n*a*B
    flux1=0
    Energy=n*flux0/t

    print(flux0,Energy)
# practice_29_2()
def char_28_5():
    I=29
    angle=pi/4
    dl=2/1000
    r=sqrt((1.5/100)**2*2)
    print(r)
    print(2*Permeability_air/4/pi*I*dl/r**2*sin(angle))
# char_28_5()
def char_28_4():
    dl=1.1/1000
    I=10
    r=5/100
    angle=atan(5/14)
    rb=r/sin(angle)
    print(rb)

    print(Permeability_air/4/pi*I*dl/r**2)
    print(Permeability_air / 4 / pi * I * dl / rb ** 2*sin(angle))
# char_28_4()
def char_28_3():
    v=3*10**5
    print(epsilon_0*Permeability_air*v**2)
# char_28_3()
def char_28_2():
    v=0.15*SPEED_OF_LIGHT
    d=2/10**6
    A_angle=30
    C_angle=90
    B_angle=150
    D_angle=180
    print(Q_magnetic_field(e,v,d,A_angle))
    print(Q_magnetic_field(e, v, d, B_angle))
    print(Q_magnetic_field(e, v, d, C_angle))
    print(Q_magnetic_field(e, v, d, D_angle))
# char_28_2()
def char_28_1():
    r=5.3/10**11
    v=2.2*10**6
    F= v**2/r*eMass
    E=Ke*e/r**2
    print(E)
    B=Permeability_air/4/pi*e*v/r**2
    print(B)
# char_28_1()
def char_27_13():

    a=6/100
    b=8/100
    m = 0.16 / 1000 *(a+b)*2*100
    print(m,0.12*2*14)
    ma=0.16 / 1000 *a*100
    mb=0.16 / 1000 *b*100
    I=9
    phi=a*b*I
    torque=m*g/2*b*sin(1/6*pi)
    B=torque/phi/sin(1/3*pi)
    print(B)
# char_27_13()
def char_27_11():
    PR=0.8
    B=3
    I=5
    F=B*I*PR
    print(F)
# char_27_11()

def char_27_10():
    print(2*0.835*1.5)
# char_27_10()
def char_27_70():
    weight=2.6
    l=1.5
    R=10
    B=1.6
    V=120
    R1=25
    R2=10
    R0=R1+R_parellel(R1,R)
    I=V/R0
    I1=I/2
    F=B*l*I
    m=weight/g
    a=F/m
    print(a)
# char_27_70()
def char_27_9():
    a=1
    b=0.5
    m=0.21
    I=2.7
    Imoment=(1/12*b**2*m+1/12/6*m*b**2)*2
    B=3.4
    A=a*b
    torque=B*A*I
    acc=torque/Imoment
    print(acc)

# char_27_9()
def char_27_8():
    print(45/100*6*0.666)
# char_27_8()
def char_27_7():
    l=3-0.95*2
    I=3.4
    B=2.2
    Fx1=I*B*l
    Fx2=I*B*3
    print(Fx2)
# char_27_7()

def char_27_6():
    m=5.6/1000
    q=2300/10**6
    v0=11
    E=25
    F=q*E+m*g
    B=F/q/v0
    print(B)
# char_27_6()
def char_27_5():
    A=28.5/10**6
    V=120
    d=8/1000
    q=2*e
    m=6.64/10**27
    V0=1600
    v0=sqrt(2*V0*q/m)
    C=epsilon_0*A/d
    E=V/d
    F=E*q
    B=F/q/v0
    print(B)
# char_27_5()
def char_27_4():
    B=0.42
    R=29
    l=0.49
    m=0.7
    R=29
    F=m*g
    I=F/B/l
    V=I*R
    print(V)
    R1=2.2
    I1=V/R1
    F1=I1*l*B
    F_all=F1-F
    a=F_all/m
    print(a)
# char_27_4()

def char_27_2():
    f=4.1/10**15
    angle=radians(63)
    B=4/10**3
    v=f/B/e/sin(angle)
    print((v))
# char_27_2()
def char_27_1():
    q=-1.5/10**8
    Vx=4.5
    Vy=-3.8
    V=sqrt(Vx**2+Vy**2)
    B=1.7
    print(-B*Vy*q)
    print(-B*Vx*q)
    print(-B * Vy * q)
# char_27_1()

def char_27_46():
    a=22/100
    b=35/100
    B=1.5
    I=1.4
    A = a * b * I
    torque=A*B*sin(radians(30))
    F=a*I*sin(radians(30))
    print(torque,F)
# char_27_46()
def char27_15():
    v=1.6*10**6
    r=10/2/100
    B=v*eMass/r/e
    F=v*e*B
    t=pi*eMass/B/e
    t2=pi*r/v

    print(B,t,t2)
# char27_15()
def lab_7():
    # print(Permeability_air,Permeability)
    def returnRawData(DList,currentList,voltage):
        rawdata={'r':[d/2/100 for d in DList],"I":currentList,"V":voltage}
        # print(rawdata)
        return rawdata
    raw_data=[returnRawData([11,9.5,9,8,7,6,5],[0.94,1.04,1.24,1.42,1.68,2.00,2.4],200),returnRawData([11,10,9,8,7,6,5],[1.3,1.48,1.68,1.9,2.16,2.52,3],300),returnRawData([11,10,9,8,7,6],[1.6,1.8,2,2.24,2.54,2.96],400)]
    theoretical=e/eMass
    a=0.33
    N=74
    for data in raw_data:
        # print(raw_data)
        # print(data)
        V=data['V']
        result_list=[]
        for r,I in zip(data['r'],data['I']):
            result=(3.91*a**2/Permeability_air**2/N**2)*(V/I**2/r**2)
            result_list.append(result)
            print(r,I,V,result,theoretical)
        print(sum(result_list)/len(result_list),theoretical)
# lab_7()
def char26_9():
    M=860
    N=15
    P=33.48
    X=M*P/N
    print(X)
# char26_9()
def char26_8():
    I1 = Symbol("I1")
    I2 = Symbol("I2")
    I3 = Symbol("I3")
    print(solve((-I1+9-(I1+I3)*8,-I2+12+(I3-I2)*5,-10*I3-9+I1-I2+12),I1,I2,I3))

# char26_8()
def char26_7():
    R0=R_parellel(50,25)+75+15
    E=100
    I0=E/R0
    R1=25+75+15+25+25
    I1=E/R1
    print(I0,I1,R1,R0)
# char26_7()
def char26_6():
    C1=15/10**6
    C2=20/10**6
    C=C1+C2
    R=80
    print(log(2.7))
    print(log(50/15)*R*C)
    print(15/R)
# char26_6()

def char26_5():
    Im=2.5/1000
    E=1.52
    R0=E/Im
    R=R0-65
    Ra=65
    I1=E/(R0+200)
    I2=Im/4*3
    R2=E/I2
    Rx=R2-R0
    print((R),R0,I1,Rx)
# char26_5()
def char26_4():
    I=40/175
    Vb=15-75*I
    I1=15/75
    I2=25/100
    Ia=I1-I2
    print(Vb,Ia,Ia)
# char26_4()
def char26_2():
    R=4.26
    E=8.93
    I0=E/(R+R_parellel(R,R))
    Va=E-I0*R
    I2=Va/R
    P2=I2**2*R
    P1=I0**2*R
    print(I0)
    print(I2,P1,P2)


# char26_2()
def char25_10():
    copper_resi = 1.72 / 10 ** 8
    T0 = 20
    l = 3
    l1 = 1.2
    r1 = 0.7 / 1000
    l2 = 1.8
    r2 = 0.35 / 1000
    I = 3 / 1000
    # T0 = 20
    # l = 3
    # l1 = 1.2
    # r1 = 0.8 / 1000
    # l2 = 1.8
    # r2 = 0.4 / 1000
    # I = 2.5 / 1000
    S1=pi*r1**2
    S2=pi*r2**2
    R1=l1*copper_resi/pi/r1**2
    R2 = l2 * copper_resi / pi / r2 ** 2
    v1=I*R1
    E1=v1/l1
    v2=I*R2
    E2=v2/l2
    print(E1,E2)
    print(I*copper_resi/S1)
    print(I*copper_resi/S2)
    print(I,copper_resi,r1,l1)
    print((R1+R2)*I)
# char25_10()
def char25_9():
    r1=1
    r2=8
    r=10
    V1=12
    V2=8
    I=(V1-V2)/r
    P=I**2*r
    P1=V1*I
    P2=V2*I
    print(P1,P2)

    print(I,P)
# char25_9()


def char25_8():
    R1=1
    R2=5
    R=R1+R2
    V=12
    P=V**2/R
    I=V/R
    Pr=I**2*R1
    Pe=I**2*R2
    print(P,Pr,Pe)
# char25_8()
def char25_7():
    V=120
    P=60
    R=V**2/P
    print(R)
    I=V/R
    print(I)
# char25_7()
def char25_5():
    v=3
    r=19
    print(v**2/r)
    print(v ** 2 / r*5.1*3600)
    print(19*(sqrt(2)-1))
# char24_6()
def char25_3():
    copper_resi=1.72/10**8
    d=12.5/100
    I=115
    l=100*1000
    S=(d/2)**2*pi
    R=l*copper_resi/S
    V=I*R
    print(V)
    print(I**2*R*3600)
    print(10/6*4)
# char25_3()
def char25_2():
    t=7.9
    Q=55*t-0.65*t**3/3
    print(Q,Q/t)
# char25_2()
def char25_1():
    print(4*2.7*3600)
# char25_1()
def char24_12():
    A=46/10**6
    d0=0.68/10**3
    C0=A*epsilon_0/d0
    C1=C0+0.3/10**12
    d1=A*epsilon_0/C1
    print((d0-d1)*1000)

    # print(C0,C1)

# char24_12()
def char24_11():
    E1=3.2*10**5
    E2=2.2*10**5
    K=E1/E2
    sur_density0=E1*epsilon_0
    sur_density1=sur_density0*(1-1/K)
    print(sur_density1)
    print(E1/E2)
# char24_11()
def char24_10():
    c1=3/10**6
    c2=6/10**6
    c3=6/10**6
    c4=c1
    v=210
    c12=C_series(c1,c2)
    q=v*c12
    v1=q/c1
    v3=q/c2
    q13=105*(c1+c3)
    q1=q13/3
    q3=q13/3*2
    print(q,q1,q3)
    print(v3-v1)
    print(v3,v1)
    print(q3-q1)

# char24_10()

def char_24_9():
    c1=8.6
    c2=4.8
    c3=6.2
    c4=11.8
    c5=3.5
    c33=c5+C_series(c3,c4)
    c=C_series(C_series(c1,c2),c5+C_series(c3,c4))
    v=12
    U=c*v**2/2
    q0=c*v
    v1=q0/c1
    v3=q0/c33
    v2=v-v1-v3
    U2=c2/2*v2**2
    print(U)
    print(U2)
# char_24_9()
def char_24_8():
    c=(35+75)/10**9
    c0=35/10**9
    c1=75/10**9
    V=220
    q=c*V
    print(q*10**6)
    q0=35*V/1000
    q1=75*V/1000
    print(q0,q1)
    U=c*V**2/2
    print(U*1000)
    print(U*1000/c*c0,U*1000/c*c1)
# char_24_8()
def char_24_7():
    c=920/10**12
    q=2.55/10**6
    v=q/c
    print(v)
    U0=q*v/2
    U1=q*v*2/2
    print(U1-U0)
# char_24_7()


def char_24_6():
    C=C_series(10,13,9)
    print(C)
    V=50
    U=V**2*C/2
    Q=C*V
    print(Q)
# char_24_6()
def char24_5():
    C=C_series(20,15)
    print(C)
# char24_5()

def char24_4():
    r=5/100
    A=4*pi*r**2
    V=220
    q=3.5/10**9
    C=q/V
    Rreci=1/r+4*pi*epsilon_0/C
    r_inner=1/Rreci
    E=Ke*q/r_inner**2

    print(C,r_inner,E)
    print(E)
# char24_4()

def char24_2():
    S=(0.03)**2
    d=0.005
    C=S*epsilon_0/d
    print(C)
# char24_2()
def char24_1():
    v=30
    c=7.9/10**6
    print(c*v)

# char24_1()
def char23_10():
    r1=25/100
    r2=60/100
    q=150/10**6
    print(potential(q,r2))
# char23_10()


def char23_10():
    v0=6.9*10**6
    E=1250
    F=E*e
    a=F/eMass
    l0=6/100
    t0=l0/v0
    r0=a*t0**2/2
    angle=degrees(atan(a*t0/v0))
    l1=12/100
    t1=l1/v0
    r1=a*t0*t1
    print(r1+r0)
    print(F,a,r0,angle)
# char23_10()
def char23_9():
    print(6*2/sqrt(2)-12-4/sqrt(3))
# char23_9()
def char23_8():
    q=7.2/10**9
    We=-3.35/10**5
    V0=We/q
    E=-V0/8.5*100
    print(V0,E)
# char23_8()

def char23_7():
    R=6.2/100
    lam=9/10**6
    h=4.8/100
    d=R+h
    print(log(d/R)*lam/2/pi/epsilon_0)
# char23_7()
def char34_6():
    q=20.5/10**9
    d=32.5/100
    r=13/100
    l=sqrt(d**2+r**2)
    P0=potential(q,l)
    P1=potential(q,r)
    v=sqrt(2*(P1-P0)*e/eMass)
    print(v)
# char34_6()

def cha323_5():
    V=350
    x0=0.65
    x1=0.9
    E=V/(x1-x0)
    print(E)
    q=-0.2/10**6
    U=V*q
    print(U)
# cha323_5()
def char23_4():
    q0=27/10**9
    E=3.2*10**4
    d=0.46
    d1=0.71
    d2=2.5
    print(E*q0*d1)
    print(E*q0*d2*cos(radians(90+45)))
# char23_4()

def char23_2():
    q1=-2.6/10**6
    q2=-7.3/10**6
    m2=1.7/1000
    v2=22
    d0=0.8
    d1=0.44
    E0=v2**2*m2/2
    U0=potential(q1,d0)
    U1=potential(q1,d1)
    V=U1*q2-U0*q2
    v1=sqrt(2/m2*(E0-V))
    print(V,E0,E0-V,v1**2*m2/2)
    print(v1)
    d1=Ke*q1*q2/(E0+U0*q2)
    print(d1)
# char23_2()

def char23_1():
    q1=3.5/10**6
    q2=-4.9/10**6
    x0 = 0.17
    x1 = 0.29
    y1 = 0.29
    U1=Ke*q1/x0

    U2=Ke*q1/sqrt(x1**2+y1**2)
    V=(U2-U1)*q2
    print(-V)

# char23_1()


def lab_1():
    experimental=[166.66,1041.67,5555,56]
    theoretical=[158,1000,5280]
    distance = [6,4.8,3.6]
    timePerCm=[0.001,0.2/1000,50/1000000]
    for x,y,grid,base in zip(experimental, theoretical,distance,timePerCm):
        t=base*grid
        t_theo=1/y
        friquency=1/t
        error=t*0.2/grid
        friMin=1/(t_theo+error)
        friMax=1/(t_theo-error)
        # print(t,error)
        print((x-y)/y,friMin,y,friMax,friquency)
        print("voltage", 7.6*2)
        print("theoretical frequency",y)
        print("experimental value", x)
        print("percentage error",(x-y)/y*100)

# lab_1()
def char22_5():
    charge_density_sheet=-2.5/10**9
    q_ball=4.8/10**8
    m=3/10**6
    e=charge_density_sheet/2/epsilon_0
    f=e*q_ball
    w=m*g
    print(f,w,degrees(atan(-f/w)))

# char22_5()
def char22_4():
    angle=radians(30)
    l=0.05
    h=l*sin(angle)
    width=0.06
    s0=width*h
    e1=2.5*10**4
    e2=8.2*10**4
    phi1=e1*s0
    phi2=-e2*s0
    phi=phi1+phi2
    q=phi*epsilon_0
    print(q)


# char22_4()
def char22_3():
    x=1
    y=3
    z=2
    E0=125
    phi1=y*z*E0
    print(phi1)
    q=-24/10**9
    phi12=q/epsilon_0
    print(phi12)
    phi2=phi12-phi1
    E2=phi2/y/z
    print(E2)
# char22_3()

def char22_2():
    d=0.12
    r=d/2
    q=-13/10**6
    S=4*pi*r**2
    S1 = 4 * pi * (r+0.06) ** 2

    print(q/epsilon_0/S)
    print(q / epsilon_0 / S1)
# char22_2()
def char22_1():
    E=19
    S=0.39
    angle=radians(69)
    print(E*S*cos(angle))

# char22_1()











def char21_10():
    q=2.5/10**6
    L=1.2
    print(2/0.6*q/L*sqrt(2)/2*Ke)
    print(2 / 0.6 * q / L * sqrt(2) / 2 * Ke*sqrt(2))
    print(Ke*q/0.3/L)
    print(Ke * q / 0.3 / L*e)


# char21_10()
def char21_6():
    print(EforceNano(6,10**9,0.15)+EforceNano(6,10**9,0.45))
    a=-0.15
    b=0.15
    x=0.15
    y=-0.4

    print(EforceNano(6, 10 ** 9, sqrt(0.3**2+y**2))/5*3)
    print(EforceNano(6, 10 ** 9, sqrt(0.3 ** 2 + y ** 2)) / 5 * 4+EforceNano(6, 10 ** 9, y))
    xcom=EforceNano(6, 10 ** 9, sqrt(0.3**2+y**2))/5*3
    ycom=EforceNano(6, 10 ** 9, sqrt(0.3 ** 2 + y ** 2)) / 5 * 4+EforceNano(6, 10 ** 9, y)
    print(sqrt(xcom**2+ycom**2))
    print(degrees(atan(-ycom/xcom))+360)
    print(EforceNano(6, 10 ** 9, sqrt(0.15 ** 2 + 0.2 ** 2)) / 5 * 4*2)

# char21_6()
def char21_5():
    E = 1.00 * 10**4
    print(e*E)
    print(eMass)
    print(e*E/eMass/g)
    print(EforceNano(12,12,0.1)/eMass/g)
    print(e*E/g)
    print(e * E / g/eMass)
# char21_5()

def char21_4():
    v0=10**6
    t=0.02/v0
    a=0.01/t**2
    E=a*eMass/e
    print(E)
    aProton=E*e/protonMass
    displacement=aProton*t**2/2
    print(displacement)
# char21_4()


def char21_3():
    print(solve(2.5/x**2-3.5/(x+0.6)**2),x)

def chap21_2():
    f1=Eforce(3/10**6,5/10**6,0.2)
    f=f1+7
    print(f1,f,solve(Ke*3/10**6*8/10**6/x**2-f,x))

def chap21_1():
    print(e)
    q=3.5/10**9
    Enumber=q/e
    print("e number",q/e)
    molarMass = 207
    m=7.6
    print("Nled",(7.6/molarMass*Avogadro))
    print("Nled",Enumber/(7.6/molarMass*Avogadro))















# w=205*2*pi/60
# o0=5.45*pi*2
# ao=w**2/2/o0
# F=1375*cos(radians(8.5))

# x=tan(radians(1.24))*0.0925



