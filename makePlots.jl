using Revise
using BenchmarkTools
using Plots
using IntervalSets
using LinearAlgebra
using DifferentialEquations
using StaticArrays
using LaTeXStrings
using NonlinearSolve
using Flux


#plotting stuff
plot_font="Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false,legendfontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12)
#formater

push!(LOAD_PATH,pwd()*"/EquationofState")
using EquationsOfStates
include("hubble_equations.jl")

HRGLow=HadronResonaceGas(Maxmass=0.5,condition=waleckacondition)
HRG=HadronResonaceGas()
LQCD=LatticeQCD()
Walecka2=WaleckaModel2()
fmGeV= 1/0.1973261


#gluing functions
function Ttrans2(mu)
    0.1+0.28*mu-0.2*mu^2#0.1+0.8*mu-0.5*mu^2
end

function fTrans2(T,mu)#,t)
    tanh((T-Ttrans2(mu))/(0.1*Ttrans2(0)))    
end

transferFunction2=Gluing(fTrans2)

highEOS=1/2*(1-transferFunction2)*HRG+1/2*(1+transferFunction2)*LQCD
fullEOS=fmGeV^3*(highEOS+Walecka2)
fluidproperties=FluidProperties(fullEOS
,EquationsOfStates.ZeroViscosity()
,EquationsOfStates.ZeroBulkViscosity()
,EquationsOfStates.ZeroDiffusion())

gammaA=2960
hrate=get_hubble_rate(gamma=gammaA)

function expansion_rate(time)
    H0=hrate(time) #hrate function gives the Hubble rate from the Landau matching
   return H0
end
tl=0.015:0.00001:0.035
hplot=plot(tl,expansion_rate.(tl),xlabel=L"$t$ [fm/c]",ylabel=L"$H$ [c/fm]")

function maxHrate(gamma)
    hrate=get_hubble_rate(;gamma)
    tlist=collect(0.015:0.000001:0.035)
    return maximum(hrate.(tlist))
end

T0=0.003
mu0=.92
check_for_transition(T0,mu0)
u0=[T0,mu0,0.0,.0]
tspan=(0.015,.04)

f(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,0.2,gammaA,1/50)
problem = ODEProblem(f, u0, tspan)
solution1 = solve(problem,AutoTsit5(Rosenbrock23(autodiff=false)),dtmax=0.001*(tspan[2]-tspan[1]))
plot(solution1)
#functions for entropy etc 

function maxT(zetaMax)
    T0=0.005
    mu0=.92
    u0=[T0,mu0,0.0,0.0]
    fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,zetaMax,1060)
    problemT = ODEProblem(fT, u0, tspan)
    solutionT = solve(problemT,AutoTsit5(Rosenbrock23()),dtmax=0.01*tspan[2])
    return maximum(solutionT[1,:])
end

function maxTH(gammaE,mexzeta,beta1)
    T0=0.005
    mu0=.92
    u0=[T0,mu0,0.0,0.0]
    fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,mexzeta,gammaE,beta1)
    problemT = ODEProblem(fT, u0, tspan)
    solutionT = solve(problemT,AutoTsit5(Rosenbrock23()),dtmax=0.01*tspan[2])
    return maximum(solutionT[1,:])
end

function finalTH(gammaE,mexzeta,beta1)
    T0=0.005
    mu0=.92
    u0=[T0,mu0,0.0,0.0]
    fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,mexzeta,gammaE,beta1)
    problemT = ODEProblem(fT, u0, tspan)
    solutionT = solve(problemT,AutoTsit5(Rosenbrock23()),dtmax=0.01*tspan[2])
    return last(solutionT[1,:])
end

function entropyProduction(sol,t,gammaE,beta1,zetaMax)
    T=sol(t)[1]
    mu=sol(t)[2]
    dtp=pressure_derivative(T,mu,Val(1),Val(0),fullEOS)
    dT=sol(t,Val{1})[1]
    dPi=sol(t,Val{1})[3]
    piB=sol(t)[3]
    hrate=get_hubble_rate(gamma=gammaE)
    Hrat=hrate(t)
    zetaParam=dtp/(1+((sqrt(T^2+0.188172^2*mu^2)-0.175)/0.024)^2)
    zeta=zetaMax*zetaParam
    a=piB^2/(T*zeta)
    b=-piB*(dPi/5+3/10*piB*Hrat+3*Hrat)#3*Hrat*piB+3*Hrat*T*s
    #return piB^2/(T*zeta)#-piB *(3* Hrat + 3*beta1 * Hrat * piB  -beta1*dT/T*piB + 2*beta1*dPi) / T
    #return -piB*(dPi/5+3/10*piB*Hrat+3*Hrat)#3*Hrat*piB+3*Hrat*T*s
    return a
end

function entropyProductionb(sol,t,gammaE,beta1,zetaMax)
    T=sol(t)[1]
    mu=sol(t)[2]
    dtp=pressure_derivative(T,mu,Val(1),Val(0),fullEOS)
    dT=sol(t,Val{1})[1]
    dPi=sol(t,Val{1})[3]
    piB=sol(t)[3]
    hrate=get_hubble_rate(gamma=gammaE)
    Hrat=hrate(t)
    zetaParam=dtp/(1+((sqrt(T^2+0.188172^2*mu^2)-0.175)/0.024)^2)
    zeta=zetaMax*zetaParam
    a=piB^2/(T*zeta)
    b=-piB *(3* Hrat + 3*beta1 * Hrat * piB  -beta1*dT/T*piB + 2*beta1*dPi) / T#3*Hrat*piB+3*Hrat*T*s
    #return piB^2/(T*zeta)#-piB *(3* Hrat + 3*beta1 * Hrat * piB  -beta1*dT/T*piB + 2*beta1*dPi) / T
    #return -piB*(dPi/5+3/10*piB*Hrat+3*Hrat)#3*Hrat*piB+3*Hrat*T*s
    return b
end
    
    function entroMax(gammaE,mexzeta,beta1)
        T0=0.005
        mu0=.92
        u0=[T0,mu0,0.0,0.0]
        fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,mexzeta,gammaE,beta1)
        problemT = ODEProblem(fT, u0, tspan)
        solutionT = solve(problemT,AutoTsit5(Rosenbrock23()),dtmax=0.01*tspan[2])
        entroT=entropyProduction.(Ref(solutionT),solutionT.t,Ref(gammaE),Ref(beta1),Ref(mexzeta))
        return maximum(entroT)
    end
    
    function entroFull(gammaE,mexzeta,beta1)
        tspan=(0.01,0.04)
        T0=0.005
        mu0=.92
        u0=[T0,mu0,0.0,0.0]
        fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,mexzeta,gammaE,beta1)
        problemT = ODEProblem(fT, u0, tspan)
        solutionT = solve(problemT,AutoTsit5(Rosenbrock23(autodiff=false)),dtmax=0.01*tspan[2])
        entroT=entropyProduction.(Ref(solutionT),solutionT.t,Ref(gammaE),Ref(beta1),Ref(mexzeta))
        return entroT, solutionT
    end

    #make the viscosity plots
    viscL=collect(0.01:0.5:1.5)

    entroLBV1=entro.(Ref(2960),viscL,Ref(1/2))
entroLBV2=entro.(Ref(2960),viscL,Ref(1/5))
entroLBV3=entro.(Ref(2960),viscL,Ref(1/10))

entroViscPlot=plot(viscL,entroLBV1,label="1/2",xlabel= L"(\zeta/s)_{max}",ylabel=ylabel=L"\int \mathrm{d}t \; \nabla_\mu S^\mu \;\; [\mathrm{fm}^{-3}]",legendtitle =L"\beta_0 \;\; [\mathrm{fm}/\mathrm{c}]")
plot!(viscL,entroLBV2,label="1/5")
plot!(viscL,entroLBV3,label="1/10")
#savefig(entroViscPlot,"PlotsThesis/entroViscTauPlot2.pdf")