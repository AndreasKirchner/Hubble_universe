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
fullEOS=fmGeV^3*(highEOS)#+Walecka2)
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
tspan=(0.02,.035)

f(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,0.2,gammaA,10)
problem = ODEProblem(f, u0, tspan)

#@time solution1 = solve(problem,AutoTsit5(Rosenbrock23(autodiff=false)),dtmax=0.01*(tspan[2]-tspan[1])) #557s
#@time solution2 = solve(problem,AutoTsit5(Rosenbrock23(autodiff=false)))
#@time solution3 = solve(problem,AutoTsit5(QNDF(autodiff=false)))
#@time solution4 = solve(problem,AutoTsit5(Rodas5P(autodiff=false)))
@time solution1 = solve(problem,AutoTsit5(Rodas4(autodiff=false)))
plotSol(solution1)

p=pressure.(solution1[1,:],solution1[2,:],Ref(fullEOS))
n=pressure_derivative.(solution1[1,:],abs.(solution1[2,:]),Val(0),Val(1),Ref(fullEOS))

plot(solution1.t,n)


plot(solution1.t,abs.(solution1[2,:]))

@time 1+1

plot(solution1)
plot(solution2)
plot(solution4)

p1=plot(solution1.t,solution1[1,:],label="T [GeV]",xlabel="t [fm/c]",ylabel="fields")
plot!(solution1.t,solution1[2,:],label="μ [GeV]")
plot!(solution1.t,1/gammaA .*expansion_rate.(solution1.t),label=L"H/$\gamma$ [c/fm]")
plot!(solution1.t, 1/1000 .*solution1[3,:],label=L"Π  [TeV/$\mathrm{fm}^3$]")
#savefig(p1,"PlotsPaper/OneEvent.pdf")

function plotSol(sol)
    p1=plot(sol.t,sol[1,:],label="T")
    plot!(sol.t,sol[2,:],label="μ")
    plot!(sol.t,1/gammaA .*expansion_rate.(sol.t),label="H/1000")
    plot!(sol.t, 1/1000 .*sol[3,:],label="Π")
    
end

plotSol(solution1)
plotSol(solution2)
plotSol(solution4)
plotSol(solution5)



pdplot=plot(abs.(solution1[2,:]),solution1[1,:],yaxis=:log,xlabel=L"\mu \;\; [\mathrm{GeV}]",ylabel=L"\mathrm{ln}(T)\;\; [\mathrm{ln(GeV)}]")
plot!((last(solution1[2,:]),last(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
plot!((first(solution1[2,:]),first(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
#savefig(pdplot,"PlotsPaper/phaseDiagramPoints.pdf")
(last(solution1[2,:]),last(solution1[1,:]))
minimum(solution1[2,:])

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
    solutionT = solve(problemT,AutoTsit5(Rodas4()),dtmax=0.01*(tspan[2]-tspan[1]))
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

    function trapezoidal_rule(x, fx)
        n = length(x)
        if n != length(fx)
        error("Length of x and fx should be the same")
        end
        
        integral = 0.0
        for i in 1:n-1
        h = x[i+1] - x[i]
        integral += (fx[i] + fx[i+1]) * h / 2
        end
        
        return integral
        end
        
    
    function entro(gammaE,mexzeta,beta1)
        T0=0.005
        mu0=.92
        u0=[T0,mu0,0.0,0.0]
        fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,mexzeta,gammaE,beta1)
        problemT = ODEProblem(fT, u0, tspan)
        solutionT = solve(problemT,AutoTsit5(Rodas5P()),dtmax=0.001*(tspan[2]-tspan[1]))
        entroT=entropyProduction.(Ref(solutionT),solutionT.t,Ref(gammaE),Ref(beta1),Ref(mexzeta))
        return trapezoidal_rule(solutionT.t,entroT)
    end

entroStuff=entropyProduction.(Ref(solution1),solution1.t,Ref(gammaA),Ref(2),Ref(1))
entroPlot1=plot(solution1.t,entroStuff,xlabel=L"$t$ [fm/c]",label=L"$\nabla_\mu S^\mu \;\; [\mathrm{fm}^{-4}]$")
plot!(solution1.t,expansion_rate.(solution1.t) .*gammaA,label=L"$\gamma H$ [c/fm]")
#savefig(entroPlot,"PlotsPaper/entropyProductionOneEvent.pdf")

entroOneEvent=plot(solution1.t,entroStuff,xlabel=L"$t$ [fm/c]",label=L"$\nabla_\mu S^\mu \;\; [\mathrm{fm}^{-4}]$")

entroStuff1=entropyProduction.(Ref(solution1),solution1.t,Ref(gammaList[1]),Ref(10),Ref(0.2))
entroPlot=plot(solution1.t,entroStuff,xlabel=L"$t$ [fm/c]",label=L"$\nabla_\mu S^\mu \;\; [\mathrm{fm}^{-4}]$")
plot!(solution1.t,expansion_rate.(solution1.t) .*gammaA,label=L"$\gamma H$ [c/fm]")
#savefig(entroPlot,"PlotsPaper/entropyProductionOneEvent.pdf")


entro(2000,1,1)
    #make the viscosity plots
    viscL=collect(0.01:0.05:1.5)

viscL

entroLBV1=entro.(Ref(2960),viscL,Ref(2))
entroLBV2=entro.(Ref(2960),viscL,Ref(5))
entroLBV3=entro.(Ref(2960),viscL,Ref(10))

finalTV1=finalTH.(Ref(2960),viscL,Ref(2))
finalTV2=finalTH.(Ref(2960),viscL,Ref(5))
finalTV3=finalTH.(Ref(2960),viscL,Ref(10))
finalTV4=finalTH.(Ref(2960),viscL,Ref(20))

finalTPlot=plot(viscL,finalTV1,label="2",xlabel= L"(\zeta/s)_{max}",ylabel=ylabel=L"T(t=\infty) \;\; [\mathrm{GeV}]",legendtitle =L"\beta_0 \;\; [\mathrm{fm}/\mathrm{c}]")
plot!(viscL,finalTV2,label="5")
plot!(viscL,finalTV3,label="10")
plot!(viscL,finalTV4,label="20")
#savefig(finalTPlot,"PlotsPaper/entroViscTempPlot.pdf")

entroViscPlot=plot(viscL,entroLBV1,label="2",xlabel= L"(\zeta/s)_{max}",ylabel=L"\int \mathrm{d}t \; \nabla_\mu S^\mu \;\; [\mathrm{fm}^{-3}]",legendtitle =L"\beta_0")
plot!(viscL,entroLBV2,label="5")
plot!(viscL,entroLBV3,label="10")
#savefig(entroViscPlot,"PlotsPaper/entroViscTauPlot2.pdf")

gammaList=collect(1500:100:3500)
maxHList=maxHrate.(gammaList)


Tmaxtf1=finalTH.(gammaList,Ref(1.0),Ref(2))
Tmaxtf2=finalTH.(gammaList,Ref(1.0),Ref(10))
Tmaxtf3=finalTH.(gammaList,Ref(1.0),Ref(20))
Tmaxtf4=finalTH.(gammaList,Ref(1.0),Ref(5))

entroGamma1=entro.(gammaList,Ref(1.0),Ref(2))
entroGamma2=entro.(gammaList,Ref(1.0),Ref(10))
entroGamma3=entro.(gammaList,Ref(1.0),Ref(20))
entroGamma4=entro.(gammaList,Ref(1.0),Ref(50))

entro(gammaList[5],1.0,10)
gammaList[4]


entroGammaPlot=plot(gammaList,entroGamma1,label="2",xlabel=L"$\gamma$",ylabel=L"\int \mathrm{d}t \; \nabla_\mu S^\mu \;\; [\mathrm{fm}^{-3}]",legendtitle = L"\beta_0")
plot!(gammaList,entroGamma2,label=L"10")
plot!(gammaList,entroGamma3,label=L"20")
plot!(gammaList,entroGamma4,label=L"50")
#savefig(entroGammaPlot,"PlotsPaper/totalEntroGamma.pdf")

solt=getSol(gammaList[1],1,2)
plotSol(solt)

plot(solt.t,solt[1,:])
plot!(solt.t,solt[2,:])
plot!(solt.t,solt[3,:])
plot!(solt.t,expansion_rate.(solt.t) ./ gammaList[1])
entroGamma1

#using Plots
TfinalPlot=plot(maxHList,Tmaxtf1,label="2",xlabel=L"$\mathrm{max} (H) \;\; [\mathrm{c}/\mathrm{fm}] $",ylabel=L"T(t=\infty) \;\; [\mathrm{GeV}]",legendtitle = L"\beta_0")
plot!(maxHList,Tmaxtf2,label=L"10")
plot!(maxHList,Tmaxtf3,label=L"20")
plot!(maxHList,Tmaxtf4,label=L"50")
#savefig(TfinalPlot2,"PlotsPaper/TfinalHRateMaxPlot.pdf")

TfinalPlot2=plot(gammaList,Tmaxtf1,label="2",xlabel=L"$\gamma $",ylabel=L"T(t=\infty) \;\; [\mathrm{GeV}]",legendtitle = L"\beta_0")
plot!(gammaList,Tmaxtf2,label=L"10")
plot!(gammaList,Tmaxtf3,label=L"20")
plot!(gammaList,Tmaxtf4,label=L"50")


function getSol(gammaE,mexzeta,beta1)
    tspan=(0.01,0.04)
    T0=0.002
    mu0=.92
    u0=[T0,mu0,0.0,0.0]
    fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,mexzeta,gammaE,beta1)
    problemT = ODEProblem(fT, u0, tspan)
    solutionT =  solve(problemT,AutoTsit5(Rosenbrock23(autodiff=false)),dtmax=0.01*tspan[2])
    return solutionT
end

sol1=getSol(gammaA,1,10)
sol2=getSol(gammaA,1,20)
sol3=getSol(gammaA,1,50)


multiPhaseD=plot(sol1[2,:],sol1[1,:],yaxis=:log,xlabel=L"$\mu$",ylabel=L"$T$",label="10",legendtitle = L"\beta_0")
plot!(sol2[2,:],sol2[1,:],label="20")
plot!(sol3[2,:],sol3[1,:],label="50")
#savefig(multiPhaseD,"PlotsPaper/multiBetaPhaseDiag.pdf")


last(sol1[1,:])
last(sol2[1,:])
last(sol3[1,:])



hrate=get_hubble_rate()
gamma=2960
v=sqrt(1-1/gamma^2)
rList=collect(-0.03:0.0005:0.03)



nlist=zeros(length(rList))
elist=zeros(length(rList))
vlist=zeros(length(rList))
nulist=zeros(length(rList))

t=0.03
for i in eachindex(rList)
    nlist[i]=doLandauMatchingBig(t,rList[i],v)[1]
    vlist[i]=doLandauMatchingBig(t,rList[i],v)[3]
    elist[i]=doLandauMatchingBig(t,rList[i],v)[2]
    nulist[i]=doLandauMatchingBig(t,rList[i],v)[4]
end

plot(nulist)
plot!(30 .*nlist)
    plot(vlist)
plot(elist)