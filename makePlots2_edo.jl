
using Plots
using IntervalSets
using LinearAlgebra
using DifferentialEquations
using StaticArrays
using LaTeXStrings
using NonlinearSolve
using CairoMakie 
xsize=300
ysize=xsize*3/4


pwd()
#plotting stuff
plot_font="Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false,legendfontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12)
#formater

push!(LOAD_PATH,pwd()*"/EquationofState")
using EquationsOfStates
include("hubble_equation.jl")

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
fullEOS=highEOS#fmGeV^3*(highEOS)#+Walecka2)
fluidproperties=FluidProperties(fullEOS
,EquationsOfStates.ZeroViscosity()
,EquationsOfStates.ZeroBulkViscosity()
,EquationsOfStates.ZeroDiffusion())

Tlist=collect(0:0.01:3)
muOverT1=0
muList1=muOverT1 .* Tlist
pList1=zeros(length(Tlist))
muOverT2=1
muList2=muOverT2 .* Tlist
pList2=zeros(length(Tlist))
muOverT3=2
muList3=muOverT3 .* Tlist
pList3=zeros(length(Tlist))

pList4=zeros(length(Tlist))
muList4=collect(0.75:0.0001:0.95)
pList4 = pressure.(Ref(0.0),muList4,Ref(fmGeV^3*(Walecka2+HRGLow)))
pList5=zeros(length(Tlist))
pList5 = pressure.(Ref(0.01),muList4,Ref(fmGeV^3*(Walecka2+HRGLow)))
pList6=zeros(length(Tlist))
pList6 = pressure.(Ref(0.025),muList4,Ref(fmGeV^3*(Walecka2+HRGLow)))
plot(muList4,pList4 )
plot!(muList4,pList5 )
plot!(muList4,pList6 )

pressure(0.0,0.9,HRGLow)


fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"\mu \;[\mathrm{GeV}]",#,ylabel = L"The y label"
    ylabel=L"p \; [\mathrm{MeV}/\mathrm{fm}^3]"
    ,xgridvisible = false,
        ygridvisible = false
)
CairoMakie.ylims!(ax,-0.01,3)
CairoMakie.xlims!(ax,0.75,0.94)
    lines!(ax, muList4,10^3 .*pList4,label=L"0 ",color=Makie.wong_colors()[1])
    lines!(ax, muList4, 10^3 .*pList5 ,label=L"10 ",color=Makie.wong_colors()[2])
    lines!(ax, muList4, 10^3 .*pList6 ,label=L"25 ",color=Makie.wong_colors()[4])
    axislegend(L"T\; [\mathrm{MeV}]", framevisible=false,position=:lt)
    
    #text!(-0.015,0.3, text=L"\frac{\nu}{\gamma n_0}")
    resize_to_layout!(fig)
    fig
end
save("LowEOSPlot.pdf",fig)



for i in eachindex(Tlist)
    pList1[i]=pressure(Tlist[i],muList1[i],fullEOS)/Tlist[i]^4
    pList2[i]=pressure(Tlist[i],muList2[i],fullEOS)/Tlist[i]^4
    pList3[i]=pressure(Tlist[i],muList3[i],fullEOS)/Tlist[i]^4
end

plot(Tlist,pList1)
plot!(Tlist,pList2)
plot!(Tlist,pList3)



fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"T \;[\mathrm{GeV}]",#,ylabel = L"The y label"
    ylabel=L"p/T^4"
    ,xgridvisible = false,
        ygridvisible = false
)
    CairoMakie.ylims!(ax,0.0,6)
    CairoMakie.xlims!(ax,0.0,3)
    lines!(ax, Tlist, pList1,label=L"0 ",color=Makie.wong_colors()[1])
    lines!(ax, Tlist, pList2,label=L"1 ",color=Makie.wong_colors()[2])
    lines!(ax, Tlist, pList3 ,label=L"2",color=Makie.wong_colors()[3])
    axislegend(L"\mu_\text/T", framevisible=false, position=:rb)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end

save("pOverT4.pdf",fig)

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

f(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,2,gammaA,7)
problem = ODEProblem(f, u0, tspan)

#@time solution1 = solve(problem,AutoTsit5(Rosenbrock23(autodiff=false)),dtmax=0.01*(tspan[2]-tspan[1])) #557s
#@time solution2 = solve(problem,AutoTsit5(Rosenbrock23(autodiff=false)))
#@time solution3 = solve(problem,AutoTsit5(QNDF(autodiff=false)))
#@time solution4 = solve(problem,AutoTsit5(Rodas5P(autodiff=false)))
@time solution1 = solve(problem,AutoTsit5(Rodas4(autodiff=false)), reltol = 1e-15)
solution2 = solve(problem,AutoTsit5(Rodas5P(autodiff=false)), reltol = 1e-15)
#plotSol(solution1)
ind=3
plot(solution1.t,solution1[ind,:])
plot!(solution2.t,solution2[ind,:])

p1=plot(solution1.t,solution1[1,:],label="T [GeV]",xlabel="t [fm/c]",ylabel="fields")
plot!(solution1.t,solution1[2,:],label="μ [GeV]")
plot!(solution1.t,1/gammaA .*expansion_rate.(solution1.t),label=L"H/$\gamma$ [c/fm]")
plot!(solution1.t, 1/1000 .* solution1[3,:],label=L"Π  [TeV/$\mathrm{fm}^3$]")

fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"t \;[\mathrm{fm}/\mathrm{c}]",#,ylabel = L"The y label"
    ylabel=L"\mathrm{fields}"
    ,xgridvisible = false,
        ygridvisible = false
)
    CairoMakie.ylims!(ax,-3.,3.5)
    CairoMakie.xlims!(ax,0.02,0.03)
    lines!(ax, solution1.t,solution1[1,:],label="T [GeV]",color=Makie.wong_colors()[1])
    lines!(ax, solution1.t,solution1[2,:],label="μ [GeV]",color=Makie.wong_colors()[2])
    lines!(ax, solution1.t,1/gammaA .*expansion_rate.(solution1.t),label=L"H/$\gamma$ [c/fm]",color=Makie.wong_colors()[3])
    lines!(ax, solution1.t, 1 ./1000 .* solution1[3,:],label=L"\pi_{\mathrm{bulk}} [TeV/$\mathrm{fm}^3$]",color=Makie.wong_colors()[4])
    axislegend(framevisible=false, position=:rb,orientation = :horizontal,nbanks = 2)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end

save("OneEvent.pdf",fig)

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

using LaTeXStrings
using CairoMakie 
xsize=300
ysize=xsize*3/4
fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],yscale = log10,
    xlabel=L"\mu \; [\mathrm{GeV}]",#,ylabel = L"The y label"
    ylabel=L"\log[T/ (1\; \mathrm{GeV})]"
    ,xgridvisible = false,
        ygridvisible = false
)
    CairoMakie.ylims!(ax,2.5*10^(-3),10.0)
    CairoMakie.xlims!(ax,-0.02,0.96)
    lines!(abs.(solution1[2,:]),solution1[1,:],label=L"0 ",color=Makie.wong_colors()[1])
    CairoMakie.scatter!(last(solution1[2,:]),last(solution1[1,:]),marker=:circle ,color=:black,markersize = 8)
    CairoMakie.scatter!(first(solution1[2,:]),first(solution1[1,:]),marker=:circle ,color=:black,markersize = 8)
    text!(0.94, 0.0045, text = "I", align = (:center, :center))
    text!(0.038, 0.95, text = "II", align = (:center, :center))
    xs = LinRange(0.7, 0.7, 1)
    ys = LinRange(0.04, 0.04, 1)
    us = [-0.2 for x in xs, y in ys]
    vs = [0.03 for x in xs, y in ys]
    arrows!(xs, ys, us, vs, arrowsize = 13)
    #point((last(solution1[2,:]),last(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    #plot!((first(solution1[2,:]),first(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    #lines!(ax, Tlist, pList2,label=L"1 ",color=Makie.wong_colors()[2])
    #lines!(ax, Tlist, pList3 ,label=L"2",color=Makie.wong_colors()[3])
    #axislegend(L"\mu_\text{B}/T", framevisible=false, position=:rb)

    xs = LinRange(0.05, 0.05, 1)
    ys = LinRange(6, 6, 1)
    us = [0 for x in xs, y in ys]
    vs = [-4 for x in xs, y in ys]
    arrows!(xs, ys, us, vs, arrowsize = 13)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end

save("phaseDiagFinal.pdf",fig)
available_marker_symbols()


#functions for entropy etc 


function maxT(zetaMax)
    T0=0.005
    mu0=.92
    u0=[T0,mu0,0.0,0.0]
    fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,zetaMax,1060)
    problemT = ODEProblem(fT, u0, tspan)
    solutionT = solve(problemT,AutoTsit5(Rodas4(autodiff=false)),dtmax=0.01*tspan[2])
    return maximum(solutionT[1,:])
end

function maxTH(gammaE,mexzeta,beta1)
    T0=0.005
    mu0=.92
    u0=[T0,mu0,0.0,0.0]
    fT(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,mexzeta,gammaE,beta1)
    problemT = ODEProblem(fT, u0, tspan)
    solutionT = solve(problemT,AutoTsit5(Rodas4(autodiff=false)),dtmax=0.01*tspan[2])
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
    #b=-piB*(dPi/5+3/10*piB*Hrat+3*Hrat)#3*Hrat*piB+3*Hrat*T*s
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
        solutionT = solve(problemT,AutoTsit5(Rodas4(autodiff=false)),dtmax=0.01*tspan[2])
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
        solutionT = solve(problemT,AutoTsit5(Rodas4(autodiff=false)),dtmax=0.01*tspan[2])
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
        solutionT = solve(problemT,AutoTsit5(Rodas5P()),dtmax=0.01*(tspan[2]-tspan[1]))#, reltol = 1e-15)
        entroT=entropyProduction.(Ref(solutionT),solutionT.t,Ref(gammaE),Ref(beta1),Ref(mexzeta))
        return trapezoidal_rule(solutionT.t,entroT)
    end

entroStuff=entropyProduction.(Ref(solution1),solution1.t,Ref(gammaA),Ref(7),Ref(2))
entroStuff2=entropyProduction.(Ref(solution2),solution2.t,Ref(gammaA),Ref(7),Ref(2))
entroPlot1=plot(solution1.t,entroStuff,xlabel=L"$t$ [fm/c]",label=L"$\nabla_\mu S^\mu \;\; [\mathrm{fm}^{-4}]$")
Plots.plot!(solution1.t,entroStuff)
Plots.plot(solution1.t,expansion_rate.(solution1.t) .*gammaA,label=L"$\gamma H$ [c/fm]")


fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"t \;[\mathrm{fm}/\mathrm{c}]"#,ylabel = L"The y label"
    ,xgridvisible = false,
        ygridvisible = false
)
    #CairoMakie.ylims!(ax,2.5*10^(-3),4.0)
    #CairoMakie.xlims!(ax,-0.02,0.96)
    #lines!(abs.(solution1[2,:]),solution1[1,:],label=L"0 ",color=Makie.wong_colors()[1])
    #CairoMakie.scatter!(last(solution1[2,:]),last(solution1[1,:]),marker=:circ,color=:black,markersize = 8)
    #CairoMakie.scatter!(first(solution1[2,:]),first(solution1[1,:]),marker=:circ,color=:black,markersize = 8)
    #text!(0.94, 0.0045, text = "I", align = (:center, :center))
    #text!(0.038, 0.95, text = "II", align = (:center, :center))
    #point((last(solution1[2,:]),last(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    #plot!((first(solution1[2,:]),first(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    lines!(ax, solution1.t,expansion_rate.(solution1.t)./gammaA ,label=L" \frac{H}{\gamma}\; [\mathrm{c}/\mathrm{fm}]",color=Makie.wong_colors()[1])
    lines!(ax, solution1.t,entroStuff./gammaA^(3/2),label=L"\frac{\nabla_\mu S^\mu}{\gamma^{3/2} }\; [\mathrm{fm}^{-4}]",color=Makie.wong_colors()[6])
    axislegend(framevisible=false, position=:rt)#orientation = :horizontal)

    #lines!(ax, Tlist, pList3 ,label=L"2",color=Makie.wong_colors()[3])
    #axislegend(L"\mu_\text{B}/T", framevisible=false, position=:rb)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end


save("entropyProductionOneEvent.pdf",fig)
trapezoidal_rule(solution1.t,entroStuff)
trapezoidal_rule(solution2.t,entroStuff2)

entroOneEvent=plot(solution1.t,entroStuff,xlabel=L"$t$ [fm/c]",label=L"$\nabla_\mu S^\mu \;\; [\mathrm{fm}^{-4}]$")

entroStuff1=entropyProduction.(Ref(solution1),solution1.t,Ref(gammaList[1]),Ref(10),Ref(0.2))
entroPlot=plot(solution1.t,entroStuff,xlabel=L"$t$ [fm/c]",label=L"$\nabla_\mu S^\mu \;\; [\mathrm{fm}^{-4}]$")
plot!(solution1.t,expansion_rate.(solution1.t) .*gammaA,label=L"$\gamma H$ [c/fm]")
#savefig(entroPlot,"PlotsPaper/entropyProductionOneEvent.pdf")

gammaA
@time entro(3500,1.5,20)
    #make the viscosity plots
viscL=collect(0.01:0.05:1.55)

viscL

entroLBV1=entro.(Ref(2960),viscL,Ref(2))
entroLBV2=entro.(Ref(2960),viscL,Ref(5))
entroLBV3=entro.(Ref(2960),viscL,Ref(10))

using JLD2
save_object("visc_entro_beta0_2_5_10.jdl2",exp1)
exp1=(viscL,entroLBV1,entroLBV2,entroLBV3)


(viscL,entroLBV1,entroLBV2,entroLBV3)=load("visc_entro_beta0_2_5_10.jdl2")["single_stored_object"]

aa=entro(gammaA,5,10)
aa
plot!(solution1.t,entroStuff,xlabel=L"$t$ [fm/c]",label=L"$\nabla_\mu S^\mu \;\; [\mathrm{fm}^{-4}]$")

1592/60

finalTV1=finalTH.(Ref(2960),viscL,Ref(2))
finalTV2=finalTH.(Ref(2960),viscL,Ref(5))
finalTV3=finalTH.(Ref(2960),viscL,Ref(10))
finalTV4=finalTH.(Ref(2960),viscL,Ref(20))

exp2=(viscL,finalTV1,finalTV2,finalTV3,finalTV4)
save_object("visc_finalT_beta0_2_5_10_20.jdl2",exp2)
(viscL,finalTV1,finalTV2,finalTV3,finalTV4)=load("visc_finalT_beta0_2_5_10_20.jdl2")["single_stored_object"]

finalTPlot=Plots.plot(viscL,finalTV1,label="2",xlabel= L"(\zeta/s)_{max}",ylabel=L"T(t=\infty) \;\; [\mathrm{GeV}]",legendtitle =L"\beta_0 \;\; [\mathrm{fm}/\mathrm{c}]")
Plots.plot!(viscL,finalTV2,label="5")
Plots.plot!(viscL,finalTV3,label="10")
Plots.plot!(viscL,finalTV4,label="20")

fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel= L"(\zeta/s)_{\mathrm{max}}",
    ylabel=L"T(t=\infty) \;[\mathrm{GeV}]"#,ylabel = L"The y label"
    ,xgridvisible = false,
        ygridvisible = false
)
    #CairoMakie.ylims!(ax,0,850.0)
    #CairoMakie.xlims!(ax,-0.02,0.96)
    #lines!(abs.(solution1[2,:]),solution1[1,:],label=L"0 ",color=Makie.wong_colors()[1])
    #CairoMakie.scatter!(last(solution1[2,:]),last(solution1[1,:]),marker=:circ,color=:black,markersize = 8)
    #CairoMakie.scatter!(first(solution1[2,:]),first(solution1[1,:]),marker=:circ,color=:black,markersize = 8)
    #text!(0.94, 0.0045, text = "I", align = (:center, :center))
    #text!(0.038, 0.95, text = "II", align = (:center, :center))
    #point((last(solution1[2,:]),last(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    #plot!((first(solution1[2,:]),first(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    lines!(ax,viscL,finalTV1,label="2",color=Makie.wong_colors()[1])
    lines!(ax,viscL,finalTV2,label="5",color=Makie.wong_colors()[6])
    lines!(ax,viscL,finalTV3,label="10",color=Makie.wong_colors()[3])
    lines!(ax,viscL,finalTV4,label="20",color=Makie.wong_colors()[4])
    #axislegend(framevisible=false, position=:lt)#orientation = :horizontal)

    #lines!(ax, Tlist, pList3 ,label=L"2",color=Makie.wong_colors()[3])
    axislegend(L"2\zeta /[\tau_{\text{bulk}}(e+p)]", framevisible=false, position=:rb,hline=:left,nbanks=2)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end

save("entroViscTempPlot.pdf",fig)

entroViscPlot=Plots.plot(viscL,entroLBV1,label="2",xlabel= L"(\zeta/s)_{max}",ylabel=L"\int \mathrm{d}t \; \nabla_\mu S^\mu \;\; [\mathrm{fm}^{-3}]",legendtitle =L"\beta_0")
Plots.plot!(viscL,entroLBV2,label="5")
Plots.plot!(viscL,entroLBV3,label="10")

fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel= L"(\zeta/s)_{\mathrm{max}}",
    ylabel=L"\int \mathrm{d}t \; \nabla_\mu S^\mu \; [\mathrm{fm}^{-3}]"#,ylabel = L"The y label"
    ,xgridvisible = false,
        ygridvisible = false
)
    CairoMakie.ylims!(ax,0,850.0)
    #CairoMakie.xlims!(ax,-0.02,0.96)
    #lines!(abs.(solution1[2,:]),solution1[1,:],label=L"0 ",color=Makie.wong_colors()[1])
    #CairoMakie.scatter!(last(solution1[2,:]),last(solution1[1,:]),marker=:circ,color=:black,markersize = 8)
    #CairoMakie.scatter!(first(solution1[2,:]),first(solution1[1,:]),marker=:circ,color=:black,markersize = 8)
    #text!(0.94, 0.0045, text = "I", align = (:center, :center))
    #text!(0.038, 0.95, text = "II", align = (:center, :center))
    #point((last(solution1[2,:]),last(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    #plot!((first(solution1[2,:]),first(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    lines!(ax,viscL,entroLBV1,label="2",color=Makie.wong_colors()[1])
    lines!(ax,viscL,entroLBV2,label="5",color=Makie.wong_colors()[6])
    lines!(ax,viscL,entroLBV3,label="10",color=Makie.wong_colors()[3])
    #axislegend(framevisible=false, position=:lt)#orientation = :horizontal)

    #lines!(ax, Tlist, pList3 ,label=L"2",color=Makie.wong_colors()[3])
    axislegend(L"2\zeta /[\tau_{\text{bulk}}(e+p)]", framevisible=false, position=:lt,hline=:left)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end


save("entroViscTauPlot2.pdf",fig)

gammaList=collect(1500:100:3500)
maxHList=maxHrate.(gammaList)


Tmaxtf1=finalTH.(gammaList,Ref(1.0),Ref(2))
Tmaxtf2=finalTH.(gammaList,Ref(1.0),Ref(10))
Tmaxtf3=finalTH.(gammaList,Ref(1.0),Ref(20))
Tmaxtf4=finalTH.(gammaList,Ref(1.0),Ref(5))

exp3=(gammaList,Tmaxtf1,Tmaxtf2,Tmaxtf3,Tmaxtf4)
save_object("gamma_finalT_beta0_2_10_20_5.jdl2",exp3)
(gammaList,Tmaxtf1,Tmaxtf2,Tmaxtf3,Tmaxtf4)=load("gamma_finalT_beta0_2_10_20_5.jdl2")["single_stored_object"]
entroGamma1=entro.(gammaList,Ref(1.0),Ref(2))
entroGamma2=entro.(gammaList,Ref(1.0),Ref(10))
entroGamma3=entro.(gammaList,Ref(1.0),Ref(20))
entroGamma4=entro.(gammaList,Ref(1.0),Ref(50))
exp4=(gammaList,entroGamma1,entroGamma2,entroGamma3,entroGamma4)
save_object("gamma_entro_beta0_2_10_20_50.jdl2",exp4)
(gammaList,entroGamma1,entroGamma2,entroGamma3,entroGamma4)=load("gamma_entro_beta0_2_10_20_50 (1).jdl2")["single_stored_object"]
entroGamma2=load("entroGamma2_correct.jdl2")["single_stored_object"]
entro(gammaList[5],1.0,10)
gammaList[4]


entroGammaPlot=plot(gammaList,entroGamma1,label="2",xlabel=L"$\gamma$",ylabel=L"\int \mathrm{d}t \; \nabla_\mu S^\mu \;\; [\mathrm{fm}^{-3}]",legendtitle = L"\beta_0")
plot!(gammaList,entroGamma2,label=L"10")
plot!(gammaList,entroGamma3,label=L"20")
plot!(gammaList,entroGamma4,label=L"50")

fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"\gamma",
    ylabel=L"\int \mathrm{d}t \; \nabla_\mu S^\mu \; [\mathrm{fm}^{-3}]"#,ylabel = L"The y label"
    ,xgridvisible = false,
        ygridvisible = false
)
    CairoMakie.ylims!(ax,0,1000.0)
    #CairoMakie.xlims!(ax,-0.02,0.96)
    #lines!(abs.(solution1[2,:]),solution1[1,:],label=L"0 ",color=Makie.wong_colors()[1])
    #CairoMakie.scatter!(last(solution1[2,:]),last(solution1[1,:]),marker=:circ,color=:black,markersize = 8)
    #CairoMakie.scatter!(first(solution1[2,:]),first(solution1[1,:]),marker=:circ,color=:black,markersize = 8)
    #text!(0.94, 0.0045, text = "I", align = (:center, :center))
    #text!(0.038, 0.95, text = "II", align = (:center, :center))
    #point((last(solution1[2,:]),last(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    #plot!((first(solution1[2,:]),first(solution1[1,:])),marker=:circ,mc=:black,markersize = 6)
    lines!(ax,gammaList,entroGamma1,label="2",color=Makie.wong_colors()[1])
    lines!(ax,gammaList,entroGamma2,label="10",color=Makie.wong_colors()[6])
    lines!(ax,gammaList,entroGamma3,label="20",color=Makie.wong_colors()[3])
    lines!(ax,gammaList,entroGamma4,label="50",color=Makie.wong_colors()[4])
    #axislegend(framevisible=false, position=:lt)#orientation = :horizontal)

    #lines!(ax, Tlist, pList3 ,label=L"2",color=Makie.wong_colors()[3])
    axislegend(L"2\zeta /[\tau_{\text{bulk}}(e+p)]", framevisible=false, position=:lt,hline=:left,orientation = :horizontal,nbanks=2)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end

save("totalEntroGamma.pdf",fig)

solt=getSol(gammaList[20],1,2)

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




function rfun(r0)
    T0=0.002
    mu0=.925
    check_for_transition(T0,mu0)
    u0=[T0,mu0,0.0,r0]
    f(du,u,p,t)=get_source(du,u,t,fullEOS,Walecka2,3,gammaA,7)
    problem = ODEProblem(f, u0, tspan)
    solutionT = solve(problem,AutoTsit5(Rosenbrock23(autodiff=false)),dtmax=0.01*(tspan[2]-tspan[1]))

    return solutionT
end

sol1=rfun(0.1)
sol3=rfun(0.25)
sol4=rfun(0.5)
sol7=rfun(0.6)
sol5=rfun(0.75)
sol6=rfun(0.9)

rplot=plot(sol1.t,sol1[4,:],xlabel=L"$t$ [fm/c]",ylabel=L"r",label=L"0.1",legendtitle=L"r(t=0)",xlim=(0.017,0.022))
plot!(sol3.t,sol3[4,:],label=L"0.25")
plot!(sol4.t,sol4[4,:],label=L"0.5")
plot!(sol5.t,sol5[4,:],label=L"0.75")
plot!(sol6.t,sol6[4,:],label=L"0.9")
#savefig(rplot,"PlotsThesis/rplot.pdf")


fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"$t$ [fm/c]",ylabel=L"r"
    ,xgridvisible = false,
        ygridvisible = false
)
    #CairoMakie.ylims!(ax,-3.,3.5)
    CairoMakie.xlims!(ax,0.02,0.023)
    lines!(ax, sol1.t,sol1[4,:],label=L"0.1",color=Makie.wong_colors()[1])
    lines!(ax, sol3.t,sol3[4,:],label=L"0.25",color=Makie.wong_colors()[2])
    lines!(ax, sol4.t,sol4[4,:],label=L"0.5",color=Makie.wong_colors()[3])
    lines!(ax,sol7.t,sol7[4,:],label=L"0.6",color=Makie.wong_colors()[4])
    #lines!(ax,sol6.t,sol6[4,:],label=L"0.9",color=Makie.wong_colors()[5])
    axislegend(L"r(t=0)",framevisible=false, position=:rb)#orientation = :horizontal)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end

save("rplot.pdf",fig)