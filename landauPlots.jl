using Revise
using BenchmarkTools
using Plots
using IntervalSets
using LinearAlgebra
using DifferentialEquations
using StaticArrays
using LaTeXStrings
using NonlinearSolve


4/3*pi*6.6^3
#plotting stuff
plot_font="Computer Modern"
default(fontfamily=plot_font,
        linewidth=4, framestyle=:box, label=nothing, grid=false,legendfontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12)



include("hubble_equations.jl")


hrate=get_hubble_rate()
gamma=1500#2960
v=sqrt(1-1/gamma^2)
rList=collect(-0.03:0.0001:0.03)

rList=collect(-0.025:0.00018:0.025)
nlist=zeros(length(rList))
elist=zeros(length(rList))
vlist=zeros(length(rList))
nulist=zeros(length(rList))

t=0.03
ind=1
for i in eachindex(rList)
    nlist[i]=doLandauMatchingBig(0.01,rList[i],v)[ind]
    vlist[i]=doLandauMatchingBig(0.02,rList[i],v)[ind]
    elist[i]=doLandauMatchingBig(0.025,rList[i],v)[ind]
    nulist[i]=doLandauMatchingBig(0.03,rList[i],v)[ind]
end

plot(rList,abs.(nlist),yaxis=:log)
plot!(rList,abs.(vlist))
plot!(rList,abs.(elist ).+0.00001)
plot!(rList,abs.(nulist))

maximum(abs.(elist))
n0=maximum(nlist)

d=2*6.6/1500

mypalette=[:lightgreen,:tomato,:orange,:blue]#palette(:turbo,4)

nPlot=plot(rList ,nlist ./n0,xlabel="z [fm]",ylabel=L"n/n_0",label="0.01",legendtitle="t [fm/c]",color=mypalette[1],size=(400,400*3/4),ylims=(0.0,2.2))
plot!([],label="0.02",color=mypalette[2])
plot!(rList,vlist ./(n0*gamma) ,linestyle=:dash,color=mypalette[2])
plot!([],label="0.025",color=mypalette[3])
plot!(rList,elist/(n0*gamma) ,linestyle=:dash,color=mypalette[3])
plot!([],label="0.03",color=mypalette[4])
fontsize=16
annotate!( -0.0225,  1.85, text(L"a)", fontsize))
annotate!(-0.0075,1.5, text(L"\frac{n}{\gamma n_0}",14))
#plot!(rList,nulist,label="0.03")
#savefig(nuplot,"PlotsPaper/nPlot.pdf")
baseEnergy=stack(doLandauMatchingBig.(Ref(0.0),rList,Ref(v)),dims=1)[:,2]
baseNumber=stack(doLandauMatchingBig.(Ref(0.0),rList,Ref(v)),dims=1)[:,1]
plot(rList,baseEnergy)
plot!(rList,baseNumber)
n0=maximum(baseNumber)
e0=maximum(baseEnergy)

lan0=stack(doLandauMatchingBig.(Ref(0.0),rList,Ref(v)),dims=1)[:,1]

nu1=stack(doLandauMatchingBig.(Ref(0.01),rList,Ref(v)),dims=1)[:,4]
nu2=stack(doLandauMatchingBig.(Ref(0.02),rList,Ref(v)),dims=1)[:,4]
nu3=stack(doLandauMatchingBig.(Ref(0.03),rList,Ref(v)),dims=1)[:,4]

nuPlot=plot(rList,nu1 ./(gamma*n0),linestyle=:dash,ylabel=L"\nu/n_0",xlabel="z [fm]",color=mypalette[1])
plot!(rList, nu2 ./(n0*gamma),linestyle=:dash,color=mypalette[2])
plot!(rList, nu3 ./(n0*gamma),linestyle=:dash,color=mypalette[4])
annotate!(-0.0225,0.4,text(L"b)", fontsize))
annotate!(0.0125,0.3, text(L"\frac{\nu}{\gamma n_0}",14))

piz1=stack(doLandauMatchingBig.(Ref(0.01),rList,Ref(v)),dims=1)[:,5]
piz2=stack(doLandauMatchingBig.(Ref(0.02),rList,Ref(v)),dims=1)[:,5]
piz3=stack(doLandauMatchingBig.(Ref(0.025),rList,Ref(v)),dims=1)[:,5]

pizPlot=plot(rList,piz1 ./e0,color=mypalette[1],ylabel=L"\pi^{zz}/\epsilon_0",xlabel="z [fm]")
plot!(rList, piz2 ./(e0),color=mypalette[2])
plot!(rList, piz3 ./(e0*gamma^2),linestyle=:dot,color=mypalette[3])
annotate!(-0.0225,1.5*10^6,text(L"c)", fontsize))
annotate!(0.003,1.2*10^6, text(L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}",14))


pib1=stack(doLandauMatchingBig.(Ref(0.01),rList,Ref(v)),dims=1)[:,6]
pib2=stack(doLandauMatchingBig.(Ref(0.02),rList,Ref(v)),dims=1)[:,6]
pib3=stack(doLandauMatchingBig.(Ref(0.025),rList,Ref(v)),dims=1)[:,6]

pibPlot=plot(rList,pib1 ./e0,color=mypalette[1],ylabel=L"\Pi/\epsilon_0",xlabel="z [fm]",xlimits=(-0.0075,0.0075))
plot!(rList, pib2 ./(e0),color=mypalette[2])
plot!(rList, pib3 ./(e0*gamma^2),color=mypalette[3],linestyle=:dot)

finalPlot=plot(nPlot,nuPlot,pizPlot,pibPlot,layout=4)
