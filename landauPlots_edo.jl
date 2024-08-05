

using IntervalSets
using LinearAlgebra
using DifferentialEquations
using StaticArrays
using LaTeXStrings
using NonlinearSolve
using CairoMakie 



include("hubble_equation.jl")


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


maximum(abs.(elist))
n0=maximum(nlist)

d=2*6.6/1500


xsize=300
ysize=xsize*3/4

fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"z [\mathrm{fm}]",#,ylabel = L"The y label"
    ylabel=L"n/n_0"
    ,xgridvisible = false,
        ygridvisible = false
)
CairoMakie.ylims!(ax,0.,5)
    lines!(ax, rList, nlist ./(n0),label=L"0.01 \frac{n}{ n_0} ",color=Makie.wong_colors()[1])
    lines!(ax, rList,vlist ./(n0*gamma) ,label=L"0.02 \frac{n}{\gamma n_0}",color=Makie.wong_colors()[2])
    lines!(ax, rList,elist/(n0*gamma) ,label=L"0.025 \frac{n}{\gamma n_0} ",color=Makie.wong_colors()[3])
    axislegend(L"t\; [\mathrm{fm/c}]", framevisible=false)
    
   # text!(-0.013,1.3, text=L"\frac{n}{\gamma n_0}")
    resize_to_layout!(fig)
    fig
end
save("density.pdf",fig)







baseEnergy=stack(doLandauMatchingBig.(Ref(0.0),rList,Ref(v)),dims=1)[:,2]
baseNumber=stack(doLandauMatchingBig.(Ref(0.0),rList,Ref(v)),dims=1)[:,1]

n0=maximum(baseNumber)
e0=maximum(baseEnergy)

lan0=stack(doLandauMatchingBig.(Ref(0.0),rList,Ref(v)),dims=1)[:,1]

nu1=stack(doLandauMatchingBig.(Ref(0.01),rList,Ref(v)),dims=1)[:,4]
nu2=stack(doLandauMatchingBig.(Ref(0.02),rList,Ref(v)),dims=1)[:,4]
nu3=stack(doLandauMatchingBig.(Ref(0.03),rList,Ref(v)),dims=1)[:,4]

fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"z [\mathrm{fm}]",#,ylabel = L"The y label"
    ylabel=L"\nu/(n_0\gamma)"
    ,xgridvisible = false,
        ygridvisible = false
)
CairoMakie.ylims!(ax,-0.6,1.2)
    lines!(ax, rList,nu1 ./(gamma*n0),label=L"0.01 ",color=Makie.wong_colors()[1])
    lines!(ax, rList, nu2 ./(n0*gamma) ,label=L"0.02 ",color=Makie.wong_colors()[2])
    lines!(ax, rList, nu3 ./(n0*gamma) ,label=L"0.03 ",color=Makie.wong_colors()[4])
    axislegend(L"t\; [\mathrm{fm/c}]", framevisible=false)
    
    #text!(-0.015,0.3, text=L"\frac{\nu}{\gamma n_0}")
    resize_to_layout!(fig)
    fig
end
save("nu.pdf",fig)


piz1=stack(doLandauMatchingBig.(Ref(0.01),rList,Ref(v)),dims=1)[:,5]
piz2=stack(doLandauMatchingBig.(Ref(0.02),rList,Ref(v)),dims=1)[:,5]
piz3=stack(doLandauMatchingBig.(Ref(0.025),rList,Ref(v)),dims=1)[:,5]


fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"z [\mathrm{fm}]",#,ylabel = L"The y label"
    ylabel=L"\pi^{zz}/(\epsilon_0\gamma^2)"
    ,xgridvisible = false,
        ygridvisible = false
)
    #CairoMakie.ylims!(ax,-0.1,0.5)
    lines!(ax, rList, piz1 ./(e0*gamma^2),label=L"0.01 ",color=Makie.wong_colors()[1])
    lines!(ax, rList, piz2 ./(e0*gamma^2),label=L"0.02 ",color=Makie.wong_colors()[2])
    lines!(ax, rList, piz3 ./(e0*gamma^2)/1e6 ,label=L"0.025  / 10^{6}",color=Makie.wong_colors()[3])
    axislegend(L"t\; [\mathrm{fm/c}]", framevisible=false)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end
save("piz.pdf",fig)





rList=collect(-0.025:0.00018:0.04)
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
pib1=stack(doLandauMatchingBig.(Ref(0.01),rList,Ref(v)),dims=1)[:,6]
pib2=stack(doLandauMatchingBig.(Ref(0.02),rList,Ref(v)),dims=1)[:,6]
pib3=stack(doLandauMatchingBig.(Ref(0.025),rList,Ref(v)),dims=1)[:,6]

fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"z [\mathrm{fm}]",#,ylabel = L"The y label"
    ylabel=L"\;\Pi/\epsilon_0"
    ,xgridvisible = false,
        ygridvisible = false
)
    #CairoMakie.ylims!(ax,-0.1,noth)
    #CairoMakie.xlims!(ax,nothing,0.04)
    lines!(ax, rList, pib1 ./(e0*gamma^2),label=L"0.01 ",color=Makie.wong_colors()[1])
    lines!(ax, rList, pib2 ./(e0)/1e5,label=L"0.02 \;\frac{\Pi}{(\epsilon_0 )10^{5}} ",color=Makie.wong_colors()[2])
    lines!(ax, rList, pib3 ./(e0*gamma^2)/1e5 ,label=L"0.025 \;\frac{\Pi}{(\epsilon_0 \gamma^2)10^{5}} ",color=Makie.wong_colors()[3])
    axislegend(L"t\; [\mathrm{fm/c}]", framevisible=false)
    
    #text!(-0.015,0.3, text=L"\frac{\pi^{zz}}{\gamma^2 \epsilon_0}")
    resize_to_layout!(fig)
    fig
end
save("pib.pdf",fig)



