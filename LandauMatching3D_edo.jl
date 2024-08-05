
using Tullio
using PolyLog
using Interpolations
using DataInterpolations
using Plots
using LinearAlgebra

@inline function WoodSaxonProfile3D(x,y,z,R,a,A)
    r=sqrt(x^2+y^2+z^2)
    A/(-8*pi*a^3*real(li3(-exp(R/a))))*1/(1+exp((r-R)/a))    

end

@inline function BoostedWoodSaxonProfile3D(x,y,z,z0,t,vel,R,a,A)

 gammaFactor=1/sqrt(1-vel^2)
 #r=sqrt(x^2+y^2+gammaFactor*(z-z0-vel*t)^2)
 return 0.938*WoodSaxonProfile3D(x,y,sqrt(gammaFactor*(z-z0-vel*t)^2),R,a,A)
end

WoodSaxonProfile3D(1,1,1,1,1,1)

rList=collect(-10:0.01:10)

a=BoostedWoodSaxonProfile3D.(0,0,rList,1,0,0.9999,3,0.1,208)
b=BoostedWoodSaxonProfile3D.(0,0,rList,-1,0,0.9999,3,0.1,208)
c=WoodSaxonProfile3D.(2.5,rList,1,3,0.1,208)

plot(rList,a)
plot!(rList,b)
plot!(c)

function findEnergy(eigensystem,metrik)   
    eigenVectors=eigensystem.vectors
    eigenValues=eigensystem.values
    energyIndex=0

    for i in 1:4
        if transpose((eigenVectors[:,i]/eigenVectors[1,i]))*metrik*(eigenVectors[:,i]/eigenVectors[1,i])<0
       # if transpose(eigenVectors[:,i])*metrik*eigenVectors[:,i]<0
            energyIndex=i
        end
        if energyIndex==0
            print("ERROR: No time-like eigenvector found")
        end
    end
    energy=-eigenValues[energyIndex]#not sure about the - in front of the fluidvelocity
    fluidVelocity=(eigenVectors[:,energyIndex]/eigenVectors[1,energyIndex])/sqrt(-transpose((eigenVectors[:,energyIndex]/eigenVectors[1,energyIndex]))*metrik*(eigenVectors[:,energyIndex]/eigenVectors[1,energyIndex]))
    velocity=fluidVelocity[4]/sqrt(1+fluidVelocity[4]^2)

    return (energy,velocity)
end

function doLandauMatchingBig3D(t,x,y,z,vel)

    #define R a A for Pb here
    R=6.6
    a=0.5
    A=208
    z0=0.5

    g=[[-1 0 0 0];[0 1 0 0];[0 0 1 0];[0 0 0 1]]
    uLeft=[1,0,0,-vel]/sqrt(1-vel^2)
    uRight=[1,0,0,vel]/sqrt(1-vel^2)
    vLeft=-vel
    vRight=vel
    gammaLeft=1/sqrt(1-vel^2)
    gammaRight=1/sqrt(1-vel^2)
    eLeft=BoostedWoodSaxonProfile3D(x,y,z,z0,t,-vel,R,a,A)
    eRight=BoostedWoodSaxonProfile3D(x,y,z,-z0,t,vel,R,a,A)
    TLeft=zeros((4,4))
    TRight=zeros((4,4))
    @tullio TLeft[i,j]=eLeft*uLeft[i]*uLeft[j]
    @tullio TRight[i,j]=eRight*uRight[i]*uRight[j]
    TColl=zeros((4,4))
    #@tullio TColl[i,j]=(TLeft[i,k]+0*TRight[i,k])*g[k,j]
    TColl=(TLeft+TRight)*g

    #return TColl
    EigenSys=eigen(TColl)

    (energy,velocity)=findEnergy(EigenSys,g)
    gamma=1/sqrt(1-velocity^2)

    #we assume that we are at T=0, so we find the initial number density by n=e/μ
    μ0=0.924#0.938
    nLeft=eLeft/μ0
    nRight=eRight/μ0


    diffusionCurrent=(nLeft*(velocity-vLeft)*gammaLeft+nRight*(velocity-vRight)*gammaRight)/(-1+velocity^2)#gamma^2*(gammaLeft*vLeft*nLeft+gammaRight*vRight*nRight-velocity*(gammaLeft*nLeft+gammaRight*nRight))
    numberDensity=(nRight*(velocity*vRight-1)*gammaRight+nLeft*(velocity*vLeft-1)*gammaLeft)/(gamma*(-1+velocity^2))#1/gamma*(gammaLeft*nLeft+gammaRight*nRight-velocity*diffusionCurrent)

    pizz=2/(3*velocity^2)*((eLeft+eRight)/(1-vRight^2)-energy/(1-velocity^2))
    piBulk=((1-velocity^2)*(eLeft+eRight)-(1-vRight^2)*energy)/(3*velocity^2*(1-vRight^2))

    #return [numberDensity,energy,velocity,diffusionCurrent,pizz,piBulk]
    return [numberDensity,energy]
end

xList=collect(-2:0.001:2)

tt=1
nx=doLandauMatchingBig3D.(tt,0,0,xList,0.999999)
nz1=doLandauMatchingBig3D.(tt,2.5,0,xList,0.999999)
#nz2=doLandauMatchingBig3D.(tt,3,0,xList,0.999999)
nz3=doLandauMatchingBig3D.(tt,5,0,xList,0.999999)
nz5=doLandauMatchingBig3D.(tt,7.5,0,xList,0.999999)
maximum(nx)/gamma^2
30.6/0.1533
e0=497

plot(xList,nx ./(gamma^2))
plot!(xList,nz1 ./(gamma^2))
#plot!(xList,nz2 ./(gamma^2))
plot!(xList,nz3 ./(gamma^2))
plot!(xList,nz5 ./(gamma^2))


plot(xList,nx)
gamma=1/sqrt(1-0.999999^2)

tlist=collect(0.5:0.001:1.5)
gamma=2960
v=sqrt(1-1/gamma^2)

nList1=zeros(length(tlist))
eList1=zeros(length(tlist))
nList2=zeros(length(tlist))
eList2=zeros(length(tlist))
nList3=zeros(length(tlist))
eList3=zeros(length(tlist))
nList4=zeros(length(tlist))
eList4=zeros(length(tlist))

z0=0.025

for i in eachindex(tlist)
    nList1[i]=doLandauMatchingBig3D(tlist[i],0,0,z0-v*tlist[i],v)[1]
    eList1[i]=doLandauMatchingBig3D(tlist[i],0,0,z0-v*tlist[i],v)[2]
    nList2[i]=doLandauMatchingBig3D(tlist[i],2.5,0,z0-v*tlist[i],v)[1]
    eList2[i]=doLandauMatchingBig3D(tlist[i],2.5,0,z0-v*tlist[i],v)[2]
    nList3[i]=doLandauMatchingBig3D(tlist[i],5,0,z0-v*tlist[i],v)[1]
    eList3[i]=doLandauMatchingBig3D(tlist[i],5,0,z0-v*tlist[i],v)[2]
    nList4[i]=doLandauMatchingBig3D(tlist[i],7.5,0,z0-v*tlist[i],v)[1]
    eList4[i]=doLandauMatchingBig3D(tlist[i],7.5,0,z0-v*tlist[i],v)[2]
end

plot(tlist,nList2)
plot!(tlist,eList2)

function get_T_mu_Landau3D(t,x,y,z,vel,eos)
    num, epsd =doLandauMatchingBig3D(t,x,y,z,vel)
    initialStateFlag=check_for_intial_state(epsd,num)
    if initialStateFlag
        return [0.0,0.938] #return T=0 and mu = mu_crit
    else
        T,mu = invertEOS(epsd,num,eos)
        return [T,mu]
    end
end

function check_for_intial_state(energyLandau,numberLandau)
    if energyLandau / numberLandau < 0.96 #check if e = μ_crit * n with some error
        return true#[0.0,0.938] #return temperature=0 and mu=mu_crit
    else
        return false
    end
end


function invertEOS(energyLandau,numberLandau,eos)
    f(u,p)=SysOfEqs(u,energyLandau,numberLandau,eos)
    problem=NonlinearProblem{false}(f,SVector{2}(.1,.5))
    #problem=NonlinearProblem{false}(f,SVector{2}(.01,1.0))
    sol=solve(problem,NewtonRaphson())#,tol = 1e-6)
    temperature = sol[1]
    chemicalPotential = sol[2]
    return [temperature,chemicalPotential]
end

function energy(T,μ,EOS)
    energy_density(T,μ,EOS)
end

function number_density(T,μ,EOS)
    pressure_derivative(T,μ,Val(0),Val(1),EOS)
end

function energy_condition(T,μ,eVal,EOS)
    energy(T,μ,EOS)-eVal
end

function number_condition(T,μ,nVal,EOS)
    number_density(T,μ,EOS)-nVal
end

function SysOfEqs(u::M,eVal,nVal,eos) where {M<:AbstractArray}
    eCond=energy_condition(u[1],u[2],eVal,eos)
    nCond=number_condition(u[1],u[2],nVal,eos)
    return SVector{2}(eCond,nCond)
end

tlist=collect(0.75:0.001:1.5)
tList1=zeros(length(tlist))
muList1=zeros(length(tlist))
tList2=zeros(length(tlist))
muList2=zeros(length(tlist))
tList3=zeros(length(tlist))
muList3=zeros(length(tlist))
tList4=zeros(length(tlist))
muList4=zeros(length(tlist))

push!(LOAD_PATH,pwd()*"/EquationofState")
using EquationsOfStates
include("hubble_equations.jl")
using NonlinearSolve
using StaticArrays


#trang=collect(0:0.0001:0.05)
#hrateplot=plot(trang, a.(trang),xlabel=L"$t$ [fm/c]",ylabel=L"$H$ [c/fm]",title=L"Hubble rate, $\sqrt{s}= 2.76$ TeV")
#savefig(hrateplot,"newhrateplot.png")
###### Construct the eos
fmGeV= 1/0.1973261 
#set up the individual eos
HRGLow=fmGeV^3*HadronResonaceGas(Maxmass=0.5,condition=waleckacondition)
HRG=fmGeV^3*HadronResonaceGas()
LQCD=fmGeV^3*LatticeQCD()
Walecka2=fmGeV^3*WaleckaModel2()


#gluing functions
function Ttrans(mu)
    0.166-0.4*(0.139*mu^2+0.053*mu^4)
end

function fTrans(T,mu)#,t)
    tanh((T-Ttrans(mu))/(0.1*Ttrans(0)))    
end

function Ttrans2(mu)
    0.1+0.28*mu-0.2*mu^2#0.1+0.8*mu-0.5*mu^2
end

function fTrans2(T,mu)#,t)
    tanh((T-Ttrans2(mu))/(0.1*Ttrans2(0)))    
end

transferFunction=Gluing(fTrans)
transferFunction2=Gluing(fTrans2)
#Combine the eos
highEOS=1/2*(1-transferFunction2)*HRG+1/2*(1+transferFunction2)*LQCD
lowEOS=HRGLow+Walecka2
fullEOS=(1/2*(1-transferFunction)*lowEOS+1/2*(1+transferFunction)*highEOS)

v
z0=0.5
tlist=collect(0.2:0.001:0.8)
tList1=zeros(length(tlist))
muList1=zeros(length(tlist))
tList2=zeros(length(tlist))
muList2=zeros(length(tlist))
tList3=zeros(length(tlist))
muList3=zeros(length(tlist))
tList4=zeros(length(tlist))
muList4=zeros(length(tlist))
for i in eachindex(tlist)
    tList1[i]=get_T_mu_Landau3D(tlist[i],0.0,0.0,z0-v*tlist[i],v,highEOS)[1]
    muList1[i]=get_T_mu_Landau3D(tlist[i],0,0,z0-v*tlist[i],v,highEOS)[2]
    tList2[i]=get_T_mu_Landau3D(tlist[i],2.5,0,z0-v*tlist[i],v,highEOS)[1]
    muList2[i]=get_T_mu_Landau3D(tlist[i],2.5,0,z0-v*tlist[i],v,highEOS)[2]
    tList3[i]=get_T_mu_Landau3D(tlist[i],5,0,z0-v*tlist[i],v,highEOS)[1]
    muList3[i]=get_T_mu_Landau3D(tlist[i],5,0,z0-v*tlist[i],v,highEOS)[2]
    tList4[i]=get_T_mu_Landau3D(tlist[i],7.5,0,z0-v*tlist[i],v,highEOS)[1]
    muList4[i]=get_T_mu_Landau3D(tlist[i],7.5,0,z0-v*tlist[i],v,highEOS)[2]
end

radius=0:0.01:22
times= range(0.2,0.5,4)
temperature_vs_r=map(t->map(r->get_T_mu_Landau3D(t,r,0,z0-v*t,v,highEOS),radius),times)


temperature_vs_r[2]

using LaTeXStrings
using CairoMakie 

plot(muList1,tList1)
plot!(muList2,tList2)
plot!(muList3,tList3)
plot!(muList4,tList4)

plot(tlist,tList1,label=L"r_\perp=0.0")
plot!(tlist,tList2,label=L"r_\perp=2.5")
plot!(tlist,tList3,label=L"r_\perp=5.0")
plot!(tlist,tList4,label=L"r_\perp=7.5")

plot(tlist,muList1)
plot!(tlist,muList2)
plot!(tlist,muList3)
plot!(tlist,muList4)

plot(collect(radius),getindex.(temperature_vs_r[4],1))
plot!(collect(radius),getindex.(temperature_vs_r[3],1))
plot!(collect(radius),getindex.(temperature_vs_r[2],1))
plot!(collect(radius),getindex.(temperature_vs_r[1],1))


xsize=300
ysize=xsize*3/4

fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"r \;\; [\mathrm{fm}/c]",#,ylabel = L"The y label"
    ylabel=L"T \;\; [\mathrm{GeV}]"
    ,xgridvisible = false,
        ygridvisible = false,   
)
#CairoMakie.ylims!(ax,-0.1,7)
CairoMakie.xlims!(ax,0,20)
    for (i,t) in enumerate(times)
    t_label=round(t;sigdigits=2)
    lines!(ax, collect(radius),getindex.(temperature_vs_r[i],1),label=L"%$t_label",color=Makie.wong_colors()[i])
    end 
    axislegend(L"t\; [\mathrm{fm}]", framevisible=false)
    
   # text!(-0.013,1.3, text=L"\frac{n}{\gamma n_0}")
    resize_to_layout!(fig)
    fig
end

save("3DTemperatureLandau.pdf",fig)


fig=with_theme(theme_latexfonts()) do
    fig = Figure(size = (xsize, ysize))
    ax = Axis(fig[1, 1],
    xlabel=L"r \;\; [\mathrm{fm}/c]",#,ylabel = L"The y label"
    ylabel=L"\mu \;\; [\mathrm{GeV}]"
    ,xgridvisible = false,
        ygridvisible = false
)
CairoMakie.xlims!(ax,0,10)
CairoMakie.ylims!(ax,0,1.5)
for (i,t) in enumerate(times)
    t_label=round(t;sigdigits=2)
    lines!(ax, collect(radius),getindex.(temperature_vs_r[i],2),label=L"%$t_label",color=Makie.wong_colors()[i])
    end 
    axislegend(L"t\; [\mathrm{fm}]", framevisible=false,position=:ct,orientation = :horizontal)
    #axislegend(L"r_\perp\; [\mathrm{fm}]", framevisible=false)
    
   # text!(-0.013,1.3, text=L"\frac{n}{\gamma n_0}")
    resize_to_layout!(fig)
    fig
end

save("3DChemicalPotLandau.pdf",fig)