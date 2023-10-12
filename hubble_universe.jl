using Revise
using BenchmarkTools
using Plots
using IntervalSets
using LinearAlgebra
using DifferentialEquations

push!(LOAD_PATH,pwd()*"/EquationofState")
using EquationsOfStates

HRGLow=HadronResonaceGas(Maxmass=0.5,condition=waleckacondition)
HRG=HadronResonaceGas()
LQCD=LatticeQCD()
WaleckaG=WaleckaModel1()

function Ttrans(mu)
    0.166-0.4*(0.139*mu^2+0.053*mu^4)
end

function fTrans(T,mu)
    tanh((T-Ttrans(mu))/(0.1*Ttrans(0)))    
end

transferFunction=Gluing(fTrans)

lowEOS=HRGLow+WaleckaG
#this part is for getting rid of the pole in the lattice
function highInterpol(T,μ)
    if 0.06 < T <= 0.08
        return 50*T-3
    elseif T > 0.08
        return 1 
    else
        return 0
    end
end

highInterpolFun=Gluing(highInterpol)

highEOS=LQCD*highInterpolFun+(1-highInterpolFun)*HRG 

fullEOS=1/2*(1-transferFunction)*lowEOS+1/2*(1+transferFunction)*highEOS

#idealQCD=IdealQCD()
#fullEOS=idealQCD

fluidproperties=FluidProperties(fullEOS
,EquationsOfStates.ZeroViscosity()
,EquationsOfStates.ZeroBulkViscosity()
,EquationsOfStates.ZeroDiffusion())

#these are the thermal fisher metric and its inverse, needed for the eom
function Thermal_Fisher_metric(T,mu,eos)
    p=pressure(T,mu,eos)
    dTp=pressure_derivative(T,mu,Val(1),Val(0),eos)
    dmup=pressure_derivative(T,mu,Val(0),Val(1),eos)
    dTdTp=pressure_derivative(T,mu,Val(2),Val(0),eos)
    dTdmup=pressure_derivative(T,mu,Val(1),Val(1),eos)
    dmudmup=pressure_derivative(T,mu,Val(0),Val(2),eos)

    return [T^3*dTdTp+2*mu*T^2*dTdmup+mu^2*T*dmudmup T*(T*dTdmup+mu*dmudmup); T*(T*dTdmup+mu*dmudmup) T*dmudmup]
end

function inverse_Thermal_Fisher_metric(T,mu,eos)
    thermal_fisher=Thermal_Fisher_metric(T,mu,eos)
    return inv(thermal_fisher)
end

#now we get the matrices, in this case only the source temperature
function get_source(source,u,time,x::FluidProperties)
    T = 1/u[1]
    mu = u[2]/u[1]
    piBulk = u[3]
    @show(T,mu)
    thermo = thermodynamic(T,mu,x)
    p = thermo.pressure
    dtp=thermo.pressure_derivative[1]
    dmp=thermo.pressure_derivative[2]
    dtdtp=thermo.pressure_hessian[1]
    dtdmp=thermo.pressure_hessian[2]
    energy = T*dtp + mu*dmp - p

    tauB=τ_bulk(T,mu,thermo,x.bulk)
    zeta=1.0#bulk_viscosity(T,mu,thermo,x.bulk)
    @show(tauB,zeta)
    H=expansion_rate(time)

    inverse_fisher_metric=inverse_Thermal_Fisher_metric(T,mu,x)
    thermodynamic_Source =[ -3*H*(energy+p+piBulk), 3*H*dmp]
    full_thermo_source = inverse_fisher_metric*thermodynamic_Source



    source[1] = full_thermo_source[1]
    source[2] = full_thermo_source[2]
    source[3] = 1/tauB*(piBulk+3*H*zeta)
    return -source
end

function expansion_rate(time)
    H0=0.1
  #H0 = (time -2.5)/10
  #H0 = 0.1*sin(time*pi)
   #end
    #H0=1.
   return H0
end

source=[0.0, 0.0, 0.0]

#get_source(source,(10,1,1),1,fluidproperties)

#now we define the problem and try to solve it
T0=.8
mu0=2
u0=[1/T0,mu0/T0,0.0]
u0m=[1/T0,-mu0/T0,0.0]
tspan=(0.0,25.)

f(u,p,t)=get_source(source,u,t,fluidproperties)
problem = ODEProblem(f, u0, tspan)
solution= solve(problem,Tsit5(),dtmax=0.1)
problemm = ODEProblem(f, u0m, tspan)
solutionm= solve(problemm,Tsit5(),dtmax=0.1)
plot(solution)
plot(solutionm)
#plot( solution[2] ./solution[1])
#plot!( 1 ./ solution[1])

#hubbleplot=plot(solution.t , 1 ./solution[1,:],xlabel="time",ylabel="fields [GeV]",label="T",title="Contracting universe, H=5, ζ=0",yaxis=:log)
#plot!(solution.t, solution[2,:] ./ solution[1,:],label="μ")
#plot!(solution.t, solution[3,:] ,label="πB")
#plot!(solution.t, expansion_rate.(solution.t),label="H")
#savefig(hubbleplot,"contraction_ideal_long.png")


hubbleplot=plot(solution.t , 1 ./solution[1,:],xlabel="time",ylabel="fields [GeV]",label="T, +μ sol",title="Expanding fullEOS universe, ζ=0,±μ0",linewidth=3)
plot!(solution.t,solution[2,:] ./ solution[1,:],label="μ, +μ sol",linewidth=3)
plot!(solutionm.t,-solutionm[2,:] ./ solutionm[1,:],label="-μ, -μ sol",linestyle=:dashdot,linewidth=3)
plot!(solutionm.t , 1 ./solutionm[1,:],label="T, -μ sol",linestyle=:dashdot,linewidth=3)
plot!(solution.t, solution[3,:] ,label="πB")
plot!(solutionm.t, solutionm[3,:] ,label="πB")


 #hline!([T0],label="initial value")
#hline!([0.08])
#savefig(hubbleplot,"contracting_pm_mu.png")

ndplot=plot(solution.t , pressure.(1 ./solution[1,:], solution[2,:] ./ solution[1,:],Ref(fullEOS) ),yaxis=:log,xlabel="time",label="number density",title="Contracting universe, H=5, ζ=0")
plot!(solution.t , 1 ./solution[1,:] .* pressure_derivative.(1 ./solution[1,:], solution[2,:] ./ solution[1,:],Val(1),Val(0),Ref(fullEOS)) .+solution[2,:] ./ solution[1,:] .* pressure_derivative.(1 ./solution[1,:], solution[2,:] ./ solution[1,:],Val(0),Val(1),Ref(fullEOS)) .-pressure.(1 ./solution[1,:], solution[2,:] ./ solution[1,:],Ref(fullEOS)),yaxis=:log,xlabel="time",label="energy density")
#savefig(ndplot,"contraction_ideal_nd_long.png")
