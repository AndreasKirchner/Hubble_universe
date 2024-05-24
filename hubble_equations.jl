### Landau matching
using Tullio
using PolyLog
using Interpolations
using DataInterpolations

#include("EquationofState/Polyloghack/Polylogarithmshack.jl")

#using .PolylogarithmsHack: li3

function doLandauMatching(t,z,vel)

    #define R a A for Pb here
    R=6.6
    a=0.5
    A=208
    z0=0.025

    g=[[-1 0 0 0];[0 1 0 0];[0 0 1 0];[0 0 0 1]]
    uLeft=[1,0,0,-vel]/sqrt(1-vel^2)
    uRight=[1,0,0,vel]/sqrt(1-vel^2)
    vLeft=-vel
    vRight=vel
    gammaLeft=1/sqrt(1-vel^2)
    gammaRight=1/sqrt(1-vel^2)
    eLeft=BoostedWoodSaxonProfile(z,z0,t,-vel,R,a,A)
    eRight=BoostedWoodSaxonProfile(z,-z0,t,vel,R,a,A)
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
    μ0=0.938
    nLeft=eLeft/μ0
    nRight=eRight/μ0


    diffusionCurrent=gamma^2*(gammaLeft*vLeft*nLeft+gammaRight*vRight*nRight-velocity*(gammaLeft*nLeft+gammaRight*nRight))
    numberDensity=1/gamma*(gammaLeft*nLeft+gammaRight*nRight-velocity*diffusionCurrent)

    return numberDensity
    #return velocity
end

function doLandauMatchingBig(t,z,vel)

    #define R a A for Pb here
    R=6.6
    a=0.5
    A=208
    z0=0.025

    g=[[-1 0 0 0];[0 1 0 0];[0 0 1 0];[0 0 0 1]]
    uLeft=[1,0,0,-vel]/sqrt(1-vel^2)
    uRight=[1,0,0,vel]/sqrt(1-vel^2)
    vLeft=-vel
    vRight=vel
    gammaLeft=1/sqrt(1-vel^2)
    gammaRight=1/sqrt(1-vel^2)
    eLeft=BoostedWoodSaxonProfile(z,z0,t,-vel,R,a,A)
    eRight=BoostedWoodSaxonProfile(z,-z0,t,vel,R,a,A)
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


    diffusionCurrent=gamma^2*(gammaLeft*vLeft*nLeft+gammaRight*vRight*nRight-velocity*(gammaLeft*nLeft+gammaRight*nRight))
    numberDensity=1/gamma*(gammaLeft*nLeft+gammaRight*nRight-velocity*diffusionCurrent)

    return [numberDensity,energy,velocity,diffusionCurrent]
end

function doLandauMatchingTMu(t,z,vel)

    #define R a A for Pb here
    R=6.6
    a=0.5
    A=208
    z0=0.025

    g=[[-1 0 0 0];[0 1 0 0];[0 0 1 0];[0 0 0 1]]
    uLeft=[1,0,0,-vel]/sqrt(1-vel^2)
    uRight=[1,0,0,vel]/sqrt(1-vel^2)
    vLeft=-vel
    vRight=vel
    gammaLeft=1/sqrt(1-vel^2)
    gammaRight=1/sqrt(1-vel^2)
    eLeft=BoostedWoodSaxonProfile(z,z0,t,-vel,R,a,A)
    eRight=BoostedWoodSaxonProfile(z,-z0,t,vel,R,a,A)
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
    μ0=0.938
    nLeft=eLeft/μ0
    nRight=eRight/μ0


    diffusionCurrent=gamma^2*(gammaLeft*vLeft*nLeft+gammaRight*vRight*nRight-velocity*(gammaLeft*nLeft+gammaRight*nRight))
    numberDensity=1/gamma*(gammaLeft*nLeft+gammaRight*nRight-velocity*diffusionCurrent)

    #return [numberDensity,energy,velocity]
end

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

@inline function WoodSaxonProfile(r,R,a,A)

    A/(-8*pi*a^3*real(li3(-exp(R/a))))*1/(1+exp((r-R)/a))    

end

@inline function BoostedWoodSaxonProfile(r,z0,t,vel,R,a,A)
 gammaFactor=1/sqrt(1-vel^2)
 return 0.938*WoodSaxonProfile(abs(gammaFactor*((r-z0)-vel*t)),R,a,A)
end


function get_hubble_rate(;gamma=2960)
    v=sqrt(1-1/gamma^2)
    tList=range(0.0,0.05,1000)
    nList=doLandauMatching.(tList,0.025 .- v .*tList,v)
    nfunction=AkimaInterpolation(nList,tList,extrapolate=true)
    dtn=DataInterpolations.derivative.(Ref(nfunction),tList)
    hubble=-dtn ./ ( 3 .*nList)
    return CubicSpline(hubble,tList,extrapolate=true)
end


### Functions for the eom
function speed_of_sound(T,mu,x)
    thermo=thermodynamic(T,mu,x)
    p = thermo.pressure
    dtp=thermo.pressure_derivative[1]
    dmp=thermo.pressure_derivative[2]
    dtdtp=thermo.pressure_hessian[1]
    dtdmp=thermo.pressure_hessian[2]
    dmdmp=thermo.pressure_hessian[3]
    energy = T*dtp + mu*dmp - p
    return (dmp^2*dtdtp-2*dtp*dmp*dtdmp+dtp^2*dmdmp)/((energy+p)*(dtdtp*dmdmp-dtdmp^2))
end


#these are the thermal fisher metric and its inverse, needed for the eom
function Thermal_Fisher_metric(T,mu,t,eos)
    thermoOne=thermodynamicPhaseOne(T,mu,eos)
    thermoTwo=thermodynamicPhaseTwo(T,mu,eos)

    dTdTpOne= thermoOne.pressure_hessian[1]
    dTdmupOne= thermoOne.pressure_hessian[2]
    dmudmupOne= thermoOne.pressure_hessian[3]

    
    dTdTpTwo= thermoTwo.pressure_hessian[1]
    dTdmupTwo= thermoTwo.pressure_hessian[2]
    dmudmupTwo= thermoTwo.pressure_hessian[3]

    fisherOne=@SMatrix [-(T^3*dTdTpOne+2*mu*T^2*dTdmupOne+mu^2*T*dmudmupOne) T*(T*dTdmupOne+mu*dmudmupOne); -T*(T*dTdmupOne+mu*dmudmupOne) T*dmudmupOne]
    fisherTwo=@SMatrix[-(T^3*dTdTpTwo+2*mu*T^2*dTdmupTwo+mu^2*T*dmudmupTwo) T*(T*dTdmupTwo+mu*dmudmupTwo); -T*(T*dTdmupTwo+mu*dmudmupTwo) T*dmudmupTwo]
    #thermo=thermodynamic(T,mu,t,eos)
   #p = thermo.pressure
   #dTp=thermo.pressure_derivative[1]
   #dmup=thermo.pressure_derivative[2]
   #dTdTp=thermo.pressure_hessian[1]
   #dTdmup=thermo.pressure_hessian[2]
   #dmudmup=thermo.pressure_hessian[3]
    #@show T
    #@show mu
    #@show dTdTp
    #return @SMatrix [T^3*dTdTp+2*mu*T^2*dTdmup+mu^2*T*dmudmup T*(T*dTdmup+mu*dmudmup); T*(T*dTdmup+mu*dmudmup) T*dmudmup]
    #return @SMatrix [-(T^3*dTdTp+2*mu*T^2*dTdmup+mu^2*T*dmudmup) T*(T*dTdmup+mu*dmudmup); -T*(T*dTdmup+mu*dmudmup) T*dmudmup]
    return (1-t)*fisherOne + t*fisherTwo
    #return t*fisherOne + (1-t)*fisherTwo
end

function Thermal_Fisher_metric(T,mu,eos)
    thermo=thermodynamic(T,mu,eos)
    p = thermo.pressure
    dTp=thermo.pressure_derivative[1]
    dmup=thermo.pressure_derivative[2]
    dTdTp=thermo.pressure_hessian[1]
    dTdmup=thermo.pressure_hessian[2]
    dmudmup=thermo.pressure_hessian[3]
    #@show T
    #@show mu
    #@show dTdTp
    #return @SMatrix [T^3*dTdTp+2*mu*T^2*dTdmup+mu^2*T*dmudmup T*(T*dTdmup+mu*dmudmup); T*(T*dTdmup+mu*dmudmup) T*dmudmup]
    return @SMatrix [(T^3*dTdTp+2*mu*T^2*dTdmup+mu^2*T*dmudmup) T*(T*dTdmup+mu*dmudmup); T*(T*dTdmup+mu*dmudmup) T*dmudmup]
end


function inverse_Thermal_Fisher_metric(T,mu,t,eos)
    thermal_fisher=Thermal_Fisher_metric(T,mu,t,eos)
    #@show det(thermal_fisher)
    return inv(thermal_fisher)
end

function inverse_Thermal_Fisher_metric(T,mu,eos)
    thermal_fisher=Thermal_Fisher_metric(T,mu,eos)
    #@show det(thermal_fisher)
    return inv(thermal_fisher)
end

#using Polylogarithms
#using DifferentialEquations
#using ForwardDiff:Dual,value,partials
#using GSL

#GSL.sf_fermi_dirac_3half(d::Dual{T,V,N}) where {T,V,N} = Dual{T}(GSL.sf_fermi_dirac_3half(value(d)),GSL.sf_fermi_dirac_half(value(d)) * partials(d))
#GSL.sf_fermi_dirac_half(d::Dual{T,V,N}) where {T,V,N} = Dual{T}(GSL.sf_fermi_dirac_half(value(d)),GSL.sf_fermi_dirac_mhalf(value(d)) * partials(d))
#GSL.sf_fermi_dirac_mhalf(d::Dual{T,V,N}) where {T,V,N} = Dual{T}(GSL.sf_fermi_dirac_mhalf(value(d)),real(-Polylogarithms.polylog(-1/2,-exp(value(d)))) * partials(d))


#now we get the matrices, in this case only the source temperature
#function get_source(source,u,time,x::FluidProperties,zetaM)
#function get_source(source,u,time,x::EquationOfState,zetaM)
function get_source_no_transition(du,u,time,x,zetaMax,gammaE,beta1)
    #T = 1/u[1]
    #mu = u[2]/u[1]
    T= max(u[1],0.0001)
    mu=max(u[2],0.0001)
    piBulk = u[3]
    mixingParam = u[4]
    fmGeV= 1/0.1973261 
    thermo = thermodynamic(T,mu,x)
    p = thermo.pressure
    dtp=thermo.pressure_derivative[1]
    dmp=thermo.pressure_derivative[2]
    dtdtp=thermo.pressure_hessian[1]
    dtdmp=thermo.pressure_hessian[2]
    dmdmp=thermo.pressure_hessian[3]
    energy = T*dtp + mu*dmp - p
    #viscosity
    zetaParam=dtp/(1+((sqrt(T^2+0.188172^2*mu^2)-0.175)/0.024)^2)
    zeta=zetaMax*zetaParam#dtp*zetaMax/(1+((sqrt(T^2+0.188172^2*mu^2)-0.175)/0.024)^2)
    tauB=2*beta1*zetaParam#zeta/(5)
    #hubble rate
    #hrate=10*sin(100*pi*time)#get_hubble_rate(gamma=gammaE)
    hrate=get_hubble_rate(gamma=gammaE)
    #H=100*sin(100*pi*time)#hrate(time)
    H=hrate(time)
    #H=expansion_rate(time)

    inverse_fisher_metric=inverse_Thermal_Fisher_metric(T,mu,x)
    #@show inverse_fisher_metric
    thermodynamic_Source =@SVector [ 3*H*(energy+p+piBulk), 3*H*dmp]
    VarChangeMatrix = @SMatrix [T^2 0 ; mu*T T]
    full_thermo_source = VarChangeMatrix*inverse_fisher_metric*thermodynamic_Source
    #@show time
    #source =SVector{3}(full_thermo_source[1],full_thermo_source[2],1/tauB*(piBulk+3*H*zeta))
    du[1] = -full_thermo_source[1]#-3*H*(-dmdmp*(piBulk+T*dtp)+T*dmp*dtdmp)/(T*(dtdmp^2-dtdtp*dmdmp))#full_thermo_source[1]
    du[2] = -full_thermo_source[2]#-3*H*((piBulk+T*dtp)*dtdmp-T*dmp*dtdtp)/(T*(dtdmp^2-dtdtp*dmdmp))#full_thermo_source[2]
    du[3] = (-piBulk -3/2*piBulk*H*tauB -3*H*zeta+ piBulk*tauB*du[1]/(2*T))/tauB  #-1/tauB*(3*H*zeta+(1+3*H*tauB/2-tauB/2*du[1]/T)*piBulk)#-1/tauB*(piBulk+3*H*zeta+1/2*tauB*3*H*piBulk)
    #du[3] = -1/tauB*(3*H*zeta+(1+3*H*tauB/2-tauB/2*du[1]/T)*piBulk)
    du[4] = zero(typeof(u[4])) #this is away from the transition
    #return -source
    #return nothing
    #return @SVector [-full_thermo_source[1],full_thermo_source[2],-1/tauB*(piBulk+3*H*zeta+1/2*tauB*3*H*piBulk)]
    #return SVector{3}(-full_thermo_source[1],full_thermo_source[2],-1/tauB*(piBulk+3*H*zeta+1/2*tauB*3*H*piBulk))
end

function conserved_charges_difference(T,mu,t,eos)
    thermoOne=thermodynamicPhaseOne(T,mu,eos)
    thermoTwo=thermodynamicPhaseTwo(T,mu,eos)
    numberDensityOne=thermoOne.pressure_derivative[2]
    numberDensityTwo=thermoTwo.pressure_derivative[2]
    numberDensity=numberDensityOne-numberDensityTwo
    energyDensityOne=T*thermoOne.pressure_derivative[1]+mu*thermoOne.pressure_derivative[2]-thermoOne.pressure
    energyDensityTwo=T*thermoTwo.pressure_derivative[1]+mu*thermoTwo.pressure_derivative[2]-thermoTwo.pressure
    energyDensity=energyDensityOne-energyDensityTwo
    return @SVector [-energyDensity,-numberDensity]
end

function normal_vector_transition_line(T)
    #mu0=0.92425
    mu1=-41.2324123
    mu2=-53938.75555664
    #return @SVector [-1, 2*T*(mu1+2*mu2*T^2)]    
    return @SVector [ 2*T*(mu1+2*mu2*T^2),-1]    
end


function get_source_transition(du,u,time,x,zetaMax,gammaE,beta1)
    #T = 1/u[1]
    #mu = u[2]/u[1]
    T= max(u[1],0.0001)
    mu=max(u[2],0.0001)
    piBulk = u[3]
    mixingParam=u[4]
    #@show(T,mu)
    thermo = thermodynamic(T,mu,mixingParam,x)#put x as walecka model here for phase transition
    #thermo = thermodynamic(T,mu,mixingParam,x)
    p = thermo.pressure
    dtp=thermo.pressure_derivative[1]
    dmp=thermo.pressure_derivative[2]
    dtdtp=thermo.pressure_hessian[1]
    dtdmp=thermo.pressure_hessian[2]
    energy = T*dtp + mu*dmp - p

    #viscosity
    zetaParam=dtp/(1+((sqrt(T^2+0.188172^2*mu^2)-0.175)/0.024)^2)
    zeta=zetaMax*zetaParam
    tauB=2*beta1*zetaParam
    #hubble rate
    hrate=get_hubble_rate(gamma=gammaE)
    H=hrate(time)
    normalVec=normal_vector_transition_line(T)
    #@show normalVec
    
    thermodynamic_Source =@SVector [ 3*H*(energy+p+piBulk), 3*H*dmp]
    VarChangeMatrix = @SMatrix [T^2 0 ; mu*T T]
    

    inverseFisher=inverse_Thermal_Fisher_metric(T,mu,mixingParam,x)
    #thermodynamic_Source =@SVector [ 3*H*(energy+p+piBulk), 3*H*dmp]
    #inverse_fisher_metric=inverse_Thermal_Fisher_metric(T,mu,mixingParam,x)
    conservedChargeDiff=conserved_charges_difference(T,mu,mixingParam,x)
    @show conservedChargeDiff
    tDer=-dot(normalVec,inverseFisher,thermodynamic_Source)/dot(normalVec,inverseFisher,conservedChargeDiff)
    #@show inverseFisher
    @show tDer
    #full_thermo_source=inverseFisher*(thermodynamic_Source+conservedChargeDiff*tDer)
    full_thermo_source = VarChangeMatrix*(inverseFisher*(thermodynamic_Source+conservedChargeDiff*tDer))
    du[1] = -full_thermo_source[1]
    @show du[1]
    du[2] = -full_thermo_source[2]
    @show du[2]
    du[3] =  (-piBulk -3/2*piBulk*H*tauB -3*H*zeta+ piBulk*tauB*du[1]/(2*T))/tauB
    du[4] = tDer
    #@show tDer
    #@show mixingParam
    

    
end


#define phase transition line
function transitionMu(T)
    if 0.0<=T<0.021
        return 0.92424-42.12324123*T^2-53938.75555664*T^4
    else
        return -one(T)
    end
end


#define area around the transition and check if T,mu is at the transition area
function check_for_transition(T,mu)
    eps=0.001 #size of transition area GeV
    transMu=transitionMu(T)
    if transMu<0
        return false
    end
    if transMu-eps < mu < transMu+eps
        return true
    else
        return false
    end
end

function get_source(source,u,time,x,wal,zetaMax,gammaE,beta1)
    #check if we are on transition line
    #T = 1/u[1]
    #μ = u[2]/u[1]
    T=u[1]
    μ=u[2]
    transitionFlag=check_for_transition(T,μ)
    #@show T
    #@show μ
    r=u[4]
   #if 0<r<1
   # if transitionFlag && 0<r<1
    #    println("Transition activated")
     #   return get_source_transition(source,u,time,wal,zetaMax,gammaE,beta1)
    #else
        return get_source_no_transition(source,u,time,x,zetaMax,gammaE,beta1)
    #end
end