#include("thermodynamic.jl") 
#include("Polyloghack/Polylogarithmshack.jl")
#using .PolylogarithmsHack: polylog
#using .PolylogarithmsHack: polylog_exp_norm
using StaticArrays
using StaticNumbers
#using ForwardDiff:
using NonlinearSolve
using SimpleNonlinearSolve
using GSL
using NLsolve

struct WaleckaModel1{T} <:EquationOfState
    m_nucleon::T 
    m_sigma::T 
    m_omega::T 
    g_sigma::T 
    g_omega::T 
    sigma_0::T 
    omega_0::T 

    function WaleckaModel1() 
        
        m_nucleon=0.939
        m_sigma=0.5
        m_omega=0.783
        g_sigma=9.69719
        g_omega=13.2025
        sigma_0=0.0444258
        omega_0= -0.026538

        new{typeof(m_nucleon)}(m_nucleon,m_sigma,m_omega,g_sigma,g_omega,sigma_0,omega_0)
    end
end

function Base.show(io::IO, z::WaleckaModel1{T}) where T 
    print(io,"Non Relativistic Walecka Model: σ=",z.sigma_0," GeV, ω₀=",z.omega_0," GeV"  )
end


#function Base.show(io::IO, ::MIME"text/plain", z::HadronResonaceGas) where{T}
#    min ,max = extrema(z.particle_list.Mass)
#    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ [",min,",",max,"] GeV\n " )
#    for part in z.particle_list
#        print(io,part.Name," mass=",part.Mass,"\n")
#    end 
#end



@inline function system_of_equation(u::M,T,μ,wal::WaleckaModel1{N}) where {M<:AbstractArray,N}
    #m_nucleon,m_sigma,m_omega,g_sigma,g_omega=p
 
    mustar=μ +wal.g_omega*u[2]
    mstar=wal.m_nucleon-wal.g_sigma*u[1]
    #@show u[2]
    #@show u[1]
    #println(mstar)
    nump=mustar-mstar
    numm=-mustar-mstar
   # @show nump
    #zplus=-exp((mustar-mstar)/T)
    #zminus=-exp((-mustar-mstar)/T)
    
    #mtfactor=(mstar*T/(2*pi))^(3/2)

    #println(mstar)
    mfactor=(mstar/(2*pi))^(3/2)

    #println(nump/T)
    n_lie32p=-sf_fermi_dirac_half(nump/T)#polylog_exp_norm(3/2,nump,T)
    n_lie52p=-sf_fermi_dirac_3half(nump/T)#polylog_exp_norm(5/2,nump,T)
    #@show nump/T    

    n_lie32m=-sf_fermi_dirac_half(numm/T)#polylog_exp_norm(3/2,numm,T)
    n_lie52m=-sf_fermi_dirac_3half(numm/T)#polylog_exp_norm(5/2,numm,T)


    #lie32p=polylog_exp_norm(3/2,zplus)
    #lie52p=polylog_exp_norm(5/2,zplus)
    

    #lie32m=polylog_exp_norm(3/2,zminus)
    #lie52m=polylog_exp_norm(5/2,zminus)
    
    #a=mtfactor*(lie32p+lie32m -3*T/(2*mstar)*( lie52p+lie52m))+ wal.m_sigma^2/(4*wal.g_sigma)*u[1]
    #b= -mtfactor*(lie32p-lie32m) + wal.m_omega^2/(4*wal.g_omega)*u[2]

    a=mfactor*(n_lie32p+n_lie32m -3/(2*mstar)*( n_lie52p+n_lie52m))+ wal.m_sigma^2/(4*wal.g_sigma)*u[1]
    b= -mfactor*(n_lie32p-n_lie32m) + wal.m_omega^2/(4*wal.g_omega)*u[2]
    #println(a)
    #a=pressure_derivativeFermiDirac(T,mustar,mstar,Val(0),Val(0),Val(1))+pressure_derivativeFermiDirac(T,-mustar,mstar,Val(0),Val(0),Val(1)) + wal.m_sigma^2/(4*wal.g_sigma)*u[1]
    
    
    #b=pressure_derivativeFermiDirac(T,mustar,mstar,Val(0),Val(1),Val(0))-pressure_derivativeFermiDirac(T,-mustar,mstar,Val(0),Val(1),Val(0)) + wal.m_omega^2/(4*wal.g_omega)*u[2]

    return SVector{2}(a,b)

end

function system_of_equation3(u::M,T,mu,wal::WaleckaModel1{N}) where {M<:AbstractArray,N}
    #Get the fields and T and mu
    sigma = u[1]
    omega = u[2]
    #T = 1/β
    #mu = α/β
    #Define the shifted mass and chemical potential
    mNStar = wal.m_nucleon - wal.g_sigma * sigma 
    muStar= mu +wal.g_omega*omega
    #check for negative mass from solver
    if mNStar<=0
        sigmaEq = wal.m_sigma^2* sigma /(4*wal.g_sigma)
        omegaEq = wal.m_omega^2 * omega /(4*wal.g_omega)
        return SVector{2}(sigmaEq,omegaEq)
    end
    if T<0.001 #check for T=0 case
        #theta function at low T
        if muStar>mNStar
            sigmaEq = sqrt(2*mNStar)*(-8*mNStar+3*muStar)*(-mNStar+muStar)^(3/2)/(15*pi^2) + wal.m_sigma^2* sigma /(4*wal.g_sigma)
            omegaEq = sqrt(2*mNStar^3)*(-mNStar+muStar)^(3/2)/(3*pi^2) + wal.m_omega^2 * omega /(4*wal.g_omega)
            return SVector{2}(sigmaEq,omegaEq)
        else 
            sigmaEq = wal.m_sigma^2* sigma /(4*wal.g_sigma)
            omegaEq = wal.m_omega^2 * omega /(4*wal.g_omega)
            return SVector{2}(sigmaEq,omegaEq)
        end
    end
    #These are the arguments for the Fermi integrations
    negArg=(-muStar-mNStar)/T
    posArg=(muStar-mNStar)/T

    if T<0.007
        if mNStar<=0 #negative mass from wrong solution
            println("mass negative")
            return SVector{2}(wal.m_sigma^2* sigma /(4*wal.g_sigma),wal.m_omega^2 * omega /(4*wal.g_omega))
        end
        if posArg<0 #theta function at low T
            return SVector{2}(wal.m_sigma^2* sigma /(4*wal.g_sigma),wal.m_omega^2 * omega /(4*wal.g_omega))
        end

        FD12pos = sf_fermi_dirac_half(posArg)
        FD32pos = sf_fermi_dirac_3half(posArg)

        if negArg<-500 #antiparticle contribution can get to small
        FD12neg = zero(negArg)
        FD32neg = zero(negArg)
        else
        FD12neg = sf_fermi_dirac_half(negArg)
        FD32neg = sf_fermi_dirac_3half(negArg)
        end
        sigmaEq = T*sqrt(mNStar*T)/(4*sqrt(2)*pi^(3/2))*(3*T*(FD32pos+FD32neg)-2*mNStar*(FD12pos+FD12neg)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
        omegaEq = (mNStar*T)^(3/2)*(FD12pos-FD12neg)/((2*pi)^(3/2)) + wal.m_omega^2/(4*wal.g_omega)*omega
        #@show negArg, posArg
        return SVector{2}(sigmaEq,omegaEq)
    else
        #These are the arguments for the Fermi integrations
        negArg=(-muStar-mNStar)/T
        posArg=(muStar-mNStar)/T
        if posArg<-500 #antiparticle contribution can get to small
            FD12pos = zero(posArg)
            FD32pos = zero(posArg)
         else
            FD12pos = sf_fermi_dirac_half(posArg)
            FD32pos = sf_fermi_dirac_3half(posArg)
         end
        #FD12pos = sf_fermi_dirac_half(posArg)
        #FD32pos = sf_fermi_dirac_3half(posArg)
        if negArg<-500 #antiparticle contribution can get to small
            FD12neg = zero(negArg)
            FD32neg = zero(negArg)
         else
            FD12neg = sf_fermi_dirac_half(negArg)
            FD32neg = sf_fermi_dirac_3half(negArg)
         end
        #FD12neg = sf_fermi_dirac_half(negArg)
        #FD32neg = sf_fermi_dirac_3half(negArg)

        sigmaEq = T*sqrt(mNStar*T)/(4*sqrt(2)*pi^(3/2))*(3*T*(FD32pos+FD32neg)-2*mNStar*(FD12pos+FD12neg)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
        omegaEq = (mNStar*T)^(3/2)*(FD12pos-FD12neg)/((2*pi)^(3/2)) + wal.m_omega^2/(4*wal.g_omega)*omega
       # @show negArg,posArg
        return SVector{2}(sigmaEq,omegaEq)
    end
end



function  gapequation3(T,mu,wal::WaleckaModel1{N}) where {N}
    _f(u,p)=system_of_equation3(u,T,mu,wal)
    #if mu>0.9
        u0=SVector{2,N}(wal.sigma_0,wal.omega_0)
    #else
      #  u0=SVector{2,N}(0.0*wal.sigma_0,0.0*wal.omega_0)
    #end
    problm=NonlinearProblem{false}(_f, u0) 
    solve(problm, NewtonRaphson(), abstol = 1e-20, reltol = 1e-20)
end

T=0.009
a=gapequation3(1/T,0.91/T,eos)
eos=WaleckaModel1()

gapequation3(0.06,0.9,eos)
gapequation3(0.059,0.9,eos)

T=0.1
murange=collect(0.:0.001:2.)
u2=gapequation3.(Ref(T),murange,Ref(eos))

sig=zeros(length(u2))
omg=zeros(length(u2))
sig2=zeros(length(u2))
omg2=zeros(length(u2))

for ind in range(1,length(murange))
    sig[ind]=u2[ind][1]
    omg[ind]=u2[ind][2]
    sig2[ind]=u2[ind][1]*u2[ind][1]
    omg2[ind]=u2[ind][2]*u2[ind][2]
end


pp=plot(murange,sig,label="sigma, T=0.039")
plot!(murange,omg,label="omega, T=0.039")
plot!(murange,sig,label="sigma, T=0.04")
plot!(murange,omg,label="omega, T=0.04")


function get_transition_line(step)
    Tmin=0.0
    Tmax=0.03
    Mumin=0.8
    Mumax=1.1
    Trange=collect(Tmin:step:Tmax)
    murange=collect(Mumin:step:Mumax)
    Tlen=length(Trange)
    mulen=length(murange)
    temp=zeros(Tlen,mulen)
    for t in 1:1:Tlen
        for mu in 1:1:mulen
            tempsig,tempomg=gapequation3(Trange[t],murange[mu],eos)
            temp[t,mu]=tempsig
        end
    end
    finalTList=[]
    finalMuList=[]
       #return temp[1,:]
       for t in 1:1:Tlen
        tempList=temp[t,:]
        indList=[]
            for i=1:length(tempList)
                if tempList[i] !=0#>0.01#!= 0
                    push!(indList,i)
                end
            end
        indMu=first(indList)
        push!(finalTList,Trange[t])
        push!(finalMuList,murange[indMu])
       end

    return finalTList, finalMuList

end

Tlist, mulist =get_transition_line(0.0005)

plot(Tlist,mulist,ylabel="μ [GeV]",xlabel="T [GeV]")
transplot=plot(mulist,Tlist,ylabel="T [GeV]",xlabel="μ [GeV]")
#savefig(transplot,"transition_line_jumps.png")
#mulist
mumodel(T,p)=p[1] .- p[2]*T.^2 .-p[3].*T.^4 
p0=[1.,1.,1.]
fit= curve_fit(mumodel,Tlist,mulist,p0)

params = fit.param
using LsqFit


model(t, p) = p[1] * exp.(-p[2] * t)
tdata = Base._linspace(0.,10.,20)
ydata = model(tdata, [1.0 2.0]) + 0.01*randn(length(tdata))
p0 = [0.5, 0.5]
fit = curve_fit(model, tdata, ydata, p0)
param = fit.param


function system_of_equation43(u::M,T,mu,wal::WaleckaModel1{N}) where {M<:AbstractArray,N}
   
    omega = u[2]
    sigma = u[1]
 
     mNStar = wal.m_nucleon - wal.g_sigma * sigma 
     #αStar = α + β*wal.g_omega * omega 
     muStar= mu +wal.g_omega*omega
 
     #negArg = -αStar-mNStar*β
     #posArg = αStar-mNStar*β
     negArg=(-muStar-mNStar)/T
     posArg=(muStar-mNStar)/T
 
 #if T>0.03
     #put all
 
     if mNStar<=0
         return SVector{2}(wal.m_sigma^2* sigma /(4*wal.g_sigma),wal.m_omega^2 * omega /(4*wal.g_omega))
     end
     if posArg<0
         return SVector{2}(wal.m_sigma^2* sigma /(4*wal.g_sigma),wal.m_omega^2 * omega /(4*wal.g_omega))
     end
 
     FD12pos = sf_fermi_dirac_half(posArg)
     FD32pos = sf_fermi_dirac_3half(posArg)
 
     if negArg<-500    
     FD12neg = zero(negArg)#sf_fermi_dirac_half(negArg)
     FD32neg = zero(negArg)#sf_fermi_dirac_3half(negArg)
     else
     FD12neg = sf_fermi_dirac_half(negArg)
     FD32neg = sf_fermi_dirac_3half(negArg)
     end
 
    # @show sigma
     #sigmaEq = 1/β^3 * (mNStar*β/(2*pi))^(3/2)*(-FD12pos-FD12neg+3/(2*mNStar*β)*(FD32pos+FD32neg)) + wal.m_sigma^2* sigma /(4*wal.g_sigma)
     #omegaEq = 1/β^4 * (mNStar*β/(2*pi))^(3/2)*(FD12pos-FD12neg) + wal.m_omega^2 * omega /(4*wal.g_omega*β)
     #these are the direct replacements
     #sigmaEq = 1/(4*β^2*sqrt(2*pi^3))*sqrt(mNStar/β)*(3*FD32neg+3FD32pos-2*mNStar*β*(FD12neg+FD12pos)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
     #omegaEq = (mNStar/(β*2*pi))^(3/2)*(FD12pos-FD12neg)+ wal.m_omega^2/(4*wal.g_omega)*omega
     #these are with T and mu
     sigmaEq = T*sqrt(mNStar*T)/(4*sqrt(2)*pi^(3/2))*(3*T*(FD32pos+FD32neg)-2*mNStar*(FD12pos+FD12neg)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
     omegaEq = (mNStar*T)^(3/2)*(FD12pos-FD12neg)/((2*pi)^(3/2)) + wal.m_omega^2/(4*wal.g_omega)*omega
 
     return SVector{2}(sigmaEq,omegaEq)
 
 end
 

function system_of_equation_zero_T(u,μ)

    
    m_nucleon=0.939
    m_sigma=0.5
    m_omega=0.783
    g_sigma=9.69719
    g_omega=13.2025
    sigma_0=0.0444258
    omega_0= -0.026538

    omega = u[2]
    sigma = u[1]
 
     mNStar = m_nucleon - g_sigma * sigma 
     μStar = μ + g_omega * omega 
 #@show μ
    #@show -mNStar*(mNStar-μStar)

    if mNStar<=0
        return SVector{2}(0.0,0.0)
    end
 if μStar>mNStar
    # @show sigma
    sigmaEq = sqrt(2*mNStar)*(-8*mNStar+3*μStar)*(-mNStar+μStar)^(3/2)/(15*pi^2) + m_sigma^2* sigma /(4*g_sigma)
     omegaEq = sqrt(2*mNStar^3)*(-mNStar+μStar)^(3/2)/(3*pi^2) + m_omega^2 * omega /(4*g_omega)
     #these are the direct replacements
     #sigmaEq = 1/(4*β^2*sqrt(2*pi^3))*sqrt(mNStar/β)*(3*FD32neg+3FD32pos-2*mNStar*β*(FD12neg+FD12pos)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
     #omegaEq = (mNStar/(β*2*pi))^(3/2)*(FD12pos-FD12neg)+ wal.m_omega^2/(4*wal.g_omega)*omega
     #these are with T and mu
     #sigmaEq = T*sqrt(mNStar*T)/(4*sqrt(2)*pi^(3/2))*(3*T*(FD32pos+FD32neg)-2*mNStar*(FD12pos+FD12neg)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
     #omegaEq = (mNStar*T)^(3/2)*(FD12pos-FD12neg)/((2*pi)^(3/2)) + wal.m_omega^2/(4*wal.g_omega)*omega
 
     return SVector{2}(sigmaEq,omegaEq)
 else 
   sigmaEq =    + m_sigma^2* sigma /(4*g_sigma)
   omegaEq =  + m_omega^2 * omega /(4*g_omega)
   return SVector{2}(sigmaEq,omegaEq)
 end
end



 function system_of_equation_small_T(u,T,μ)
   
    m_nucleon=0.939
    m_sigma=0.5
    m_omega=0.783
    g_sigma=9.69719
    g_omega=13.2025
    sigma_0=0.0444258
    omega_0= -0.026538

    omega = u[2]
    sigma = u[1]
 
     mNStar = m_nucleon - g_sigma * sigma 
     μStar = μ + g_omega * omega 
if mNStar<=0

    sigmaEq = m_sigma^2* sigma /(4*g_sigma)
    omegaEq = m_omega^2 * omega /(4*g_omega)
    return SVector{2}(sigmaEq,omegaEq)
    #is wrong I think; m only makes the first terms zero
#    return SVector{2}(0.,0.)
end
#@show mNStar,μ
 if μStar>mNStar
    sigmaEq = sqrt(2*mNStar)*(-8*mNStar+3*μStar)*(-mNStar+μStar)^(3/2)/(15*pi^2) + (sqrt(2)*T^2*(3*μStar-4*mNStar))/(12*sqrt(2)*sqrt(μStar-mNStar)) + m_sigma^2* sigma /(4*g_sigma)
    omegaEq = sqrt(2*mNStar^3)*(-mNStar+μStar)^(3/2)/(3*pi^2) + mNStar^(3/2)*T^2/(12*sqrt(2)*sqrt(μStar-mNStar)) + m_omega^2 * omega /(4*g_omega)
    #these are the direct replacements
    #sigmaEq = 1/(4*β^2*sqrt(2*pi^3))*sqrt(mNStar/β)*(3*FD32neg+3FD32pos-2*mNStar*β*(FD12neg+FD12pos)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
    #omegaEq = (mNStar/(β*2*pi))^(3/2)*(FD12pos-FD12neg)+ wal.m_omega^2/(4*wal.g_omega)*omega
    #these are with T and mu
    #sigmaEq = T*sqrt(mNStar*T)/(4*sqrt(2)*pi^(3/2))*(3*T*(FD32pos+FD32neg)-2*mNStar*(FD12pos+FD12neg)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
    #omegaEq = (mNStar*T)^(3/2)*(FD12pos-FD12neg)/((2*pi)^(3/2)) + wal.m_omega^2/(4*wal.g_omega)*omega
    return SVector{2}(sigmaEq,omegaEq)
 else 
   sigmaEq = m_sigma^2* sigma /(4*g_sigma)
   omegaEq = m_omega^2 * omega /(4*g_omega)
   return SVector{2}(sigmaEq,omegaEq)
 end

end

system_of_equation_small_T((0.03,-0.01),0.1,1)

function system_of_equation0(u,T,mu)
    if T<1e-5
        system_of_equation_zero_T(u,mu)
    else
        system_of_equation_small_T(u,T,mu)
    end
end

system_of_equation0((0.03,-0.02),2*1e-7,1)

using StaticArrays
system_of_equation0((0,0),4)

function gapequation0(T,mu)
    _f(u,p)=system_of_equation0(u,T,mu)
   # if mu>0.9
    u0=SVector{2}(0.0444258,-0.026538)
    #else
#u0=SVector{2}(0,0)
   # end
        problm=NonlinearProblem{false}(_f, u0)
    solve(problm, NewtonRaphson(), reltol = 1e-10) 
end

function gapequation00(mu)
    _f(u,p)=system_of_equation_zero_T(u,mu)
   # if mu>0.9
    u0=SVector{2}(0.0444258,-0.026538)
    #else
#u0=SVector{2}(0,0)
 #   end
        problm=NonlinearProblem{false}(_f, u0)
    solve(problm, NewtonRaphson()) 
end

function sigma_coefficients(index)
    internal_list=(0.009552652804179272,0.05892556509887895,-0.01696251715190005,-0.06619659887773394,-1.0538966296794874,-37.782944508552646,-2410.394222997936,-240480.4700457461,
    -3.457064380249338e7,-6.767280831712292e9,-1.7307370204456975e12)
    return internal_list[index]    
end

function omega_coefficients(index)
    internal_list=(0.047763264020896354,0.05892556509887895,0.05088755145570016,0.46337619214413756,11.59286292647436,566.7441676282896,45797.49023696079,5.53105081105216e6,9.334073826673213e8,
    2.0978570578308105e11,6.057579571559942e13) 
    return internal_list[index]    
end

function do_omega_equation(T,μ,u,nmax)
    m_nucleon=0.939
    m_sigma=0.5
    m_omega=0.783
    g_sigma=9.69719
    g_omega=13.2025
    sigma_0=0.0444258
    omega_0= -0.026538

    omega = u[2]
    sigma = u[1]
 
     mNStar = m_nucleon - g_sigma * sigma 
     μStar = μ + g_omega * omega 

     temp=0
    if mNStar<0
        return m_omega^2 * omega /(4*g_omega)
    end

    if μStar<mNStar
        return m_omega^2 * omega /(4*g_omega)
    end
     for n in 0:1:nmax
        temp += omega_coefficients(n+1)*mNStar^(3/2)*T^(2*n)*(μStar-mNStar)^(3/2-2n)
     end
return temp + m_omega^2 * omega /(4*g_omega)
end

function do_sigma_equation(T,μ,u,nmax)
    m_nucleon=0.939
    m_sigma=0.5
    m_omega=0.783
    g_sigma=9.69719
    g_omega=13.2025
    sigma_0=0.0444258
    omega_0= -0.026538

    omega = u[2]
    sigma = u[1]
 
     mNStar = m_nucleon - g_sigma * sigma 
     μStar = μ + g_omega * omega 

     temp=0
    if mNStar<0
        return m_sigma^2 * sigma /(4*g_sigma)
    end

    if μStar<mNStar
        return m_sigma^2 * sigma /(4*g_sigma)
    end
     for n in 0:1:nmax
        temp += sigma_coefficients(n+1)*mNStar^(1/2)*T^(2*n)*(μStar-mNStar)^(3/2-2n)*(3μStar+4*mNStar*(n-2))
     end
return temp + m_sigma^2 * sigma /(4*g_sigma)

end




function system_of_equation_expansion(u,T,mu,nmax)
    sigmaEq = do_sigma_equation(T,mu,u,nmax)
    omegaEq = do_omega_equation(T,mu,u,nmax)
    return SVector{2}(sigmaEq,omegaEq)
end

system_of_equation_expansion((0.02,-0.04),0.002,0.99,3)

function gapequation_expansion(T,mu,nmax)
    _f(u,p)=system_of_equation_expansion(u,T,mu,nmax)
    #if mu>0.95
        u0=SVector{2}(0.0444258,-0.026538)
    #else
       # u0=SVector{2}(0.0,0.0)
    #end
    problm=NonlinearProblem{false}(_f, u0)
    solve(problm, NewtonRaphson(), reltol = 1e-10, abstol=1e-15) 
end


function gapequation_expansion_dynamic(T,mu,nmax)

    m_nucleon=0.939
    m_sigma=0.5
    m_omega=0.783
    g_sigma=9.69719
    g_omega=13.2025
    sigma_0=0.0444258
    omega_0= -0.026538

    _f0(u,p)=system_of_equation_expansion(u,T,mu,nmax)
    u0=SVector{2}(0.0444258*0,-0.026538*0)
    problm0=NonlinearProblem{false}(_f0, u0)
    sigmaSeed, omegaSeed =solve(problm0, NewtonRaphson(), reltol = 1e-10) 
    u1=SVector{2}(0.0444258,-0.026538)
    _f(u,p)=system_of_equation_expansion(u,T,mu,nmax)
    problm=NonlinearProblem{false}(_f, u1)
    sigma, omega= solve(problm, NewtonRaphson(), abstol=1e-15 , reltol = 1e-15)
    pot0=-0.5*m_sigma^2*sigmaSeed^2+0.5*m_omega^2*omegaSeed^2
    pot=-0.5*m_sigma^2*sigma^2+0.5*m_omega^2*omega^2
    if pot0<pot
        return (sigmaSeed,omegaSeed)
    else
        return (sigma,omega)
    end
end


using SimpleNonlinearSolve

function gapequation_expansion_braket(T,mu,nmax)

    m_nucleon=0.939
    m_sigma=0.5
    m_omega=0.783
    g_sigma=9.69719
    g_omega=13.2025
    sigma_0=0.0444258
    omega_0= -0.026538

    _f0(u,p)=system_of_equation_expansion(u,T,mu,nmax)
    u0=SVector{2}(0.0444258,-0.026538)
    uspan=(SVector{2}(0.,0.),SVector{2}(sigma_0,omega_0))
    problm0=IntervalNonlinearProblem(_f0,uspan)
    #problm0=NonlinearProblem{false}(_f0, u0)
    #sigmaSeed, omegaSeed =solve(problm0, NewtonRaphson(), reltol = 1e-10) 
    solve(problm0, ITP())
end

T=0.001
alpharange=murange/T
T=0.02
murange=collect(0.9:0.001:1)
u2=gapequation3.(Ref(T),murange,Ref(eos))
#u2=gapequation_expansion_braket.(Ref(0.007),murange,Ref(5))
#u2=gapequation_expansion_dynamic.(Ref(0.005),murange,Ref(10))

#gapequation3(1/T,0.9/T,eos)

#gapequation_expansion_dynamic(0.01,1,1)

sig=zeros(length(u2))
omg=zeros(length(u2))
sig2=zeros(length(u2))
omg2=zeros(length(u2))

for ind in range(1,length(murange))
    sig[ind]=u2[ind][1]
    omg[ind]=u2[ind][2]
    sig2[ind]=u2[ind][1]*u2[ind][1]
    omg2[ind]=u2[ind][2]*u2[ind][2]
end


plot(murange,sig,label="sigma, T=0")
plot!(murange,omg,label="omega, T=0")
plot!(murange,sig,label="sigma, T=0.015")
plot!(murange,omg,label="omega, T=0.015")

sigma,omega =gapequation3.(Ref(T),murange,Ref(eos))

a=collect(0.01:0.01:0.3)
a[3]

function get_transition_line(step)
    Tmin=0.001
    Tmax=0.030
    Mumin=0.8
    Mumax=1.1
    Trange=collect(Tmin:step:Tmax)
    murange=collect(Mumin:step:Mumax)
    Tlen=length(Trange)
    mulen=length(murange)
    temp=zeros(Tlen,mulen)
    for t in 1:1:Tlen
        for mu in 1:1:mulen
            tempsig,tempomg=gapequation3(Trange[t],murange[mu],eos)
            temp[t,mu]=tempsig
        end
    end
    finalTList=[]
    finalMuList=[]
       #return temp[1,:]
       for t in 1:1:Tlen
        tempList=temp[t,:]
        indList=[]
        for i=1:length(tempList)
            if tempList[i] != 0
                push!(indList,i)
            end
        end
        indMu=first(indList)
        push!(finalTList,Trange[t])
        push!(finalMuList,murange[indMu])
       end

    return finalTList, finalMuList

end

Tlist, mulist =get_transition_line(0.0005)

plot(Tlist,mulist,ylabel="μ [GeV]",xlabel="T [GeV]")
plot(mulist,Tlist,ylabel="T [GeV]",xlabel="μ [GeV]")

A = [1,2,0,0,4,0]
B = []


for i=1:length(test)
    if test[i] != 0
        push!(B,i)
    end
end
B
Tlist, mulist =get_transition_line()

plot(Tlist,mulist)

plot(test)
test[1]
test[18]
sigma,omega =gapequation3(T,0.96,eos)
sigma
plot(sigma)
sigma
rand()

using Plots

plot(do_sigma_equation.(Ref(0.),Ref(1.9),Ref((0.001,-0.001)),collect(0:1:10)))
plot!(do_sigma_equation.(Ref(0.001),Ref(1.9),Ref((0.001,-0.001)),collect(0:1:10)))
plot!(do_sigma_equation.(Ref(0.002),Ref(1.9),Ref((0.001,-0.001)),collect(0:1:10)))
plot!(do_sigma_equation.(Ref(0.12),Ref(1.9),Ref((0.001,-0.001)),collect(0:1:10)))
plot(collect(0.93:0.001:1),do_omega_equation.(Ref(0.),collect(0.93:0.001:1),Ref((0.001,-0.001)),Ref(5)))

sigma_coefficients(1)
omega_coefficients(11)
#gapequation0(0,1)

#gapequation0.(collect(0.8:0.01:1))
murange=collect(0.85:0.01:1.0)



u2=gapequation0.(Ref(0.01),murange)
sig=zeros(length(u2))
omg=zeros(length(u2))
sig2=zeros(length(u2))
omg2=zeros(length(u2))


for ind in range(1,length(murange))
    sig[ind]=u2[ind][1]
    omg[ind]=u2[ind][2]
    sig2[ind]=u2[ind][1]*u2[ind][1]
    omg2[ind]=u2[ind][2]*u2[ind][2]
end
#using Plots
zerop=plot(murange,sig,label="sigma")
plot!(murange,omg,label="omega")
plot!(murange,sig,label="sigma2")
plot!(murange,omg,label="omega2")

#savefig(zerop,"zero_T_fields_bigR.png")

function gapequation(T,mu,wal::WaleckaModel1{N}) where {N}
    _f(u,p)=system_of_equation(u,T,mu,wal)
    problm=NonlinearProblem{false}(_f, SVector{2,N}(wal.sigma_0,wal.omega_0))
    #problm=NonlinearProblem{false}(_f, SVector{2,N}(0.0,0))
    #solve(problm, NewtonRaphson(autodiff = Val{false}()), tol = 1e-6) 
    #solve(problm, NewtonRaphson(), tol = 1e-6) 
    solve(problm, NewtonRaphson(), abstol = 1e-20) 
    #solve(problm,Broyden(), tol = 1e-6)
    #solve(problm,Bisection(), tol = 1e-6)
end

#solve(problm,NewtonRaphson(), tol = 1e-6)  
#struct Thermodyanmic{T,N} 
#    pressure::T
#    pressure_derivative::SVector{N,T}
#    pressure_hessian::SHermitianCompact{N,T}
#end

# Import the package and define the problem to optimize
using Optimization
rosenbrock(u, p) = (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
u0 = zeros(2)
p = [1.0, 100.0]

prob = OptimizationProblem(rosenbrock, u0, p)

using ForwardDiff
# Import a solver package and solve the optimization problem
using OptimizationOptimJL
sol = solve(prob, NelderMead())

function pressure(u,p)

    m_nucleon=0.939
    m_sigma=0.5
    m_omega=0.783
    g_sigma=9.69719
    g_omega=13.2025
    sigma_0=0.0444258
    omega_0= -0.026538

    omega = u[2]
    sigma = u[1]

    β=p[1]
    α=p[2]

    mNStar = m_nucleon - g_sigma * sigma 
    αStar = α + β*g_omega * omega 

    negArg = -αStar-mNStar*β
    posArg = αStar-mNStar*β

    FD12pos = sf_fermi_dirac_half(posArg)
    FD12neg = sf_fermi_dirac_half(negArg)
    FD32pos = sf_fermi_dirac_3half(posArg)
    FD32neg = sf_fermi_dirac_3half(negArg)
    return -(4/β^4*(mNStar*β/(2*pi))^(3/2)*(FD32neg+FD32pos) -1/2*m_sigma^2*sigma^2 + 1/2*m_omega^2*omega^2 - g_sigma * sigma^4 - g_omega*omega^4 )
end

u=zeros(2)
pressure((0,0),(1,1))

#problem=OptimizationProblem(pressure,u,(1,1),lb=[-Inf,-Inf],ub=[0.939/9.69719,Inf])
problem=OptimizationProblem(pressure,u,(1,1))

cons(res, x, p) = (res .= [x[1]])

using Optimization, OptimizationMOI, OptimizationOptimJL#, Ipopt
using ForwardDiff, ModelingToolkit

optprob = OptimizationFunction(pressure, Optimization.AutoForwardDiff(), cons = cons)
prob = OptimizationProblem(optprob, u, (1,1), lcons = [-Inf], ucons = [0.939/9.69719])
sol = solve(prob, BFGS())

optf = OptimizationFunction(pressure)
prob = OptimizationProblem(optf, u, (1,1))
sol = solve(prob, NelderMead())
using Optimization, OptimizationBBO

u=ones(2)
f = OptimizationFunction(pressure)
prob = Optimization.OptimizationProblem(f, u, (1,1), lb = [-0.1, -0.1], ub = [0.939/9.69719, 0.1])
#sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 10,maxtime = 10.0)
sol = solve(prob,BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 1000, maxtime = 100.0)
0.939/9.69719
1
#sol = solve(problem, NelderMead())
using Plots
sigmarange=collect(-1:0.01:0.939/9.69719)
urange=(sigmarange,ones(length(sigmarange)))
omegarange=collect(-1:0.001:1)
plist=zeros(length(sigmarange))
plist2=zeros(length(omegarange))
for i in range(1,length(sigmarange))
   plist[i]= pressure((sigmarange[i],1),(10,10))

end

for i in range(1,length(omegarange))
    plist2[i]= pressure((0.03,omegarange[i]),(1,1))
 
 end

plot(sigmarange,plist)
plot(omegarange,plist2)
urange[1]
#Thermodyanmic(a::T,b::NTuple{N,T},c::NTuple{M,T}) where {T,N,M}= Thermodyanmic{T,N}(a,SVector{N,T}(b),SHermitianCompact{N,T}(SVector{M}(c))) 

#Thermodyanmic(a::T,b::SVector{N,T},c::SVector{M,T}) where {T,N,M}= Thermodyanmic{T,N}(a,b,SHermitianCompact{N,T}(c)) 


plot(collect(-10:0.1:10),sf_fermi_dirac_3half.(collect(-10:0.1:10)))

function zeroTpressure(u,p)

    mu=p[1]
    m_nucleon=0.939
    m_sigma=0.5
    m_omega=0.783
    g_sigma=9.69719
    g_omega=13.2025
    sigma_0=0.0444258
    omega_0= -0.026538

    omega = u[2]
    sigma = u[1]

    mNStar = m_nucleon - g_sigma * sigma 
    muStar = mu + g_omega * omega 
if -mNStar>muStar
    return -(4*(2*mNStar)^(5/2)/(60*pi^2*mNStar)*((mNStar-muStar)^(5/2)) -1/2*m_sigma^2*sigma^2 + 1/2*m_omega^2*omega^2)
elseif muStar>mNStar
    return - ( 4*(2*mNStar)^(5/2)/(60*pi^2*mNStar)*((muStar-mNStar)^(5/2)) -1/2*m_sigma^2*sigma^2 + 1/2*m_omega^2*omega^2)
else
    return 0
end
end

function getp(mu)
u0=zeros(2)
prob = OptimizationProblem(zeroTpressure, u0, (mu))
sol = solve(prob, NelderMead())
zeroTpressure(sol,(mu))
end

getp(1)

plot(collect(0.8:0.01:1.1),zeroTpressure.(collect(0.8:0.01:1.1),Ref((0.05,-0.07))))

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))

heaviside(-1.0)

function thermodynamic(T::A,mu::B,wal::WaleckaModel1{L}) where {A,B,L}

    u=gapequation(T,mu,wal::WaleckaModel1{L})

    mustar=mu +wal.g_omega*u[2]
    mstar=wal.m_nucleon-wal.g_sigma*u[1]
    mSigma2=wal.m_sigma^2
    gSigma2=wal.g_sigma^2
    mOmega2=wal.m_omega^2
    gOmega2=wal.g_omega^2
    nump=mustar-mstar
    numm=-mustar-mstar

    #zplus=-exp((mustar-mstar)/T)
    #zminus=-exp((-mustar-mstar)/T)
    
    #mtfactor=(mstar*T/(2*pi))^(3/2)

    mfactor=(mstar/(2*pi))^(3/2)
    
    n_lie32p=-BigFloat(sf_fermi_dirac_half(nump/T))#polylog_exp_norm(3/2,nump,T)
    n_lie52p=-BigFloat(sf_fermi_dirac_3half(nump/T))#polylog_exp_norm(5/2,nump,T)
    n_lie12p=-BigFloat(sf_fermi_dirac_mhalf(nump/T))#polylog_exp_norm(1/2,nump,T)
    

    n_lie32m=-BigFloat(sf_fermi_dirac_half(numm/T))#polylog_exp_norm(3/2,numm,T)
    n_lie52m=-BigFloat(sf_fermi_dirac_half(numm/T))#polylog_exp_norm(5/2,numm,T)
    n_lie12m=-BigFloat(sf_fermi_dirac_half(numm/T))#polylog_exp_norm(1/2,numm,T)

    # derivative of the fermi dirac distro
    p100p= mfactor*(2*(-mstar+mustar)*n_lie32p-5*n_lie52p)/(2*T)
    p100m= mfactor*(2*(-mstar-mustar)*n_lie32p-5*n_lie52m)/(2*T)
    p010p= -mfactor*n_lie32p
    p010m= -mfactor*n_lie32m

    p200p=-mfactor/(4*T^2)*(4*(mstar-mustar)*((mstar-mustar)*n_lie12p+3*n_lie32p) +15*n_lie52p)
    p110p=-mfactor/(2*T)*(2*(mstar-mustar)*n_lie12p+3*n_lie32p)
    p101p=mfactor/(4*mstar*T)*(4*mstar*(mstar-mustar)*n_lie12p +(6*mustar*n_lie32p- 15*n_lie52p ) )
    p011p=mfactor/(2*mstar)*(2*mstar*n_lie12p -3*n_lie32p)
    p020p= -mfactor*n_lie12p
    p002p=mfactor/(4*mstar)*( -4*mstar*n_lie12p+ 3*(   4*mstar*n_lie32p-n_lie52p))
    
    p200m=-mfactor/(4*T^2)*(4*(mstar+mustar)*((mstar+mustar)*n_lie12m+3*n_lie32m) +15*n_lie52m)
    p110m=-mfactor/(2*T)*(2*(mstar+mustar)*n_lie12m+3*n_lie32m)
    p101m=mfactor/(4*mstar*T)*(4*mstar*(mstar+mustar)*n_lie12m +(-6*mustar*n_lie32p- 15*n_lie52m ) )
    p011m=mfactor/(2*mstar)*(2*mstar*n_lie12m -3*n_lie32m)
    p020m= -mfactor*n_lie12m
    p002m=mfactor/(4*mstar)*( -4*mstar*n_lie12m+ 3*(   4*mstar*n_lie32m-n_lie52m))
    




    pressure=-4*mfactor*(n_lie52m+n_lie52p)  
    -0.5*mSigma2*u[1]^2 +0.5*mOmega2*u[2]^2
    #first derivative
    p10= 4 *(p100m + p100p)
    p01=-4 *p010m + 4*p010p 


    #second derivative
    p20=4*((4*(gSigma2*mOmega2*^(p101m + p101p,2) + 
    4*gOmega2*gSigma2*(p020m*^(p101m + p101p,2) + p020p*^(p101m + 
    p101p,2) + (-2*(p011m - p011p)*(p101m + p101p) + (p002m + 
    p002p)*(p110m - p110p))*(p110m - p110p)) - gOmega2*mSigma2*^(p110m - 
    p110p,2)))/(mOmega2*(mSigma2 - 4*gSigma2*(p002m + p002p)) + 
    4*gOmega2*(mSigma2*(p020m + p020p) + 4*gSigma2*(^(p011m - p011p,2) - 
    (p002m + p002p)*(p020m + p020p)))) + p200m + p200p)

    p11=(4*mOmega2*(4*gSigma2*(-(p011m*(p101m + p101p)) + p011p*(p101m + 
    p101p) + (p002m + p002p)*(p110m - p110p)) + mSigma2*(-p110m + 
    p110p)))/(mOmega2*(mSigma2 - 4*gSigma2*(p002m + p002p)) + 
    4*gOmega2*(mSigma2*(p020m + p020p) + 4*gSigma2*(^(p011m - p011p,2) - 
    (p002m + p002p)*(p020m + p020p))))

    p02=(4*mOmega2*(mSigma2*(p020m + p020p) + 4*gSigma2*(^(p011m - p011p,2) 
    - (p002m + p002p)*(p020m + p020p))))/(mOmega2*(mSigma2 - 
    4*gSigma2*(p002m + p002p)) + 4*gOmega2*(mSigma2*(p020m + p020p) + 
    4*gSigma2*(^(p011m - p011p,2) - (p002m + p002p)*(p020m + p020p))))


    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
    

end 

function thermodynamic2(T::A,mu::B,wal::WaleckaModel1{L}) where {A,B,L}

    u=gapequation2(T,mu,wal::WaleckaModel1{L})

    omega = u[2]
    sigma = u[1]
 
     mNStar = wal.m_nucleon - wal.g_sigma * sigma 
     muStar = μ + wal.g_omega * omega 
 
     negArg = (-muStar-mNStar)/T
     posArg = (muStar-mNStar)/T
 
     FD12pos = sf_fermi_dirac_half(posArg)
     FD12neg = sf_fermi_dirac_half(negArg)
     FD32pos = sf_fermi_dirac_3half(posArg)
     FD32neg = sf_fermi_dirac_3half(negArg)

    pressure = sqrt(2)*T*(mNStar*T)^(3/2)*(FD32pos+FD32neg)/(pi^(3/2))-1/2*wal.m_sigma^2*sigma^2+1/2*wal.m_omega^2*omega^2

end
using Plots

eos=WaleckaModel1()

T=0.2
mu=1

gapequation3(1/T,mu/T,eos)

T=0.01
1/T
murange=collect(0.8:0.005:0.99)
alpharange=murange/T
u2=gapequation3.(Ref(1/T),alpharange,Ref(eos))
sig=zeros(length(u2))
omg=zeros(length(u2))
sig2=zeros(length(u2))
omg2=zeros(length(u2))

#gapequation3(1/T,0.91/T,eos)

#0.939/9.69

for ind in range(1,length(murange))
    sig[ind]=u2[ind][1]
    omg[ind]=u2[ind][2]
    sig2[ind]=u2[ind][1]*u2[ind][1]
    omg2[ind]=u2[ind][2]*u2[ind][2]
end



plot(murange,sig)
plot!(murange,omg)

using LinearAlgebra
mul!
0.01
 #pressure(T,mu,wal::WaleckaModel1{L}) where {L}= thermodynamics(T,mu,wal).pressure
#
 #pressure_derivative(T,mu,::Val{1},::Val{0},wal::WaleckaModel1{L}) where {L}= thermodynamics(T,mu,wal).pressure_derivative[1]
 #pressure_derivative(T,mu,::Val{0},::Val{1},wal::WaleckaModel1{L}) where {L}= thermodynamics(T,mu,wal).pressure_derivative[2]
 #pressure_derivative(T,mu,::Val{2},::Val{0},wal::WaleckaModel1{L}) where {L}= thermodynamics(T,mu,wal).pressure_hessian[1]
 #pressure_derivative(T,mu,::Val{1},::Val{1},wal::WaleckaModel1{L}) where {L}= thermodynamics(T,mu,wal).pressure_hessian[2]
 #pressure_derivative(T,mu,::Val{0},::Val{2},wal::WaleckaModel1{L}) where {L}= thermodynamics(T,mu,wal).pressure_hessian[3]
    
#=
using Plots

eos=WaleckaModel1()

murange=collect(0.855:0.009:0.95)
u=gapequation.(Ref(0.005),murange,Ref(eos))

0.009
murange=collect(0.855:0.005:0.95)
u2=gapequation2.(Ref(0.007),murange,Ref(eos))

using SteadyStateDiffEq
import Pkg; Pkg.add("SteadyStateDiffEq")


gapequation2(0.003,0.9,eos)
gapequation2(0.003,0.95,eos)
gapequation(0.005,0.9,eos)
gapequation(0.005,0.95,eos)


sig=zeros(length(u2))
omg=zeros(length(u2))
sig2=zeros(length(u2))
omg2=zeros(length(u2))


a=zeros(length(u))
b=zeros(length(u))
a2=zeros(length(u))
b2=zeros(length(u))

for ind in range(1,21)
    a[ind]=u[ind][1]
    b[ind]=u[ind][2]
    a2[ind]=u[ind][1]*u[ind][1]
    b2[ind]=u[ind][2]*u[ind][2]
    println(u[ind][1])
end

murange=collect(0.851:0.01:0.95)
u2=gapequation2.(Ref(0.015),murange,Ref(eos))
sig=zeros(length(u2))
omg=zeros(length(u2))
sig2=zeros(length(u2))
omg2=zeros(length(u2))
for ind in range(1,length(murange))
    sig[ind]=u2[ind][1]
    omg[ind]=u2[ind][2]
    sig2[ind]=u2[ind][1]*u2[ind][1]
    omg2[ind]=u2[ind][2]*u2[ind][2]
end

gapequation2(0.015,0.92,eos)

pp=plot(murange,sig,label="sigma",xlabel="μ [GeV]",ylabel="fields",title="Phase transition, T=7 MeV")
plot!(murange,omg,label="omega")

savefig(pp,"t07_fields.png")
show(pp)

plot(murange,sig./omg,label="sigma")


plot(-eos.m_omega^2*b2+eos.m_sigma^2*a2)
potplot=plot(murange,-0.5*(-eos.m_omega^2*omg2+eos.m_sigma^2*sig2),xlabel="μ [GeV]",ylabel="potential",title="-0.5*mσ^2*σ+0.5*mω^2*ω, T=7 MeV")
savefig(potplot,"potential_7MeV.png")


murange2=collect(0.8:0.01:0.95)
nlist=pressure_derivative.(Ref(0.007),murange2,Val(0),Val(1),Ref(eos))

plot(murange2,nlist)

plot(murange,a)
plot(murange,b)

plot(a)
=#
#=
#####
TODO

-solver with domain condition to make sure mN does not get negative 
-chenge implementation of thermodynamic function to also have gsl
-code breaks in two situations: 1. underflow from gsl, negative number arg/T gives e^-100
2. domain error for mass, when sigma gets to large and makes mN negative-> need solver with constraints

#####
=#
#= using Plots
using LaTeXStrings
using Measures
import Pkg; Pkg.add("Measures")
pressure(0.01,5.0,eos)

x=collect(0.01:0.01:1.5)


[ i for i in 0:0.001:0.003]
[ L"%$i"  for i in 0:0.001:0.003]
p1=scatter(x,pressure.(Ref(0.01),x,Ref(eos)),
    
    xlabel=L"\mu\quad \mathrm{GeV}",ylabel=L"p(\mu,T)\quad \mathrm{GeV}^4",
    yticks=([ 0, 0.0001, 0.0002,0.0003],[  L"0", L"0.0001", L"0.0002",L"0.0003"]),
    xticks=([ 0, 0.5, 1,1.5],[  L"0", L"0.5", L"1",L"1.5"]),
    label=L"T=0.01 \mathrm{GeV}",
    bottom_margin=8mm,
    left_margin=1mm,
    top_margin=5mm,
    xtickfont=font(10),
    ytickfont=font(10),
    legendfont=font(12),
    )


    savefig(p,"pressure.pdf")



p=scatter(x,pressure_derivative.(Ref(0.01),x,Val(0),Val(1),Ref(eos)),
    
    xlabel=L"\mu\quad \mathrm{GeV}",
    ylabel=L"$\frac{\partial p(\mu,T)}{\partial \mu}\quad \mathrm{GeV}^3$",
    yticks=([0,0.0005,0.0010,0.0015,0.002,0.0025,0.003],[L"0",L"0.0005 ",L"0.0010 ",L"0.0015 ",L"0.002 ",L"0.0025 ",L"0.003 "]),
    xticks=([0,0.5,1,1.5],[L"0",L"0.5",L"1",L"1.5"]),
    label=L"T=0.01 \mathrm{GeV}",
    bottom_margin=8mm,
    left_margin=5mm,
    top_margin=5mm,
    xtickfont=font(10),
    ytickfont=font(10),
    legendfont=font(12),
    )

    savefig(p,"second_derivative.pdf")



    p=scatter(x,pressure_derivative.(Ref(0.01),x,Val(0),Val(2),Ref(eos)),
    
    xlabel=L"\mu\quad \mathrm{GeV}",
    ylabel=L"$\frac{\partial^2 p(\mu,T)}{\partial \mu^2}\quad \mathrm{GeV}^2$",
    yticks=([-0.001,-0.0005,0,0.0005,0.0010,0.0015,0.002,0.0025,0.003],[L"-0.001",L"-0.0005",L"0",L"0.0005 ",L"0.0010 ",L"0.0015 ",L"0.002 ",L"0.0025 ",L"0.003 "]),
    xticks=([0,0.5,1,1.5],[L"0",L"0.5",L"1",L"1.5"]),
    label=L"T=0.01 \mathrm{GeV}",
    bottom_margin=8mm,
    left_margin=5mm,
    top_margin=5mm,
    xtickfont=font(10),
    ytickfont=font(10),
    legendfont=font(12),
    )    

scatter(x,pressure_derivative.(Ref(0.01),x,Ref(Val(0)),Ref(Val(1)),Ref(eos)))



p=plot(x,pressure_derivative.(Ref(0.01),x,Val(0),Val(1),Ref(eos)),
    
    xlabel=L"\mu\quad \mathrm{GeV}",
    ylabel=L"$\frac{\partial p(\mu,T)}{\partial \mu}\quad \mathrm{GeV}^3$",
    yticks=([0,2*10^(-7),4*10^(-7),6*10^(-7),8*10^(-7),1*10^(-6),1.2*10^(-6)],[L"0",L"2\times 10^{-7} ",L"4\times 10^{-7} ",L"6\times 10^{-7} ",L"8\times 10^{-7} ",L"1\times 10^{-6}",L"1.2\times 10^{-6}"]),
    xticks=([0,0.5,1,1.5],[L"0",L"0.5",L"1",L"1.5"]),
    label=L"T=0.01 \mathrm{GeV}",
    bottom_margin=8mm,
    left_margin=5mm,
    top_margin=5mm,
    xtickfont=font(10),
    ytickfont=font(10),
    legendfont=font(12),
    )


    

using BenchmarkTools
    @benchmark thermodynamics(.01,0.1,$eos)

    @code_llvm thermodynamics(.01,0. 1,eos)=#