#include("thermodynamic.jl") 
#include("Polyloghack/Polylogarithmshack.jl")

using .PolylogarithmsHack: polylog
using .PolylogarithmsHack: polylog_exp_norm
using StaticArrays
#using StaticNumbers
#using ForwardDiff:
using NonlinearSolve
using SimpleNonlinearSolve
using GSL




#polylog_exp_norm(x::Number,d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(polylog_exp_norm(x,ForwardDiff.value(d)), polylog_exp_norm(x-1,ForwardDiff.value(d))/ForwardDiff.value(d) * ForwardDiff.partials(d))

#@inline pressureFermiDirac(T,mu,mass)=-T*(mass* T /(2*pi))^(3/2)  * polylog_exp_norm(5/2, -exp((mu-mass)/T))

#@inline pressure_derivativeFermiDirac(T,mu,mass,::Val{0},::Val{1},::Val{0})=-(mass* T /(2*pi))^(3/2)  * polylog_exp_norm(3/2, -exp((mu-mass)/T))

#@inline pressure_derivativeFermiDirac(T,mu,mass,::Val{0},::Val{2},::Val{0})=-(mass* T /(2*pi))^(3/2) /T  * polylog_exp_norm(1/2, -exp((mu-mass)/T))

#@inline function pressure_derivativeFermiDirac(T::Real,mu::Real,mass::Real,::Val{0},::Val{0},::Val{1})
#    (mass* T /(2*pi))^(3/2)  * (
#   polylog_exp_norm(3/2, -exp((mu-mass)/T))-
#    3*T/(2*mass)*polylog_exp_norm(5/2, -exp((mu-mass)/T)))
#end

#@inline function pressure_derivativeFermiDirac(T,mu,mass,::Val{0},::Val{0},::Val{2})
# T(4*mass^2 * polylog_exp_norm(1/2, -exp((mu-mass)/T) ) +3*T*(
#   -4*mass* polylog_exp_norm(3/2, -exp((mu-mass)/T) )
#   + T* polylog_exp_norm(5/2, -exp((mu-mass)/T) )
# ))/(4*pi^2* sqrt(mass*T/(2*pi)))
#end

#@inline function pressure_derivativeFermiDirac(T,mu,mass,::Val{0},::Val{1},::Val{1})
#    sqrt(m* T /(2*pi))*(2*mass * polylog_exp_norm(1/2, -exp((mu-mass)/T) ) -3*T*
#    polylog_exp_norm(3/2, -exp((mu-mass)/T)) )/(2*pi)
#end



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


function system_of_equation2(u::M,T,μ,wal::WaleckaModel1{N}) where {M<:AbstractArray,N}
   
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

    sigmaEq = T*sqrt(mNStar*T)/(4*sqrt(2)*pi^(3/2))*(3*T*(FD32pos+FD32neg)-2*mNStar*(FD12pos+FD12neg)) + wal.m_sigma^2/(4*wal.g_sigma)*sigma
    omegaEq = (mNStar*T)^(3/2)*(FD12pos-FD12neg)/((2*pi)^(3/2)) + wal.m_omega^2/(4*wal.g_omega)*omega

    return SVector{2}(sigmaEq,omegaEq)

end


function  gapequation2(T,mu,wal::WaleckaModel1{N}) where {N}
    _f(u,p)=system_of_equation2(u,T,mu,wal)
    if mu>1
        u0=SVector{2}(0.0444258,-0.026538)
        else
    u0=SVector{2}(0,0)
        end

    problm=NonlinearProblem{false}(_f, u0) 
    solve(problm, NewtonRaphson(), abstol = 1e-10)
    #solve(problm, DynamicSS(), abstol = 1e-10)
end



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




#Thermodyanmic(a::T,b::NTuple{N,T},c::NTuple{M,T}) where {T,N,M}= Thermodyanmic{T,N}(a,SVector{N,T}(b),SHermitianCompact{N,T}(SVector{M}(c))) 

#Thermodyanmic(a::T,b::SVector{N,T},c::SVector{M,T}) where {T,N,M}= Thermodyanmic{T,N}(a,b,SHermitianCompact{N,T}(c)) 





 
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


#gapequation2(0.2,1,eos)



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