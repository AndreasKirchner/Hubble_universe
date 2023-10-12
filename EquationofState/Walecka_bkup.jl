#include("thermodynamic.jl") 
include("Polyloghack/Polylogarithmshack.jl")

using .PolylogarithmsHack: polylog
using .PolylogarithmsHack: polylog_exp_norm
using StaticArrays
#using StaticNumbers
#using ForwardDiff:
using NonlinearSolve




#polylog(x::Number,d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(polylog(x,ForwardDiff.value(d)), polylog(x-1,ForwardDiff.value(d))/ForwardDiff.value(d) * ForwardDiff.partials(d))

#@inline pressureFermiDirac(T,mu,mass)=-T*(mass* T /(2*pi))^(3/2)  * polylog(5/2, -exp((mu-mass)/T))

#@inline pressure_derivativeFermiDirac(T,mu,mass,::Val{0},::Val{1},::Val{0})=-(mass* T /(2*pi))^(3/2)  * polylog(3/2, -exp((mu-mass)/T))

#@inline pressure_derivativeFermiDirac(T,mu,mass,::Val{0},::Val{2},::Val{0})=-(mass* T /(2*pi))^(3/2) /T  * polylog(1/2, -exp((mu-mass)/T))

#@inline function pressure_derivativeFermiDirac(T::Real,mu::Real,mass::Real,::Val{0},::Val{0},::Val{1})
#    (mass* T /(2*pi))^(3/2)  * (
#   polylog(3/2, -exp((mu-mass)/T))-
#    3*T/(2*mass)*polylog(5/2, -exp((mu-mass)/T)))
#end

#@inline function pressure_derivativeFermiDirac(T,mu,mass,::Val{0},::Val{0},::Val{2})
# T(4*mass^2 * polylog(1/2, -exp((mu-mass)/T) ) +3*T*(
#   -4*mass* polylog(3/2, -exp((mu-mass)/T) )
#   + T* polylog(5/2, -exp((mu-mass)/T) )
# ))/(4*pi^2* sqrt(mass*T/(2*pi)))
#end

#@inline function pressure_derivativeFermiDirac(T,mu,mass,::Val{0},::Val{1},::Val{1})
#    sqrt(m* T /(2*pi))*(2*mass * polylog(1/2, -exp((mu-mass)/T) ) -3*T*
#    polylog(3/2, -exp((mu-mass)/T)) )/(2*pi)
#end



struct WaleckaModel{T} <:EquationOfState
    m_nucleon::T 
    m_sigma::T 
    m_omega::T 
    g_sigma::T 
    g_omega::T 
    sigma_0::T 
    omega_0::T 

    function WaleckaModel() 
        
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

function Base.show(io::IO, z::WaleckaModel{T}) where T 
    print(io,"Non Relativistic Walecka Model: σ=",z.sigma_0," GeV, ω₀=",z.omega_0," GeV"  )
end


#function Base.show(io::IO, ::MIME"text/plain", z::HadronResonaceGas) where{T}
#    min ,max = extrema(z.particle_list.Mass)
#    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ [",min,",",max,"] GeV\n " )
#    for part in z.particle_list
#        print(io,part.Name," mass=",part.Mass,"\n")
#    end 
#end



@inline function system_of_equation(u::M,T,μ,wal::WaleckaModel{N}) where {M<:AbstractArray,N}

    #m_nucleon,m_sigma,m_omega,g_sigma,g_omega=p
    
    mustar=μ +wal.g_omega*u[2]
    mstar=wal.m_nucleon-wal.g_sigma*u[1]


    
    zplus=-exp((mustar-mstar)/T)
    zminus=-exp((-mustar-mstar)/T)
    
    mtfactor=(mstar*T/(2*pi))^(3/2)
    
    lie32p=polylog(3/2,zplus)
    lie52p=polylog(5/2,zplus)
    

    lie32m=polylog(3/2,zminus)
    lie52m=polylog(5/2,zminus)
    
    a=mtfactor*(lie32p+lie32m -3*T/(2*mstar)*( lie52p+lie52m))+ wal.m_sigma^2/(4*wal.g_sigma)*u[1]
    b= -mtfactor*(lie32p-lie32m) + wal.m_omega^2/(4*wal.g_omega)*u[2]
    
    #a=pressure_derivativeFermiDirac(T,mustar,mstar,Val(0),Val(0),Val(1))+pressure_derivativeFermiDirac(T,-mustar,mstar,Val(0),Val(0),Val(1)) + wal.m_sigma^2/(4*wal.g_sigma)*u[1]
    
    
    #b=pressure_derivativeFermiDirac(T,mustar,mstar,Val(0),Val(1),Val(0))-pressure_derivativeFermiDirac(T,-mustar,mstar,Val(0),Val(1),Val(0)) + wal.m_omega^2/(4*wal.g_omega)*u[2]

    return SVector{2}(a,b)

end



function gapequation(T,mu,wal::WaleckaModel{N}) where {N}
    _f(u,p)=system_of_equation(u,T,mu,wal)
    problm=NonlinearProblem{false}(_f, SVector{2,N}(wal.sigma_0,wal.omega_0))
    solve(problm,NewtonRaphson(), tol = 1e-6)  
end


#struct Thermodyanmic{T,N} 
#    pressure::T
#    pressure_derivative::SVector{N,T}
#    pressure_hessian::SHermitianCompact{N,T}
#end




#Thermodyanmic(a::T,b::NTuple{N,T},c::NTuple{M,T}) where {T,N,M}= Thermodyanmic{T,N}(a,SVector{N,T}(b),SHermitianCompact{N,T}(SVector{M}(c))) 

#Thermodyanmic(a::T,b::SVector{N,T},c::SVector{M,T}) where {T,N,M}= Thermodyanmic{T,N}(a,b,SHermitianCompact{N,T}(c)) 





 
function thermodynamic(T::A,mu::B,wal::WaleckaModel{L}) where {A,B,L}

    u=gapequation(T,mu,wal::WaleckaModel{L})

    mustar=mu +wal.g_omega*u[2]
    mstar=wal.m_nucleon-wal.g_sigma*u[1]
    mSigma2=wal.m_sigma^2
    gSigma2=wal.g_sigma^2
    mOmega2=wal.m_omega^2
    gOmega2=wal.g_omega^2
    zplus=-exp((mustar-mstar)/T)
    zminus=-exp((-mustar-mstar)/T)
    
    mtfactor=(mstar*T/(2*pi))^(3/2)

    lie32p=polylog(3/2,zplus)
    lie52p=polylog(5/2,zplus)
    lie12p=polylog(1/2,zplus)

    lie32m=polylog(3/2,zminus)
    lie52m=polylog(5/2,zminus)
    lie12m=polylog(1/2,zminus)

    # derivative of the fermi dirac distro
    p100p= mtfactor*(2*(-mstar+mustar)*lie32p-5*T*lie52p)/(2*T)
    p100m= mtfactor*(2*(-mstar-mustar)*lie32p-5*T*lie52m)/(2*T)
    p010p= -mtfactor*lie32p
    p010m= -mtfactor*lie32m

    p200p=-mtfactor/(4*T^3)*(4*(mstar-mustar)*((mstar-mustar)*lie12p+3*T*lie32p) +15*T^2*lie52p)
    p110p=-mtfactor/(2*T^2)*(2*(mstar-mustar)*lie12p+3*T*lie32p)
    p101p=mtfactor/(4*mstar*T^2)*(4*mstar*(mstar-mustar)*lie12p +T*(6*mustar*lie32p- 15*T*lie52p ) )
    p011p=mtfactor/(2*mstar*T)*(2*mstar*lie12p -3*T*lie32p)
    p020p= -mtfactor/T*lie12p
    p002p=mtfactor/(4*mstar*T)*( -4*mstar*lie12p+ 3*T*(   4*mstar*lie32p-T*lie52p))
    
    p200m=-mtfactor/(4*T^3)*(4*(mstar+mustar)*((mstar+mustar)*lie12m+3*T*lie32m) +15*T^2*lie52m)
    p110m=-mtfactor/(2*T^2)*(2*(mstar+mustar)*lie12m+3*T*lie32m)
    p101m=mtfactor/(4*mstar*T^2)*(4*mstar*(mstar+mustar)*lie12m +T*(-6*mustar*lie32p- 15*T*lie52m ) )
    p011m=mtfactor/(2*mstar*T)*(2*mstar*lie12m -3*T*lie32m)
    p020m= -mtfactor/T*lie12m
    p002m=mtfactor/(4*mstar*T)*( -4*mstar*lie12m+ 3*T*(   4*mstar*lie32m-T*lie52m))
    




    pressure=-4*T*mtfactor*(lie52m+lie52p)  
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

 #pressure(T,mu,wal::WaleckaModel{L}) where {L}= thermodynamics(T,mu,wal).pressure
#
 #pressure_derivative(T,mu,::Val{1},::Val{0},wal::WaleckaModel{L}) where {L}= thermodynamics(T,mu,wal).pressure_derivative[1]
 #pressure_derivative(T,mu,::Val{0},::Val{1},wal::WaleckaModel{L}) where {L}= thermodynamics(T,mu,wal).pressure_derivative[2]
 #pressure_derivative(T,mu,::Val{2},::Val{0},wal::WaleckaModel{L}) where {L}= thermodynamics(T,mu,wal).pressure_hessian[1]
 #pressure_derivative(T,mu,::Val{1},::Val{1},wal::WaleckaModel{L}) where {L}= thermodynamics(T,mu,wal).pressure_hessian[2]
 #pressure_derivative(T,mu,::Val{0},::Val{2},wal::WaleckaModel{L}) where {L}= thermodynamics(T,mu,wal).pressure_hessian[3]
    


eos=WaleckaModel()





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