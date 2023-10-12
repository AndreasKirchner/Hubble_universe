using StaticArrays
using StaticNumbers
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


function system_of_equation(u::M,T,mu,wal::WaleckaModel1{N}) where {M<:AbstractArray,N}
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



function  gapequation(T,mu,wal::WaleckaModel1{N}) where {N}
    _f(u,p)=system_of_equation(u,T,mu,wal)
    #if mu>0.9
        u0=SVector{2,N}(wal.sigma_0,wal.omega_0)
    #else
      #  u0=SVector{2,N}(0.0*wal.sigma_0,0.0*wal.omega_0)
    #end
    problm=NonlinearProblem{false}(_f, u0) 
    solve(problm, NewtonRaphson(), abstol = 1e-20, reltol = 1e-20)
end


function thermodynamic(T::A,mu::B,wal::WaleckaModel1{L}) where {A,B,L}

    u=gapequation(T,mu,wal::WaleckaModel1{L})

    sigma = u[1]
    omega = u[2]
 
    mNStar = wal.m_nucleon - wal.g_sigma * sigma 
    muStar = mu + wal.g_omega * omega 
 
    negArg = (-muStar-mNStar)/T
    posArg = (muStar-mNStar)/T
 
    if posArg<-500 #antiparticle contribution can get to small
        FD12pos = zero(posArg)
        FD32pos = zero(posArg)
        FDm12pos = zero(posArg)
     else
        FD12pos = sf_fermi_dirac_half(posArg)
        FD32pos = sf_fermi_dirac_3half(posArg)
        FDm12pos = sf_fermi_dirac_mhalf(posArg)
     end
    #FD12pos = sf_fermi_dirac_half(posArg)
    #FD32pos = sf_fermi_dirac_3half(posArg)
    if negArg<-500 #antiparticle contribution can get to small
        FD12neg = zero(negArg)
        FD32neg = zero(negArg)
        FDm12neg = zero(negArg)
     else
        FD12neg = sf_fermi_dirac_half(negArg)
        FD32neg = sf_fermi_dirac_3half(negArg)
        FDm12neg = sf_fermi_dirac_mhalf(negArg)
     end
    #FD12pos = sf_fermi_dirac_half(posArg)
    #FD12neg = sf_fermi_dirac_half(negArg)
    #FD32pos = sf_fermi_dirac_3half(posArg)
    #FD32neg = sf_fermi_dirac_3half(negArg)
    #FDm12pos = sf_fermi_dirac_mhalf(posArg)
    #FDm12neg = sf_fermi_dirac_mhalf(negArg)
     if mNStar<=0
    pressure = -1/2*wal.m_sigma^2*sigma^2+1/2*wal.m_omega^2*omega^2
    #first derivatives
    p10 = 0.0
    p01 = 0.0
    #second derivatives
    p20 = 0.0
    p11 = 0.0
    p02 = 0.0

    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
        
     end
    #this is the pressure
    pressure = sqrt(2)*T*(mNStar*T)^(3/2)*(FD32pos+FD32neg)/(pi^(3/2))-1/2*wal.m_sigma^2*sigma^2+1/2*wal.m_omega^2*omega^2

    #first derivatives
    p10 = mNStar*sqrt(mNStar*T)*(5*T*(FD32pos+FD32neg)+2*(mNStar-muStar)*FD12pos+2*(mNStar+muStar)*FD12neg)/(sqrt(2)*pi^(3/2))
    p01 = sqrt(2)*(mNStar*T)^(3/2)*(FD12pos+FD12neg)/(pi^(3/2))
    #second derivatives
    p20 = (mNStar^3*(15*T^2*FD32pos+15*T^2*FD32neg+4*(mNStar-muStar)*(3*T*FD12pos+(mNStar-muStar)*FDm12pos)+4*(mNStar+muStar)*(3*T*FD12neg+(mNStar+muStar)*FDm12neg)))/(2*sqrt(2)*pi^(3/2)*(mNStar*T)^(3/2))
    p11 = (mNStar^2*(3*T*(FD12pos+FD12neg)+2*(mNStar-muStar)*FDm12pos+2*(mNStar+muStar)*FDm12neg))/(sqrt(2)*pi^(3/2)*sqrt(mNStar*T))
    p02 = sqrt(2)*mNStar*sqrt(mNStar*T)*(FDm12pos+FDm12neg)/(pi^(3/2))

    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
end
