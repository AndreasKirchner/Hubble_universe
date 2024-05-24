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
        
        m_nucleon=0.939#0.939
        m_sigma=0.640#0.5#0.55#0.5
        m_omega=0.783
        g_sigma=10.#10.1#9.69719#8.25284#9.69719
        g_omega=9.5#9.25544#13.2025#13.5992#13.2025
        sigma_0=0.0697#0.0444258#0.0698#0.0444258
        omega_0=0.0#-0.018#0.0#-0.026538#-0.018#-0.026538

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
    muStar= mu + wal.g_omega*omega
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
            omegaEq = sqrt(2)*mNStar^(3/2)*(-mNStar+muStar)^(3/2)/(3*pi^2) + wal.m_omega^2 * omega /(4*wal.g_omega) #was +
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

        sigmaEq = T*sqrt(mNStar*T)/(4*sqrt(2)*pi^(3/2))*(3*T*(FD32pos+FD32neg)-2*mNStar*(FD12pos+FD12neg)) - wal.m_sigma^2/(4*wal.g_sigma)*sigma
        omegaEq = (mNStar*T)^(3/2)*(FD12pos-FD12neg)/((2*pi)^(3/2)) - wal.m_omega^2/(4*wal.g_omega)*omega
       # @show negArg,posArg
        return SVector{2}(sigmaEq,omegaEq)
    end
end



function  gapequation(T,mu,wal::WaleckaModel1{N}) where {N}
    _f(u,p)=system_of_equation(u,T,mu,wal)
    #if mu>0.9
#        u0=SVector{2,N}(0.0*wal.sigma_0,0.0*wal.omega_0)
    #else
      #  u0=SVector{2,N}(0.0*wal.sigma_0,0.0*wal.omega_0)
    #end
   # u0=SVector{2,N}(1.,-1.)
    u0=SVector{2,N}(wal.sigma_0,wal.omega_0)
    problm=NonlinearProblem{false}(_f, u0) 
    solve(problm, NewtonRaphson(), abstol = 1e-20, reltol = 1e-20)
end

"""
using Plots
eos=WaleckaModel1()
murang=collect(0.8:0.001:1.2)

sigma=zeros(length(murang))
omega=zeros(length(murang))

for i in collect(1:1:length(murang))
    sigma[i], omega[i] = gapequation(0.0,murang[i],eos)
end


plot(murang,sigma,label="sigma")
plot!(murang,omega,label="omega")

sigma, mu = gapequation.(Ref(0.0),murang,Ref(eos))

sigma
"""
function thermodynamic(T::A,mu::B,wal::WaleckaModel1{L}) where {A,B,L}

    u=gapequation(T,mu,wal::WaleckaModel1{L})

    sigma = u[1]
    omega = u[2]
 
    mNStar = wal.m_nucleon - wal.g_sigma * sigma 
    muStar = mu + wal.g_omega * omega 
 
    
    if T<0.001 #check for T=0 case
        if muStar>mNStar
            pressure = (2*mNStar)^(3/2)*(muStar-mNStar)^(5/2)/(15*pi^2) -1/2*wal.m_sigma^2*sigma^2+1/2*wal.m_omega^2*omega^2
            p10 = zero(pressure)
            p01 = sqrt(2)*mNStar^(3/2)*(muStar-mNStar)^(3/2)/(3*pi^2)
            p20 = zero(pressure)
            p11 = zero(pressure)
            p02 = mNStar^(3/2)*sqrt(muStar-mNStar)/(sqrt(2)*pi^2)
            return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
            #println("hallo1")
            #@show (2*mNStar)^(3/2)*(muStar-mNStar)^(5/2)/(15*pi^2)
            #@show Thermodynamic(1.,(1.,1.),(1.,1.,1.))
            #return Thermodynamic(1.,(1.,1.),(1.,1.,1.))
        else
            #println(-1/2*wal.m_sigma^2*sigma^2)
            #println(+ 1/2*wal.m_omega^2*omega^2)
            pressure = -1/2*wal.m_sigma^2*sigma^2 + 1/2*wal.m_omega^2*omega^2
            p10 = zero(pressure)
            p01 = zero(pressure)
            p20 = zero(pressure)
            p11 = zero(pressure)
            p02 = zero(pressure)
            #@show -1/2*wal.m_sigma^2*sigma^2+1/2*wal.m_omega^2*omega^2
            #println("hallo2")
            return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
            #return Thermodynamic(1.,(1.,1.),(1.,1.,1.))
        end
        
    end
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
    pressure = +1/2*wal.m_sigma^2*sigma^2-1/2*wal.m_omega^2*omega^2
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
    pressure = sqrt(2)*T*(mNStar*T)^(3/2)*(FD32pos+FD32neg)/(pi^(3/2))+1/2*wal.m_sigma^2*sigma^2-1/2*wal.m_omega^2*omega^2

    #first derivatives
    p10 = mNStar*sqrt(mNStar*T)*(5*T*(FD32pos+FD32neg)+2*(mNStar-muStar)*FD12pos+2*(mNStar+muStar)*FD12neg)/(sqrt(2)*pi^(3/2))
    p01 = sqrt(2)*(mNStar*T)^(3/2)*(FD12pos+FD12neg)/(pi^(3/2))
    #second derivatives
    p20 = (mNStar^3*(15*T^2*FD32pos+15*T^2*FD32neg+4*(mNStar-muStar)*(3*T*FD12pos+(mNStar-muStar)*FDm12pos)+4*(mNStar+muStar)*(3*T*FD12neg+(mNStar+muStar)*FDm12neg)))/(2*sqrt(2)*pi^(3/2)*(mNStar*T)^(3/2))
    p11 = (mNStar^2*(3*T*(FD12pos+FD12neg)+2*(mNStar-muStar)*FDm12pos+2*(mNStar+muStar)*FDm12neg))/(sqrt(2)*pi^(3/2)*sqrt(mNStar*T))
    p02 = sqrt(2)*mNStar*sqrt(mNStar*T)*(FDm12pos+FDm12neg)/(pi^(3/2))
    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
end
"""

#here I want to get the potential. for this I need a fixed T and mu, together with a fixed value for one field and then solve for the other
function system_of_equation_fixed_sigma(sigma,omega,T,mu,wal::WaleckaModel1{N}) where {N}
    #Get the fields and T and mu
    #sigma = u[1]
    #omega = u[2]
    #T = 1/β
    #mu = α/β
    #Define the shifted mass and chemical potential
    mNStar = wal.m_nucleon - wal.g_sigma * sigma 
    muStar= mu +wal.g_omega*omega
    @show muStar
    @show mNStar
    #check for negative mass from solver
    if mNStar<=0
        sigmaEq = wal.m_sigma^2* sigma /(4*wal.g_sigma)
        omegaEq = wal.m_omega^2 * omega /(4*wal.g_omega)
        #return SVector{2}(sigmaEq,omegaEq)
        return omegaEq
    end
    if T<0.001 #check for T=0 case
        #theta function at low T
        if muStar>mNStar
            #sigmaEq = sqrt(2*mNStar)*(-8*mNStar+3*muStar)*(-mNStar+muStar)^(3/2)/(15*pi^2) + wal.m_sigma^2* sigma /(4*wal.g_sigma)
            omegaEq = sqrt(2)*mNStar^(3/2)*((-mNStar+muStar)^(3/2))/(3*pi^2) + wal.m_omega^2 * omega /(4*wal.g_omega)
           # sigmaEq = sqrt(2*mNStar)*(-8*mNStar+3*muStar)*(-mNStar+muStar)^(3/2)/(15*pi^2) + wal.m_sigma^2* sigma /(4*wal.g_sigma)
            #omegaEq = sqrt(2*mNStar^3)*(-mNStar+muStar)^(3/2)/(3*pi^2) + wal.m_omega^2 * omega /(4*wal.g_omega)
            #return SVector{2}(sigmaEq,omegaEq)
            #print(omegaEq)
            @show wal.m_omega^2/(4*wal.g_omega)
            return omegaEq
            print("1")
        else 
         #   sigmaEq = wal.m_sigma^2* sigma /(4*wal.g_sigma)
            omegaEq = wal.m_omega^2 * omega /(4*wal.g_omega)
            print("2")
            #return SVector{2}(sigmaEq,omegaEq)
            return omegaEq
        end
    end
    #These are the arguments for the Fermi integrations
    negArg=(-muStar-mNStar)/T
    posArg=(muStar-mNStar)/T

    if T<0.007
        if mNStar<=0 #negative mass from wrong solution
            println("mass negative")
            #return SVector{2}(wal.m_sigma^2* sigma /(4*wal.g_sigma),wal.m_omega^2 * omega /(4*wal.g_omega))
            return wal.m_omega^2 * omega /(4*wal.g_omega)
        end
        if posArg<0 #theta function at low T
            #return SVector{2}(wal.m_sigma^2* sigma /(4*wal.g_sigma),wal.m_omega^2 * omega /(4*wal.g_omega))
            return wal.m_omega^2 * omega /(4*wal.g_omega)
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
        #return SVector{2}(sigmaEq,omegaEq)
        return omegaEq
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
        #return SVector{2}(sigmaEq,omegaEq)
        return omegaEq
    end
end


function  gapequation_fixed_sigma(sigma,T,mu,wal::WaleckaModel1{N}) where {N}
    _f(omega,p)=system_of_equation_fixed_sigma(sigma,omega,T,mu,wal)
    #if mu>0.9
       # u0=SVector{2,N}(wal.sigma_0,wal.omega_0)
       omega0= wal.omega_0
    #else
      #  u0=SVector{2,N}(0.0*wal.sigma_0,0.0*wal.omega_0)
    #end
    problm=NonlinearProblem{false}(_f, omega0) 
    sol= solve(problm, NewtonRaphson(), abstol = 1e-20, reltol = 1e-20)
    return sol[1]
end

gapequation_fixed_sigma(0.,0.0,0.99,eos)
system_of_equation_fixed_sigma(0.0,0.0,0.0,0.0,eos)
sigmaRange = collect(-0.01:0.001:0.05)
T=0.0
#mu=0.9
function Tline(T)
    return 0.91684 - 11.6948*T^2 - 12096.1799*T^4
end
mu=0.9
#mu=Tline(T)
omega = gapequation_fixed_sigma.(sigmaRange,Ref(T),Ref(mu),Ref(eos))
omega[60]
plot(sigmaRange,omega)

gapequation_fixed_sigma(0.01,0.0,0.9,eos)

system_of_equation_fixed_sigma(0,0.01,0,0.9,eos)

using Plots
m_nucleon=0.939#0.939
m_sigma=0.640#0.5#0.55#0.5
m_omega=0.783#0.783
g_sigma=10.1#9.69719#8.25284#9.69719
g_omega=9.25544#13.2025#13.5992#13.2025
sigma_0=0.0697#0.0444258#0.0698#0.0444258
omega_0=0.0#-0.018#0.0#-0.026538#-0.018#-0.026538

WaleckaG=WaleckaModel1()
sigmaRange = collect(-0.2:0.001:0.2)
T=0.0
#mu=0.9
function Tline(T)
    return 0.91684 - 11.6948*T^2 - 12096.1799*T^4
end
mu=1.2
#mu=Tline(T)
omegaList=zeros(length(sigmaRange))
mNStar=zeros(length(sigmaRange))
muStar=zeros(length(sigmaRange))
omega = gapequation_fixed_sigma.(sigmaRange,Ref(T),Ref(mu),Ref(WaleckaG))
mNStar .= m_nucleon .- sigmaRange .* g_sigma
muStar .= mu .+ g_omega .* omega
mNStarp=check_neg_mass(mNStar)
press=zeros(length(sigmaRange))
diffp=zeros(length(sigmaRange))
diffp=check_neg_mass(muStar-mNStarp)

press .= (2 .* mNStarp).^(3 ./ 2) .*(diffp) .^(5 ./2) ./(15 .*pi .^2) .-1 ./2* m_sigma .^2 .* sigmaRange .^2 .+1 ./2 .*m_omega .^2 .*omega .^2
potPlot=plot!(sigmaRange,-press,label="μ=1.2 GeV",xlabel="σ [MeV]",ylabel="U [MeV^4]")
savefig(potPlot,"pot.png")
plot(sigmaRange,omega)

omega[1]

1 .+ sigmaRange

omegaList=zeros(length(sigmaRange))
omegaList2=zeros(length(sigmaRange))
potential=zeros(length(sigmaRange))
mNStar=zeros(length(sigmaRange))
muStar=zeros(length(sigmaRange))
posArg=zeros(length(sigmaRange))
negArg=zeros(length(sigmaRange))
FD32pos=zeros(length(sigmaRange))
FD32neg=zeros(length(sigmaRange))
press=zeros(length(sigmaRange))

 mNStar .= m_nucleon .- sigmaRange .* g_sigma
muStar .= mu .+ g_omega .* omega
function check_neg_mass(list)
    for i in collect(1:1:length(list))
        if list[i]<0
            list[i]=0
        end
    end
    return list
end
mNStarp=check_neg_mass(mNStar)

    muStar = mu + wal.g_omega * omega 
    @. negArg = (-muStar-mNStar)/T
    @. posArg = (muStar-mNStar)/T
plot(posArg)
posArg[1]
 #  @. FD12pos = sf_fermi_dirac_half(posArg)
     @.       FD32pos = sf_fermi_dirac_3half(posArg)
     @.       FD32neg = sf_fermi_dirac_3half(negArg)

   @.  press = sqrt(2)*T*(mNStarp*T)^(3/2)*(FD32pos+0*FD32neg)/(pi^(3/2))-1/2*m_sigma^2*sigmaRange^2+1/2*m_omega^2*omega^2
for i in collect(1:1:length(sigmaRange))
    omegaList[i]=omega[i,1]
    omegaList2[i]=omega[i,1]^2
    potential[i]=0.5*m_omega*omega[i,1]^2+0.5*m_sigma*sigmaRange[i]^2+ sqrt(2)*mNStar*sqrt(mNStar*T)*(FDm12pos+FDm12neg)/(pi^(3/2))
end
#plot(mNStar)
#plot(FD32pos)
#plot(sigmaRange,potential)
#plot(press)
#plot(sigmaRange)
#plot!(sigmaRange,potential)
a=plot(sigmaRange,-press,ylabel="-pressure",xlabel="sigma")
savefig(a,"potential.png")
"""