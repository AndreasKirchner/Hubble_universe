using StaticArrays
using StaticNumbers
using NonlinearSolve
using SimpleNonlinearSolve
using GSL
using NLsolve
using Profile
using PProf
using InteractiveUtils
using BenchmarkTools


struct WaleckaModel2{T} <:EquationOfState
    fπ::T
    mπ::T
    m_omega::T
    h::T
    g::T
    lambda::T
    gamma3::T
    gamma4::T
    omega0::T
    sigma0::T

    function WaleckaModel2() 
        
        fπ=0.093
        mπ=0.135
        m_omega=0.783
        h=10.096774193548388
        g=9.5
        lambda=50.
        gamma3=2.8062
        gamma4=50.41027
        omega0=-0.018
        sigma0=0.0698

        new{typeof(fπ)}(fπ,mπ,m_omega,h,g,lambda,gamma3,gamma4,omega0,sigma0)
    end
end

function Base.show(io::IO, z::WaleckaModel2{T}) where T 
    print(io,"Non Relativistic Walecka Model: σ=",z.sigma0," GeV, ω₀=",z.omega0," GeV"  )
end


function system_of_equation(u::M,T,mu,wal::WaleckaModel2{N}) where {M<:AbstractArray,N}
    sigma = exp(u[1])
    omega = u[2]
    shiftedSigma2 = sigma^2 -wal.fπ^2
    muStar = mu + wal.g*omega
    if muStar>wal.h*sigma
        sigmaEq = -wal.fπ*wal.mπ^2 + wal.mπ^2*sigma + 1/2 *wal.lambda*sigma*shiftedSigma2 + 2*wal.gamma3*sigma*shiftedSigma2^2/wal.fπ^2 + 2*wal.gamma4*sigma*shiftedSigma2^3/wal.fπ^4 +(4*sqrt(2)*wal.h*(wal.h*sigma)^(3/2)*(muStar-wal.h*sigma)^(3/2))/(3*pi^2)-(4*sqrt(2)*wal.h*sqrt(wal.h*sigma)*(muStar-wal.h*sigma)^(5/2))/(5*pi^2)
        omegaEq = -wal.m_omega^2*omega - 4*wal.g *(sqrt(2)*(wal.h*sigma)^(3/2)*(muStar-wal.h*sigma)^(3/2))/(3*pi^2)
    else
        sigmaEq = -wal.fπ*wal.mπ^2 + wal.mπ^2*sigma + 1/2 *wal.lambda*sigma*shiftedSigma2 + 2*wal.gamma3*sigma*shiftedSigma2^2/wal.fπ^2 + 2*wal.gamma4*sigma*shiftedSigma2^3/wal.fπ^4 
        omegaEq = -wal.m_omega^2*omega
    end
    return SVector{2}(sigmaEq,omegaEq)
end

function  gapequation(T,mu,wal::WaleckaModel2{N}) where {N}
    _f(u,p)=system_of_equation(u,T,mu,wal)
   u0=SVector{2,N}(log(0.0689),-0.018)
    problm=NonlinearProblem{false}(_f, u0) 
    solve(problm, NewtonRaphson(), abstol = 1e-20, reltol = 1e-20)
end


#Fermi pressure
function pFG(T,mu,m)
    if m<0
        m=zero(m)
    end
    if T<0.001 # T=0 case
        if mu<m
            return zero(mu)
        else
            return (2* m)^(3/2)*(mu - m)^(5/2)/(15*pi^2)
        end
    end
    if (mu-m)/T<-600
        return zero(T)
    end
    return T*(m*T/(2*pi))^(3/2)*sf_fermi_dirac_3half((mu-m)/T)
end


function vacuumPotential(sigma,omega)
    fπ=0.093
        mπ=0.135
        m_omega=0.783
        h=10.096774193548388
        g=9.5
        lambda=50.
        gamma3=2.8062
        gamma4=50.41027
        omega0=-0.018
        sigma0=0.0698
   shiftedSigma2 = sigma^2 - fπ^2
   return 1/2*mπ^2*shiftedSigma2 + 1/8*lambda*shiftedSigma2^2+1/3*gamma3*shiftedSigma2^3/fπ^2+1/4*gamma4*shiftedSigma2^4/fπ^4-mπ^2*fπ*(sigma-fπ)-1/2*m_omega^2*omega^2

end

function effective_potential(T,mu,sigma,omega)#,wal::WaleckaModel2{N}) where N
    fπ=0.093
        mπ=0.135
        m_omega=0.783
        h=10.096774193548388
        g=9.5
        lambda=50.
        gamma3=2.8062
        gamma4=50.41027
        omega0=-0.018
        sigma0=0.0698
        mStar=if sigma>0 #hotfix
            h*sigma
        else
            zero(sigma)
        end
        if (-mu-g*omega-h*sigma)/T>-600
            antiParticleContribution = 4*pFG(T,-(mu+g*omega),h*sigma)
        else
            antiParticleContribution = zero(T)
        end
    vacuumPotential(sigma,omega)-4*pFG(T,mu+g*omega,h*sigma)-antiParticleContribution
end

## this sovles for the omega ##

function pterm_omega(T,muStar,mStar)
    h=10.096774193548388
    g=9.5
    if T<0.001
        if muStar>mStar
            return (4*sqrt(2)*g*(mStar)^(3/2)*(muStar-mStar)^(3/2))/(3*pi^2)
        else
            return zero(T)
        end
    end
    posArg=(muStar-mStar)/T
    negArg=(-muStar-mStar)/T
    if posArg>-600
        fd12posArg=sf_fermi_dirac_half(posArg)
    else
        fd12posArg=zero(T)
    end
    if negArg>-600
        fd12negArg=sf_fermi_dirac_half(negArg)
    else
        fd12negArg=zero(T)
    end
    return sqrt(2)*g*(T*mStar)^(3/2)*(fd12posArg-fd12negArg)/(pi^(3/2))
end


function omega_min_fun(omega,p)
    fπ=0.093
    mπ=0.135
    m_omega=0.783
    h=10.096774193548388
    g=9.5
    T=p[1]
    mu=p[2]
    sigma=p[3]
    muStar=mu +g*omega
    mStar=if sigma>0 #hotfix
        h*sigma
    else
        zero(sigma)
    end
    pterm=pterm_omega(T,muStar,mStar)
    return m_omega^2*omega+ pterm
end


function find_omega(T,mu,sigma)#,wal::WaleckaModel2{N}) where N    
    p = (T,mu,sigma)#,wal)
        #_f(omega, p) = omega_min_fun(omega, p)
    omega0=10.0
    #println("find_omega T: $(T)")
    problem=NonlinearProblem{false}(omega_min_fun,omega0, p)
    sol=solve(problem,FastShortcutNonlinearPolyalg())#SimpleNewtonRaphson())#autodiff=false))
    return sol[1]
end

"""
Benchmark stuff and check for allocations

@time find_omega(0.0,0.98,0.07)#,eos)
@benchmark find_omega(0.0,0.0,0.0)#,eos)
@InteractiveUtils.code_warntype find_omega(0.0,0.0,0.0)#,eos)
Profile.Allocs.clear()
@time Profile.Allocs.@profile sample_rate=1 find_omega(0.0,0.0,0.0)#,eos)
PProf.Allocs.pprof(from_c=false)
@code_llvm find_omega(0.0,0.0,0.0,eos)
"""

## solve for sigma ##
function find_sigma(T,mu,sigma0)
    p=(T,mu)
    problem=NonlinearProblem{false}(sigma_min_fun,sigma0,p)
    sol=solve(problem,FastShortcutNonlinearPolyalg())#SimpleNewtonRaphson())#autodiff=false))
    return sol[1]
end



function sigma_min_fun(sigma,p)#,T,mu,wal::WaleckaModel2{N}) where N
    T=p[1]
    mu=p[2]
    #wal=p[3]
    fπ=0.093
    mπ=0.135
    m_omega=0.783
    h=10.096774193548388
    g=9.5
    lambda=50.
    gamma3=2.8062
    gamma4=50.41027
    omega0=-0.018
    sigma0=0.0698
    omega=find_omega(T,mu,sigma)#,wal)
    muStar=mu +g*omega
    mStar=if sigma>0
        h*sigma
    else
        zero(sigma)
    end
    shiftedSigma2 = sigma^2-fπ^2
    pterm = pterm_sigma(T,muStar,mStar)
    -fπ*mπ^2 + mπ^2*sigma +1/2*lambda*sigma*shiftedSigma2+2*gamma3*sigma*shiftedSigma2^2/fπ^2+2*gamma4*sigma*shiftedSigma2^3/fπ^4+pterm
end

function pterm_sigma(T,muStar,mStar)
    h=10.096774193548388
    g=9.5
    if T<0.001
        if muStar>mStar
            return (4*sqrt(2)*h*(mStar)^(3/2)*(muStar-mStar)^(3/2))/(3*pi^2)-(4*sqrt(2)*h*sqrt(mStar)*(muStar-mStar)^(5/2))/(5*pi^2)
        else
            return zero(T)
        end
    end
    posArg=(muStar-mStar)/T
    negArg=(-muStar-mStar)/T
    

    if posArg>-600
        fd12posArg=sf_fermi_dirac_half(posArg)
        fd32posArg=sf_fermi_dirac_3half(posArg)
    else
        fd12posArg=zero(T)
        fd32posArg=zero(T)
    end

    if negArg>-600
        fd12negArg=sf_fermi_dirac_half(negArg)
        fd32negArg=sf_fermi_dirac_3half(negArg)
    else
        fd12negArg=zero(T)
        fd32negArg=zero(T)
    end
    return sqrt(2)*h*T*sqrt(T*mStar)/(2*pi^(3/2))*(-3*T*(fd32posArg+fd32negArg)+2*mStar*(fd12posArg+fd12negArg))
end

"""
More benchmarks
@code_warntype sigma_min_fun(0.0,(0.0,0.0))
@benchmark sigma_min_fun(0.0,(0.0,0.0))
@time find_sigma(0.0,0.0,0.09)#,eos)
@benchmark find_sigma(0.0,0.0,0.09)
find_sigma(0.0,0.0,0.09)
@code_warntype find_sigma(0.0,0.0,0.09)
find_sigma(0.0,0.0,0.09)
Profile.Allocs.clear()
@time Profile.Allocs.@profile sample_rate=1 find_sigma(0.0,0.0,0.09)
PProf.Allocs.pprof(from_c=false)
PProf.Allocs.pprof()
"""

function get_potential_one_val(T,mu,sigmaSeed)
    sigmaVal=find_sigma(T,mu,sigmaSeed)
    omegaVal=find_omega(T,mu,sigmaVal)
    potVal=round(effective_potential(T,mu,sigmaVal,omegaVal),digits=8)
    return [sigmaVal, omegaVal, potVal]
end


function get_potential_one_val2(T,mu,sigmaSeed)
    sigmaVal=find_sigma(T,mu,sigmaSeed)
    omegaVal=find_omega(T,mu,sigmaVal)
    potVal=effective_potential(T,mu,sigmaVal,omegaVal)
    return @SVector [sigmaVal, omegaVal, potVal]
end

function get_minima(potVals)
    eps=0.01
    sortIndex=sortperm(potVals)
    sort!(potVals)
    jumpIndexList=[]
    front(itr, n=1) = Iterators.take(itr, length(itr)-n)
    for i in front(eachindex(potVals))
        delta=potVals[i+1]-potVals[i]
        if delta>eps
            push!(jumpIndexList,i)
        end
    end
    finalList=[]
    prepend!(jumpIndexList,1)
    append!(jumpIndexList,length(potVals))
    print(jumpIndexList)
    for i in front(eachindex(jumpIndexList))
        push!(finalList,mean(potVals[jumpIndexList[i]:jumpIndexList[i+1]]))
    end
    
    return finalList
end

function get_minima2(potVals)
    front(itr, n=1) = Iterators.take(itr, length(itr)-n)
    eps=0.5
    sortIndex=sortperm(potVals)
    sort!(potVals)
    finalList=[]
    push!(finalList,potVals[1])
    seedVal=potVals[1]
    for i in front(eachindex(potVals))
        if seedVal + eps < potVals[i+1]
            push!(finalList,potVals[i+1])
            seedVal = potVals[i+1]
        end
    end
return finalList

end


function filter_approx!(cur_vector; atol=1e-6)
    delete_indices = []
  
    for (cur_index, cur_value) in enumerate(cur_vector)
      is_duplicate = any(
        tmp_value -> isapprox(tmp_value, cur_value, atol=atol),
        cur_vector[1:cur_index-1]
      )
  
      is_duplicate || continue
      push!(delete_indices, cur_index)
    end
  
    deleteat!(cur_vector, delete_indices)
  
    return (cur_vector,delete_indices) 
  end





function get_fields2(T,mu,flag=false)
   # @show T
    sigmaSeedRange=0.0:0.01:0.1
    ResMatrixRaw=reduce(vcat,transpose.(get_potential_one_val2.(Ref(T),Ref(mu),sigmaSeedRange)))
    #@show ResMatrixRaw
    #this is to remove possible nans from the solver
    ResMatrix=mapreduce(permutedims, vcat, filter(x->all(isfinite.(x)),eachrow(ResMatrixRaw)))
   # @show ResMatrix
    sigmaVals=ResMatrix[:,1]
    omegaVals=ResMatrix[:,2]
    potVals=ResMatrix[:,3]

    sortIndex=sortperm(potVals)
    sortedSigmaBig=sigmaVals[sortIndex]
    sortedOmegaBig=omegaVals[sortIndex]
    sortedPotBig=potVals[sortIndex]
    #@show potVals
    #@show omegaVals

    reducedPot , deleteIndex = filter_approx!(sortedPotBig)
    #i = unique(i -> potVals[i], eachindex(potVals))

    reducedSigma=deleteat!(sortedSigmaBig,deleteIndex)#sigmaVals[i]
    reducedOmega=deleteat!(sortedOmegaBig,deleteIndex)#omegaVals[i]
    #reducedPot=#potVals[i]
    PotIndices=sortperm(reducedPot)
    sortedSigma=reducedSigma[PotIndices]
    sortedOmega=reducedOmega[PotIndices]
    sortedPot=reducedPot[PotIndices]
    if flag
        return [[sortedSigma[1],sortedOmega[1],sortedPot[1]],[sortedSigma[2],sortedOmega[2],sortedPot[2]]]
    else
        return @SVector [[sortedSigma[1],sortedOmega[1],sortedPot[1]]]
    end
end

#@benchmark get_fields2(0.1,0.9)

function transitionMu(T)
    if 0.0<=T<0.021
        return 0.92424-42.12324123*T^2-53938.75555664*T^4
    else
        return -one(T)
    end
end


#define area around the transition and check if T,mu is at the transition area
function check_for_transition(T,mu)
    eps=0.005 #size of transition area GeV
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

function thermodynamicPhaseOne(T,mu,wal::WaleckaModel2{L}) where {L}
    fields=get_fields2(T,mu,check_for_transition(T,mu))
    sigma=fields[1][1]
    #@show sigma
    omega=fields[1][2]
    pressure=-fields[1][3]
    mStar = wal.h*sigma
    muStar = mu + wal.g*omega
    if mStar<0
        mStar=0
    end
    if T<0.001
        p10=zero(pressure)
        if muStar> mStar
            p01 = 4*sqrt(2)*(mStar)^(3/2)*(muStar-wal.h*sigma)^(3/2)/(3*pi^2)
            p02=8*sqrt(2)*15/4*(mStar)^(3/2)*sqrt(muStar-mStar)/(15*pi^2)
        else
            p01 = zero(pressure)
            p02 = zero(pressure)
        end
        p20=zero(pressure)
        p11=zero(pressure)
    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
    end
    posArg=(muStar-mStar)/T
    negArg=(-muStar-mStar)/T
    if posArg>-600
        fd32posArg=sf_fermi_dirac_3half(posArg)
        fd12posArg=sf_fermi_dirac_half(posArg)
        fdm12posArg=sf_fermi_dirac_mhalf(posArg)
    else
        fd32posArg=zero(T)
        fd12posArg=zero(T)
        fdm12posArg=zero(T)
    end
    if negArg>-600
        fd32negArg=sf_fermi_dirac_3half(negArg)
        fd12negArg=sf_fermi_dirac_half(negArg)
        fdm12negArg=sf_fermi_dirac_mhalf(negArg)
    else
        fd32negArg=zero(T)
        fd12negArg=zero(T)
        fdm12negArg=zero(T)
    end

    if mStar<0
        mStar=0
    end

    p10=mStar*sqrt(mStar*T)/(sqrt(2)*pi^(3/2))*(5T*(fd32posArg+fd32negArg)-2*(muStar-mStar)*fd12posArg+2*(muStar+mStar)*fd12negArg)
    p01=sqrt(2)*(mStar*T/pi)^(3/2)*(fd12posArg-fd12negArg)

    p20=-mStar^3/(2*sqrt(2)*(pi*T*mStar)^(3/2))*(-15*T^2*fd32posArg-15*T^2*fd32negArg+12*T*(muStar-mStar)*fd12posArg-12*T*(muStar+mStar)*fd12negArg-4*(muStar-mStar)^2*fdm12posArg-4*(muStar+mStar)^2*fdm12negArg)
    p11=-(T*mStar)^(3/2)/(sqrt(2)*T^2*pi^(3/2))*(-3*T*fd12posArg+3*T*fd12negArg+2*(muStar-mStar)*fdm12posArg+2*(muStar+mStar)*fdm12negArg)
    p02=sqrt(2*T*mStar)*mStar/(pi^(3/2))*(fdm12posArg+fdm12negArg)
    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
end

function thermodynamicPhaseTwo(T,mu,wal::WaleckaModel2{L}) where {L}
    fields=get_fields2(T,mu,check_for_transition(T,mu))
    #@show fields
    print(check_for_transition(T,mu))
    sigma=fields[2][1]
    #@show sigma
    omega=fields[2][2]
    pressure=-fields[2][3]
    mStar = wal.h*sigma
    muStar = mu + wal.g*omega
    if mStar<0
        mStar=0.00001
    end
    if T<0.001
        p10=zero(pressure)
        if muStar> mStar
            p01 = 4*sqrt(2)*(mStar)^(3/2)*(muStar-wal.h*sigma)^(3/2)/(3*pi^2)
            p02=8*sqrt(2)*15/4*(mStar)^(3/2)*sqrt(muStar-mStar)/(15*pi^2)
        else
            p01 = zero(pressure)
            p02 = zero(pressure)
        end
        p20=zero(pressure)
        p11=zero(pressure)
   # return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
    end
    posArg=(muStar-mStar)/T
    negArg=(-muStar-mStar)/T
    if posArg>-600
        fd32posArg=sf_fermi_dirac_3half(posArg)
        fd12posArg=sf_fermi_dirac_half(posArg)
        fdm12posArg=sf_fermi_dirac_mhalf(posArg)
    else
        fd32posArg=zero(T)
        fd12posArg=zero(T)
        fdm12posArg=zero(T)
    end
    if negArg>-600
        fd32negArg=sf_fermi_dirac_3half(negArg)
        fd12negArg=sf_fermi_dirac_half(negArg)
        fdm12negArg=sf_fermi_dirac_mhalf(negArg)
    else
        fd32negArg=zero(T)
        fd12negArg=zero(T)
        fdm12negArg=zero(T)
    end

    if mStar<0
        mStar=0.00001
    end

    p10=mStar*sqrt(mStar*T)/(sqrt(2)*pi^(3/2))*(5T*(fd32posArg+fd32negArg)-2*(muStar-mStar)*fd12posArg+2*(muStar+mStar)*fd12negArg)
    p01=sqrt(2)*(mStar*T/pi)^(3/2)*(fd12posArg-fd12negArg)

    p20=-mStar^3/(2*sqrt(2)*(pi*T*mStar)^(3/2))*(-15*T^2*fd32posArg-15*T^2*fd32negArg+12*T*(muStar-mStar)*fd12posArg-12*T*(muStar+mStar)*fd12negArg-4*(muStar-mStar)^2*fdm12posArg-4*(muStar+mStar)^2*fdm12negArg)
    p11=-(T*mStar)^(3/2)/(sqrt(2)*T^2*pi^(3/2))*(-3*T*fd12posArg+3*T*fd12negArg+2*(muStar-mStar)*fdm12posArg+2*(muStar+mStar)*fdm12negArg)
    p02=sqrt(2*T*mStar)*mStar/(pi^(3/2))*(fdm12posArg+fdm12negArg)
    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
end


function thermodynamic(T,mu,wal::WaleckaModel2{L}) where {L}
    if T>1.0
        return Thermodynamic(zero(T),(zero(T),zero(T)),(zero(T),zero(T),zero(T)))
    end
    return thermodynamicPhaseOne(T,mu,wal)
end

function thermodynamic(T,mu,t,wal::WaleckaModel2{L}) where {L}
    if T>1.0
        return Thermodynamic(zero(T),(zero(T),zero(T)),(zero(T),zero(T),zero(T)))
    end
    if check_for_transition(T,mu)
        thermoPhaseOne=thermodynamicPhaseOne(T,mu,wal)
        thermoPhaseTwo=thermodynamicPhaseTwo(T,mu,wal)
        #@show thermoPhaseOne
        #@show thermoPhaseTwo
        return (1-t)*thermoPhaseOne+t*thermoPhaseTwo
        #return t*thermoPhaseOne+(1-t)*thermoPhaseTwo
    else
        return thermodynamicPhaseOne(T,mu,wal)
    end

end



"""
function thermodynamic(T,mu,wal::WaleckaModel2{L}) where {L}
    fields=get_fields(T,mu)
    sigma=fields[1]
    #@show sigma
    omega=fields[2]
    pressure=-fields[3]
    mStar = wal.h*sigma
    muStar = mu + wal.g*omega
    if T<0.001
        p10=zero(pressure)
        if muStar> mStar
            p01 = 4*sqrt(2)*(mStar)^(3/2)*(muStar-wal.h*sigma)^(3/2)/(3*pi^2)
            p02=8*sqrt(2)*15/4*(mStar)^(3/2)*sqrt(muStar-mStar)/(15*pi^2)
        else
            p01 = zero(pressure)
            p02 = zero(pressure)
        end
        p20=zero(pressure)
        p11=zero(pressure)
   # return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
    end
    posArg=(muStar-mStar)/T
    negArg=(-muStar-mStar)/T
    if posArg>-600
        fd32posArg=sf_fermi_dirac_3half(posArg)
        fd12posArg=sf_fermi_dirac_half(posArg)
        fdm12posArg=sf_fermi_dirac_mhalf(posArg)
    else
        fd32posArg=zero(T)
        fd12posArg=zero(T)
        fdm12posArg=zero(T)
    end
    if negArg>-600
        fd32negArg=sf_fermi_dirac_3half(negArg)
        fd12negArg=sf_fermi_dirac_half(negArg)
        fdm12negArg=sf_fermi_dirac_mhalf(negArg)
    else
        fd32negArg=zero(T)
        fd12negArg=zero(T)
        fdm12negArg=zero(T)
    end

    if mStar<0
        mStar=0
    end

    p10=mStar*sqrt(mStar*T)/(sqrt(2)*pi^(3/2))*(5T*(fd32posArg+fd32negArg)-2*(muStar-mStar)*fd12posArg+2*(muStar+mStar)*fd12negArg)
    p01=sqrt(2)*(mStar*T/pi)^(3/2)*(fd12posArg-fd12negArg)

    p20=-mStar^3/(2*sqrt(2)*(pi*T*mStar)^(3/2))*(-15*T^2*fd32posArg-15*T^2*fd32negArg+12*T*(muStar-mStar)*fd12posArg-12*T*(muStar+mStar)*fd12negArg-4*(muStar-mStar)^2*fdm12posArg-4*(muStar+mStar)^2*fdm12negArg)
    p11=-(T*mStar)^(3/2)/(sqrt(2)*T^2*pi^(3/2))*(-3*T*fd12posArg+3*T*fd12negArg+2*(muStar-mStar)*fdm12posArg+2*(muStar+mStar)*fdm12negArg)
    p02=sqrt(2*T*mStar)*mStar/(pi^(3/2))*(fdm12posArg+fdm12negArg)
    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
end

using Plots
using LaTeXStrings
plot_font="Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false,legendfontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12)


muRang=collect(0.8:0.005:1.3)
sigmaL=zeros(length(muRang))
omegaL=zeros(length(muRang))
potL=zeros(length(muRang))
#@benchmark get_fields(0.0,0.926)
#Profile.Allocs.clear()
#@time Profile.Allocs.@profile sample_rate=1 get_fields(0.0,0.926)
#PProf.Allocs.pprof(from_c=false)
T=0.01
for i in eachindex(muRang)
    #sigmaL[i]=get_fields(T,muRang[i])[1]
    omegaL[i]=get_fields(T,muRang[i])[2]
    #potL[i]=get_fields(0.0,muRang[i])[3]
end
"""
#plot!(muRang,omegaL,label=L"$\omega$")
#fieldP=plot(muRang,sigmaL,xlabel=L"$\mu$ [GeV]",ylabel=L"Field $[\mathrm{GeV}]$",label=L"$\sigma$",legendtitle=L"Field $[\mathrm{GeV}]$")#,legend=:bottomright)
#plot!(muRang,sigmaL,label="$(T)")


#pPlot=plot(muRang,-potL,xlabel=L"$\mu$ [GeV]",ylabel=L"$p\;\; [\mathrm{GeV}^4]$")
#fieldP=plot(muRang,sigmaL,xlabel=L"$\mu$ [GeV]",ylabel=L"Field $[\mathrm{GeV}]$",label=L"$\sigma$",legendtitle=L"Field $[\mathrm{GeV}]$")#,legend=:bottomright)
#plot(muRang,omegaL,label=L"$\omega$")

#savefig(fieldP,"fieldsT0.png")
#savefig(pPlot,"pPlotT0.png")

"""

eos=WaleckaModel2()

omega_min_fun(-0.0,[0.0,0.98,0.01,eos])
T=0.0
find_omega(0.0,1.0,0.04,eos)
sigmaList=collect(0.04:0.001:0.11)
omgList=find_omega.(Ref(0.0),Ref(0.93),sigmaList)#,Ref(eos))
plot(sigmaList,omgList)
T=0.015
omgList=find_omega.(Ref(T),Ref(0.93),sigmaList)#,Ref(eos))
effPot=effective_potential.(Ref(T),Ref(0.93),sigmaList,omgList)#,Ref(eos))
omgList2=find_omega.(Ref(T),Ref(0.924),sigmaList)#,Ref(eos))
effPot2=effective_potential.(Ref(T),Ref(0.924),sigmaList,omgList2)#,Ref(eos))
omgList3=find_omega.(Ref(T),Ref(0.91),sigmaList)#,Ref(eos))
effPot3=effective_potential.(Ref(T),Ref(0.91),sigmaList,omgList3)#,Ref(eos))
"""
#effpotplot=plot(sigmaList,effPot3,xlabel=L"$\sigma$ [GeV]",ylabel=L"$U_{eff}\;\; [\mathrm{GeV}^4]$",label=L"$0.915$",legendtitle=L"$\mu\;\;[\mathrm{GeV}]$",legend=:top)

#plot!(sigmaList,effPot2,label=L"$0.924$")
#plot!(sigmaList,effPot,label=L"$0.930$")
#savefig(effpotplot,"eff_pot_plotT0.png")

"""
function get_transition_line(step)
    Tmin=0.0
    Tmax=0.0207
    Mumin=0.8
    Mumax=1.1
    Trange=collect(Tmin:step:Tmax)
    murange=collect(Mumin:step:Mumax)
    Tlen=length(Trange)
    mulen=length(murange)
    temp=zeros(Tlen,mulen)
    for t in 1:1:Tlen
        for mu in 1:1:mulen
            fields=get_fields(Trange[t],murange[mu])
            tempsig=fields[1]
            #@show sigma
            #omega=fields[2]
            #pressure=-fields[3] tempsig,tempomg=gapequation3(Trange[t],murange[mu],eos)
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
            if tempList[i] < 0.09
                push!(indList,i)
            end
        end
        indMu=first(indList)
        push!(finalTList,Trange[t])
        push!(finalMuList,murange[indMu])
       end

    return finalTList, finalMuList

end
using Plots
Tlist, mulist =get_transition_line(0.001)
plot(mulist,Tlist,ylabel="T [GeV]",xlabel="μ [GeV]")

using LsqFit

m(t,p) = p[1] .+ p[2] .* t.^2 + p[3] .* t.^4
p0=[0.92,1.,1.]
fit=curve_fit(m,Tlist,mulist,p0)
fit.param
plot(Tlist,mulist)
xrang=range(0.,0.020,length=100)
y= fit.param[1] .+fit.param[2] .* xrang.^2 .+ fit.param[3] .*xrang .^4  
plot!(xrang,y)
using LaTeXStrings
plot_font="Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false,legendfontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12)
"""
#transitionLinePlot=plot(mulist,Tlist,label="Phase transition",xlabel=L"$\mu$ [GeV]",ylabel=L"$T$ [GeV]")
#plot!(y,xrang,label="Fit function")
#savefig(transitionLinePlot,"phase_transition_line.png")
