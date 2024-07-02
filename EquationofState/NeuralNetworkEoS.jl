using BSON
using BSON: @load
using Zygote
using Flux
#using Plots

struct NeuralNet{T,S} <:EquationOfState
    normList::T
    model::S

    function NeuralNet()

        #a1=19

        #load the model
        #@load "OneOutDerTraining.bson" model2
        #@load "hpnorms.bson" normList
        @load "pressure_fitLongCorErr.json" model
        @load "normListP.json" normList

        new{typeof(normList),typeof(model)}(normList,model)
    end
end

#function NeuralNet(a1)
#    NeuralNet(a1)
#end

#NeuralNet()=NeuralNet(19.0)

#normFac=19

idQCD=IdealQCD()

function pressure(T,mu,x::NeuralNet{L}) where {L}
    return x.model([T,mu])[1]*x.normList[1]#*pressure(T,mu,idQCD)#*x.a1*(T^4+mu^4)
end

function pressure(T,mu)
    return pressure(T,mu,NeuralNet())
end

function pressure_derivative(T,mu,::Val{1},::Val{0},x::NeuralNet{L}) where {L}
    return max(0.0,gradient(pressure,T,mu)[1])
end

function pressure_derivative(T,mu,::Val{0},::Val{1},x::NeuralNet{L}) where {L}
    return max(0.0,gradient(pressure,T,mu)[2])
end


function pressure_derivative(T,mu,::Val{2},::Val{0},x::NeuralNet{L}) where {L}
    return max(0.0,hessian(a-> pressure(a[1],a[2]),[T,mu])[1,1])
end


function pressure_derivative(T,mu,::Val{1},::Val{1},x::NeuralNet{L}) where {L}
    return max(0.0,hessian(a-> pressure(a[1],a[2]),[T,mu])[1,2])
end

function pressure_derivative(T,mu,::Val{0},::Val{2},x::NeuralNet{L}) where {L}
    return max(0.0,hessian(a-> pressure(a[1],a[2]),[T,mu])[2,2])
end

function thermodynamic(T::N,mu::S,x::NeuralNet{L}) where {N,S,L}
    Thermodynamic{promote_type(N,S),2,3}(pressure(T,mu,x),
    (pressure_derivative(T,mu,Val{1}(),Val{0}(),x),pressure_derivative(T,mu,Val{0}(),Val{1}(),x)), 
    (pressure_derivative(T,mu,Val{2}(),Val{0}(),x),pressure_derivative(T,mu,Val{1}(),Val{1}(),x),pressure_derivative(T,mu,Val{0}(),Val{2}(),x) ))
end