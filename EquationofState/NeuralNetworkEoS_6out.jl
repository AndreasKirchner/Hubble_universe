using BSON
using BSON: @load
using Zygote
using Flux
#using Plots

s(x) = vcat(softplus(x[[1],:]),softplus(x[[2],:]),softplus(x[[3],:]),sigmoid(x[[4],:]),sigmoid(x[[5],:]),sigmoid(x[[6],:]))

struct NeuralNet6{T,S} <:EquationOfState
    #norm1::T
    #norm2::T
    #norm3::T
    #norm4::T
    #norm5::T
    #norm6::T
    normList::T
    model::S

    function NeuralNet6()

        

        #load the model
        @load "whileFit.bson" model#"firstTest2.bson" model#"6out_mod_T45.bson" model #model
        @load "6norms.bson" normList

        new{typeof(normList),typeof(model)}(normList,model)
    end
end

#function NeuralNet(a1)
#    NeuralNet(a1)
#end

#NeuralNet()=NeuralNet(19.0)

#normFac=19

idQCD=IdealQCD()

function pressure(T,mu,x::NeuralNet6{L}) where {L}
    #return x.model([T,mu])[1]*x.normList[1]*(T^4+mu^4)
    return x.model([T,mu])[1]*x.normList[1]*pressure(T,mu,idQCD)
end

#function pressure(T,mu)
#    return pressure(T,mu,NeuralNet6())
#end

function pressure_derivative(T,mu,::Val{1},::Val{0},x::NeuralNet6{L}) where {L}
    #return gradient(pressure,T,mu)[1]
    #return (x.model([T,mu])[2]*x*normList[2]*(T^4+mu^4)^2+4*T^3*pressure(T,mu,x))/(T^4+mu^4)
    #return x.model([T,mu])[2]*x.normList[2]*(T^3+mu^3)
    return x.model([T,mu])[2]*x.normList[2]*pressure_derivative(T,mu,Val(1),Val(0),idQCD)
end

function pressure_derivative(T,mu,::Val{0},::Val{1},x::NeuralNet6{L}) where {L}
    #return x.model([T,mu])[3]*x.normList[3]*(T^3+mu^3)
    return x.model([T,mu])[3]*x.normList[3]*pressure_derivative(T,mu,Val(0),Val(1),idQCD)
end

function pressure_derivative(T,mu,::Val{2},::Val{0},x::NeuralNet6{L}) where {L}
    #return x.model([T,mu])[4]*x.normList[4]*(T^2+mu^2)
    return x.model([T,mu])[4]*x.normList[4]*pressure_derivative(T,mu,Val(2),Val(0),idQCD)
end

function pressure_derivative(T,mu,::Val{1},::Val{1},x::NeuralNet6{L}) where {L}
    #return x.model([T,mu])[5]*x.normList[5]*(T^2+mu^2)
    return x.model([T,mu])[5]*x.normList[5]*pressure_derivative(T,mu,Val(1),Val(1),idQCD)
end

function pressure_derivative(T,mu,::Val{0},::Val{2},x::NeuralNet6{L}) where {L}
    #return x.model([T,mu])[6]*x.normList[6]*(T^2+mu^2)
    return x.model([T,mu])[6]*x.normList[6]*pressure_derivative(T,mu,Val(0),Val(2),idQCD)
end

function thermodynamic(T::N,mu::S,x::NeuralNet6{L}) where {N,S,L}
    Thermodynamic{promote_type(N,S),2,3}(pressure(T,mu,x),
    (pressure_derivative(T,mu,Val{1}(),Val{0}(),x),pressure_derivative(T,mu,Val{0}(),Val{1}(),x)), 
    (pressure_derivative(T,mu,Val{2}(),Val{0}(),x),pressure_derivative(T,mu,Val{1}(),Val{1}(),x),pressure_derivative(T,mu,Val{0}(),Val{2}(),x) ))
end

#=
function pressure_derivative(T,mu,::Val{1},::Val{0},x::NeuralNet{L}) where {L}
    return gradient(pressure,T,mu)[2]
end


function pressure_derivative(T,mu,::Val{2},::Val{0},x::NeuralNet{L}) where {L}
    return hessian(a-> pressure(a[1],a[2]),[T,mu])[1,1]
end


function pressure_derivative(T,mu,::Val{1},::Val{1},x::NeuralNet{L}) where {L}
    return hessian(a-> pressure(a[1],a[2]),[T,mu])[1,2]
end

function pressure_derivative(T,mu,::Val{0},::Val{2},x::NeuralNet{L}) where {L}
    return hessian(a-> pressure(a[1],a[2]),[T,mu])[2,2]
end


function thermodynamic(T::N,mu::S,x::NeuralNet{L}) where {N,S,L}
    Thermodynamic{promote_type(N,S),2,3}(pressure(T,mu,x),
    (pressure_derivative(T,mu,Val{1}(),Val{0}(),x),pressure_derivative(T,mu,Val{0}(),Val{1}(),x)), 
    (pressure_derivative(T,mu,Val{2}(),Val{0}(),x),pressure_derivative(T,mu,Val{1}(),Val{1}(),x),pressure_derivative(T,mu,Val{0}(),Val{2}(),x) ))
end
=#