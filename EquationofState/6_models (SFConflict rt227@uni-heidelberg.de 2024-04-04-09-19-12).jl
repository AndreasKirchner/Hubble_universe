using BSON
using BSON: @load
using Flux

struct NeuralNet6M{T,R,S} <:EquationOfState
    model1::T
    model2::T
    model3::T
    model4::R
    model5::R
    model6::R
    normList::S


    function NeuralNet6M()
        #load the models
        @load "6Net_1R.bson" model1
        @load "6Net_2R.bson" model2
        @load "6Net_3R.bson" model3
        @load "6Net_4R.bson" model4
        @load "6Net_5R.bson" model5
        @load "6Net_6R.bson" model6
        @load "hpnorms6Rat.bson" normList
    
        #new{typeof(normList),typeof(model1),typeof(model2),typeof(model3),typeof(model4),typeof(model5),typeof(model6)}(normList,model1,model2,model3,model4,model5,model6)
        #new{typeof(normList),typeof(model1,model2)}(normList,model1,model2,model3,model4,model5,model6)
        new{typeof(model1),typeof(model4),typeof(normList)}(model1,model2,model3,model4,model5,model6,normList)
    end

end




idQCD=IdealQCD()

function pressure(T,mu,x::NeuralNet6M{L}) where {L}
    return x.model1([T,mu])[1]*x.normList[1]*pressure(T,mu,idQCD)
end

function pressure_derivative(T,mu,::Val{1},::Val{0},x::NeuralNet6M{L}) where {L}
    return x.normList[2]*x.model2([T,mu])[1]*pressure_derivative(T,mu,Val(1),Val(0),idQCD)#(x.normList[1]*x.model2([T,mu])[1]*pressure(T,mu,idQCD)^2+pressure(T,mu,NeuralNet6M())*pressure_derivative(T,mu,Val(1),Val(0),idQCD))/pressure(T,mu,idQCD)
end

function pressure_derivative(T,mu,::Val{0},::Val{1},x::NeuralNet6M{L}) where {L}
    return x.normList[3]*x.model3([T,mu])[1]*pressure_derivative(T,mu,Val(0),Val(1),idQCD)#(x.normList[1]*x.model3([T,mu])[1]*pressure(T,mu,idQCD)^2+x.normList[1]*x.model1([T,mu])[1]*pressure(T,mu,idQCD)*pressure_derivative(T,mu,Val(0),Val(1),idQCD))/pressure(T,mu,idQCD)
end

function pressure_derivative(T,mu,::Val{1},::Val{1},x::NeuralNet6M{L}) where {L}
    return x.normList[4]*x.model4([T,mu])[1]*pressure_derivative(T,mu,Val(1),Val(1),idQCD)#(x.normList[1]*x.model4([T,mu])[1]*pressure(T,mu,idQCD)^3-2*x.normList[1]*x.model1([T,mu])[1]*pressure(T,mu,idQCD)*pressure_derivative(T,mu,Val(1),Val(0),idQCD)*pressure_derivative(T,mu,Val(0),Val(1),idQCD)+pressure(T,mu,idQCD)*(pressure_derivative(T,mu,Val(1),Val(0),NeuralNet6M())*pressure_derivative(T,mu,Val(0),Val(1),idQCD)+pressure_derivative(T,mu,Val(1),Val(0),idQCD)+pressure_derivative(T,mu,Val(1),Val(1),idQCD)))/(pressure(T,mu,idQCD)^2)
end

function pressure_derivative(T,mu,::Val{2},::Val{0},x::NeuralNet6M{L}) where {L}
    return x.normList[5]*x.model5([T,mu])[1]*pressure_derivative(T,mu,Val(2),Val(0),idQCD)#1.0
end

function pressure_derivative(T,mu,::Val{0},::Val{2},x::NeuralNet6M{L}) where {L}
    return x.normList[6]*x.model6([T,mu])[1]*pressure_derivative(T,mu,Val(0),Val(2),idQCD)#1.0
end

function thermodynamic(T::N,mu::S,x::NeuralNet6M{L}) where {N,S,L}
    Thermodynamic{promote_type(N,S),2,3}(pressure(T,mu,x),
    (pressure_derivative(T,mu,Val{1}(),Val{0}(),x),pressure_derivative(T,mu,Val{0}(),Val{1}(),x)), 
    (pressure_derivative(T,mu,Val{2}(),Val{0}(),x),pressure_derivative(T,mu,Val{1}(),Val{1}(),x),pressure_derivative(T,mu,Val{0}(),Val{2}(),x) ))
end

