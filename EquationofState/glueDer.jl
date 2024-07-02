struct glueDer{T} <:EquationOfState
    a1::T
    a2::T
    a3::T
end

function glueDer(a1,a2,a3)
    glueDer(a1,a2,a3)
end

glueDer()=glueDer(0.166,-0.0556,-0.0212)


function glueDer(T,mu,x::glueDer)
    return tanh((10. *(-x.a1 - x.a2*^(mu,2) - x.a3*^(mu,4) + T)) / x.a1)
end

function glueDer_derivative(T,mu,::Val{1},::Val{0},x::glueDer)
    return (10. *^(sech((10. *(-x.a1 - x.a2*^(mu,2) - x.a3*^(mu,4) +     T))/x.a1),2))/x.a1
end

function glueDer_derivative(T,mu,::Val{0},::Val{1},x::glueDer)
    return (-20. *(x.a2*mu + 2. *x.a3*^(mu,3))*^(sech((10. *(-x.a1 -     x.a2*^(mu,2) - x.a3*^(mu,4) + T))/x.a1),2))/x.a1
end

function glueDer_derivative(T,mu,::Val{2},::Val{0},x::glueDer)
    return (-200. *^(sech((10. *(-x.a1 - x.a2*^(mu,2) - x.a3*^(mu,4) +     T))/x.a1),2)*tanh((10. *(-x.a1 - x.a2*^(mu,2) - x.a3*^(mu,4) +     T))/x.a1))/^(x.a1,2)
end

function glueDer_derivative(T,mu,::Val{1},::Val{1},x::glueDer)
    return (400. *(x.a2*mu + 2. *x.a3*^(mu,3))*^(sech((10. *(-x.a1 -     x.a2*^(mu,2) - x.a3*^(mu,4) + T))/x.a1),2)*tanh((10. *(-x.a1 -     x.a2*^(mu,2) - x.a3*^(mu,4) + T))/x.a1))/^(x.a1,2)
end

function glueDer_derivative(T,mu,::Val{0},::Val{2},x::glueDer)
    return (-20. *^(sech((10. *(-x.a1 - x.a2*^(mu,2) - x.a3*^(mu,4) +  T))/x.a1),2)*(x.a1*(1. *x.a2 + 6. *x.a3*^(mu,2)) + (40.     *^(x.a2,2)*^(mu,2) + 160. *x.a2*x.a3*^(mu,4) + 160.     *^(x.a3,2)*^(mu,6))*tanh((10. *(-x.a1 - x.a2*^(mu,2) - x.a3*^(mu,4) + T))/x.a1)))/^(x.a1,2)
end

function thermodynamic(T::N,mu::S,x::glueDer) where {N,S}
    Thermodynamic{promote_type(N,S),2,3}(glueDer(T,mu,x),
    (glueDer_derivative(T,mu,Val{1}(),Val{0}(),x),glueDer_derivative(T,mu,Val{0}(),Val{1}(),x)),
    (glueDer_derivative(T,mu,Val{2}(),Val{0}(),x),glueDer_derivative(T,mu,Val{1}(),Val{1}(),x),glueDer_derivative(T,mu,Val{0}(),Val{2}(),x)))
end