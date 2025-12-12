


"""
    ̇x = A x + B u
"""
mutable struct FirstOrder
    A;
    B;
end


"""
    M ̈x + D ̇x + S x = G u
"""
mutable struct SecondOrder
    M;
    D;
    S;
    G;
end


# This is NEW ------------- MUST BE TESTED
#=
function create_sos(Y::Vector{T1}, X1::Vector{T2}, X2::Vector{T3}, ident)

    N = length(R)

    tri1_m = buildUpperTriangular(X1, Y)
    tri2_m = buildUpperTriangular(X2, Y)
    diff_m = diagm(-1 => -1*ones(N-1), 0 => ones(N))

    G = vcat(1,zeros(N-1))


    if ident == :RLC
        # Connector: Capacitor
        M,D,S,G = (tri2_m, tri1_m, diff_m, G)
    else if ident == :LCR
        # Connector: Resistor
        M,D,S,G = (tri1_m,diff_m,tri2_m,G)
    else if ident == :CRL
        # Connector: Inductor
        M,D,S,G = (diff_m,tri2_m,tri1_m,G)
    end

    return (M,D,S,G)
end
=#

function createSecondOrder(R::Vector{T1},L::Vector{T2},C::Vector{T3}, circuit) where {T1 <: Real, T2 <: Real, T3 <: Real}

    N = length(R)

    if circuit == :RLC
        Cinv = inv.(C)

        M = buildUpperTriangular(L, Cinv)
        D = buildUpperTriangular(R, Cinv)
        S = diagm(-1 => -1*ones(N-1), 0 => ones(N))
    end

    if circuit == :LCR
        Cinv = inv.(C)

        M = buildUpperTriangular(L, R)
        D = diagm(-1 => -1*ones(N-1), 0 => ones(N))
        S = buildUpperTriangular(Cinv, R)
    end

    if circuit == :CRL
        Cinv = inv.(C)

        M = diagm(-1 => -1*ones(N-1), 0 => ones(N))
        D = buildUpperTriangular(R, L)
        S = buildUpperTriangular(Cinv, L)
    end

    G = vcat(1,zeros(N-1))

    return (M,D,S,G)
end


function createSubFirstOrder(R::Vector{T1},L::Vector{T2},C::Vector{T3}, circuit) where {T1 <: Real, T2 <: Real, T3 <: Real}

    N = length(R)
    Cinv = inv.(C)

    if circuit == :RLC        
        A1 = buildDiffusion(L, Cinv)
        A2 = buildReaction(L,R,Cinv)
        B = Cinv[1]/L[1] * vcat(1,zeros(N-1))
    elseif circuit == :LCR   
        A1 = buildReaction(L,Cinv,R)
        A2 = buildDiffusion(L, R)
        B = R[1]/L[1] * vcat(1,zeros(N-1))
    elseif circuit == :CRL
        A1 = buildDoubleIntegration(Cinv, L)
        A2 = buildDoubleIntegration(R, L)
        B = ones(N)
    end

    return (A1,A2,B)
end

function createFirstOrder_new(R::Vector{T1},L::Vector{T2},C::Vector{T3}, circuit) where {T1 <: Real, T2 <: Real, T3 <: Real}
    N = length(R)
    (A1,A2,B) = createSubFirstOrder(R,L,C, circuit)

    Afo = vcat(hcat(zeros(N,N), I),hcat(A1,A2))
    Bfo = vcat(zeros(N),B)

    out = hcat(zeros(1,N-1),1,zeros(1,N))

    return (Afo,Bfo,out)
end



function createFirstOrder(R::Vector{T1},L::Vector{T2},C::Vector{T3}, circuit) where {T1 <: Real, T2 <: Real, T3 <: Real}

    N = length(R)

    if circuit == :RLC
        Cinv = inv.(C)
        A1 = buildDiffusion(L, Cinv)
        A2 = buildReaction(L,R,Cinv)
        A = vcat(hcat(zeros(N,N),I(N)),hcat(A1, A2))
        B = Cinv[1]/L[1] * vcat(zeros(N),1,zeros(N-1))
    elseif circuit == :LCR   
        Cinv = inv.(C)
        A1 = buildReaction(L,Cinv,R)
        A2 = buildDiffusion(L, R)
        A = vcat(hcat(zeros(N,N),I(N)),hcat(A1, A2))
        B = R[1]/L[1] * vcat(zeros(N),1,zeros(N-1))
    elseif circuit == :CRL
        Cinv = inv.(C)
        A1 = buildDoubleIntegration(Cinv, L)
        A2 = buildDoubleIntegration(R, L)
        A = vcat(hcat(zeros(N,N),I(N)),hcat(A1, A2))

        #dm = diagm(-1 => -1*ones(N-1), 0 => ones(N))
        #B = inv(dm) * vcat(1,zeros(N-1))
        B = vcat(zeros(N),ones(N))
    end

    out = hcat(zeros(1,N-1),1,zeros(1,N))

    return (A,B,out)
end