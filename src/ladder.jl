abstract type LadderObject end
abstract type SingleLadderObject <: LadderObject end
# abstract type MultiLadderObject <: LadderObject end

struct SingleLadderSimple{T}
    n_circuits :: Int
    circuitType :: Symbol
    R :: T
    C :: T
    L :: T
end

struct SingleLadderVariable{T}
    n_circuits :: Int
    circuitType :: Symbol
    R::Vector{T}
    C::Vector{T}
    L::Vector{T}
end

"""
    Equation
    M x'' + D x' + S x = G μ
"""
struct SecondOrderSystem{T}
    M :: Matrix{T}
    D :: Matrix{T}
    S :: Matrix{T}
    G :: Matrix{T}
    in_deriv_order :: Int # Order of derivative of input signal
end

"""
    Equation
    [x' ]   [0  I]  [x ]   [0 ] 
    [   ] = [    ]  [  ] + [  ] μ
    [x'']   [P  Q]  [x']   [Bg]
"""
struct FirstOrderSystem{T}
    A :: Matrix{T}
    B :: Matrix{T}
    C :: Matrix{T}
    in_deriv_order :: Int # Order of derivative of input signal
end


function SingleLadder(n_circuits::Int, circuitType::Symbol, Rval::T, Cval::T, Lval::T) where {T <: Real}
    N = n_circuits
    if  circuitType in (:R_con_C, :C_con_R) 
        return SingleLadderSimple(N, circuitType, Rval, Cval, 0 )

    elseif circuitType in (:L_con_C, :C_con_L) 
        return SingleLadderSimple(N, circuitType, 0, Cval, Lval )

    elseif circuitType in (:L_con_R, :R_con_L) 
        return SingleLadderSimple(N, circuitType, Rval, 0, Lval )

    elseif circuitType in (:RL_con_C, :RC_con_L, :LC_con_R)
        return SingleLadderSimple(N, circuitType, Rval, Cval, Lval )   
    end
end

function SingleLadder(n_circuits::Int, circuitType::Symbol, Rvec::Vector{T}, Cvec::Vector{T}, Lvec::Vector{T}) where {T <: Real}
    N = n_circuits
    @assert n_circuits > 0 "At least one circuit is required."
    @assert circuitType in getCircuitTypes() "Type of circuit must exist."
    if  circuitType == :R_con_C ||
        circuitType == :R_con_L ||
        circuitType == :L_con_R || 
        circuitType == :C_con_R || 
        circuitType == :RL_con_C || 
        circuitType == :RC_con_L || 
        circuitType == :LC_con_R 
        
        len_R = length(Rvec)
        @assert len_R == n_circuits "Number of resistors and circuits must be equal: length(R)= $len_R != $n_circuits =n_circuits."
        @assert all(Rvec .>= 0) "Resistor values cannot be negative."
    end
    if  circuitType == :R_con_C ||
        circuitType == :L_con_C ||
        circuitType == :C_con_L || 
        circuitType == :C_con_R || 
        circuitType == :RL_con_C || 
        circuitType == :RC_con_L || 
        circuitType == :LC_con_R 
        
        len_C = length(Cvec)
        @assert len_C == n_circuits "Number of capacitors and circuits must be equal: length(C)= $len_C != $n_circuits =n_circuits."
        @assert all(Cvec .>= 0) "Capacitor values cannot be negative."
    end
    if  circuitType == :R_con_L ||
        circuitType == :L_con_C ||
        circuitType == :C_con_L || 
        circuitType == :L_con_R || 
        circuitType == :RL_con_C || 
        circuitType == :RC_con_L || 
        circuitType == :LC_con_R 
        
        len_L = length(Lvec)
        @assert len_L == n_circuits "Number of inductors and circuits must be equal: length(L)= $len_L != $n_circuits =n_circuits."
        @assert all(Lvec .>= 0) "Capacitor values cannot be negative."
    end

    return SingleLadderVariable(n_circuits::Int, circuitType::Symbol, Rvec, Cvec, Lvec)

end

function getLadderType(ladder)
    return typeof(ladder)
end

function getCircuitTypes()
    return [
        :R_con_C,   # Series: Resistor  ; Shunt: Capacitor
        :L_con_C,   # Series: Inductor  ; Shunt: Capacitor

        :C_con_R,   # Series: Capacitor ; Shunt: Resistor
        :L_con_R,   # Series: Inductor  ; Shunt: Resistor

        :R_con_L,   # Series: Resistor  ; Shunt: Inductor
        :C_con_L,   # Series: Capacitor ; Shunt: Inductor

        :RL_con_C,  # Series: Resistor & Inductor ; Shunt: Capacitor
        :LC_con_R,  # Series: Inductor &  Capacitor; Shunt: Resistor
        :RC_con_L,  # Series: Resistor & Capacitor ; Shunt: Inductor
    ]
end


function _getDerivOrderInput(ladder)
    ctype = ladder.circuitType

    if ctype in (:R_con_C, :L_con_C,:RL_con_C)
        return 0
    elseif ctype in (:C_con_R, :L_con_R, :LC_con_R)
        return 1
    elseif ctype in (:R_con_L, :C_con_L, :RC_con_L)
        return 2
    else
        return -1
    end

end


function secondOrder(ladder :: SingleLadderSimple)
   
    ctype = ladder.circuitType
    N = ladder.n_circuits
    R = ladder.R
    C = ladder.C
    L = ladder.L
    T = typeof(R)

    M, D, S = [], [], []
    in_deriv_order = _getDerivOrderInput(ladder);  # Order of differentiation of input signal   
    # G₁, G₂, G₃ = [], [], []

    # Upper triangular matrix
    M_ut = mapreduce( i-> vcat(ones(T, i),zeros(T, N-i)), hcat, 1:N )
    
    # Finite difference Stencil
    M_fd = LinearAlgebra.diagm(-1 => -1*ones(T, N-1), 0 => ones(T, N))
    
    # RHS input weighing vector
    V_g = vcat(T(1), zeros(T, N-1,1))

    if  ctype == :R_con_C       
        # Series: Resistor | Shunt: Capacitor
        M = zeros(T,1,1)
        D = R*C*M_ut 
        S = M_fd 
        # G₁ = V_g

    elseif ctype == :L_con_C    
        # Series: Inductor | Shunt: Capacitor
        M = L*C*M_ut  # buildUpperTriangular(L, Cinv)
        D = zeros(T,1,1)
        S = M_fd
        # G₁ = V_g

    elseif ctype == :C_con_R    
        # Series: Capacitor | Shunt: Resistor
        M = zeros(T,1,1)
        S = M_ut // (R*C) 
        D = M_fd
#        G₂ = V_g

    elseif ctype == :L_con_R    
        # Series: Inductor  | Shunt: Resistor
        M = (L//R) * M_ut
        D = M_fd
        S = zeros(T,1,1)
#        G₂ = V_g       

    elseif ctype == :R_con_L    
        # Series: Resistor  | Shunt: Inductor
        M = M_fd
        D = (R//L) * M_ut
        S = zeros(T,1,1)
        G₃ = V_g

    elseif ctype == :C_con_L    
        # Series: Capacitor | Shunt: Inductor
        M = M_fd 
        D = zeros(T,1,1)
        S = M_ut // (L*C)
        G₃ = V_g
 
    elseif ctype == :RL_con_C   
        # Series: Resistor & Inductor | Shunt: Capacitor
        M = L*C*M_ut 
        D = R*C*M_ut
        S = M_fd
#        G₁ = V_g
        
    elseif ctype == :LC_con_R   
        # Series: Inductor & Capacitor | Shunt: Resistor
        M = (L//R) * M_ut 
        D = M_fd 
        S = M_ut // (R*C)
#        G₂ = V_g          
                
    elseif ctype == :RC_con_L   
        # Series: Resistor & Capacitor | Shunt: Inductor
        M = M_fd 
        D = (R//L) * M_ut 
        S = M_ut // (L*C) 
#        G₃ = V_g 
    end

    return SecondOrderSystem(M, D, S, V_g, in_deriv_order) #(M, D, S, G₁, G₂, G₃)
end

#=
    M x'' + D x' + S x = G₁ v + G₂ v' + G₃ v''
=#
function secondOrder(ladder :: SingleLadderVariable)
   
    ctype = ladder.circuitType
    N = ladder.n_circuits
    R = ladder.R
    C = ladder.C
    L = ladder.L
    T = typeof(R[end])

    M, D, S = [], [], []
    in_deriv_order = _getDerivOrderInput(ladder);  # Order of differentiation of input signal   
    # G₁, G₂, G₃ = [], [], []

    # Finite difference Stencil
    M_fd = LinearAlgebra.diagm(-1 => -1*ones(T, N-1), 0 => ones(T, N))
    
    # RHS input weighing vector
    V_g = vcat(T(1), zeros(T, N-1,1))

    if  ctype == :R_con_C       
        # Series: Resistor | Shunt: Capacitor
        M = zeros(T,1,1)
        D = upperTriangular(R, C) # buildUpperTriangular(R, Cinv)
        S = M_fd
#        G₁ = V_g

    elseif ctype == :L_con_C    
        # Series: Inductor | Shunt: Capacitor
        M = upperTriangular(L, C) # buildUpperTriangular(L, Cinv)
        D = zeros(T,1,1)
        S = M_fd
#        G₁ = V_g

    elseif ctype == :C_con_R    
        # Series: Capacitor | Shunt: Resistor
        M = zeros(T,1,1)
        S = upperTriangular(C, R, is_inverted=true) # buildUpperTriangular(Cinv, R)
        D = M_fd
#        G₂ = V_g

    elseif ctype == :L_con_R    
        # Series: Inductor  | Shunt: Resistor
        M = upperTriangular(L, (1 ./ R)) # buildUpperTriangular(L, R)
        D = M_fd
        S = zeros(T,1,1)
#        G₂ = V_g     

    elseif ctype == :R_con_L    
        # Series: Resistor  | Shunt: Inductor
        M = M_fd
        D = upperTriangular(R, (1 ./ L)) # buildUpperTriangular(R, L)
        S = zeros(T,1,1)
#        G₃ = V_g

    elseif ctype == :C_con_L    
        # Series: Capacitor | Shunt: Inductor
        # Cinv = inv.(C);
        M = M_fd
        D = zeros(T,1,1)
        S = upperTriangular(C, L, is_inverted=true) # buildUpperTriangular(Cinv, L)
#        G₃ = V_g     

    elseif ctype == :RL_con_C   
        # Series: Resistor & Inductor | Shunt: Capacitor
        # Cinv = inv.(C);
        M = upperTriangular(L, C) # buildUpperTriangular(L, Cinv)
        D = upperTriangular(R, C) # buildUpperTriangular(R, Cinv)
        S = M_fd
#        G₁ = V_g
        
    elseif ctype == :LC_con_R   
        # Series: Inductor &  Capacitor | Shunt: Resistor
        # Cinv = inv.(C);
        M = upperTriangular(L, (1 ./ R)) # buildUpperTriangular(L, R)
        D = M_fd
        S = upperTriangular(C, R, is_inverted=true) # buildUpperTriangular(Cinv, R)
#        G₂ = V_g            
                
    elseif ctype == :RC_con_L   
        # Series: Resistor & Capacitor | Shunt: Inductor
        # Cinv = inv.(C);
        M = M_fd
        D = upperTriangular(R, (1 ./ L)) # buildUpperTriangular(R, L)
        S = upperTriangular(C, L, is_inverted=true) # buildUpperTriangular(Cinv, L)
#        G₃ = V_g
    end

    return SecondOrderSystem(M, D, S, V_g, in_deriv_order) # (M, D, S, G₁, G₂, G₃)
end


function firstOrderSubmatrices(ladder :: SingleLadderSimple)
    ctype = ladder.circuitType
    N = ladder.n_circuits
    R = ladder.R
    C = ladder.C
    L = ladder.L
    T = typeof(R)

    A₁, A₂, B = [], [], [];

    # Upper triangular matrix
    # M_ut = mapreduce( i-> vcat(ones(T, i),zeros(T, N-i)), hcat, 1:N )
    
    # Finite difference Stencil
    # M_fd = LinearAlgebra.diagm(-1 => -1*ones(T, N-1), 0 => ones(T, N))
    
    # Diffusion Stencil
    M_di = LinearAlgebra.diagm(-1 => ones(T, N-1), 0 => -2ones(T, N), 1 => ones(T, N-1))
    M_di[end,end] = T(-1)

    # Inverse Diffusion Matrix
    M_id = T(-1)*_build_InvDiffusion(T,N)

    # RHS input weighing vector
    V_g = vcat(T(1), zeros(T, N-1))

    if ctype == :L_con_C        
        # Series: Inductor | Shunt: Capacitor
        # Cinv = inv.(C);
        A₁ = M_di // (L*C) # buildDiffusion(L, Cinv)
        A₂ = zeros(T, N, N)
        B  = V_g // (L*C)
    
    elseif ctype == :C_con_L    
        # Series: Capacitor | Shunt: Inductor
        # Cinv = inv.(C);
        A₁ = M_id // (L*C) # buildDoubleIntegration(Cinv, L)
        A₂ = zeros(T, N, N)
        B  = ones(T, N)

    elseif  ctype == :RL_con_C   
        # Series: Resistor & Inductor | Shunt: Capacitor
        # Cinv = inv.(C);
        A₁ = M_di // (L*C)  # buildDiffusion(L, Cinv)
        A₂ = -(R//L)*LinearAlgebra.I      # buildReaction(L,R,Cinv)
        B  = V_g // (L*C)   # Cinv[1]/L[1] * vcat(1,zeros(N-1))
        
    elseif ctype == :LC_con_R   
        # Series: Inductor & Capacitor | Shunt: Resistor
        # Cinv = inv.(C);
        A₁ = -LinearAlgebra.I // (L*C) # buildReaction(L,Cinv,R)
        A₂ = (R//L)*M_di # buildDiffusion(L, R)
        B  = (R//L)*V_g  # vcat(R//L, zeros(T,N-1))             
                
    elseif ctype == :RC_con_L   
        # Series: Resistor & Capacitor | Shunt: Inductor
        # Cinv = inv.(C);
        A₁ = M_id // (L*C) # buildDoubleIntegration(Cinv, L)
        A₂ = (R//L) * M_id # buildDoubleIntegration(R, L)
        B  = ones(T, N)
    end

    return (A₁, A₂, B)
end

#=
    A = [0 , I;
         A₁, A₂]
    B
=#
function buildFOSubmatrices(ladder)
    ctype = ladder.circuitType
    N = ladder.n_circuits
    R = ladder.R
    C = ladder.C
    L = ladder.L

    A₁, A₂, B = [], [], [];

    if ctype == :L_con_C        # Series: Inductor ; Shunt: Capacitor
        Cinv = inv.(C);
        A₁ = buildDiffusion(L, Cinv)
        A₂ = zeros(N,N)
        B = Cinv[1]/L[1] * vcat(1,zeros(N-1))
    
    elseif ctype == :C_con_L    # Series: Capacitor ; Shunt: Inductor
        Cinv = inv.(C);
        A₁ = buildDoubleIntegration(Cinv, L)
        A₂ = zeros(N,N)
        B = ones(N)

    elseif  ctype == :RL_con_C   # Series: Resistor & Inductor ; Shunt: Capacitor
        Cinv = inv.(C);
        A₁ = buildDiffusion(L, Cinv)
        A₂ = buildReaction(L,R,Cinv)
        B = Cinv[1]/L[1] * vcat(1,zeros(N-1))
        
    elseif ctype == :LC_con_R   # Series: Inductor &  Capacitor; Shunt: Resistor
        Cinv = inv.(C);
        A₁ = buildReaction(L,Cinv,R)
        A₂ = buildDiffusion(L, R)
        B = R[1]/L[1] * vcat(1,zeros(N-1))             
                
    elseif ctype == :RC_con_L   # Series: Resistor & Capacitor ; Shunt: Inductor
        Cinv = inv.(C);
        A₁ = buildDoubleIntegration(Cinv, L)
        A₂ = buildDoubleIntegration(R, L)
        B = ones(N)
    end

    return (A₁, A₂, B)
end

#=
    x' = A x + B u
    y = C x + D u

    ---------------
    TODO: update this
    ---------------
    If aug = true
    u' = ..
    u'' = v

=#
function firstOrder(ladder :: SingleLadderSimple)
    ctype = ladder.circuitType
    N = ladder.n_circuits
    R = ladder.R
    C = ladder.C
    L = ladder.L

    A, B, Cout = [], [], [];
    T = typeof(R)

    in_deriv_order = _getDerivOrderInput(ladder);  # Order of differentiation of input signal   

    # Diffusion Stencil
    M_di = LinearAlgebra.diagm(-1 => ones(T, N-1), 0 => -2ones(T, N), 1 => ones(T, N-1))
    M_di[end,end] = T(-1)

    # Inverse Diffusion Matrix
    M_id = T(-1)*_build_InvDiffusion(T,N)

    # RHS input weighing vector
    V_g = vcat(T(1), zeros(T, N-1,1))

    if ctype in (:L_con_C, :C_con_L, :RL_con_C, :LC_con_R, :RC_con_L)
        A₁, A₂, B̃ = firstOrderSubmatrices(ladder)
        A = vcat( hcat(zeros(T, N, N), LinearAlgebra.I), hcat(A₁, A₂) )
        B = vcat(zeros(T, N, 1), B̃)
        Cout = hcat(zeros(T, 1, N-1),1,zeros(T, 1, N))

    elseif ctype == :R_con_C
        # Series: Resistor | Shunt: Capacitor 
        # Cinv = inv.(C);
        A = M_di // (R*C) # buildDiffusion(R, Cinv)
        B = V_g  // (R*C) # Cinv[1]/R[1] * vcat(1,zeros(N-1))
        Cout = hcat(zeros(T, 1, N-1), T(1))

    elseif ctype == :C_con_R    
        # Series: Capacitor | Shunt: Resistor 
        # Cinv = inv.(C);
        A = M_id // (R*C) # buildDoubleIntegration(Cinv, R)
        B = ones(T, N, 1)
        Cout = hcat(zeros(T, 1, N-1), T(1))

    elseif ctype == :L_con_R    
        # Series: Inductor | Shunt: Resistor
        A = (R//L) * M_di # buildDiffusion(L, R)
        B = (R//L) * V_g  # R[1]/L[1] * vcat(1,zeros(N-1)) 
        Cout = hcat(zeros(T, 1, N-1), T(1))

    elseif ctype == :R_con_L    
        # Series: Resistor | Shunt: Inductor
        A = (R//L) * M_id # buildDoubleIntegration(R, L)
        B = ones(T, N, 1)
        Cout = hcat(zeros(T, 1, N-1), T(1))

    end

    return FirstOrderSystem(A, B, Cout, in_deriv_order)
end

function buildFirstOrderSystem(ladder; aug=false)
    ctype = ladder.circuitType
    N = ladder.n_circuits
    R = ladder.R
    C = ladder.C
    L = ladder.L

    A, B, Cout = [], [], [];

    if ctype in (:L_con_C, :C_con_L, :RL_con_C, :LC_con_R, :RC_con_L)
        A₁, A₂, B̃ = buildFOSubmatrices(ladder)
        A = vcat(hcat(zeros(N,N), LinearAlgebra.I),hcat(A₁,A₂))
        B = vcat(zeros(N),B̃)
        Cout = hcat(zeros(1,N-1),1,zeros(1,N))

    elseif ctype == :R_con_C
        Cinv = inv.(C);
        A = buildDiffusion(R, Cinv)
        B = Cinv[1]/R[1] * vcat(1,zeros(N-1))
        Cout = hcat(zeros(1,N-1),1)

    elseif ctype == :C_con_R    # Series: Capacitor ; Shunt: Resistor 
        Cinv = inv.(C);
        A = buildDoubleIntegration(Cinv, R)
        B = ones(N)
        Cout = hcat(zeros(1,N-1),1)

    elseif ctype == :L_con_R    # Series: Inductor  ; Shunt: Resistor
        A = buildDiffusion(L, R)
        B = R[1]/L[1] * vcat(1,zeros(N-1)) 
        Cout = hcat(zeros(1,N-1),1)

    elseif ctype == :R_con_L    # Series: Resistor  ; Shunt: Inductor
        A = buildDoubleIntegration(R, L)
        B = ones(N)
        Cout = hcat(zeros(1,N-1),1)

    end

    return (A, B, Cout)
end

