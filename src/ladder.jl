

struct SingleLadder{T}
    n_circuits :: Int
    circuitType :: Symbol
    R::Vector{T}
    C::Vector{T}
    L::Vector{T}
    function SingleLadder(
        n_circuits::Int=1,
        circuitType::Symbol=:R_con_C,
        R::Vector{T}=ones(n_circuits),
        C::Vector{T}=ones(n_circuits),    
        L::Vector{T}=zeros(n_circuits),
    ) where {T}
    @assert n_circuits > 0 "At least one circuit is required."
    @assert circuitType in getCircuitTypes() "Type of circuit must exist."

    if  circuitType == :R_con_C ||
        circuitType == :R_con_L ||
        circuitType == :L_con_R || 
        circuitType == :C_con_R || 
        circuitType == :RL_con_C || 
        circuitType == :RC_con_L || 
        circuitType == :LC_con_R 
        
        len_R = length(R)
        @assert len_R == n_circuits "Number of resistors and circuits must be equal: length(R)= $len_R != $n_circuits =n_circuits."
        @assert all(R .>= 0) "Resistor values cannot be negative."
    end
    if  circuitType == :R_con_C ||
        circuitType == :L_con_C ||
        circuitType == :C_con_L || 
        circuitType == :C_con_R || 
        circuitType == :RL_con_C || 
        circuitType == :RC_con_L || 
        circuitType == :LC_con_R 
        
        len_C = length(C)
        @assert len_C == n_circuits "Number of capacitors and circuits must be equal: length(C)= $len_C != $n_circuits =n_circuits."
        @assert all(C .>= 0) "Capacitor values cannot be negative."
    end
    if  circuitType == :R_con_L ||
        circuitType == :L_con_C ||
        circuitType == :C_con_L || 
        circuitType == :L_con_R || 
        circuitType == :RL_con_C || 
        circuitType == :RC_con_L || 
        circuitType == :LC_con_R 
        
        len_L = length(L)
        @assert len_L == n_circuits "Number of inductors and circuits must be equal: length(L)= $len_L != $n_circuits =n_circuits."
        @assert all(L .>= 0) "Capacitor values cannot be negative."
    end

    return new{T}(n_circuits, circuitType, R, C, L)
    end
end

function SingleLadder(n_circuits::Int, circuitType::Symbol)
    N = n_circuits
    if  circuitType in (:R_con_C, :C_con_R) 
        return SingleLadder(N, circuitType, ones(N), ones(N), zeros(N) )

    elseif circuitType in (:L_con_C, :C_con_L) 
        return SingleLadder(N, circuitType, zeros(N), ones(N), ones(N) )

    elseif circuitType in (:L_con_R || :R_con_L) 
        return SingleLadder(N, circuitType, ones(N), zeros(N), ones(N) )

    elseif circuitType in (:RL_con_C, :RC_con_L, :LC_con_R)
        return SingleLadder(N, circuitType, ones(N), ones(N), ones(N) )   
    end

end


function getCircuitTypes()
    circuitTypes = [
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
    return circuitTypes
end

#=
    M x'' + D x' + S x = G₁ v + G₂ v' G₃ v''
=#
function buildSecondOrderSystem(ladder)
   
    ctype = ladder.circuitType
    N = ladder.n_circuits
    R = ladder.R
    C = ladder.C
    L = ladder.L

    M, D, S = [], [], []
    G₁, G₂, G₃ = [], [], []


    if  ctype == :R_con_C       # Series: Resistor  ; Shunt: Capacitor
        Cinv = inv.(C);
        D = buildUpperTriangular(R, Cinv)
        S = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        G₁ = vcat(1, zeros(N-1))

    elseif ctype == :L_con_C    # Series: Inductor  ; Shunt: Capacitor
        Cinv = inv.(C);
        M = buildUpperTriangular(L, Cinv)
        S = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        G₁ = vcat(1, zeros(N-1))

    elseif ctype == :C_con_R    # Series: Capacitor ; Shunt: Resistor
        Cinv = inv.(C);
        S = buildUpperTriangular(Cinv, R)
        D = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        G₂ = vcat(1, zeros(N-1))     

    elseif ctype == :L_con_R    # Series: Inductor  ; Shunt: Resistor
        M = buildUpperTriangular(L, R)
        D = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        G₂ = vcat(1, zeros(N-1))       

    elseif ctype == :R_con_L    # Series: Resistor  ; Shunt: Inductor
        M = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        D = buildUpperTriangular(R, L)
        G₃ = vcat(1, zeros(N-1)) 

    elseif ctype == :C_con_L    # Series: Capacitor ; Shunt: Inductor
        Cinv = inv.(C);
        M = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        S = buildUpperTriangular(Cinv, L)
        G₃ = vcat(1, zeros(N-1))         

    elseif ctype == :RL_con_C   # Series: Resistor & Inductor ; Shunt: Capacitor
        Cinv = inv.(C);
        M = buildUpperTriangular(L, Cinv)
        D = buildUpperTriangular(R, Cinv)
        S = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        G₁ = vcat(1, zeros(N-1))
        
    elseif ctype == :LC_con_R   # Series: Inductor &  Capacitor; Shunt: Resistor
        Cinv = inv.(C);
        M = buildUpperTriangular(L, R)
        D = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        S = buildUpperTriangular(Cinv, R)
        G₂ = vcat(1, zeros(N-1))                
                
    elseif ctype == :RC_con_L   # Series: Resistor & Capacitor ; Shunt: Inductor
        Cinv = inv.(C);
        M = LinearAlgebra.diagm(-1 => -1*ones(N-1), 0 => ones(N))
        D = buildUpperTriangular(R, L)
        S = buildUpperTriangular(Cinv, L)
        G₃ = vcat(1, zeros(N-1)) 
    end

    return (M, D, S, G₁, G₂, G₃)
end