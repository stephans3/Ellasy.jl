
"""
Builds matrix of shape
[ 1  0  0  0... ]
[-1  1  0  0... ]
[ 0 -1  1  0... ]
[ 0  0 -1  1 ... ]
"""
function buildFiniteDiff(dim :: T) where T <: Integer
    M = LinearAlgebra.diagm(ones(T,dim))    
    for i = 2 : dim
        M[i,i-1] = T(-1)
    end
    return M
end



"""
    Builds matrix of shape

    [1 1 1 ... ]
    [0 1 1 ... ]
    [0 0 1 ... ]
"""

function buildUpperTriangular(dim :: T where T <: Integer)
    M = zeros(typeof(dim), dim,dim);
    for i = 1 : dim, j = i : dim
        if j >= i 
            M[i,j] = 1
        end
    end
    return M
end

"""
    Builds matrix of shape

    [n[1]/d[1] n[1]/d[2] n[1]/d[3] ... ]
    [          n[2]/d[2] n[2]/d[3] ... ]
    [                    n[3]/d[3] ... ]
"""
function buildUpperTriangular(num :: Vector{T1}, den :: Vector{T2} ) where {T1 <: Real, T2 <: Real}
    N = length(num) 

    if N != length(den)
        error("Both vectors must have the same length! Length Num=", N, " Length Den=", length(den))
    end

    M = zeros(N, N);
    for i = 1 : N, j = i : N
        if j >= i 
            M[i,j] = num[i] / den[j]
        end
    end
    return M
end


# Calculates -∫ dx
function buildDoubleIntegration(num :: Vector{T1}, den :: Vector{T2} ) where {T1 <: Real, T2 <: Real}
    N = length(num) 

    if N != length(den)
        error("Both vectors must have the same length! Length Num=", N, " Length Den=", length(den))
    end

    M = zeros(N, N);
    for i = 1 : N, j = 1 : N
        M[i,j] = -sum(num[1:min(i,j)]) / den[j]    
    end
    return M
end



"""
    nocon: not connector elements
    con: connector elements
"""
function buildDiffusion(nocon :: Vector{T1}, con :: Vector{T2} ) where {T1 <: Real, T2 <: Real}
    N = length(nocon) 

    if N != length(con)
        error("Both vectors must have the same length! Length nocon=", N, " Length con=", length(con))
    end

    M = zeros(N, N);
   
    for i = 1 : N-1   
        M[i,i] = -(con[i]/nocon[i] + con[i]/nocon[i+1])
        M[i,i+1] = con[i]/nocon[i+1]
    
        if i > 1
            M[i,i-1] = con[i] / nocon[i]
        end

    end

    M[N,N-1] = con[N] / (nocon[N])
    M[N,N] = -con[N]/ (nocon[N])

    return M
end

function buildReaction(nocon_a :: Vector{T1}, nocon_b :: Vector{T2}, con :: Vector{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}
    N = length(nocon_a) 

    if N != length(con) || N != length(nocon_b) 
        error("Both vectors must have the same length! Length nocon_a=", N, " Length con=", length(nocon_b), " Length con=", length(con))
    end

    M = zeros(N, N);
   
    for i = 1 : N-1    
        M[i,i] = -nocon_b[i] / nocon_a[i]
        for j = i+1 : N
            M[i,j] = -(nocon_b[i]*con[i] / (nocon_a[i]*con[j]) - nocon_b[i+1]*con[i]/(nocon_a[i+1]*con[j]))
        end
    end
 
    M[N,N] = -nocon_b[N] / nocon_a[N]

    return M
end