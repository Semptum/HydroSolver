function invert_2(matrix)
    a,c,b,d=matrix
    return [d  -b 
             -c  a]/(a*d-b*c)
end

function +(trX::AffineTransform, trY::AffineTransform)
    A = trX.A + trY.A
    B = trX.B + trY.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function +(trX::AffineTransform, M::AbstractMatrix)
    A = trX.A + M
    B = trX.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function +(M::AbstractMatrix, trX::AffineTransform)
    trX+M
end

function -(M::AbstractMatrix, trX::AffineTransform)
    -trX+M
end

function -(trX::AffineTransform, M::AbstractMatrix)
    trX + -1*M
end

function -(trX::AffineTransform, trY::AffineTransform)
    trX + -1*trY
end

function -(trX::AffineTransform)
    trX *-1
end

function *(trX::AffineTransform, trY::AffineTransform)
    A = trX.A * trY.A
    B = trX.A * trY.B + trX.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function *(M::AbstractMatrix, trX::AffineTransform)
    A = M * trX.A
    B = M * trX.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function *(trX::AffineTransform, M::AbstractMatrix)
    A = trX.A * M
    B = trX.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function *(trX::AffineTransform, X::AbstractVector)
    return trX.A * X + trX.B
end


function *(trX::AffineTransform{MT,VT}, X::Number) where {VT,MT}
    return AffineTransform{MT,VT}(trX.A*X,trX.B*X)
end

function *(X::Number, trX::AffineTransform{MT,VT}) where {VT,MT}
    return AffineTransform{MT,VT}(trX.A*X,trX.B*X)
end


function *(x::SVector{2,T}, y::SVector{2,T}) where T
    x⋅y
end


function *(A::AbstractMatrix{SVector{Dim,T}}, B::AbstractMatrix{SVector{Dim,T}}) where {T,Dim}
    a = split(A)
    b = split(B)
    sum(a[i]*b[i] for i in 1:Dim)
end

function *(A::AbstractMatrix{SVector{Dim,T}}, B::AbstractMatrix{T}) where {T,Dim}
    a = split(A)
    sum(a[i]*B for i in 1:Dim)
end

function *(A::AbstractMatrix{T}, B::AbstractMatrix{SVector{Dim,T}}) where {T,Dim}
    b = split(B)
    sum(A*b[i] for i in 1:Dim)
end

function split(A)
    Dim = length(first(A))
    [getindex.(A,i) for i in 1:Dim]
end

function split(A::AffineTransform)
    Dim = length(first(A.B))
    [AffineTransform(getindex.(A.A,i),getindex.(A.B,i)) for i in 1:Dim]
end



function reunite(VA)
    T=eltype(VA[1])
    Dim=length(VA)
    PointT = SVector{Dim,T}
    base = LinearAlgebra.I(Dim)
    B = [PointT(base[:,k]) for k in 1:Dim]
    sum(VA[i].*[B[i]] for i in 1:Dim)
end