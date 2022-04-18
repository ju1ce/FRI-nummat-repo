using LinearAlgebra

""" 
Podatkovni tip za pasovno matriko.
Hrani vektor za diagonalo, vektor vektorjev za neničelne vrstice pod diagonalo in vektor vektorjev za neničelne vrstice nad diagonalo

A = PasovnaMatrika([1,1,1],[[2,2],[3]],[[4,4]]) -> 
[1 4 0]
[2 1 4]
[3 2 1]

"""

struct PasovnaMatrika{T} <: AbstractArray{T,2}
        d::Vector{T}  
        s::Vector{Vector{T}}
        z::Vector{Vector{T}}
end

import Base:size,getindex,setindex!

"""
    size(A)

    Vrne velikost pasovne matrike

"""

function size(A::PasovnaMatrika)
    return (size(A.d)[1],size(A.d)[1])
end

"""
    getindex(A,I) ali A[i,j]

    Vrne element na mestu i,j

"""

function getindex(A::PasovnaMatrika, I::Vararg{Int,2})
    
    if I > size(A)
        return 0
    end
    
    i =  I[1] - I[2]

    if i == 0
        return A.d[I[1]]
    end
    if i > 0 && i <= size(A.s)[1]
        return A.s[i][I[1]-i]
    end
    i = -i
    if i > 0 && i <= size(A.z)[1]
        return A.z[i][I[2]-i]
    end
    return 0

end

"""
    setindex(A,v,I) ali A[i,j] = v

    Nastavi element na mestu i, j na v, ce je to mesto v pasu, drugače ne naredi ničesar.

"""

function setindex!(A::PasovnaMatrika, v, I::Vararg{Int,2})
    
    if I > size(A)
        return
    end
    
    i =  I[1] - I[2]

    if i == 0
        A.d[I[1]] = v
    end
    if i > 0 && i <= size(A.s)[1]
        A.s[i][I[1]-i] = v
    end
    i = -i
    if i > 0 && i <= size(A.z)[1]
        A.z[i][I[2]-i] = v
    end
end

""" 
Podatkovni tip za zgornje pasovno matriko.
Hrani vektor za diagonalo in vektor vektorjev za neničelne vrstice nad diagonalo

A = ZgornjePasovnaMatrika([1,1,1],[[4,4]]) -> 
[1 4 0]
[0 1 4]
[0 0 1]

"""

struct ZgornjePasovnaMatrika{T} <: AbstractArray{T,2}
        d::Vector{T}  
        z::Vector{Vector{T}}
end

import Base:size,getindex,setindex!

"""
    size(A)

    Vrne velikost pasovne matrike

"""

function size(A::ZgornjePasovnaMatrika)
    return (size(A.d)[1],size(A.d)[1])
end

"""
    getindex(A,I) ali A[i,j]

    Vrne element na mestu i,j

"""

function getindex(A::ZgornjePasovnaMatrika, I::Vararg{Int,2})
    
    if I > size(A)
        return 0
    end
    
    i =  I[1] - I[2]

    if i == 0
        return A.d[I[1]]
    end
    i = -i
    if i > 0 && i <= size(A.z)[1]
        return A.z[i][I[2]-i]
    end
    return 0

end

"""
    setindex(A,v,I) ali A[i,j] = v

    Nastavi element na mestu i, j na v, ce je to mesto v pasu, drugače ne naredi ničesar.

"""

function setindex!(A::ZgornjePasovnaMatrika, v, I::Vararg{Int,2})
    
    if I > size(A)
        return
    end
    
    i =  I[1] - I[2]

    if i == 0
        A.d[I[1]] = v
    end
    i = -i
    if i > 0 && i <= size(A.z)[1]
        A.z[i][I[2]-i] = v
    end
end

""" 
Podatkovni tip za spodnje pasovno matriko.
Hrani vektor za diagonalo in vektor vektorjev za neničelne vrstice pod diagonalo

A = PasovnaMatrika([1,1,1],[[2,2],[3]]) -> 
[1 0 0]
[2 1 0]
[3 2 1]

"""

struct SpodnjePasovnaMatrika{T} <: AbstractArray{T,2}
        d::Vector{T}  
        s::Vector{Vector{T}}
end

import Base:size,getindex,setindex!

"""
    size(A)

    Vrne velikost pasovne matrike

"""

function size(A::SpodnjePasovnaMatrika)
    return (size(A.d)[1],size(A.d)[1])
end

"""
    getindex(A,I) ali A[i,j]

    Vrne element na mestu i,j

"""

function getindex(A::SpodnjePasovnaMatrika, I::Vararg{Int,2})
    
    if I > size(A)
        return 0
    end
    
    i =  I[1] - I[2]

    if i == 0
        return A.d[I[1]]
    end
    if i > 0 && i <= size(A.s)[1]
        return A.s[i][I[1]-i]
    end
    return 0

end

"""
    setindex(A,v,I) ali A[i,j] = v

    Nastavi element na mestu i, j na v, ce je to mesto v pasu, drugače ne naredi ničesar.

"""

function setindex!(A::SpodnjePasovnaMatrika, v, I::Vararg{Int,2})
    
    if I > size(A)
        return
    end
    
    i =  I[1] - I[2]

    if i == 0
        A.d[I[1]] = v
    end
    if i > 0 && i <= size(A.s)[1]
        A.s[i][I[1]-i] = v
    end
end

import Base:*, \
import LinearAlgebra:lu

"""
    A*v
    Množenje pasovne matrike A s vektorjem v
"""

function *(A::PasovnaMatrika, v::Vector)
    y = zeros(size(A.d))
    
    for i = 1 : size(A.d)[1]
       y[i] = v[i] * A.d[i]
    end
    
    for i = 1 : size(A.s)[1]
       for j = 1 : size(A.s[i])[1]
            y[j+i] += v[j] * A.s[i][j]
        end
    end
    
    for i = 1 : size(A.z)[1]
       for j = 1 : size(A.z[i])[1]
            y[j] += v[j+i] * A.z[i][j]
        end
    end
    
    return y
end

"""
    A*v
    Množenje zgornje pasovne matrike A s vektorjem v
"""

function *(A::ZgornjePasovnaMatrika, v::Vector)
    y = zeros(size(A.d))
    
    for i = 1 : size(A.d)[1]
       y[i] = v[i] * A.d[i]
    end
    
    for i = 1 : size(A.z)[1]
       for j = 1 : size(A.z[i])[1]
            y[j] += v[j+i] * A.z[i][j]
        end
    end
    
    return y
end

"""
    A*v
    Množenje spodnje pasovne matrike A s vektorjem v
"""

function *(A::SpodnjePasovnaMatrika, v::Vector)
    y = zeros(size(A.d))
    
    for i = 1 : size(A.d)[1]
       y[i] = v[i] * A.d[i]
    end
    
    for i = 1 : size(A.s)[1]
       for j = 1 : size(A.s[i])[1]
            y[j+i] += v[j] * A.s[i][j]
        end
    end
    
    return y
end

"""
    L, U = lu(A)
    Izračuna LU razcep, če je matrika diagonalno dominantna.
    L je tipa SpodnjeTrikotnaMatrika, U pa ZgornjeTrikotnaMatrika
"""

function lu(A::PasovnaMatrika)
    
    #test whether matrix is diagonaly dominant
    for i = 1:length(A.d)   
        @assert 2*A[i,i] >= sum(A[i,:])     #matrix is not diagonaly dominant!
    end
    
    n = length(A.d)
    ns = length(A.s)
    nz = length(A.z)
    
    U = PasovnaMatrika(deepcopy(A.d),deepcopy(A.s),deepcopy(A.z))
    #U = copy(A)
    L = SpodnjePasovnaMatrika(ones(length(A.d)),deepcopy(A.s))

    #Calculate L and U through gauss elimination
    for i = 1:n-1
        for j = i+1:i+ns
            if j > n
                break
            end
            l = U[j,i]/U[i,i]
            L[j,i] = l
            U[j,:] -= U[i,:]*l
        end
    end

    U = ZgornjePasovnaMatrika(U.d,U.z)
    
    return L, U
end

"""
    x = A\b
    Izračuna sistem A*x=b
"""

function \(A::PasovnaMatrika, b::Vector)
     
    x = copy(b)
    n = length(b)
    ns = length(A.s)
    nz = length(A.z)

    #gauss elimination
    for i = 1:n-1
        for j = i+1:i+ns
            if j > n
                break
            end
            #row = A[i,:]
            l = A[j,i]/A[i,i]
            A[j,:] -= A[i,:]*l
            x[j] -= x[i]*l
        end
    end
    
    #calculate solution row by row
    for i = n:-1:1
        temp = x[i]
        for j = i+1:n

            if j > nz+i
                break
            end
            temp -= x[j]*A[i,j]

        end
        x[i] = temp/A[i,i]
    end
    
    return x
end

"""
    x = A\b
    Izračuna sistem A*x=b
"""

function \(A::ZgornjePasovnaMatrika, b::Vector)
     
    x = copy(b)
    n = length(b)
    nz = length(A.z)
    
    #calculate solution row by row
    for i = n:-1:1
        temp = x[i]
        for j = i+1:n
            if j > nz+i
                break
            end
            temp -= x[j]*A[i,j]
        end
        x[i] = temp/A[i,i]
    end
    
    return x
end

"""
    x = A\b
    Izračuna sistem A*x=b
"""

function \(A::SpodnjePasovnaMatrika, b::Vector)
     
    x = copy(b)
    n = length(b)
    ns = length(A.s)

    #gauss elimination
    for i = 1:n-1
        for j = i+1:i+ns
            if j > n
                break
            end
            l = A[j,i]/A[i,i]
            A[j,:] -= A[i,:]*l
            x[j] -= x[i]*l
        end
    end
    
    #read the solution
    for i = n:-1:1
        x[i] = x[i]/A[i,i]
    end
    
    return x
end


