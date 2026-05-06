module PythonASAP

using ASAP
using PythonCall
using Random
using TypeUtils
using TypeUtils: @public
using Base: @propagate_inbounds

@public(
    build,
    inv_R_Rt,
    inv_Rt_R,
    learn,
    mul!,
    mul,
)

struct FastNumPyArray{T,N,M} <: AbstractArray{T,N}
    ptr::Ptr{T}
    size::Dims{N}
    length::Int
    handle::Py
end
FastNumPyArray(A::FastNumPyArray) = A

function FastNumPyArray(A::PyArray{T,N,M,L,R}) where {T,N,M,L,R}
    R === T || error("parameter `R` must be the same as element type `T`")
    isbitstype(T) || error("element type `T` is not a \"plain data\" type")
    len = length(A)
    if N == 0
        return FastNumPyArray{T,0,M}(A.ptr, (), 1, A.handle)
    else
        if A.strides[1] == sizeof(T)
            column_major = true
            for d in 2:N
                if A.strides[d] != A.size[d-1]*A.strides[d-1]
                    column_major = false
                    break
                end
            end
            if column_major
                return FastNumPyArray{T,N,M}(A.ptr, A.size, len, A.handle)
            end
        elseif A.strides[N] == sizeof(T)
            row_major = true
            for d in 1:N-1
                if A.strides[d] != A.size[d+1]*A.strides[d+1]
                    row_major = false
                    break
                end
            end
            if row_major
                return FastNumPyArray{T,N,M}(A.ptr, reverse(A.size), len, A.handle)
            end
        end
    end
    error("NumPy array is neither in column-major, nor in row-major contiguous storage order")
end

# Abstract array API for fast NumPy arrays.
Base.length(A::FastNumPyArray) = getfield(A, :length)
Base.size(A::FastNumPyArray) = getfield(A, :size)
Base.IndexStyle(::Type{<:FastNumPyArray}) = IndexLinear()
Base.pointer(A::FastNumPyArray) = getfield(A, :ptr)
Base.unsafe_convert(::Type{Ptr{T}}, A::FastNumPyArray{T,N}) where {T,N} = pointer(A)
Base.elsize(::Type{FastNumPyArray{T,N}}) where {T,N} = sizeof(T)
@propagate_inbounds function Base.getindex(A::FastNumPyArray{T,N}, i::Int) where {T,N}
    return unsafe_load(pointer(A), i)
end
@propagate_inbounds function Base.setindex!(A::FastNumPyArray{T,N,true}, x, i::Int) where {T,N}
    unsafe_store!(pointer(A), x, i)
    return A
end

"""
    B = PythonASAP.fast_array(A)

Return a view `B` for fast access to the elements of `A`. `B` has linear indexing style with
1-based indices and contiguous elements in memory.

"""
fast_array(A::Array) = A
if isdefined(Base, :Memory)
    fast_array(A::Memory) = A
end
fast_array(A::PyArray) = FastNumPyArray(A)
function fast_array(A::DenseArray)
    IndexStyle(A) === IndexLinear() || throw(ArgumentError(
        "array index style is not linear indexing"))
    firstindex(A) == 1 || throw(ArgumentError("first index of array is not 1"))
    return A
end

"""
    B = PythonASAP.inv_R_Rt(R)

Return the model of a symmetric positive-definite matrix given by `inv(R*R')` with `R` an
ASAP sparse factor.

"""
inv_R_Rt(R::ASAP.SparseFactor) = inv(R*R')

"""
    B = PythonASAP.inv_Rt_R(R)

Return the model of a symmetric positive-definite matrix given by `inv(R'*R)` with `R` an
ASAP sparse factor.

"""
inv_Rt_R(R::ASAP.SparseFactor) = inv(R'*R)

build_format(::Type{F}) where {F <: ASAP.Format} = F
function build_format(fmt::AbstractString)
    if fmt == "RowWiseLower"
        return ASAP.RowWiseLower
    elseif fmt == "RowWiseUpper"
        return ASAP.RowWiseUpper
    elseif fmt == "ColumnWiseLower"
        return ASAP.ColumnWiseLower
    elseif fmt == "ColumnWiseUpper"
        return ASAP.ColumnWiseUpper
    else
        throw(ArgumentError("$(repr(fmt)) is not a known ASAP matrix format"))
    end
end

function build_perm(perm, n::Integer)
    if perm isa Nothing || perm == "none"
        return Base.OneTo{Int}(n)
    elseif perm == "random"
        return Random.randperm(n)
    elseif perm ∈ ("multiscale", "FRiM")
        return Symbol(perm)
    else
        throw(ArgumentError("argument `perm = $(repr(perm))` not supported"))
    end
end

function build(fmt::AbstractString, msk::AbstractArray{Bool}, perm, args...; kwds...)
    return build(build_format(fmt), msk, perm, args...; kwds...)
end

function build(::Type{F}, msk::AbstractArray{Bool}, perm, args...;
               kwds...) where {F<:ASAP.Format}
    return ASAP.SparseFactor{F,Bool}(fast_array(msk), build_perm(perm, length(msk)),
                                     args...; kwds...)
end

"""
    PythonASAP.learn(mdl, ref) -> res

Learn the coefficients of the ASAP model `mdl` from the the reference matrix `ref` and
return a new model `res` similar to `mdl` but whose coefficients are independent from those
of `mdl`.

"""
function learn(mdl::ASAP.Gram, ref)
    return ASAP.learn(mdl, ref)
end

"""
    PythonASAP.learn!(mdl, ref) -> mdl

Learn the coefficients of the ASAP model `mdl` from the the reference matrix `ref`. The
coefficients of `mdl` are overwritten by the learned ones.

"""
function learn!(mdl::ASAP.Gram, ref)
    return ASAP.learn(mdl, ref)
end

"""
"""
function mul(op, src::AbstractArray)
    return ASAP.mul(op, fast_array(src))
end

"""
"""
function mul!(dst::AbstractArray, op, src::AbstractArray)
    ASAP.mul!(fast_array(dst), op, fast_array(src))
    return dst
end

end
