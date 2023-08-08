#### FUNCTIONS TO HANDLE CONVOLUTIONS in CARS.JL

# Convolutions are about the most expensive functions for the library based approach,
# because for every pump region, there are three spectral components to convolve
# and we need to account for the instrumental function.
# For a dual pump spectrum, these are 7 convolutions + cross-terms if multiple species are inside the same pumping regions
# which have to be done on fairly large arrays. Thus, getting the convolution fast is really the key for fast fitting.
# 
# One further important aspect is that when using Automatic Differentiation during the fit
# instead of finite differences, JuMP introduces a Dual Number type which transports not only
# the value but also the partial derivatives of the decision variables through the function.
# This is much better than finite differences, because it requires less function calls, but
# all functions within the fit have to handle T<:Real, which includes Dual Numbers.
# Unfortunately, DSP.conv or FFTs do not support Dual Numbers out of the box, so we have to 
# use some trickery and handle the partials ourselves.
# This is what dualConv does. Essentially, it unwrapps the values and partials and convolves them
# just if they were ordinary Floats. For the underlying convolution, you can choose your favorite
# algorithm in conv() directly below this comment. Make sure, that it allows to return the same
# array size as input a. So for typical applications, the syntax is dualConv(longVector,shortKernel)
# and the returned vector is the convolution longVector ⋆ shortKernel but with length(longVector)
# The default behavior is to return the same size!

# choose your favorite convolution algorithm here

∗(a,b) = conv(a,b)

function conv(a::Vector{T}, b::Vector{T};same::Bool=true) where {T}
    # check if vector b and vector a are non-zero.
    # sufficient to check the center values
    # this is very helpful for convolving dual numbers
    # as very often, partials are actually zero, so no convolution is necessary
    center_a = Int(floor(length(a)/2))
    center_b = Int(floor(length(b)/2))
   
    if a[center_a] == 0. || b[center_b] == 0.
        return zeros(T,length(a))
    end

    #return dspconv(a,b;same=true)
    return dconv_edge(a,b)
end

function dspconv(a::Vector{T}, b::Vector{T};same::Bool=true) where {T}
    # the parameter same is similar to MATLABs behavior
    # where the size returned is the size of a. for some applications
    # (and certainly this one), this is very beneficial
    same ? retlength = length(a) : retlength = length(a)+length(b)-1

    # do your convolution
    c = DSP.conv(a,b)

    # crop array to correct size
    if same
        # return only the central part
        idx1 = Int(ceil(length(b)/2))
        idx2 = Int(floor(length(b)/2))
        return c[idx1:end-idx2]
    end
    return c
end

# FUNCTION THAT DOES FAST DIRECT CONVOLUTION USING LoopVectorization
# Returns same size array, but has edge effects. should be no problem
# for this particular case because simulated spectra are always broader than experiments
# (in terms of wavenumberrange)
function dconv_edge(signal::Vector{T},kern::Vector{T}) where {T}
    M = length(kern)
    N = length(signal)
    m = ceil(Int,M/2) # this has to be ceil in case the length of the kernel is odd.
    out = zeros(T,N)
    @avxt for j in M:N
        tmp = zero(T)
        for i in 1:M
            tmp += signal[j-i+1] * kern[i]
        end
        out[j-m+1]=tmp
    end
    return out
end


function dualConv(signal::Vector{T}, a::Vector{T}) where {T<:AbstractFloat}
    # if the input is Float64, just use DSP.conv. No need for fancy stuff.
    return conv(signal,a)
end

function dualConv(signal::Vector{V},a::AbstractVector{<:ForwardDiff.Dual{T,V,K}}) where {T,V,K}
    # convolve a vector of type V with a Dual{T,V,K}
    #Sre   = reinterpret(reshape, V, signal)
    Xre   = reinterpret(reshape, V, a)
    # convolve the signal with values and partials
    # this assumes d/dt(f⋆g) = df/dt ⋆ g
    # Note: using dualConv here allows also to convolve nested Duals,
    # which is required when the Hessian matrix is evaluated, not only the gradient
    values = dualConv(signal, Xre[1,:]) # values
    partials = [dualConv(signal, Xre[i,:]) for i in 2:size(Xre,1)] # partials
    ptl = [ForwardDiff.Partials{K,V}(NTuple{K,V}([partials[j][i] for j in 1:length(partials)])) for i in 1:length(values)]
    # return a dual type with the same type as the input - as if nothing happened.
    return [ForwardDiff.Dual{T,V,K}(values[i], ptl[i]) for i in 1:length(values)]
end

function dualConv(signal::AbstractVector{<:ForwardDiff.Dual{T,V,K}},a::AbstractVector{<:ForwardDiff.Dual{T,V,K}}) where {T,V,K}
    # function to convolve to Dual Numbers
    Sre   = reinterpret(reshape, V, signal)
    Xre   = reinterpret(reshape, V, a)
    # Note: using dualConv here allows also to convolve nested Duals,
    # which is required when the Hessian matrix is evaluated, not only the gradient
    values = dualConv(Sre[1,:], Xre[1,:]) # values
    # these are a lot of convolutions
    # at the moment, these are two for each variable that is used for the gradient
    partials = [dualConv(Sre[1,:], Xre[i,:]) + dualConv(Sre[i,:],Xre[1,:]) for i in 2:size(Xre,1)] # partials
    ptl = [ForwardDiff.Partials{K,V}(NTuple{K,V}([partials[j][i] for j in 1:length(partials)])) for i in 1:length(values)]
    return [ForwardDiff.Dual{T,V,K}(values[i], ptl[i]) for i in 1:length(values)]
end