function interp1(ω,ωres,χ,ωq;extrapolation=true)
    # if extrapolation is not needed, don't do it. it is actually rather expensive
    if extrapolation
        intfM = linear_interpolation(ω,χ,extrapolation_bc=0)
    else
        intfM = scale(interpolate(χ, BSpline(CARSopt.interp)), LinRange(ω[1],ω[end],length(ω)))
    end
    χq = intfM[ωq]

    return χq
end

"""
    findNeighbours(wexp,wlib;extrap=false)
Finds neighbours of each element in wexp in the finer array wlib.
Used for [`fastLinearInt`](@ref).
"""
function findNeighbours(wexp,wlib;extrap=false)
    r=1
    qp = Array{Int32}(undef,length(wexp))
    w = Array{Float64}(undef,length(wexp))
    dw = wlib[2]-wlib[1]
    # auto detect reversal
    if wexp[2]>wexp[1]
        iterator = 1:length(wexp)
    else
        iterator = length(wexp):-1:1
    end
    for i in iterator
        for j in r:length(wlib)
            if wexp[i]<wlib[j]
                r=j
                qp[i]=j
                w[i]=(wlib[j]-wexp[i])/dw
                break
            else
                qp[i]=length(wlib)
                w[i]=0
            end
        end
    end

    if extrap
        qp[qp.==1] .= 2
    end
    return qp,qp.-1,w,1/dw
end

function fastLinearInt(wlib::AbstractArray,data::AbstractArray{T},wexp::AbstractArray{T};extrap=false)  where {T<:AbstractFloat} 
    right,left,w = findNeighbours(wexp,wlib,extrap=extrap)
    interpolated = @avxt @. (w)*data[left]+(1-w)*data[right]
    return interpolated
end

"""
    fastLinearInt(wlib::AbstractArray{R},data::AbstractArray{<:ForwardDiff.Dual{T,V,K}},wexp::AbstractArray{<:AbstractFloat};extrap=false)  where {R,T,V,K}  
Do a fast linear interpolation of a Dual number array `data` from the original grid `wlib` to the destination grid `wexp`.
"""
function fastLinearInt(wlib::AbstractArray{R},data::AbstractArray{<:ForwardDiff.Dual{T,V,K}},wexp::AbstractArray{<:AbstractFloat};extrap=false)  where {R,T,V,K}  
   # in this function data, is always of dual type
   # wlib is dual type if wavenumber shift is fitted
   # wexp is always float here, because dual case is handled in another implementation below
   
   # reinterpret values and partials
   dataE = reinterpret(reshape, V, data)

   # preallocate output array
   out = Array{eltype(data)}(undef,length(wexp))
   # and reinterpret it here to fill with data.
   outE = reinterpret(reshape, V, out)

   if R<:ForwardDiff.Dual
        wlibE = reinterpret(reshape,V,wlib)
        # get index where wlibE has a partial derivative of 1.
        idx = findlast(wlibE[:,1].==1.)
        # remove the partials, we just need the wavenumber array for interpolation
        wlibE = wlibE[1,:]
   else
        wlibE = wlib
   end
   
   # find neighbors and weights
   right,left,w,dw = findNeighbours(wexp,wlibE,extrap=extrap)

   # interpolate values and partials
   @tturbo for i in axes(dataE,1), j in axes(outE,2)
        outE[i,j] = (w[j])*dataE[i,left[j]]+(1-w[j])*dataE[i,right[j]]
   end

   # special case: w depends on partial of wlib
   # so if wlib is a Dual (ie. WavenumberShift is a FitParameter)
   # account for the partial derivative
    if R<:ForwardDiff.Dual
        for j in axes(outE,2)
                outE[idx,j]=dw*dataE[1,left[j]]-dw*dataE[1,right[j]]
        end
    end

    # reference implementation to compare against:
    #    @show outE[:,300]
    #    outO = interp1(wlib,0.,data,wexp;extrapolation=extrap)
    #    outOE = reinterpret(reshape, V, outO)
    #    @show outOE[:,300]
    # out = interp1(wlib,0.,data,wexp;extrapolation=extrap)

   return out   # out has contents of outE!
end

function fastLinearInt(wlib::AbstractArray,data::AbstractArray{<:ForwardDiff.Dual{T,V,K}},wexp::AbstractArray{<:ForwardDiff.Dual};extrap=false)  where {T,V,K}  
    # if wavenumber polynomial is fitted, always resort to this case:
    return interp1(wlib,0.,data,wexp;extrapolation=extrap)
 end
 
