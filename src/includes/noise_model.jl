## this function covers the noise models that can be used
"""
    noiseModel!(type::Symbol,x...)
Abstraction to extract a noiseModel with various methods. Only implemented at the moment: `:argon` which calls [`argonModel!`](@ref)
"""
function noiseModel!(type::Symbol,x...)
    if type == :none
        return 1
    end

    if type == :argon
        return argonModel!(x...)
    end

end

"""
    argonModel!(preprocessed,nr,roi)
Extract a weight based on Stokes noise estimated from measurements in Argon.
"""
function argonModel!(preprocessed,nr,roi)
    # extract the non resonant signal in the region of interest:
    nr = nr.I[roi[1]:roi[end],:]
    # take sqrt if specified in options
    if CARSopt.fittype == :sqrt
        @. nr = sqrt(nr)
    end
    # normalize each spectrum to its area
    nr = nr ./ sum(nr,dims=1)
    # normalize to temporal mean = 1
    nr = nr ./ mean(nr,dims=2)
    # get variance
    nr_var = var(nr,dims=2)
    for i in eachindex(preprocessed)
        # multiply with the signal (squared, because variance, not rms)
        noise_variance = nr_var .* preprocessed[i].I.^2
        # # extract covariance matrix           # this is not robust because of an ill conditioned matrix
        # noise_cov = cov(noise,dims=2)
        # # compute the inverse
        # noise_cov_inv = inv(noise_cov)
        preprocessed[i].w = 1 ./ noise_variance
    end
    # display(heatmap(noise_cov))
end

# ## check the model
# import LinearAlgebra


# nr = load_spectra("./src/examples/data/210622A144.SPE";roi=700:1000)
# nr.I = (nr.I ./ sum(nr.I,dims=1))


# nr_mean = mean(nr.I,dims=2)
# nr_std = std(nr.I,dims=2)
# nr_cov = cov(nr.I,dims=2)


# heatmap(nr_cov)
# LinearAlgebra.cond(nr_cov)
# nr_cov_inv = inv(nr_cov)
# nr_cov_inv*nr_cov
# heatmap(nr.I)
# plot(nr_mean)
# plot!(nr_std)
# plot(nr_std./nr_mean)
# heatmap(nr_cov)
# heatmap(nr_cov_inv[50:60,50:60])

# plot(LinearAlgebra.diag(nr_cov))
# for i = 1:10
# plot!(LinearAlgebra.diag(nr_cov,i))
# end
# display(plot!())

# @time nr_cov_inv*nr_cov
# @time nr_cov_inv*nr_cov
# @time nr_cov_inv*nr_cov
# @time nr_cov_inv*nr_cov

# @time LinearAlgebra.inv!(LinearAlgebra.lu!(nr_cov))
# heatmap(nr_cov)
# heatmap(nr_cov_inv)
