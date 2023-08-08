## THIS FILE CONTAINS ALL CODE TO PREPROCESS MEASURED SPECTRA
# it uses mulitple dispatch to have a high flexibility
# if a Nx1xF or NxF Matrix is passed as an argument, it will assume N to be the spectral direction and F the number of frames
# if a NxMxF Matrix is passed (e.g. for multiple strips containing other useful data), a function has to be implemented here
# this could be useful if e.g. the background subtraction is to be performed in situ or a single-shot NR signal is captured at the same time on a different strip
# these methods all then call preprocess_spectrum to do the actual preprocessing
#
#   preprocessing includes:
#   eliminating samples based on the following criteria: saturation, signal too low, breakdown detection
#   background subtraction of the resonant signal
#   background correction of the non-resonant signal
#   correction of the resonant signal to the non-resonant signal
#   crop to area of interest


Base.@kwdef mutable struct ErrorCriteria
    max_signal::Float64 = 60000.   # If any countlevel in the spectrum is higher than this, it will be excluded from the fit
    min_signal::Float64 = 500.     # If the peak of a spectrum after(!) background subtraction is below this limit, it will be excluded
    breakdown_threshold::Float64 = 200  # if the signal is more than 100 counts above background at the specified location, it is assumed to be a breakdown
    breakdown_location::Int32 = -1  # Specify pixel for breakdown detection. Will be ignored if == -1.
end

# this function does the preprocessing of all spectra in a RawSpectra struct
"""
    preprocess_spectra(signal::RawSpectra,nr_signal::RawSpectra,roi;errorCriteria::ErrorCriteria=ErrorCriteria())
Preprocesses spectra passed to this function as signal.
It uses the some criteria to exclude spectra from fitting, e.g. when the signal is too high or too low or if a breakdown occured.
Check out the example ```fit_experimental_data.jl``` for how to adapt the ErrorCriteria to your needs.
Signals are normalized to the supplied non-resonant signal nr_signal of type [`RawSpectra`](@ref) and, if specified in [`CARSOptions`](@ref) the square root is taken.
"""
function preprocess_spectra(signal::RawSpectra,nr_signal::RawSpectra,roi;errorCriteria::ErrorCriteria=ErrorCriteria())

    # assume valid spectra in the beginning. will be overwritten by errors.
    errorcode = zeros(size(signal.I,2))

    # eliminate samples. preprocessing will continue, but samples with a non-zero errorcode will be omitted in fitting
    errorcode[vec(maximum(signal.I,dims=1).>errorCriteria.max_signal)] .= -1


    # check the non-resonant signal for saturation
    if sum(nr_signal.I .> errorCriteria.max_signal) > 0
        warn("NR Signal has saturated samples. Removing them before normalization.")
        warn("You might want to check if everything seems reaonsable.")
        idx_to_remove = vec(maximum(nr_signal.I;dims=1) .> errorCriteria.max_signal)
        idx_to_keep = vec(maximum(nr_signal.I;dims=1) .< errorCriteria.max_signal)
        warn("Removing samples $(findall(idx_to_remove)) from NR signal")
        nr_signal.I = nr_signal.I[:,idx_to_keep]
    end

    # background subtractions
    spectrum_bgcorr = signal.I .- signal.BG
    nr_bgcorr = nr_signal.I .- nr_signal.BG

    # eliminate samples with too weak signal
    errorcode[vec(maximum(signal.I,dims=1).<errorCriteria.min_signal)] .= -2

    # eliminate breakdowns
    if errorCriteria.breakdown_location != -1 
        errorcode[vec(signal.I[errorCriteria.breakdown_location,:].>errorCriteria.breakdown_threshold)] .= -3
    end

    # normalize to Argon
    spectrum_norm = spectrum_bgcorr ./ mean(nr_bgcorr,dims=2)

    # take square root if desired
    if CARSopt.fittype == :sqrt
        spectrum_norm = @. sign(spectrum_norm)*sqrt(abs(spectrum_norm))
    end

    # crop to desired region of interest
    spectrum_norm = spectrum_norm[roi[1]:roi[2],:]

    # normalize to area
    A = sum(spectrum_norm,dims=1)
    spectrum_norm = spectrum_norm ./ A

    # store in an array of structs,return mean optionally
    return [PreProcSpectrum(spectrum_norm[:,i],[1.],A[i],errorcode[i]) for i in axes(spectrum_norm,2)]

end

# extend statistics' mean function to handle arrays of PreProcSpectrum easily
function Statistics.mean(preprocessed_signals::AbstractArray{PreProcSpectrum})
    # get all spectra with a zero errorcode
    idxs = [preprocessed_signals[i].ErrorCode for i in eachindex(preprocessed_signals)] .== 0
    # return the mean spectrum of all valid samples  
    meanspec = mean([preprocessed_signals[i].I for i in findall(idxs)])
    return PreProcSpectrum(meanspec,[1.],1.,0)
end


function showErrorCodes(spectra::AbstractArray{PreProcSpectrum})

    errorcodes = [spectra[i].ErrorCode for i in eachindex(spectra)]
    errorcount = sum( errorcodes.!= 0)
    totalsamples = length(spectra)
    println("Spectra contain $(totalsamples-errorcount) valid samples out of $(totalsamples) in total.")
    println("Signal too high: $(sum(errorcodes .== -1))")
    println("Signal too low: $(sum(errorcodes .== -2))")
    println("Breakdowns: $(sum(errorcodes .== -3))")


end