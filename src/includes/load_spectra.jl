## THIS FILE CAN BE USED TO IMPLEMENT DIFFERENT FUNCTIONS TO LOAD MEASURED CARS DATA
# Some structs
"""
    mutable struct RawSpectra
        I::AbstractArray{Float64}      # signal
        BG::AbstractArray{Float64}     # background
        SingleShotNR::AbstractArray{Float64}     # preparation for single shot NR correction. not yet implemented.
    end
Struct containing raw spectra after loading from a file.
"""
mutable struct RawSpectra
    I::AbstractArray{Float64}      # signal
    BG::AbstractArray{Float64}     # background
    SingleShotNR::AbstractArray{Float64}     # preparation for single shot NR correction. not yet implemented, see below        
end
# use this to automatically load data with the correct method based on the extension
function load_spectra(filename::String;roi=-1:1:1)

    if contains(lowercase(filename),".spe")
        spectra = load_spe(filename)
    end
    # add new if clauses here for new filetypes and implement the correct method somewhere in this file
    # they all shoud return a NxMxF matrix, where N is the spectral direction, M the number of strips and F the number of frames
    # if M is equal to one, it assumes a rolling background was used (i.e. camera records twice as fast as the lasers are running)
    # if M is not equal to one, a method has to be implemented to handle this case
    # this may be helpful if e.g. the background is extracted above or below the cars signal in a different strip
    # also, this could be used for single-shot NR correction on another strip. be creative and implement the methods.

    # separate the signal from the background. depending on the size of spectra,
    # this uses multiple dispatch to do the separation.
    # in case multiple strips are used, a method stub is implemented below which has to be implemented for the situation at hand
    if size(spectra,2) == 1
        signal,background = separateRollingBackground(spectra)
        nrsignal = zeros(size(signal));
    else
        signal,background,nrsignal = separateBackground(spectra)
    end

    if roi[1] != -1
        return RawSpectra(signal[roi,:],background[roi,:],nrsignal[roi,:])
    else
        return RawSpectra(signal,background,nrsignal)
    end
end

function separateBackground(spectra)
    error("No method implemented for loaded spectra with size $(size(spectra)) of type $(typeof(spectra)).")
end

function separateRollingBackground(spectra)
    # remove singleton dimensions
    spectra = dropdims(spectra;dims=2)

    # detect the background. sum the frames to detect it
    specsum = vec(sum(spectra;dims=1))

    # detect the signal
    if sum(specsum[1:2:end]) > sum(specsum[2:2:end])
        signal = spectra[:,1:2:end]
        background = spectra[:,2:2:end]    
    else
        signal = spectra[:,2:2:end]
        background = spectra[:,1:2:end]    
    end

    # check the non-resonant signal for saturation
    if sum(background .> 64000) > 0
        warn("Background Signal has saturated samples. Removing them before normalization.")
        warn("You might want to check if everything seems reaonsable.")
        idx_to_remove = vec(maximum(background;dims=1) .> 64000)
        idx_to_keep = vec(maximum(background;dims=1) .< 64000)
        warn("Removing samples $(findall(idx_to_remove)) from background signal")
        background = background[:,idx_to_keep]
    end

    # rolling background is averaged
    background = mean(background;dims=2)


    return signal,background
end


################################
# *.SPE files for PIXIS camera #
# Ported from speReadHead_pm.m #
################################
function load_spe(filename)
    
    # helper function to read a specific type of certain length starting from a certain position
    function readbin(io::IO,T,n,pos)
        data = Array{T,1}(undef,n)
        seek(io,pos)
        read!(io,data)

        # convert character array to string
        if T == Char
            data = String(data)
        end

        # if it's just one element, return a scalar
        if n == 1
            data = data[1]
        end
      
        return data
    end
  
    # % open file and read data
    io = open(filename)

    # read meta information. note that not all is used
    xdim = readbin(io,UInt16,1,42)
    # avgx = readbin(io,Int16,1,2)
    datyp = readbin(io,Int16,1,108)
    # date = readbin(io,Char,10,20)
    # hour = readbin(io,Char,2,30)
    # min = readbin(io,Char,2,32)
    # nosc = readbin(io,Int16,1,34)
    # secs = readbin(io,Char,2,38)
    ydim = readbin(io,Int16,1,656)
    # lnosc = readbin(io,Int32,1,664)
    # lavg = readbin(io,Int32,1,668)
    frames = readbin(io,Int16,1,1446)
    
    DT = [Float32,Int32,Int16,UInt16]

    # read header
    # head = readbin(io,Int16,2050,0)
    
    # read the data
    data = readbin(io,DT[1+datyp],Int(xdim)*Int(ydim)*(frames),4100) # need to convert xdim to int to avoid overflow error
    data = convert(Vector{Float64},data)
    # convert to a 3D array
    data = reshape(data,(xdim,ydim,frames))
    close(io);

    return data

end
