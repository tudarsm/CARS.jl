function RecipesBase.plot(preprocessed::AbstractArray{PreProcSpectrum})
    return plot([preprocessed[i].I for i in eachindex(preprocessed)],legend=false,xlabel="px",ylabel="√I")
end

function RecipesBase.plot(preprocessed::PreProcSpectrum)
    return plot(preprocessed.I,legend=false,xlabel="px",ylabel="√I")
end

function plotSpec(ωE,IE)
    # plot the experimental data
    plot(ωE,IE,label="Experiment",legend=:topright)
    xlabel!("Raman shift in cm-1")
    ylabel!("CARS Signal")
end

function plotSpec(params)
    # plot the fit and annotate with temperature and mole Fractions
    ann = ""
    for (key,value) in params.X
        ann = ann .* "$(key) = $(round(value;sigdigits=4)), "
    end
    ann = ann .* "T=$(round(params.T;sigdigits=5))"
    # plot the experimental data
    plot(params.Spectra.omega,params.Spectra.Exp,label="Experiment",legend=:topright)
    plot!(params.Spectra.omega,params.Spectra.Fit,label="Fit: $(ann)")
    xlabel!("Raman shift in cm-1")
    ylabel!("CARS Signal")
    display(plot!(params.Spectra.omega,params.Spectra.Res,label="Residual"))
end

# new and fancy version of plotSpec including self updating scatter plots!
# function plotSpec(params)
# @async begin
#     # scatter plots with temperature and species mole fractions
#     names = reshape([key for (key,value) in params.X],1,length(params.X))
#     molefractions = reshape([value for (key,value) in params.X],1,length(params.X))
#     temperatures = ones(size(molefractions)).*params.T
#     global pa
#     if !@isdefined pa
#         pa = scatter()
#     else
#         pa = scatter(pa,temperatures,molefractions,legend=false,layout=[length(params.X)],markercolor=:black,ylabel=names)
#     end
        

#     # # plot the fit and annotate with temperature and mole Fractions
#     # ann = ""
#     # for (key,value) in params.X
#     #     ann = ann .* "$(key) = $(round(value;digits=3)), "
#     # end
#     # ann = ann .* "T=$(round(params.T;sigdigits=5))"

#     # # plot the experimental data
#     ps = plot(params.Spectra.omega,params.Spectra.Exp,label="Experiment")
#     ps = plot!(params.Spectra.omega,params.Spectra.Fit,label="Fit")
#     ps = plot!(params.Spectra.omega,params.Spectra.Res,label="Residual")
#     ps = plot!(legend=:outertop,legend_column=-1)
#     ps = xlabel!("Raman shift in cm-1")
#     ps = ylabel!("CARS Signal")
#     p = plot(pa,ps,layout=[1;1],reuse=true)
#     display(p)
# end
# end

function plotSpec(ωE,IE,params)
    # plot the experimental data
    plot(ωE,IE./sum(IE),label="Experiment",legend=:topright)
    # plot the fit and annotate with temperature and mole Fractions
    ann = ""
    for (key,value) in params.X
        ann = ann .* "$(key) = $(round(value;sigdigits=4)), "
    end
    ann = ann .* "T=$(round(params.T;sigdigits=5))"
    plot!(ωE,params.Spectra.Fit,label="Fit: $(ann)")
    xlabel!("Raman shift in cm-1")
    ylabel!("CARS Signal")
    display(plot!(ωE,params.Spectra.Res,label="Residual"))

end
