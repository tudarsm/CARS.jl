###### this file contains functions that can be used for live fitting

function quasiLiveFit(folder::String,load::Function,preprocess::Function,fit::Function;skip=1,newest_only=true)
    try
        # list content of folder
        files = readdir(folder)

        println("Waiting for new files... Press CTRL+C to stop.")
        
        initplot = true
        ps = nothing    # handle on scatter plots
        pl = nothing    # handle on line plots
        pr = nothing    # handle on spectra
        idx = 1

        # outer while loop, run until CTRL+C is pressed
        while true
            
            
            # check for new files
            newfile = watch_folder(folder)
            if newfile.second.renamed
                print("Found file $(newfile.first). Start processing...")
                if newest_only
                    unwatch_folder(folder)
                end
                try
                    # load spectra after waiting until size does not change anymore for more than 200 ms
                    initsize = sizeof(newfile.first)
                    sleep(0.2)
                    while initsize != sizeof(newfile.first)
                        initsize = sizeof(newfile.first)
                        sleep(0.2)
                    end
                    signal = load(newfile.first)
                    # preprocess
                    preprocessed = preprocess(signal)
                    # skip samples
                    preprocessed = preprocessed[1:skip:end]
                    # fit
                    results = fit(preprocessed)

                    # extract data for plot
                    X = reduce(hcat,[[value for (key,value) in results[i].X] for i in eachindex(results)])'
                    T = [results[i].T for i in eachindex(results)]
                    Tx = ones(size(X)).*T
                    names = reshape([key for (key,value) in results[1].X],1,length(results[1].X))

                    if initplot
                        ps = scatter(layout=[length(names)])
                        pl = plot(layout=[1;length(names)])
                        initplot = false
                    end

                    scatter!(ps,Tx,X,legend=false,ylabel=names,layout=[length(names)])
                    scatter!(pl,[idx],[mean(T)],ylabel="T",legend=false,subplot=1)
                    scatter!(pl,[idx],[mean(X,dims=1)],ylabel=names,legend=false,subplot=2:(1+length(names)))

                    idx=idx+1

                    display(plot(ps,pl,layout=[1;1;1],size=(1000,1000),reuse=true))

                    # show data

                    println("done!")
                    println("Waiting for new files... Press CTRL+C to stop.")
                catch e
                    warn("Something went wrong when processing $(newfile.first)...")
                    println(e)
                end
            end

        end

    # catch interrupt
    catch e 
        if e isa InterruptException
            unwatch_folder(folder) 
            println("Stopped.")
            return
        end
        rethrow(e)
    end

end

