__precompile__()

module AdaptiveLASSO

    ### exports
    export
    # main function
    LSIAF,
    # input names
    CoordinateDescent3, CoordinateDescentTM

    ### include source files

    # core algorithm
    include("utils.jl")
    include("ComposeMatrices.jl")
    include("CoordinateDescentTM.jl")
    include("CoordinateDescent3.jl")
    include("adaLASSO.jl")
    include("LSIAF.jl")

end
