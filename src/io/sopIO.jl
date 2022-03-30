"""
    readSOPFile(filePath::String)

Read a properly formatted (SOP according to TSPLIB) .sop file and return a SequencingModel

# Arguments
- `filePath::String`: The location of the sop file
"""
function readSOPFile(filePath::String)
    model = nothing
    open(filePath) do openedFile
        model = parseSOPInput!(read(openedFile, String))
    end
    return model
end

"""
    parseSOPInput!(input_data::String)

Convert the contents of a file, and return a SequencingModel

# Arguments
- `input_data::String`: The contents of file
"""
function parseSOPInput!(input_data::String)
    #split the file into lines and store each line as a node object
    lines = split(input_data, '\n')

    #values
    dimension = 0
    objectiveValues = Dict{Tuple{Int, Int}, Int}()
    precedenceP2F = Dict{Int, Vector{Int}}()#map of priors to followers
    precedenceF2P = Dict{Int, Vector{Int}}()#map of followers to priors

    #bools for reaing file
    matrixSectionHeading = false
    matrixSection = false

    #counter for reading file
    n = 1

    #main loop
    for i in 1:length(lines)
        splitLine = split(lines[i], (' ','\t','\r'))
        filter!(x -> x!="", splitLine)
        if length(splitLine)>=1
            if occursin("EDGE_WEIGHT_SECTION",splitLine[1])
                matrixSectionHeading = true
            elseif matrixSectionHeading
                matrixSectionHeading = false
                matrixSection = true
                dimension = parse(Int,splitLine[1])
            elseif matrixSection && n <= dimension
                for j in 1:length(splitLine)
                    value  = splitLine[j]
                    if length(value)>=1
                        val = parse(Int,value)
                        #if precedence constraint
                        if val == -1 && n!=1 && j!=1
                            if haskey(precedenceP2F,j)
                                push!(precedenceP2F[j],n)
                            else
                                precedenceP2F[j] = [n]
                            end
                            if haskey(precedenceF2P,n)
                                push!(precedenceF2P[n],j)
                            else
                                precedenceF2P[n] = [j]
                            end
                        elseif n != j && val!= -1 #else if it does not point to itself
                            objectiveValues[(n, j)] = val
                        end
                    end
                end
                n = n+1
            end
        end
    end
    delete!(objectiveValues, (1, dimension))
    precedenceNumbers = getPrecedenceNumbers(precedenceF2P)
    for i in 2:dimension
        if !haskey(precedenceNumbers,i)
            precedenceNumbers[i] = 1
        else
            precedenceNumbers[i] += 1
        end
    end
    return SequencingModel(
        ExtensionalArcObjective(ExtensionalArcFunction(objectiveValues),minimization),
        dimension,
        1,
        dimension,
        true,
        precedenceP2F,
        precedenceF2P,
        precedenceNumbers
    )
end

"""
    printSOPSolution(bestKnownSolution::Vector{Int}, bestKnownValue::T) where T<:Real

Print a SOP solution to the terminal

# Arguments
- `bestKnownSolution::Vector{Int}`: The solution as an ordered list of nodes
- `bestKnownValue::T`: The cost of the solution
"""
function printSOPSolution(bestKnownSolution::Vector{Int}, bestKnownValue::T) where T<:Real
    println()
    print("Sequence: ",)
    for i in 1:length(bestKnownSolution)
        print(bestKnownSolution[i], " ")
    end
    if isinteger(bestKnownValue)
        println("\nCost: ", Int(bestKnownValue))
    else
        println("\nCost: ", bestKnownValue)
    end
end

"""
    writeSOPSolution(bestKnownSolution::Vector{Int}, bestKnownValue::T, fileName::String; time::Union{U,Nothing}=nothing, timeUnit::Union{String,Nothing}=nothing) where{T<:Real, U<:Real}

Print a SOP solution to the terminal

# Arguments
- `bestKnownSolution::Vector{Int}`: The solution as an ordered list of nodes
- `bestKnownValue::T`: The cost of the solution
- `bestKnownValue::T`: The best known relaxed bound
- `fileName::String`: The file to write to
- `time::U`: The time elapsed
- `timeUnit`: The unit of time being used
"""
function writeSOPSolution(bestKnownSolution::Vector{Int}, bestKnownValue::T, bestKnownBound::T, fileName::String; time::Union{U,Nothing}=nothing, timeUnit::Union{String,Nothing}=nothing) where{T<:Real, U<:Real}
    f = open(fileName,"a")
        write(f, "Seqeunce: ")
        for i in 1:length(bestKnownSolution)
            write(f,string(bestKnownSolution[i], " "))
        end
        if isinteger(bestKnownValue)
            write(f,"\nCost: ", string(Int(bestKnownValue)))
        else
            write(f,"\nCost: ", string(bestKnownValue))
        end
        if isinteger(bestKnownBound)
            write(f,"\nRelaxed Bound: ", string(Int(bestKnownBound)))
        else
            write(f,"\nRelaxed Bound: ", string(bestKnownBound))
        end
        write(f,"\n")
        if !isnothing(time)
            write(f, "Time Elapsed: ", string(time), " $timeUnit\n")
        end
    close(f)
end
