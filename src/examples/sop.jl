"""
    solveSOPPeel(num::T,maxWidth::Int;peelSetting::PeelSetting=peellen, rrbSetting::Union{RRBSetting,Nothing}=sop,loggingOn::Bool=true, debugOn::Bool=false, timeLimit::Union{Nothing, T}=nothing, useN::Bool=false, useParallel::Bool=false,fileName::Union{Nothing,String}=nothing) where T<:Int

A method for testing various CVRP solver functions

# Arguments
- `num::T`: which problem from the set to get
- `maxWidth::Int`: The max width allowed in a decision diagram layer
- `peelSetting::PeelSetting`:How to pick the peel node
- `loggingOn::Bool`: An optional parameter that can be used to turn off logging statements in other functions
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions
- `timeLimit::Union{Nothing, U}`: An optional parameter that adds a time limit
- `useN::Bool`: An optional paramter to use a width of n
- `fileName::Union{Nothing,String}`: The file to write results to
"""
function solveSOPPeel(num::T,maxWidth::Int;peelSetting::PeelSetting=peellen, rrbSetting::Union{RRBSetting,Nothing}=sop,loggingOn::Bool=true, debugOn::Bool=false, timeLimit::Union{Nothing, T}=nothing, useN::Bool=false, useParallel::Bool=false,fileName::Union{Nothing,String}=nothing) where T<:Int
    model = readSOPFile(getSOPFilePath(num))

    if useN
        maxWidth = length(model)
    end
    start = time()
    rrbValues, rrbLastArc = precalculateSOPRRB(model)

    if !isnothing(fileName)
        f = open(fileName,"a")
            write(f, "__________________________________\n\n")
            write(f, string("Peel Setting: ", peelSetting, "\n"))
            write(f, string("Max Width: ", maxWidth, "\n"))
            if useParallel
                write(f, string("Processors: ", Threads.nthreads(), "\n"))
            end
            write(f, "Problem $num: \n")
        close(f)
    end
    if !useParallel
        bestKnownSolution, bestKnownValue, bestBound = peelAndBound(
            model,
            maxWidth,
            orderingType=soporder,
            peelSetting=peelSetting,
            rrbSetting=rrbSetting,
            rrbPacket=(rrbValues, rrbLastArc),
            loggingOn=loggingOn,
            debugOn=debugOn,
            timeElapsed=time()-start,
            timeLimit=timeLimit,
            fileName = fileName
        )
    else
        bestKnownSolution, bestKnownValue, bestBound = parallelPeelAndBound(
            model,
            maxWidth,
            orderingType=soporder,
            peelSetting=peelSetting,
            rrbSetting=rrbSetting,
            rrbPacket=(rrbValues, rrbLastArc),
            loggingOn=loggingOn,
            debugOn=debugOn,
            timeElapsed=time()-start,
            timeLimit=timeLimit,
            fileName = fileName
        )
    end
    if !isnothing(fileName)
        f = open(fileName,"a")
            writeSOPSolution(bestKnownSolution, bestKnownValue, bestBound,fileName)
            write(f, string(bestKnownSolution,",",bestKnownValue,",",bestBound, "\n"))
        close(f)
    end
    return bestKnownSolution, bestKnownValue, bestBound
end


"""
    solveSOPBnB(num::T,maxWidth::Int = 32; rrbSetting::Union{RRBSetting,Nothing}=sop,loggingOn::Bool=true, debugOn::Bool=false, timeLimit::Union{Nothing, T}=nothing, useN::Bool=false,fileName::Union{Nothing,String}=nothing) where T<:Int

A method for testing various CVRP solver functions

# Arguments
- `num::T`: which problem from the set to get
- `maxWidth::Int`: The max width allowed in a decision diagram layer
- `loggingOn::Bool`: An optional parameter that can be used to turn off logging statements in other functions
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions
- `timeLimit::Union{Nothing, U}`: An optional parameter that adds a time limit
- `useN::Bool`: An optional paramter to use a width of n
- `fileName::Union{Nothing,String}`: The file to write results to
"""
function solveSOPBnB(num::T,maxWidth::Int = 32; rrbSetting::Union{RRBSetting,Nothing}=sop,loggingOn::Bool=true, debugOn::Bool=false, timeLimit::Union{Nothing, T}=nothing, useN::Bool=false,fileName::Union{Nothing,String}=nothing) where T<:Int
    model = readSOPFile(getSOPFilePath(num))

    if useN
        maxWidth = length(model)
    end
    start = time()
    rrbValues, rrbLastArc = precalculateSOPRRB(model)

    if !isnothing(fileName)
        f = open(fileName,"a")
        write(f, "__________________________________\n\n")
            write(f, "Problem $num: \n")
            write(f, string("Max Width: ", maxWidth, "\n"))
        close(f)
    end
    bestKnownSolution, bestKnownValue, bestBound = branchAndBound(
        model,
        maxWidth,
        orderingType=soporder,
        rrbSetting=rrbSetting,
        rrbPacket=(rrbValues, rrbLastArc),
        loggingOn=loggingOn,
        debugOn=debugOn,
        timeElapsed=time()-start,
        timeLimit=timeLimit,
        fileName=fileName
    )
    if !isnothing(fileName)
        f = open(fileName,"a")
            writeSOPSolution(bestKnownSolution, bestKnownValue, bestBound,fileName)
            write(f, string(bestKnownSolution,",",bestKnownValue,",",bestBound, "\n"))
        close(f)
    end
    return bestKnownSolution, bestKnownValue, bestBound
end

function benchmarkSOPs(width::Int, fileName::String, timeLimit::Int;numStart=1,numEnd=41,usePeel::Bool=true, peelSetting::PeelSetting=peellen, useN::Bool=false, useParallel::Bool=false)
    if usePeel
        solveSOPPeel(8,64, loggingOn=true, debugOn=false)
    else
        solveSOPBnB(8,324, loggingOn=true, debugOn=false)
    end
    for i in numStart:numEnd
        if usePeel
            block = @timed solveSOPPeel(i,width,peelSetting=peelSetting, loggingOn=true, debugOn=false, timeLimit=timeLimit, useN=useN, useParallel=useParallel, fileName=fileName)
        else
            
            block = @timed solveSOPBnB(i,width, loggingOn=true, debugOn=false, timeLimit=timeLimit, fileName=fileName)
        end
        f = open(fileName,"a")
            write(f, string(block[2], " seconds \n"))
        close(f)
    end
end

"""
    getSOPFilePath(num::T;getSet::Bool=false) where T<:Int

The locations of the benchmark problems for convenience

# Arguments
- `num::T`: which problem from the selected set to get
- `getSet::Bool`: An optional parameter that can be used to get an entire set instead of just one
"""
function getSOPFilePath(num::T;getSet::Bool=false) where T<:Int
    #SOPLIB: http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/
    paths = String[
        "./examples/sop/ESC07.sop",#2125 (1)
        "./examples/sop/ESC11.sop",#2075 (2)
        "./examples/sop/ESC12.sop",#1675 (3)
        "./examples/sop/ESC25.sop",#1681 (4)
        "./examples/sop/ESC47.sop",#1288 (5)
        "./examples/sop/ESC63.sop",#62 (6)
        "./examples/sop/ESC78.sop",#18230 (7)
        "./examples/sop/br17.10.sop",#55 (8)
        "./examples/sop/br17.12.sop",#55 (9)
        "./examples/sop/ft53.1.sop",#7531 (10)
        "./examples/sop/ft53.2.sop",#8026 (11)
        "./examples/sop/ft53.3.sop",#10262 (12)
        "./examples/sop/ft53.4.sop",#14425 (13)
        "./examples/sop/ft70.1.sop",#39313 (14)
        "./examples/sop/ft70.2.sop",#[40101,40419]  (15)
        "./examples/sop/ft70.3.sop",#42535  (16)
        "./examples/sop/ft70.4.sop",#53530  (17)
        "./examples/sop/kro124p.1.sop",#[38762,39420]  (18)
        "./examples/sop/kro124p.2.sop",#[39841,41336]  (19)
        "./examples/sop/kro124p.3.sop",#[43904,49499]  (20)
        "./examples/sop/kro124p.4.sop",#[73021,76103]  (21)
        "./examples/sop/p43.1.sop",#28140  (22)
        "./examples/sop/p43.2.sop",#28480  (23)
        "./examples/sop/p43.3.sop",#28835 (24)
        "./examples/sop/p43.4.sop",#83005  (25)
        "./examples/sop/prob.42.sop",#243 (26)
        "./examples/sop/prob.100.sop",#[1045,1163]  (27)
        "./examples/sop/rbg048a.sop",#351 (28)
        "./examples/sop/rbg050c.sop",#467 (29)
        "./examples/sop/rbg109a.sop",#1038 (30)
        "./examples/sop/rbg150a.sop",#1750 (31)
        "./examples/sop/rbg174a.sop",#2033 (32)
        "./examples/sop/rbg253a.sop",#2950 (33)
        "./examples/sop/rbg323a.sop",#3140 (34)
        "./examples/sop/rbg341a.sop",#2568 (35)
        "./examples/sop/rbg358a.sop",#2545 (36)
        "./examples/sop/rbg378a.sop",#[2809,2816]  (37)
        "./examples/sop/ry48p.1.sop",#15805 (38)
        "./examples/sop/ry48p.2.sop",#[16074,16666]  (39)
        "./examples/sop/ry48p.3.sop",#[19490,19894]  (40)
        "./examples/sop/ry48p.4.sop"#31446  (41)
    ]
    if getSet
        return paths
    else
        return paths[num]
    end
end
