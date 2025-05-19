using Random, Distributions, ProgressMeter, StatsBase, PlotlyJS, BenchmarkTools, HypothesisTests

function compute_difference(x::Int, y::Int)::Vector{Float64}
    v1 = [1 / (y + 1):1 / (y + 1):y / (y + 1); ones(x)]
    v2 = [zeros(x); 1 / (y + 1):1 / (y + 1):y / (y + 1)]
    return v1 .- v2
end

function compute_difference(x::Int, y::Int)::Vector{Float64}
    v1 = [1:2:2*y; fill(2*y,x)]
    v2 = [fill(0,x); 1:2:2*y]
    return v1 .- v2
end

function random_vector_efficient(x::Int, y::Int)::Vector{Int}
    indices = randperm(x + y)  # 随机排列索引
    vec = zeros(Int, x + y)    # 初始化为全 0
    vec[indices[1:y]] .= 1     # 将前 y 个随机位置设置为 1
    return vec
end

function cal_MSE(i::Int, y::Vector)
    p = i/10
    n_max = length(y)
    n = 1:n_max
    P = (1-p)^2 .* n .* p .^(n.-1)
    MSE = mean((y .- P).^2)
    return MSE
end

function Data_Plot(Sys)
    fig = make_subplots(
        rows=3, cols=3,
        specs=[
            Spec(kind="xy", colspan=2)   missing                     Spec(kind="xy")
            # missing                                 missing                     missing
            Spec(kind="xy", colspan=2)   missing                     Spec(kind="xy")
            # missing                                 missing                     missing
            Spec(kind="xy")              Spec(kind="xy")             Spec(kind="xy")
            # missing                                 missing                     missing
        ]
    )
    color_map = Dict(
        80 => "purple",
        70 => "blue",
        60 => "cyan ",
        50 => "green",
        40 => "yellow",
        30 => "orange",
        20 => "red",
        10 => "brown" # gray pink
    )

    max_x = 3*maximum(Sys[:Series][:X_n])
    for i in 80:-10:10
        trace = scatter(x=1:max_x, y=Sys[:Critical]["$(i)%"][:Data], name="$(i)%", mode="lines", line=attr(color=color_map[i]))
        add_trace!(fig, trace, row=1, col=1)
    end
    
    k1 = 1
    k2 = 1
    Gap = 2
    Critical = ["10%","20%","30%","40%","50%","60%","70%","80%"]
    function cut(values, breaks; extend::Bool = false)
        intervals = Vector{Int}(undef, length(values))
        breaks = sort(breaks)  # Ensure breaks are sorted
        if extend
            breaks = [minimum([minimum(values); breaks]) - eps(), breaks..., maximum([maximum(values); breaks]) + eps()]
        end
        length_value = length(values)
        for i in 1:length_value
            intervals[i] = findfirst(x -> values[i] <= x, breaks) - 1
        end
        
        return intervals
    end

    for i in 80:-10:10
        Weight = Sys[:Critical]["$(i)%"][:Data]
        Max_MW = findlast(x -> x != 0, Weight)
        Breaks = 0:Gap:Max_MW
        Intervals = cut(1:Max_MW, Breaks, extend = true)
        Temp = Weight[1:Max_MW] .* (1:Max_MW)
        Weight_intervals = [sum(Temp[findall(x -> x == interval, Intervals)]) for interval in unique(Intervals)]
        Numbers = [sum(Weight[1:Max_MW][findall(x -> x == interval, Intervals)]) for interval in unique(Intervals)]

        GPC_Result = zeros(length(Numbers), 2)
        non_zero_indices = findall(x -> x != 0, Weight_intervals)
        if !isempty(non_zero_indices)
            GPC_Result[non_zero_indices, 1] .= log10.(Weight_intervals[non_zero_indices] ./ Numbers[non_zero_indices]) .* k1 .+ k2
            GPC_Result[non_zero_indices, 2] .= (-0.4228 .* GPC_Result[non_zero_indices, 1] .+ 10.38) .* Weight_intervals[non_zero_indices]
            GPC_Result = GPC_Result[findall(x -> x != 0, GPC_Result[:, 2]), :]
            trace = scatter(x = GPC_Result[:, 1], y = GPC_Result[:, 2], name="$(i)%", mode = "lines", line=attr(color=color_map[i]), showlegend=false)
            add_trace!(fig, trace, row=2, col=1)
        end
    end

    add_trace!(
        fig,
        scatter(x=Sys[:Series][:Time], y=Sys[:Series][:Radical], mode="markers", name="Radical", line=attr(color="gray")),
        row=1, col=3
    )

    add_trace!(
        fig,
        scatter(x=Sys[:Series][:Time], y=sqrt.(5.981*10^(-14).*exp.(-3.780*10^(-5)*Sys[:Series][:Time]))*Sys[:Initial][:M]/5, mode="lines", name="Theoretical", line=attr(color="red")),
        row=1, col=3
    )

    add_trace!(
        fig,
        scatter(
            x=Sys[:Series][:Time], 
            y=Sys[:Series][:P], 
            name="Time v.s. Conv.",
            mode="lines", line=attr(color="black"), 
            showlegend=false
        ),
        row=2, col=3
    )

    add_trace!(
        fig,
        scatter(
            x=Sys[:Series][:P], 
            y=Sys[:Series][:X_n], 
            mode="lines", line=attr(color="black"), 
            name="X_n",
            showlegend=false
            ),
        row=3, col=1
    )

    add_trace!(
        fig,
        scatter(
            x=Sys[:Series][:P], 
            y=Sys[:Series][:X_w], 
            name="X_n",
            mode="lines", line=attr(color="black"), 
            showlegend=false
        ),
        row=3, col=2
    )

    add_trace!(
        fig,
        scatter(x=Sys[:Series][:P], 
        y=Sys[:Series][:PDI], 
        name="X_n",
        mode="lines", line=attr(color="black"), 
        showlegend=false ),
        row=3, col=3
    )

    relayout!(fig, 
        xaxis = attr(title="Chain Length", autorange="reversed"),
        yaxis = attr(title="???"), 
        xaxis2 = attr(title="Time"),
        yaxis2 = attr(title="Amount"), 
        xaxis3 = attr(title="Chain Length"),
        yaxis3 = attr(title="???"), 
        xaxis4= attr(title="Time"),
        yaxis4 = attr(title="Conv."), 
        xaxis5= attr(title="Conv. "),
        yaxis5 = attr(title="X_n"),
        xaxis6= attr(title="Conv. "),
        yaxis6 = attr(title="X_w"),
        xaxis7= attr(title="Conv. "),
        yaxis7 = attr(title="PDI"),
        legend = attr(x = 1, y = 1, xanchor = "left", yanchor = "top"),
        template=templates.plotly_white
    )
    fig
end

function SSA_FRP(M::Int=10^9)
    # M = 10^9
    I = div(M, 100)
    # CTA = div(M, 100)

    P_lists = 0.0001:0.0001:0.1
    P_critical = 0:0.1:0.8
    Max_count = 6000
    
    k_I       = 3.780 * 10^(-5)
    k_p       = 1.039 * 10^( 3) * 5 / M 
    k_T       = 3.160 * 10^( 7) * 5 / M * 2
    k_CTAgent = 1.000 * 10^  5  * 0
    k_CTADE   = 1.000 * 10^  5  * 0

    # Current Stage
    Monomer = M
    Initiator = I
    Radical = Int[]
    # CTAgent = Int[]
    # Intermediate = Int[]
    Termination = zeros(Int, Max_count+1)
    P_Current = 0.0
    Time_Current = 0.0

    M_list           = zeros(Int, 8001)
    Time_list        = zeros(8001)
    Radical_list     = zeros(Int, 8001)
    Termination_list = zeros(Int, 8001)
    X_n_list         = ones(8001)
    X_w_list         = ones(8001)

    System = Dict{Symbol, Dict}(
        :k => Dict{Symbol,Float64}(
            :I           => k_I,
            :p           => k_p,
            :Termination => k_T,
            :CTAgent     => k_CTAgent,
            :CTADE       => k_CTADE
        ),
        :Initial => Dict(
            :M => M,
            :I => I
        ),
        :Critical => Dict{String, Dict}("$(i*10)%" => Dict{String,Vector{Int}}() for i in 1:8),
    )
    Temp_react = [k_I * Initiator, 0, 0]
    length_Radical = 0

    for i in 1:8
        println("\n$(P_critical[i]*100)% ~ $(P_critical[i+1]*100)%")

        @showprogress for j in 1:1000
            P_threshold = P_critical[i] + P_lists[j]

            while P_Current < P_threshold
                Temp_time = rand(Exponential(1 / sum(Temp_react)))
                Time_Current += Temp_time

                Temp_react_class = cumsum(Temp_react) / sum(Temp_react)
                Index_01 = rand()

                if Index_01 < Temp_react_class[1]
                    Initiator -= 1
                    push!(Radical, 0, 0)
                    length_Radical += 2
                    Temp_react = [k_I*Initiator, k_p*length_Radical*Monomer, k_T*(length_Radical)^2/2]
                elseif Index_01 < Temp_react_class[2]
                    Monomer -= 1
                    Index = rand(1:length_Radical)
                    Radical[Index] += 1
                    P_Current = 1 - Monomer / M
                    Temp_react[2] = k_p*length_Radical*Monomer
                elseif Index_01 < Temp_react_class[3] && length_Radical >= 2
                    selected = sort(sample(1:length_Radical, 2, replace=false))
                    Length = sum(Radical[selected])
                    Termination[Length + 1] += 1
                    deleteat!(Radical, selected)
                    length_Radical -= 2
                    Temp_react[[2,3]] = [k_p*length_Radical*Monomer, k_T*(length_Radical)^2/2]
                end
            end

            Index_State = 1000 * (i - 1) + j + 1
            M_list[Index_State] = Monomer
            Time_list[Index_State] = Time_Current
            Radical_list[Index_State] = length(Radical)
            Termination_list[Index_State] = sum(Termination)
            Data = copy(Termination[2:end])
            Counts = countmap(Radical)
            Counts = filter(x -> x[1] != 0, Counts)
            for (idx, count) in Counts
                if idx <= length(Data)
                    Data[idx] += count
                else
                    push!(Data, count)
                end
            end
            X_n_list[Index_State] = sum(Data .* (1:Max_count)) / sum(Data)
            X_w_list[Index_State] = sum(Data .* (1:Max_count).^2) / sum(Data .* (1:Max_count))
        end

        Data = copy(Termination[2:end])
        Counts = countmap(Radical)
        Counts = filter(x -> x[1] != 0, Counts)
        for (idx, count) in Counts
            if idx <= length(Data)
                Data[idx] += count
            else
                push!(Data, count)
            end
        end

        System[:Critical]["$(i*10)%"] = Dict(
            # :CTAgent => CTAgent,
            # :Intermediate => Intermediate,
            :Data => Data
        )
    end

    PDI_list = X_w_list ./ X_n_list
    System[:Series] = Dict{Symbol, AbstractVector}(
        :P           => 0:0.0001:0.8,
        :M           => M_list,
        :Time        => Time_list,
        :Radical     => Radical_list,
        :Termination => Termination_list,
        :X_n         => X_n_list,
        :X_w         => X_w_list,
        :PDI         => PDI_list
    )

    Data = Termination[2:end]
    Counts = countmap(Radical)
    Counts = filter(x -> x[1] != 0, Counts)
    for (idx, count) in Counts
        if idx <= length(Data)
            Data[idx] += count
        else
            push!(Data, count)
        end
    end
    System[:Bias] = Dict{Symbol,Int64}(
        :M => sum(Data .* (1:Max_count)) + Monomer             - M,
        :I => Initiator + length(Radical)/2 + sum(Termination) - I
    )

    return System
end

function Quasi_Leaping_FRP(M::Int=10^9)
    # M = 10^9
    I = M/100
    # CTA = M/100

    P_lists = 0.0001:0.0001:0.1
    P_critical = 0:0.1:0.8
    Max_count = 6000
    
    k_I       = 3.780 * 10^(-5)
    k_p       = 1.039 * 10^( 3) * 5 / M 
    k_T       = 3.160 * 10^( 7) * 5 / M * 2
    k_CTAgent = 1.000 * 10^  5  * 0
    k_CTADE   = 1.000 * 10^  5  * 0

    # Current Stage
    Monomer = M
    Initiator = I
    Radical = Int[]
    # CTAgent = Int[]
    # Intermediate = Int[]
    Termination = zeros(Int, Max_count+1)
    P_Current = 0.0
    Time_Current = 0.0

    M_list           = zeros(Int, 8001)
    Time_list        = zeros(8001)
    Radical_list     = zeros(Int, 8001)
    Termination_list = zeros(Int, 8001)
    X_n_list         = ones(8001)
    X_w_list         = ones(8001)

    System = Dict{Symbol, Dict}(
        :k => Dict{Symbol,Float64}(
            :I           => k_I,
            :p           => k_p,
            :Termination => k_T,
            :CTAgent     => k_CTAgent,
            :CTADE       => k_CTADE
        ),
        :Initial => Dict{Symbol,Int}(
            :M => M,
            :I => I
        ),
        :Critical => Dict{String, Dict}("$(i*10)%" => Dict{String,Vector{Int}}() for i in 1:8),
    )

    Temp_react = [k_I * Initiator, 0]
    length_Radical = 0

    for i in 1:8
        println("\n$(P_critical[i]*100)% ~ $(P_critical[i+1]*100)%")

        @showprogress for j in 1:1000
            P_threshold = P_critical[i] + P_lists[j]

            while P_Current < P_threshold
                Temp_time = rand(Exponential(1 / sum(Temp_react)))
                Time_Current += Temp_time

                Propagation = rand(Poisson(k_p * length_Radical * Monomer * Temp_time))
                Monomer -= Propagation
                
                Index = sample(1:length_Radical, Propagation, replace=true)
                    Counts = countmap(Index)
                    # Update Radical according to Counts
                    for (idx, count) in Counts
                        Radical[idx] += count
                    end
                P_Current = 1 - Monomer / M

                Temp_react_class = cumsum(Temp_react) / sum(Temp_react)
                Index_01 = rand()

                if Index_01 < Temp_react_class[1]
                    Initiator -= 1
                    push!(Radical, 0, 0)
                    length_Radical += 2
                    Temp_react = [k_I * Initiator,k_T * (length_Radical)^2 / 2]

                elseif Index_01 < Temp_react_class[2] && length_Radical >= 2
                    selected = sort(sample(1:length_Radical, 2, replace=false))
                    Length = sum(Radical[selected])
                    Termination[Length + 1] += 1
                    deleteat!(Radical, selected)
                    length_Radical -= 2
                    Temp_react[2] = k_T * (length_Radical)^2 / 2
                end
            end

            Index_State = 1000 * (i - 1) + j + 1
            M_list[Index_State] = Monomer
            Time_list[Index_State] = Time_Current
            Radical_list[Index_State] = length(Radical)
            Termination_list[Index_State] = sum(Termination)
            Data = copy(Termination[2:end])
            Counts = countmap(Radical)
            Counts = filter(x -> x[1] != 0, Counts)
            for (idx, count) in Counts
                if idx <= length(Data)
                    Data[idx] += count
                else
                    push!(Data, count)
                end
            end
            X_n_list[Index_State] = sum(Data .* (1:Max_count)) / sum(Data)
            X_w_list[Index_State] = sum(Data .* (1:Max_count).^2) / sum(Data .* (1:Max_count))
        end

        Data = Termination[2:end]
        Counts = countmap(Radical)
        Counts = filter(x -> x[1] != 0, Counts)
        for (idx, count) in Counts
            if idx <= length(Data)
                Data[idx] += count
            else
                push!(Data, count)
            end
        end
        #all = sum(Data)
        #Data = Data / all

        System[:Critical]["$(i*10)%"] = Dict(
            # :CTAgent => CTAgent,
            # :Intermediate => Intermediate,
            :Data => Data
        )
    end

    PDI_list = X_w_list ./ X_n_list
    System[:Series] = Dict{Symbol, AbstractVector}(
        :P           => 0:0.0001:0.8,
        :M           => M_list,
        :Time        => Time_list,
        :Radical     => Radical_list,
        :Termination => Termination_list,
        :X_n         => X_n_list,
        :X_w         => X_w_list,
        :PDI         => PDI_list
    )

    # Test Bias
    Data = Termination[2:end]
    Counts = countmap(Radical)
    Counts = filter(x -> x[1] != 0, Counts)
    for (idx, count) in Counts
        if idx <= length(Data)
            Data[idx] += count
        else
            push!(Data, count)
        end
    end
    System[:Bias] = Dict{Symbol,Int64}(
        :M => sum(Data .* (1:Max_count)) + Monomer             - M,
        :I => Initiator + length(Radical)/2 + sum(Termination) - I
    )
    
    # System[:MSE] = cal_MSE(8,Data)

    return System
end

function SSA_DT(M::Int=10^9)
    # M = 10^9
    I = div(M, 100)
    CTA = div(M, 100)

    P_lists = 0.0001:0.0001:0.1
    P_critical = 0:0.1:0.8
    Max_count = 2000

    k_I       = 3.780 * 10^(-5)
    k_p       = 1.039 * 10^( 3) * 5 / M 
    k_T       = 3.160 * 10^( 7) * 5 / M * 2
    k_CTAgent = 1.000 * 10^  5  * 5 / M
    k_CTADE   = 1.000 * 10^  5  * 0

    # Current Stage
    Monomer = M
    Initiator = I
    Radical = Int[]
    CTAgent = zeros(Int, Max_count+1)
    CTAgent[1] = CTA
    # Intermediate = Int[]
    Termination = zeros(Int, Max_count+1)
    P_Current = 0.0
    Time_Current = 0.0

    M_list           = zeros(Int, 8001)
    Time_list        = zeros(8001)
    Radical_list     = zeros(Int, 8001)
    Termination_list = zeros(Int, 8001)
    X_n_list         = ones(8001)
    X_w_list         = ones(8001)

    System = Dict{Symbol, Dict}(
        :k => Dict{Symbol,Float64}(
            :I           => k_I,
            :p           => k_p,
            :Termination => k_T,
            :CTAgent     => k_CTAgent,
            :CTADE       => k_CTADE
        ),
        :Initial => Dict{Symbol,Int}(
            :M => M,
            :I => I,
            :CTAgent   => CTA
        ),
        :Critical => Dict{String, Dict}("$(i*10)%" => Dict{String,Vector{Int}}() for i in 1:8),
    )

    Temp_react = [k_I * Initiator, 0,0,0]
    length_Radical = length(Radical)

    for i in 1:8
        println("\n$(P_critical[i]*100)% ~ $(P_critical[i+1]*100)%")

        @showprogress for j in 1:1000
            P_threshold = P_critical[i] + P_lists[j]

            while P_Current < P_threshold
                Temp_time = rand(Exponential(1 / sum(Temp_react)))
                Time_Current += Temp_time

                Temp_react_class = cumsum(Temp_react) / sum(Temp_react)
                Index_01 = rand()

                if Index_01 < Temp_react_class[1]
                    Initiator -= 1
                    push!(Radical, 0, 0)
                    length_Radical += 2
                    Temp_react = [k_I*Initiator, k_p*length_Radical*Monomer, k_CTAgent*length_Radical*CTA, k_T*(length_Radical)^2/2]
                elseif Index_01 < Temp_react_class[2]
                    Monomer -= 1
                    Index = rand(1:length_Radical)
                    Radical[Index] += 1
                    P_Current = 1 - Monomer / M
                    Temp_react[2] = k_p*length_Radical*Monomer
                elseif Index_01 < Temp_react_class[3]
                    Index = rand(1:length_Radical)
                    Selected_CTA = sample(1:Max_count+1, Weights(CTAgent))
                    CTAgent[Selected_CTA]     -= 1
                    CTAgent[Radical[Index]+1] += 1
                    Radical[Index] = Radical[end]
                    pop!(Radical)  
                    push!(Radical,Selected_CTA-1)
                elseif Index_01 < Temp_react_class[4] && length_Radical >= 2
                    selected = sort(sample(1:length_Radical, 2, replace=false))
                    Length = sum(Radical[selected])
                    Termination[Length + 1] += 1
                    deleteat!(Radical, selected)
                    length_Radical -= 2
                    Temp_react[[2,3,4]] = [k_p*length_Radical*Monomer, k_CTAgent*length_Radical*CTA, k_T*(length_Radical)^2/2]
                end
            end

            Index_State = 1000 * (i - 1) + j + 1
            M_list[Index_State] = Monomer
            Time_list[Index_State] = Time_Current
            Radical_list[Index_State] = length(Radical)
            Termination_list[Index_State] = sum(Termination)
            Data = copy(Termination[2:end])
            Counts = countmap(Radical)
            Counts = filter(x -> x[1] != 0, Counts)
            for (idx, count) in Counts
                if idx <= length(Data)
                    Data[idx] += count
                else
                    push!(Data, count)
                end
            end
            Counts = copy(CTAgent[2:end])
            Data += Counts
            X_n_list[Index_State] = sum(Data .* (1:Max_count)) / sum(Data)
            X_w_list[Index_State] = sum(Data .* (1:Max_count).^2) / sum(Data .* (1:Max_count))
        end

        Data = Termination[2:end]
        Counts = countmap(Radical)
        Counts = filter(x -> x[1] != 0, Counts)
        for (idx, count) in Counts
            if idx <= length(Data)
                Data[idx] += count
            else
                push!(Data, count)
            end
        end
        Counts = copy(CTAgent[2:end])
        Data += Counts
        
        System[:Critical]["$(i*10)%"] = Dict(
            # :CTAgent => CTAgent,
            # :Intermediate => Intermediate,
            :Data => Data
        )
    end

    PDI_list = X_w_list ./ X_n_list
    System[:Series] = Dict{Symbol, AbstractVector}(
        :P           => 0:0.0001:0.8,
        :M           => M_list,
        :Time        => Time_list,
        :Radical     => Radical_list,
        :Termination => Termination_list,
        :X_n         => X_n_list,
        :X_w         => X_w_list,
        :PDI         => PDI_list
    )

    # Test Bias
    Data = Termination[2:end]
    Counts = countmap(Radical)
    Counts = filter(x -> x[1] != 0, Counts)
    for (idx, count) in Counts
        if idx <= length(Data)
            Data[idx] += count
        else
            push!(Data, count)
        end
    end
    Counts = copy(CTAgent[2:end])
    Data += Counts
    System[:Bias] = Dict{Symbol,Int64}(
        :M => sum(Data .* (1:Max_count)) + Monomer             - M,
        :I => Initiator + length(Radical)/2 + sum(Termination) - I,
        :C => sum(CTAgent) - CTA
    )

    return System
end

#矩形
function Quasi_Leaping_DT(M::Int=10^9)
    # M = 10^9
    I = div(M, 100)
    CTA = div(M, 100)

    P_lists = 0.0001:0.0001:0.1
    P_critical = 0:0.1:0.8
    Max_count = 2000

    k_I       = 3.780 * 10^(-5)
    k_p       = 1.039 * 10^( 3) * 5 / M 
    k_T       = 3.160 * 10^( 7) * 5 / M * 2
    k_CTAgent = 1.000 * 10^  5  * 5 / M
    k_CTADE   = 1.000 * 10^  5  * 0

    # Current Stage
    Monomer = M
    Initiator = I
    Radical = Int[]
    CTAgent = zeros(Int, Max_count+1)
    CTAgent[1] = CTA
    # Intermediate = Int[]
    Termination = zeros(Int, Max_count+1)
    P_Current = 0.0
    Time_Current = 0.0

    M_list           = zeros(Int, 8001)
    Time_list        = zeros(8001)
    Radical_list     = zeros(Int, 8001)
    Termination_list = zeros(Int, 8001)
    X_n_list         = ones(8001)
    X_w_list         = ones(8001)

    System = Dict{Symbol, Dict}(
        :k => Dict{Symbol,Float64}(
            :I           => k_I,
            :p           => k_p,
            :Termination => k_T,
            :CTAgent     => k_CTAgent,
            :CTADE       => k_CTADE
        ),
        :Initial => Dict{Symbol,Int}(
            :M => M,
            :I => I,
            :CTAgent   => CTA
        ),
        :Critical => Dict{String, Dict}("$(i*10)%" => Dict{String,Vector{Int}}() for i in 1:8),
    )

    Temp_react = [k_I * Initiator, 0]
    length_Radical = length(Radical)

    for i in 1:8
        println("\n$(P_critical[i]*100)% ~ $(P_critical[i+1]*100)%")

        @showprogress for j in 1:1000
            P_threshold = P_critical[i] + P_lists[j]

            while P_Current < P_threshold
                Temp_time = rand(Exponential(1 / sum(Temp_react)))
                Time_Current += Temp_time

                Propagation = rand(Poisson(k_p * length_Radical * Monomer * Temp_time))
                Chain_transfer = rand(Poisson(k_CTAgent * length_Radical * CTA * Temp_time))
                Monomer -= Propagation
                
                Selected_CTA = sample(1:1:Max_count+1, Weights(CTAgent), Chain_transfer, replace=true)  # false
                All_elements = append!(Radical, Selected_CTA.-1)
                # All_elements = [Radical; Selected_CTA.-1]
                Selected_CTA_counts = countmap(Selected_CTA)
                for (idx, count) in Selected_CTA_counts
                    CTAgent[idx] -= count
                end

                Index = sample(1:length(All_elements), Propagation, replace=true)
                Counts = countmap(Index)
                for (idx, count) in Counts
                    All_elements[idx] += count
                end
                
                Radical = All_elements[(Chain_transfer+1):end]
                CTAgent_Back = countmap(All_elements[1:Chain_transfer].+1)
                for (idx, count) in CTAgent_Back
                    CTAgent[idx] += count
                end

                P_Current = 1 - Monomer / M

                Temp_react_class = cumsum(Temp_react) / sum(Temp_react)
                Index_01 = rand()

                if Index_01 < Temp_react_class[1]
                    Initiator -= 1
                    push!(Radical, 0, 0)
                    length_Radical += 2
                    Temp_react = [k_I * Initiator,k_T * (length_Radical)^2 / 2]
                elseif Index_01 < Temp_react_class[2] && length_Radical >= 2
                    selected = sort(sample(1:length_Radical, 2, replace=false))
                    Length = sum(Radical[selected])
                    Termination[Length + 1] += 1
                    deleteat!(Radical, selected)
                    length_Radical -= 2
                    Temp_react[2] = k_T * (length_Radical)^2 / 2
                end
            end

            Index_State = 1000 * (i - 1) + j + 1
            M_list[Index_State] = Monomer
            Time_list[Index_State] = Time_Current
            Radical_list[Index_State] = length(Radical)
            Termination_list[Index_State] = sum(Termination)
            Data = copy(Termination[2:end])
            Counts = countmap(Radical)
            Counts = filter(x -> x[1] != 0, Counts)
            for (idx, count) in Counts
                if idx <= length(Data)
                    Data[idx] += count
                else
                    push!(Data, count)
                end
            end
            Counts = copy(CTAgent[2:end])
            Data += Counts
            X_n_list[Index_State] = sum(Data .* (1:Max_count)) / sum(Data)
            X_w_list[Index_State] = sum(Data .* (1:Max_count).^2) / sum(Data .* (1:Max_count))
        end

        Data = Termination[2:end]
        Counts = countmap(Radical)
        Counts = filter(x -> x[1] != 0, Counts)
        for (idx, count) in Counts
            if idx <= length(Data)
                Data[idx] += count
            else
                push!(Data, count)
            end
        end
        Counts = copy(CTAgent[2:end])
        Data += Counts
        
        System[:Critical]["$(i*10)%"] = Dict(
            # :CTAgent => CTAgent,
            # :Intermediate => Intermediate,
            :Data => Data
        )
    end

    PDI_list = X_w_list ./ X_n_list
    System[:Series] = Dict{Symbol, AbstractVector}(
        :P           => 0:0.0001:0.8,
        :M           => M_list,
        :Time        => Time_list,
        :Radical     => Radical_list,
        :Termination => Termination_list,
        :X_n         => X_n_list,
        :X_w         => X_w_list,
        :PDI         => PDI_list
    )

    # Test Bias
    Data = Termination[2:end]
    Counts = countmap(Radical)
    Counts = filter(x -> x[1] != 0, Counts)
    for (idx, count) in Counts
        if idx <= length(Data)
            Data[idx] += count
        else
            push!(Data, count)
        end
    end
    Counts = copy(CTAgent[2:end])
    Data += Counts
    System[:Bias] = Dict{Symbol,Int64}(
        :M => sum(Data .* (1:Max_count)) + Monomer             - M,
        :I => Initiator + length(Radical)/2 + sum(Termination) - I,
        :C => sum(CTAgent) - CTA
    )

    return System
end

# 平行四边形
function Quasi_Leaping_DT1(M::Int=10^9)
    # M = 10^9
    I = div(M, 100)
    CTA = div(M, 100)

    P_lists = 0.0001:0.0001:0.1
    P_critical = 0:0.1:0.8
    Max_count = 2000

    k_I       = 3.780 * 10^(-5)
    k_p       = 1.039 * 10^( 3) * 5 / M 
    k_T       = 3.160 * 10^( 7) * 5 / M * 2
    k_CTAgent = 1.000 * 10^  5  * 5 / M
    k_CTADE   = 1.000 * 10^  5  * 0

    # Current Stage
    Monomer = M
    Initiator = I
    Radical = Int[]
    CTAgent = zeros(Int, Max_count+1)
    CTAgent[1] = CTA
    # Intermediate = Int[]
    Termination = zeros(Int, Max_count+1)
    P_Current = 0.0
    Time_Current = 0.0

    M_list           = zeros(Int, 8001)
    Time_list        = zeros(8001)
    Radical_list     = zeros(Int, 8001)
    Termination_list = zeros(Int, 8001)
    X_n_list         = ones(8001)
    X_w_list         = ones(8001)

    System = Dict{Symbol, Dict}(
        :k => Dict{Symbol,Float64}(
            :I           => k_I,
            :p           => k_p,
            :Termination => k_T,
            :CTAgent     => k_CTAgent,
            :CTADE       => k_CTADE
        ),
        :Initial => Dict{Symbol,Int}(
            :M => M,
            :I => I,
            :CTAgent   => CTA
        ),
        :Critical => Dict{String, Dict}("$(i*10)%" => Dict{String,Vector{Int}}() for i in 1:8),
    )

    Temp_react = [k_I * Initiator, 0]
    length_Radical = length(Radical)

    for i in 1:8
        println("\n$(P_critical[i]*100)% ~ $(P_critical[i+1]*100)%")

        @showprogress for j in 1:1000
            P_threshold = P_critical[i] + P_lists[j]

            while P_Current < P_threshold
                Temp_time = rand(Exponential(1 / sum(Temp_react)))
                Time_Current += Temp_time

                Propagation = rand(Poisson(k_p * length_Radical * Monomer * Temp_time))
                Chain_transfer = rand(Poisson(k_CTAgent * length_Radical * CTA * Temp_time))
                Monomer -= Propagation
                
                Selected_CTA = sample(1:1:Max_count+1, Weights(CTAgent), Chain_transfer, replace=true)  # false
                All_elements = append!(shuffle!(Radical), Selected_CTA.-1)
                # All_elements = [Radical; Selected_CTA.-1]
                Selected_CTA_counts = countmap(Selected_CTA)
                for (idx, count) in Selected_CTA_counts
                    CTAgent[idx] -= count
                end

                weights::Vector{Int} = compute_difference(length_Radical,Chain_transfer)
                Index = sample(1:length(All_elements), Weights(weights), Propagation, replace=true)
                Counts = countmap(Index)
                for (idx, count) in Counts
                    All_elements[idx] += count
                end
                
                Radical = All_elements[(Chain_transfer+1):end]
                CTAgent_Back = countmap(All_elements[1:Chain_transfer].+1)
                for (idx, count) in CTAgent_Back
                    CTAgent[idx] += count
                end

                P_Current = 1 - Monomer / M

                Temp_react_class = cumsum(Temp_react) / sum(Temp_react)
                Index_01 = rand()

                if Index_01 < Temp_react_class[1]
                    Initiator -= 1
                    push!(Radical, 0, 0)
                    length_Radical += 2
                    Temp_react = [k_I * Initiator,k_T * (length_Radical)^2 / 2]
                elseif Index_01 < Temp_react_class[2] && length_Radical >= 2
                    selected = sort(sample(1:length_Radical, 2, replace=false))
                    Length = sum(Radical[selected])
                    Termination[Length + 1] += 1
                    deleteat!(Radical, selected)
                    length_Radical -= 2
                    Temp_react[2] = k_T * (length_Radical)^2 / 2
                end
            end

            Index_State = 1000 * (i - 1) + j + 1
            M_list[Index_State] = Monomer
            Time_list[Index_State] = Time_Current
            Radical_list[Index_State] = length(Radical)
            Termination_list[Index_State] = sum(Termination)
            Data = copy(Termination[2:end])
            Counts = countmap(Radical)
            Counts = filter(x -> x[1] != 0, Counts)
            for (idx, count) in Counts
                if idx <= length(Data)
                    Data[idx] += count
                else
                    push!(Data, count)
                end
            end
            Counts = copy(CTAgent[2:end])
            Data += Counts
            X_n_list[Index_State] = sum(Data .* (1:Max_count)) / sum(Data)
            X_w_list[Index_State] = sum(Data .* (1:Max_count).^2) / sum(Data .* (1:Max_count))
        end

        Data = Termination[2:end]
        Counts = countmap(Radical)
        Counts = filter(x -> x[1] != 0, Counts)
        for (idx, count) in Counts
            if idx <= length(Data)
                Data[idx] += count
            else
                push!(Data, count)
            end
        end
        Counts = copy(CTAgent[2:end])
        Data += Counts
        
        System[:Critical]["$(i*10)%"] = Dict(
            # :CTAgent => CTAgent,
            # :Intermediate => Intermediate,
            :Data => Data
        )
    end

    PDI_list = X_w_list ./ X_n_list
    System[:Series] = Dict{Symbol, AbstractVector}(
        :P           => 0:0.0001:0.8,
        :M           => M_list,
        :Time        => Time_list,
        :Radical     => Radical_list,
        :Termination => Termination_list,
        :X_n         => X_n_list,
        :X_w         => X_w_list,
        :PDI         => PDI_list
    )

    # Test Bias
    Data = Termination[2:end]
    Counts = countmap(Radical)
    Counts = filter(x -> x[1] != 0, Counts)
    for (idx, count) in Counts
        if idx <= length(Data)
            Data[idx] += count
        else
            push!(Data, count)
        end
    end
    Counts = copy(CTAgent[2:end])
    Data += Counts
    System[:Bias] = Dict{Symbol,Int64}(
        :M => sum(Data .* (1:Max_count)) + Monomer             - M,
        :I => Initiator + length(Radical)/2 + sum(Termination) - I,
        :C => sum(CTAgent) - CTA
    )

    return System
end

# for loop
function Quasi_Leaping_DT2(M::Int=10^9)
    # M = 10^9
    I = div(M, 100)
    CTA = div(M, 100)

    P_lists = 0.0001:0.0001:0.1
    P_critical = 0:0.1:0.8
    Max_count = 2000

    k_I       = 3.780 * 10^(-5)
    k_p       = 1.039 * 10^( 3) * 5 / M 
    k_T       = 3.160 * 10^( 7) * 5 / M * 2
    k_CTAgent = 1.000 * 10^  5  * 5 / M
    k_CTADE   = 1.000 * 10^  5  * 0

    # Current Stage
    Monomer = M
    Initiator = I
    Radical = Int[]
    CTAgent = zeros(Int, Max_count+1)
    CTAgent[1] = CTA
    # Intermediate = Int[]
    Termination = zeros(Int, Max_count+1)
    P_Current = 0.0
    Time_Current = 0.0

    M_list           = zeros(Int, 8001)
    Time_list        = zeros(8001)
    Radical_list     = zeros(Int, 8001)
    Termination_list = zeros(Int, 8001)
    X_n_list         = ones(8001)
    X_w_list         = ones(8001)

    System = Dict{Symbol, Dict}(
        :k => Dict{Symbol,Float64}(
            :I           => k_I,
            :p           => k_p,
            :Termination => k_T,
            :CTAgent     => k_CTAgent,
            :CTADE       => k_CTADE
        ),
        :Initial => Dict{Symbol,Int}(
            :M => M,
            :I => I,
            :CTAgent   => CTA
        ),
        :Critical => Dict{String, Dict}("$(i*10)%" => Dict{String,Vector{Int}}() for i in 1:8),
    )

    Temp_react = [k_I * Initiator, 0]
    length_Radical = length(Radical)

    for i in 1:8
        println("\n$(P_critical[i]*100)% ~ $(P_critical[i+1]*100)%")

        @showprogress for j in 1:1000
            P_threshold = P_critical[i] + P_lists[j]

            while P_Current < P_threshold
                Temp_time = rand(Exponential(1 / sum(Temp_react)))
                Time_Current += Temp_time

                Propagation = rand(Poisson(k_p * length_Radical * Monomer * Temp_time))
                Chain_transfer = rand(Poisson(k_CTAgent * length_Radical * CTA * Temp_time))
                Monomer -= Propagation
                
                seq::Vector{Int} = random_vector_efficient(Propagation, Chain_transfer)
                for k::Int in seq
                    if k == 0
                        Index::Int = rand(1:length_Radical)
                        Radical[Index] += 1
                    else 
                        Index = rand(1:length_Radical)
                        Selected_CTA::Int = sample(1:Max_count+1, Weights(CTAgent))
                        CTAgent[Selected_CTA]     -= 1
                        CTAgent[Radical[Index]+1] += 1
                        Radical[Index] = Radical[end]
                        pop!(Radical)  
                        push!(Radical,Selected_CTA-1)
                    end
                end

                P_Current = 1 - Monomer / M

                Temp_react_class = cumsum(Temp_react) / sum(Temp_react)
                Index_01 = rand()

                if Index_01 < Temp_react_class[1]
                    Initiator -= 1
                    push!(Radical, 0, 0)
                    length_Radical += 2
                    Temp_react = [k_I * Initiator,k_T * (length_Radical)^2 / 2]
                elseif Index_01 < Temp_react_class[2] && length_Radical >= 2
                    selected = sort(sample(1:length_Radical, 2, replace=false))
                    Length = sum(Radical[selected])
                    Termination[Length + 1] += 1
                    deleteat!(Radical, selected)
                    length_Radical -= 2
                    Temp_react[2] = k_T * (length_Radical)^2 / 2
                end
            end

            Index_State = 1000 * (i - 1) + j + 1
            M_list[Index_State] = Monomer
            Time_list[Index_State] = Time_Current
            Radical_list[Index_State] = length(Radical)
            Termination_list[Index_State] = sum(Termination)
            Data = copy(Termination[2:end])
            Counts = countmap(Radical)
            Counts = filter(x -> x[1] != 0, Counts)
            for (idx, count) in Counts
                if idx <= length(Data)
                    Data[idx] += count
                else
                    push!(Data, count)
                end
            end
            Counts = copy(CTAgent[2:end])
            Data += Counts
            X_n_list[Index_State] = sum(Data .* (1:Max_count)) / sum(Data)
            X_w_list[Index_State] = sum(Data .* (1:Max_count).^2) / sum(Data .* (1:Max_count))
        end

        Data = Termination[2:end]
        Counts = countmap(Radical)
        Counts = filter(x -> x[1] != 0, Counts)
        for (idx, count) in Counts
            if idx <= length(Data)
                Data[idx] += count
            else
                push!(Data, count)
            end
        end
        Counts = copy(CTAgent[2:end])
        Data += Counts
        
        System[:Critical]["$(i*10)%"] = Dict(
            # :CTAgent => CTAgent,
            # :Intermediate => Intermediate,
            :Data => Data
        )
    end

    PDI_list = X_w_list ./ X_n_list
    System[:Series] = Dict{Symbol, AbstractVector}(
        :P           => 0:0.0001:0.8,
        :M           => M_list,
        :Time        => Time_list,
        :Radical     => Radical_list,
        :Termination => Termination_list,
        :X_n         => X_n_list,
        :X_w         => X_w_list,
        :PDI         => PDI_list
    )

    # Test Bias
    Data = Termination[2:end]
    Counts = countmap(Radical)
    Counts = filter(x -> x[1] != 0, Counts)
    for (idx, count) in Counts
        if idx <= length(Data)
            Data[idx] += count
        else
            push!(Data, count)
        end
    end
    Counts = copy(CTAgent[2:end])
    Data += Counts
    System[:Bias] = Dict{Symbol,Int64}(
        :M => sum(Data .* (1:Max_count)) + Monomer             - M,
        :I => Initiator + length(Radical)/2 + sum(Termination) - I,
        :C => sum(CTAgent) - CTA
    )

    return System
end

# matrix
function Quasi_Leaping_DT3(M::Int=10^9)
    # M = 10^9
    I = div(M, 100)
    CTA = div(M, 100)

    P_lists = 0.0001:0.0001:0.1
    P_critical = 0:0.1:0.8
    Max_count = 2000

    k_I       = 3.780 * 10^(-5)
    k_p       = 1.039 * 10^( 3) * 5 / M 
    k_T       = 3.160 * 10^( 7) * 5 / M * 2
    k_CTAgent = 1.000 * 10^  5  * 5 / M
    k_CTADE   = 1.000 * 10^  5  * 0

    # Current Stage
    Monomer = M
    Initiator = I
    Radical = Int[]
    CTAgent = zeros(Int, Max_count+1)
    CTAgent[1] = CTA
    # Intermediate = Int[]
    Termination = zeros(Int, Max_count+1)
    P_Current = 0.0
    Time_Current = 0.0

    M_list           = zeros(Int, 8001)
    Time_list        = zeros(8001)
    Radical_list     = zeros(Int, 8001)
    Termination_list = zeros(Int, 8001)
    X_n_list         = ones(8001)
    X_w_list         = ones(8001)

    System = Dict{Symbol, Dict}(
        :k => Dict{Symbol,Float64}(
            :I           => k_I,
            :p           => k_p,
            :Termination => k_T,
            :CTAgent     => k_CTAgent,
            :CTADE       => k_CTADE
        ),
        :Initial => Dict{Symbol,Int}(
            :M => M,
            :I => I,
            :CTAgent   => CTA
        ),
        :Critical => Dict{String, Dict}("$(i*10)%" => Dict{String,Vector{Int}}() for i in 1:8),
    )

    Temp_react = [k_I * Initiator, 0]
    length_Radical = length(Radical)

    for i in 1:8
        println("\n$(P_critical[i]*100)% ~ $(P_critical[i+1]*100)%")

        @showprogress for j in 1:1000
            P_threshold = P_critical[i] + P_lists[j]

            while P_Current < P_threshold
                Temp_time = rand(Exponential(1 / sum(Temp_react)))
                Time_Current += Temp_time

                Propagation = rand(Poisson(k_p * length_Radical * Monomer * Temp_time))
                Chain_transfer = rand(Poisson(k_CTAgent * length_Radical * CTA * Temp_time))
                Monomer -= Propagation
                
                Selected_CTA = sample(1:1:Max_count+1, Weights(CTAgent), Chain_transfer, replace=true)  # false
                All_elements = append!(shuffle!(Radical), Selected_CTA.-1)
                Selected_CTA_counts = countmap(Selected_CTA)
                for (idx, count) in Selected_CTA_counts
                    CTAgent[idx] -= count
                end

                seq = sample(1:length_Radical,Chain_transfer,replace=true)
                matrix=generate_matrix(length_Radical, Chain_transfer+1, seq)
                matrix[:, 1:end-1] .= ifelse.(matrix[:, 1:end-1] .< matrix[:, 2:end], 0, matrix[:, 1:end-1])
                Radicalx = extract_and_zero!(matrix)
                # matrix
                remainx = findall(x -> x == Chain_transfer+1, Radicalx)
                matrix = matrix[:, end:-1:1]
                CTAx = vec(matrix[matrix.!=0])

                Index = sample(1:length(All_elements), Weights(vcat(Radicalx,CTAx)), Propagation, replace=true)
                Counts = countmap(Index)
                for (idx, count) in Counts
                    All_elements[idx] += count
                end

                Index = [remainx; (length_Radical+1):(2*length_Radical-length(remainx))]
                Radical = All_elements[Index]
                CTAgent_Back = countmap(All_elements[setdiff(1:length(All_elements), Index)].+1)
                for (idx, count) in CTAgent_Back
                    CTAgent[idx] += count
                end
                
                P_Current = 1 - Monomer / M

                Temp_react_class = cumsum(Temp_react) / sum(Temp_react)
                Index_01 = rand()

                if Index_01 < Temp_react_class[1]
                    Initiator -= 1
                    push!(Radical, 0, 0)
                    length_Radical += 2
                    Temp_react = [k_I * Initiator,k_T * (length_Radical)^2 / 2]
                elseif Index_01 < Temp_react_class[2] && length_Radical >= 2
                    selected = sort(sample(1:length_Radical, 2, replace=false))
                    Length = sum(Radical[selected])
                    Termination[Length + 1] += 1
                    deleteat!(Radical, selected)
                    length_Radical -= 2
                    Temp_react[2] = k_T * (length_Radical)^2 / 2
                end
            end

            Index_State = 1000 * (i - 1) + j + 1
            M_list[Index_State] = Monomer
            Time_list[Index_State] = Time_Current
            Radical_list[Index_State] = length(Radical)
            Termination_list[Index_State] = sum(Termination)
            Data = copy(Termination[2:end])
            Counts = countmap(Radical)
            Counts = filter(x -> x[1] != 0, Counts)
            for (idx, count) in Counts
                if idx <= length(Data)
                    Data[idx] += count
                else
                    push!(Data, count)
                end
            end
            Counts = copy(CTAgent[2:end])
            Data += Counts
            X_n_list[Index_State] = sum(Data .* (1:Max_count)) / sum(Data)
            X_w_list[Index_State] = sum(Data .* (1:Max_count).^2) / sum(Data .* (1:Max_count))
        end

        Data = Termination[2:end]
        Counts = countmap(Radical)
        Counts = filter(x -> x[1] != 0, Counts)
        for (idx, count) in Counts
            if idx <= length(Data)
                Data[idx] += count
            else
                push!(Data, count)
            end
        end
        Counts = copy(CTAgent[2:end])
        Data += Counts
        
        System[:Critical]["$(i*10)%"] = Dict(
            # :CTAgent => CTAgent,
            # :Intermediate => Intermediate,
            :Data => Data
        )
    end

    PDI_list = X_w_list ./ X_n_list
    System[:Series] = Dict{Symbol, AbstractVector}(
        :P           => 0:0.0001:0.8,
        :M           => M_list,
        :Time        => Time_list,
        :Radical     => Radical_list,
        :Termination => Termination_list,
        :X_n         => X_n_list,
        :X_w         => X_w_list,
        :PDI         => PDI_list
    )

    # Test Bias
    Data = Termination[2:end]
    Counts = countmap(Radical)
    Counts = filter(x -> x[1] != 0, Counts)
    for (idx, count) in Counts
        if idx <= length(Data)
            Data[idx] += count
        else
            push!(Data, count)
        end
    end
    Counts = copy(CTAgent[2:end])
    Data += Counts
    System[:Bias] = Dict{Symbol,Int64}(
        :M => sum(Data .* (1:Max_count)) + Monomer             - M,
        :I => Initiator + length(Radical)/2 + sum(Termination) - I,
        :C => sum(CTAgent) - CTA
    )

    return System
end

function SSA_ATRP(M::Int=10^9)
    # M = 10^9
    I = div(M, 100)

    # CTA = M/100

    P_lists = 0.0001:0.0001:0.1
    P_critical = 0:0.1:0.8
    Max_count = 6001
    
    k_act   = 1.000 * 10^( 0) * 5 / M 
    k_deact = 1.000 * 10^( 7) * 5 / M
    k_p     = 1.039 * 10^( 3) * 5 / M 
    k_T     = 3.160 * 10^( 7) * 5 / M * 2

    # Current Stage
    Monomer = M
    CuI     = I
    CuII    = 0
    R_Xn    = I
    R_X     = [I;zeros(Int, Max_count-1)]
    Radical = Int[]
    length_Radical = 0

    Termination = zeros(Int, Max_count)
    P_Current = 0.0
    Time_Current = 0.0

    M_list           = zeros(Int, 8001)
    Time_list        = zeros(8001)
    Radical_list     = zeros(Int, 8001)
    Termination_list = zeros(Int, 8001)
    X_n_list         = ones(8001)
    X_w_list         = ones(8001)

    System = Dict{Symbol, Dict}(
        :k => Dict{Symbol,Float64}(
            :act         => k_act,
            :deact       => k_deact,
            :p           => k_p,
            :Termination => k_T,
        ),
        :Initial => Dict(
            :M => M,
            :I => I
        ),
        :Critical => Dict{String, Dict}("$(i*10)%" => Dict() for i in 1:8)
    )
    Temp_react = [k_act * CuI * R_Xn, 0, 0, 0]
    TTT = [[0.0,0.0,0.0,0.0],missing,missing,missing,missing,missing,missing,missing,missing]

    for i in 1:8
        println("\n$(P_critical[i]*100)% ~ $(P_critical[i+1]*100)%")

        @showprogress for j in 1:1000
            P_threshold = P_critical[i] + P_lists[j]

            while P_Current < P_threshold
                Temp_time = rand(Exponential(1 / sum(Temp_react)))
                Time_Current += Temp_time
                Temp_react_class = cumsum(Temp_react) / sum(Temp_react)
                Index_01 = rand()

                if Index_01 < Temp_react_class[1]
                    NewR    = sample(1:1:Max_count, Weights(R_X), 1, replace=true)[1]
                    Radical = push!(Radical,NewR-1)
                    R_X[NewR] -= 1

                    CuI            -= 1
                    CuII           += 1
                    R_Xn           -= 1
                    length_Radical += 1

                    Temp_react = [k_act*CuI*R_Xn, k_deact*CuII*length_Radical, k_p*length_Radical*Monomer, k_T*(length_Radical)^2/2]
                elseif Index_01 < Temp_react_class[2]
                    Index = rand(1:length_Radical)
                    R_X[Radical[Index]+1]   += 1
                    Radical = Radical[setdiff(1:length_Radical, Index)]

                    CuI            += 1
                    CuII           -= 1
                    R_Xn           += 1
                    length_Radical -= 1

                    Temp_react = [k_act*CuI*R_Xn, k_deact*CuII*length_Radical, k_p*length_Radical*Monomer, k_T*(length_Radical)^2/2]
                elseif Index_01 < Temp_react_class[3]
                    Monomer -= 1
                    Index = rand(1:length_Radical)
                    Radical[Index] += 1
                    P_Current = 1 - Monomer / M

                    Temp_react[3] = k_p*length_Radical*Monomer
                elseif Index_01 < Temp_react_class[4] && length_Radical >= 2
                    selected = sample(1:length_Radical, 2, replace=false)
                    Length = sum(Radical[selected]) + 1
                    Termination[Length] += 1
                    Radical = Radical[setdiff(1:length(Radical), selected)]
                    length_Radical -= 2

                    Temp_react[[2,3,4]] = [k_deact*CuII*length_Radical, k_p*length_Radical*Monomer, k_T*(length_Radical)^2/2]
                end
            end

            Index_State = 1000 * (i - 1) + j + 1
            M_list[Index_State] = Monomer
            Time_list[Index_State] = Time_Current
            Radical_list[Index_State] = length(Radical)
            Termination_list[Index_State] = sum(Termination)
            Data = copy(Termination[2:end])
            Counts = countmap(Radical)
            Counts = filter(x -> x[1] != 0, Counts)
            for (idx, count) in Counts
                if idx <= length(Data)
                    Data[idx] += count
                else
                    push!(Data, count)
                end
            end
            X_n_list[Index_State] = sum(Data .* (1:Max_count-1)) / sum(Data)
            X_w_list[Index_State] = sum(Data .* (1:Max_count-1).^2) / sum(Data .* (1:Max_count-1))
        end
        TTT[i+1] = [k_act*CuI*R_Xn, k_deact*CuII*length_Radical, k_p*length_Radical*Monomer, k_T*(length_Radical)^2/2]

        Data = copy(Termination[2:end])
        Counts = countmap(Radical)
        Counts = filter(x -> x[1] != 0, Counts)
        for (idx, count) in Counts
            if idx <= length(Data)
                Data[idx] += count
            else
                push!(Data, count)
            end
        end

        System[:Critical]["$(i*10)%"] = Dict(
            # :CTAgent => CTAgent,
            # :Intermediate => Intermediate,
            :Data => Data
        )
    end

    PDI_list = X_w_list ./ X_n_list
    System[:Series] = Dict{Symbol, AbstractVector}(
        :P           => 0:0.0001:0.8,
        :M           => M_list,
        :Time        => Time_list,
        :Radical     => Radical_list,
        :Termination => Termination_list,
        :X_n         => X_n_list,
        :X_w         => X_w_list,
        :PDI         => PDI_list
    )

    Data = copy(Termination[2:end])
    Data += R_X[2:end]
    Counts = countmap(Radical)
    Counts = filter(x -> x[1] != 0, Counts)
    for (idx, count) in Counts
        if idx <= length(Data)
            Data[idx] += count
        else
            push!(Data, count)
        end
    end

    System[:Bias] = Dict{Symbol,Int64}(
        :M => sum(Data .* (1:Max_count-1)) + Monomer - M,
        :I => R_Xn + length(Radical) + 2*sum(Termination) - I
    )

    return System, TTT
end

# @benchmark Quasi_Leaping_FRP() samples=10

@time Data_FRP1 = SSA_FRP()
@time Data_FRP2 = Quasi_Leaping_FRP()
@time Data_DT1 = SSA_DT()
@time Data_DT2 = Quasi_Leaping_DT()  #jx
@time Data_DT3 = Quasi_Leaping_DT1() #pxsbx
@time Data_DT4 = Quasi_Leaping_DT2() #for


#Data_Plot(Data_FRP1)
#Data_Plot(Data_FRP2)

# Chi
function find_cutoff(freq_vector, percentile=0.95)
    total = sum(freq_vector)
    cumsum_ratio = cumsum(freq_vector) ./ total
    findfirst(x -> x >= percentile, cumsum_ratio)
end


FRP_80_SSA = Data_FRP1[:Critical]["80%"][:Data]
FRP_80_NEW = Data_FRP2[:Critical]["80%"][:Data]
cutoff_SSA = find_cutoff(FRP_80_SSA, 0.95)
cutoff_NEW = find_cutoff(FRP_80_NEW, 0.95)
cutoff = max(cutoff_SSA,cutoff_NEW)
FRP_80_SSA_trunc = FRP_80_SSA[1:cutoff]
FRP_80_NEW_trunc = FRP_80_NEW[1:cutoff]
contingency_table = vcat(FRP_80_SSA_trunc', FRP_80_NEW_trunc')
chi2_test = ChisqTest(contingency_table)



@time Data_DT1 = SSA_DT()

@time Data_DT2 = Quasi_Leaping_DT3()


contingency_table = vcat(fx', fy')
test_result = ChisqTest(contingency_table)



@time Data_DT1 = SSA_DT()
@time Data_DT2 = Quasi_Leaping_DT2()

DT_80_SSA = Data_DT1[:Critical]["80%"][:Data]
DT_80_NEW = Data_DT2[:Critical]["80%"][:Data]

cutoff_SSA = find_cutoff(DT_80_SSA, 0.95)
cutoff_NEW = find_cutoff(DT_80_NEW, 0.95)
cutoff = max(cutoff_SSA,cutoff_NEW)
DT_80_SSA_trunc = DT_80_SSA[1:cutoff]   
DT_80_NEW_trunc = DT_80_NEW[1:cutoff]
contingency_table_1 = vcat(DT_80_SSA_trunc', DT_80_NEW_trunc')
chi2_test = ChisqTest(contingency_table_1)


@time Data_DT2 = Quasi_Leaping_DT1()
DT_80_NEW = Data_DT2[:Critical]["80%"][:Data]
cutoff_NEW = find_cutoff(DT_80_NEW, 0.5)
cutoff = max(cutoff_SSA,cutoff_NEW)
DT_80_SSA_trunc = DT_80_SSA[1:cutoff]
DT_80_NEW_trunc = DT_80_NEW[1:cutoff]
contingency_table_1 = vcat(DT_80_SSA_trunc', DT_80_NEW_trunc')
chi2_test = ChisqTest(contingency_table_1)

@time Data_DT2 = Quasi_Leaping_DT()
DT_80_NEW = Data_DT2[:Critical]["80%"][:Data]
cutoff_NEW = find_cutoff(DT_80_NEW, 0.95)
cutoff = max(cutoff_SSA,cutoff_NEW)
DT_80_SSA_trunc = DT_80_SSA[1:cutoff]
DT_80_NEW_trunc = DT_80_NEW[1:cutoff]
contingency_table_1 = vcat(DT_80_SSA_trunc', DT_80_NEW_trunc')
chi2_test = ChisqTest(contingency_table_1)

plot(FRP_80_SSA_trunc-FRP_80_NEW_trunc)
plot(DT_80_SSA-DT_80_NEW)

@time Data_DT1 = SSA_DT()
@time Data_DT2 = Quasi_Leaping_DT2()
@time Data_DT3 = Quasi_Leaping_DT3()
DT_80_SSA = Data_DT1[:Critical]["80%"][:Data]
DT_80_NEW = Data_DT2[:Critical]["80%"][:Data]

cutoff_SSA = find_cutoff(DT_80_SSA, 0.95)
cutoff_NEW = find_cutoff(DT_80_NEW, 0.95)
cutoff = max(cutoff_SSA,cutoff_NEW)
DT_80_SSA_trunc = DT_80_SSA[1:cutoff]
DT_80_NEW_trunc = DT_80_NEW[1:cutoff]
contingency_table_1 = vcat(DT_80_SSA_trunc', DT_80_NEW_trunc')
chi2_test = ChisqTest(contingency_table_1)

function Chi_Square_Test(x, y)
    Data_x = x[:Critical]["80%"][:Data]
    Data_y = y[:Critical]["80%"][:Data]
    cutoff_x = find_cutoff(Data_x, 0.95)
    cutoff_y = find_cutoff(Data_y, 0.95)
    cutoff = max(cutoff_x, cutoff_y)
    Data_x_trunc = Data_x[1:cutoff] 
    Data_y_trunc = Data_y[1:cutoff]
    contingency_table = vcat(Data_x_trunc', Data_y_trunc')
    chi2_test = ChisqTest(contingency_table)

    p = Data_x_trunc ./ sum(Data_x_trunc)
    q = Data_y_trunc ./ sum(Data_y_trunc)
    m = 0.5 .* (p .+ q)
    kl_pm = sum(p .* log.(p ./ m)) / log(2)
    kl_qm = sum(q .* log.(q ./ m)) / log(2)
    js_divergence = 0.5 * (kl_pm + kl_qm)

    #return (chi2_test=chi2_test, js_divergence=js_divergence)
    return(Data_x_trunc-Data_y_trunc)
end

