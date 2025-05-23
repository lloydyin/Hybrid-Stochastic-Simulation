# Direct Use
using Plots, JLD2
data = load("mydata.jld2")
pyplot()

function Data_Plot(Sys,xlimit)
    color_map = Dict(
        80 => :purple,
        70 => :blue,
        60 => :cyan,
        50 => :green,
        40 => :yellow,
        30 => :orange,
        20 => :red,
        10 => :brown
    )
    # 第一张大图 CLD
    p1 = plot(title="CLD", legend=:topleft)
    for i in 80:-10:10
        ydata = Sys[:Critical]["$(i)%"][:Data][1:xlimit]
        xdata = 1:xlimit
        plot!(p1, xdata, ydata, label="$(i)%", color=color_map[i])
    end
    xlabel!(p1, "Chain Length")
    ylabel!(p1, "Count")
    xflip!(p1)
    # 第二张图 GPC
    p2 = plot(title="GPC", legend=false)
    k1, k2 = 1, 1
    Gap = 2
    for i in 80:-10:10
        Weight = Sys[:Critical]["$(i)%"][:Data]
        Max_MW = findlast(x -> x != 0, Weight)
        if Max_MW === nothing
            continue
        end
        Breaks = 0:Gap:Max_MW
        function cut(values, breaks; extend=false)
            breaks = sort(breaks)
            if extend
                breaks = [minimum([minimum(values); breaks]) - eps(), breaks..., maximum([maximum(values); breaks]) + eps()]
            end
            intervals = similar(values, Int)
            for idx in eachindex(values)
                intervals[idx] = findfirst(x -> values[idx] <= x, breaks) - 1
            end
            return intervals
        end
        Intervals = cut(1:Max_MW, Breaks, extend=true)
        Temp = Weight[1:Max_MW] .* (1:Max_MW)
        Weight_intervals = [sum(Temp[findall(x->x==interval, Intervals)]) for interval in unique(Intervals)]
        Numbers = [sum(Weight[1:Max_MW][findall(x->x==interval, Intervals)]) for interval in unique(Intervals)]
        nonzero_idx = findall(x->x!=0, Weight_intervals)
        if isempty(nonzero_idx)
            continue
        end
        GPC_Result_x = log10.(Weight_intervals[nonzero_idx] ./ Numbers[nonzero_idx]) .* k1 .+ k2
        GPC_Result_y = (-0.4228 .* GPC_Result_x .+ 10.38) .* Weight_intervals[nonzero_idx]
        
        plot!(p2, GPC_Result_x, GPC_Result_y, label="$(i)%", color=color_map[i])
    end
    xlabel!(p2, "RT")
    ylabel!(p2, "DR")
    # 第三张图 PDI
    p3 = plot(title="PDI", legend=false)
    plot!(p3, Sys[:Series][:P], Sys[:Series][:PDI], color=:black, label="X_n")
    xlabel!(p3, "Conv.(%)")
    ylabel!(p3, "PDI")
    # 布局定义，第一行占2列，第二行分成2列
    layout = @layout [a{0.6h}; b c]
    plt = plot(p1, p2, p3, layout=layout, size=(800, 450))
    return plt
end

default(
    guidefont = font("Times New Roman", 13), 
    tickfont  = font("Times New Roman", 10), 
    titlefont = font("Times New Roman", 17), 
    legendfont = font("Times New Roman", 8) 
)

#p10 = Data_Plot(data["Data_FRP1"],3000)
#p11 = Data_Plot(data["Data_FRP2"],3000)
#p20 = Data_Plot(data["Data_DT1"],200)
#p21 = Data_Plot(data["Data_DT2"],200) # for loop
#p22 = Data_Plot(data["Data_DT3"],200) # Parallelogram
#p23 = Data_Plot(data["Data_DT4"],200) # Rectangle

savefig(Data_Plot(data["Data_FRP1"],3000), "Plot_FRP_SSA.pdf")
savefig(Data_Plot(data["Data_FRP2"],3000), "Plot_FRP_New.pdf")
savefig(Data_Plot(data["Data_DT1"],200), "Plot_DT_SSA.pdf")
savefig(Data_Plot(data["Data_DT2"],200), "Plot_DT_New.pdf")
savefig(Data_Plot(data["Data_DT3"],200), "Plot_DT_New1.pdf")
savefig(Data_Plot(data["Data_DT4"],200), "Plot_DT_New2.pdf")
