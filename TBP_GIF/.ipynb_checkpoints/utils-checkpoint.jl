include("TBP.jl");

function Irudikatu(z0, t0, T, dtau, par)
    u0 = HasierakoEgoera(z0)
    U0 = LCinvFcn(u0)
    tt, UU = EMI_Irudikatzeko(U0, t0, T, dtau, par)
    pl = plot(title = "Hiru gorputzen ibilbideak [0, $T] tartean",
                  xlabel = "x", ylabel = "y", aspect_ratio=1)#, xlims=(-0.01,0.01), ylims=(-0.01,0.01))
    
    solu = LCFcn.(UU)
    
    for j in 1:3
        xx = [(u[2*(mod(j-2,3)+1)-1] - u[2j-1])/3 for u in solu]
        yy = [(u[2*(mod(j-2,3)+1)] - u[2j])/3 for u in solu]
        plot!(xx,yy,legend=false)
    end
    u0 = LCFcn(U0)
    
    # Hasierako eta bukaerako triangeluak marraztu
    xx = [(u0[2*(mod(j-2,3)+1)-1] - u0[2j-1])/3 for j=[1,2,3,1]]
    yy = [(u0[2*(mod(j-2,3)+1)] - u0[2j])/3  for j=[1,2,3,1]]
    plot!(xx,yy,color="gray")
    scatter!(xx,yy,color="gray")
    Uend = UU[end]
    uend = LCFcn(Uend)
    xx = [(uend[2*(mod(j-2,3)+1)-1] - uend[2j-1])/3 for j=[1,2,3,1]]
    yy = [(uend[2*(mod(j-2,3)+1)] - uend[2j])/3 for j=[1,2,3,1]]
    scatter!(xx,yy,color="purple")
    plot!(xx,yy,color="purple")
    
    display(pl)
    return tt, UU, solu
end

function ItzuliIrudikatu(U0, t0, T, dtau, par, xlims, ylims, background_color, size, color)
    tt, UU = EMI_Irudikatzeko(U0, t0, T, dtau, par)
    pl = plot(aspect_ratio=1, xlims = xlims, ylims = ylims, background_color=background_color, framestyle=:none, size = size)    
    uu = LCFcn.(UU)
    
    for j in 1:3
        xx = [(u[2*(mod(j-2,3)+1)-1] - u[2j-1])/3 for u in uu]
        yy = [(u[2*(mod(j-2,3)+1)] - u[2j])/3 for u in uu]
        plot!(pl, xx,yy,legend=false, linewidth = 3, color = color[j])
    end
        
    return pl
end



function GifErregularizatuGabe(z0, uu, tt, dtau, par; background_color, color, elapsed_time = 0.08)
    u0 = HasierakoEgoera(z0)
    U0 = LCinvFcn(u0)
    uuabs = rel2abs.(uu)
    
    # Eskuratu balioak
    x1 = [u[1] for u in uuabs]; y1 = [u[2] for u in uuabs]
    x2 = [u[3] for u in uuabs]; y2 = [u[4] for u in uuabs]
    x3 = [u[5] for u in uuabs]; y3 = [u[6] for u in uuabs]
    real_tt = [u[13] for u in uu]  # Denbora erreala
    tt_ind = Int64[]
    t_ = real_tt[1]
    push!(tt_ind, 1)
    for i in 2:size(real_tt)[1]
        if real_tt[i]-t_>=elapsed_time
            push!(tt_ind, i)
            t_ = real_tt[i]
        end
    end
    
    # Gorputzen araberako ardatz mugak
    xlims = extrema(vcat(x1, x2, x3))
    xlims = xlims .* 1.05
    ylims = extrema(vcat(y1, y2, y3))
    ylims = ylims .* 1.05

    xrange = Int64(ceil(xlims[2] - xlims[1]))
    yrange = Int64(ceil(ylims[2] - ylims[1]))
    
    # Sortuko den irudiaren tamaina (px). Zehaztu altuera eta zabalera 1:1 ratioan lortuko da
    height = 800 
    width = Int64(ceil(height * (xrange / yrange)))
    
    # Animazioaren guztizko denbora 
    total_time = real_tt[end] - real_tt[1]
    
    # Frame kopuru totala
    total_frames = length(tt_ind)
    
    # FPS kalkulatu
    fps = total_frames / total_time
        
    # ANIMAZIOAK SORTU:
        #1. Gorputzen ibilbidea gorputza mugitu ahala sortuko da
    anim1 = @animate for i in tt_ind
        plt = ItzuliIrudikatu(U0, 0., tt[i], dtau, par, xlims, ylims,background_color, (width, height), color)
        scatter!(plt, [x1[i], x2[i], x3[i]],
                [y1[i], y2[i], y3[i]],
                legend = false, grid = false,
                aspect_ratio = 1, 
                xlims = xlims, ylims = ylims, 
                background_color= background_color, 
                framestyle=:none, 
                size = (width, height), 
                markersize = [8,8,8], color = color)
        xaxis!(plt, false)
        yaxis!(plt, false)
    end
        #1. Gorputzen ibilbidea marraztuta egongo da uneoro eta gorputzak gainean mugituko dira
    anim2 = @animate for i in tt_ind
        plt = ItzuliIrudikatu(U0, 0., tt[tt_ind[end]], dtau, par, xlims, ylims, background_color, (width, height), color)
        scatter!(plt, [x1[i], x2[i], x3[i]],
                [y1[i], y2[i], y3[i]],
                legend = false, grid = false,
                aspect_ratio = 1, 
                xlims = xlims, ylims = ylims, 
                background_color= background_color, 
                framestyle=:none, 
                size = (width, height), 
                markersize = [8,8,8], color = color)
        xaxis!(plt, false)
        yaxis!(plt, false)
    end
    return anim1, anim2, fps, total_frames
end



function CSVIrakurri(filepath)
    data = CSV.File(filepath; skipto=2) |> Tables.matrix
    
    z0_list = Vector{Vector{Float64}}()
    Zopt_list = Vector{Vector{Float64}}()
    dtau_list = Float64[]
    sigma_list = Vector{NTuple{3, Int64}}()
    
    for i in 1:size(data,1)
        row = data[i, :]
        push!(z0_list, row[1:4])
        push!(Zopt_list, row[5:10])
        push!(dtau_list, row[11])
        push!(sigma_list, (row[12], row[13], row[14]))
    end
    return z0_list, Zopt_list, dtau_list, sigma_list
end