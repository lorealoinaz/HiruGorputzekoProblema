include("TBP.jl");

"""
    displayPlotTBPState(u0; background_color = :white, color = ["#9d6dc2","#52b788", "#277da1"], markersize = [15])

Hiru gorputzen problemaren u0 egoerako bateko gorputzak eta haren abiadura bektoreak irudikatzen ditu (ardatzik gabe), background_color parametroak adierazitako atzeko kolorean, gorputz bakoitzak color taulako kolore bat izanik eta markersize tamaina.
"""
function displayPlotTBPState(u0; background_color = :white, color = ["#9d6dc2","#52b788", "#277da1"], markersize = [15])
    x1,y1,x2,y2,x3,y3 = u0
    pl = scatter(aspect_ratio = 1, background_color = background_color, xaxis = false, yaxis = false, grid = false)
    
    #Hasierako eta amaierako abiadura bektoreak marraztu
    positions0 = [(x1, y1), (x2, y2), (x3, y3)]
    velocities0 = [(u0[7], u0[8]), (u0[9], u0[10]), (u0[11], u0[12])]
    
    for i in 1:3
        x, y = positions0[i]
        vx, vy = velocities0[i]
        quiver!([x], [y], quiver=([vx], [vy]), arrow=true, color=:white, label=false)
    end;
    
    scatter!(pl, [x1],[y1], legend = false, markersize = markersize, color = color[1])
    scatter!(pl, [x2],[y2], legend = false, markersize = markersize, color = color[2])
    scatter!(pl, [x3],[y3], legend = false, markersize = markersize, color = color[3])
    return pl
end



function Irudikatu(z0, t0, T, dtau, par; yrange = nothing, triangeluak = true)
    u0 = HasierakoEgoera(z0)
    U0 = LCinvFcn(u0)
    tt, UU = EMI_Irudikatzeko(U0, t0, T, dtau, par)
    if yrange!=nothing
        if triangeluak
             pl = plot(xlabel = "x", ylabel = "y", yrange = yrange, title = "Hiru gorputzen ibilbideak [0, $T] tartean")
        else
            pl = plot(yrange = yrange)
        end
    else
        if triangeluak
            pl = plot(xlabel = "x", ylabel = "y", aspect_ratio=1, title = "Hiru gorputzen ibilbideak [0, $T] tartean")
        else
            pl = plot(aspect_ratio = 1)
        end
    end
    solu = LCFcn.(UU)

    for j in 1:3
        xx = [(u[2*(mod(j-2,3)+1)-1] - u[2j-1])/3 for u in solu]
        yy = [(u[2*(mod(j-2,3)+1)] - u[2j])/3 for u in solu]
        plot!(xx,yy,legend=false)
    end
    u0 = LCFcn(U0)
    
    if triangeluak
        # Hasierako eta bukaerako triangeluak marraztu
        xx = [(u0[2*(mod(j-2,3)+1)-1] - u0[2j-1])/3 for j=[1,2,3,1]]
        yy = [(u0[2*(mod(j-2,3)+1)] - u0[2j])/3  for j=[1,2,3,1]]
        plot!(xx,yy,color=:gray)
        scatter!(xx,yy,color=:gray)
        Uend = UU[end]
        uend = LCFcn(Uend)
        xx = [(uend[2*(mod(j-2,3)+1)-1] - uend[2j-1])/3 for j=[1,2,3,1]]
        yy = [(uend[2*(mod(j-2,3)+1)] - uend[2j])/3 for j=[1,2,3,1]]
        scatter!(xx,yy,color=:purple)
        plot!(xx,yy,color=:purple)
    
        #Hasierako eta amaierako abiadura bektoreak marraztu
        absu0 = rel2abs(solu[1])
        positions0 = [(absu0[1], absu0[2]), (absu0[3], absu0[4]), (absu0[5], absu0[6])]
        velocities0 = [(u0[7], u0[8]), (u0[9], u0[10]), (u0[11], u0[12])]
    
        absuend = rel2abs(solu[end])
        positionsend = [(absuend[1], absuend[2]), (absuend[3], absuend[4]), (absuend[5], absuend[6])]
        velocitiesend = [(absuend[7], absuend[8]), (absuend[9], absuend[10]), (absuend[11], absuend[12])]
        for i in 1:3
            x, y = positions0[i]
            vx, vy = velocities0[i]
            quiver!([x], [y], quiver=([vx], [vy]), arrow=true, color=:gray, label=false)
            x, y = positionsend[i]
            vx, vy = velocitiesend[i]
            quiver!([x], [y], quiver=([vx], [vy]), arrow=true, color=:purple, label=false)
        end
    end
    if triangeluak
        display(pl)
        return tt, UU, solu
    else
        return pl
    end
end




function IrudikatuIbilbidea(z0, t0, T, dtau, par; background = :black, color = nothing, yrange = nothing, linewidth = 3)
    u0 = HasierakoEgoera(z0)
    U0 = LCinvFcn(u0)
    tt, UU = EMI_Irudikatzeko(U0, t0, T, dtau, par)
    solu = LCFcn.(UU)

    # Oinarrizko grafikoa sortu
    if yrange != nothing
        pl = plot(; yrange = yrange)
    else
        pl = plot(; aspect_ratio = 1)
    end
    plot!(pl, legend = false, 
            framestyle = :none,
            grid = false,
            ticks = false,
            background_color_inside = background,
            dpi=300)
    
    if color ==  :white
        color = [:white, :white, :white]
    end
    for j in 1:3
        xx = [(u[2*(mod(j-2,3)+1)-1] - u[2j-1])/3 for u in solu]
        yy = [(u[2*(mod(j-2,3)+1)] - u[2j])/3 for u in solu]
        if color ==nothing
            plot!(xx, yy, linewidth = linewidth)
        else
            plot!(xx, yy, color = color[j], linewidth = linewidth)
        end
    end

    return pl
end


function IrudikatuIELog(tt, UU, uu; potentziala = false)
    d = size(UU)[1]
    II = zeros(d)
    log_rmin = zeros(d)
    logII = zeros(d)
    logratio = zeros(d)
    EE_kanpo = zeros(d)
    EE_barne = zeros(d)
    EE_dif = zeros(d)
    
    for i in eachindex(UU)
        I, rmin, r_ratio, E_barne, E_kanpo = IELog(UU[i], uu[i])
        II[i] = I
        logII[i] = log10(I)
        log_rmin[i] = log10(rmin)
        logratio[i] = log10(r_ratio)
        EE_barne[i] = E_barne
        EE_kanpo[i] = E_kanpo
        EE_dif[i] = E_barne + E_kanpo
    end

    pl1 = plot(tt, II, title = "Inertzia momentua", legend=false)

    pl2 = plot(tt, log_rmin, title = "Log rmin, I eta r_ratio", legend=false)
    plot!(pl2, tt, logII, legend=false)
    plot!(pl2, tt, logratio, legend=false)

    pl3 = plot(tt, EE_kanpo, title = "Kanpo eta barne energia", legend=false)
    plot!(pl3, tt, EE_barne)
    plot!(pl3, tt, EE_dif)

    if potentziala
        PP = [-ThreeBodyRelPotential(LCFcn(U)) for U in UU]
        pl4 = plot(tt, PP, title = "Energia potentzialaren aurkakoa", legend=false)
    
        pl = plot(pl1,pl2,pl3,pl4, layout = (2,2),size = (900, 500))
    else
        pl = plot(pl1,pl2,pl3, layout = (1,3),size = (1200, 300))
    end 
    display(pl)
end


function IrudikatuErrEH(tt, UU, solu)
    yrange = (-20, 0)
    energia_erroreakTP = [ThreeBodyRelEnergy(u)/h - 1  for u in solu]
    energia_errore_lokalakTP = energia_erroreakTP[2:end]-energia_erroreakTP[1:end-1]
    
    pl1 = plot(tt, log10.(abs.(energia_erroreakTP)), legend = false,  title = "Energia errorea TwicePrecision", ylims=yrange)
    pl2 = plot(tt[2:end], log10.(abs.(energia_errore_lokalakTP)), legend = false, title = "Energia errore lokala TwicePrecision", ylims=yrange)
    
    
    plot(pl1,pl2, size = (1000, 300))
    
    yrange = (-30,-1)
    
    Ham_erroreakTP = [ThreeBodyHamiltonianLC(U,h)  for U in UU]
    Ham_errore_lokalakTP = Ham_erroreakTP[2:end]-Ham_erroreakTP[1:end-1]
    
    pl3 = plot(tt, log10.(abs.(Ham_erroreakTP)),title = "Ham_erroreak TwicePrecision", legend=false, ylims=yrange)
    pl4 = plot(tt[2:end], log10.(abs.(Ham_errore_lokalakTP)), title ="Ham_errore lokalak TwicePrecision", legend=false, ylims=yrange)
    
    pl = plot(pl1,pl2, pl3,pl4, layout = (2,2), size = (1000, 600))
    display(pl)
end


function IrudikatuRet(tt, UU)
    ret = [norm(U[1:6]-UU[1][1:6]) for U in UU]
    pl = plot(tt, ret, title = "Uk-U0 Diferentziaren norma", legend=false)
    display(pl)
end


function IrudikatuIDeribatuak(tt, UU, uu; extra = false, n = 4)
    d = size(UU)[1]
    II = zeros(d)
    dI = zeros(d)
    d2I = zeros(d)
    d3I = zeros(d)
    if n==4
        d4I = zeros(d)
    end
    if extra
        extras = zeros(d)
    end
    for i in 1:d
        res = I_dnI_fcn(UU[i], uu[i], extra = extra, n = n)
        II[i] = res[1]
        dI[i] = res[2]
        d2I[i] = res[3]
        d3I[i] = res[4]
        if n==4
            d4I[i] = res[5]
            if extra
                extras[i] = res[6]
            end
        end
    end
    pl1 = plot(tt, II, title = "I", legend = false)
    pl2 = plot(tt, dI, title = "dI", legend = false)
    pl3 = plot(tt, d2I, title = "d2I", legend=false)
    pl4 = plot(tt, d3I, title = "d3I", legend=false)

    if n ==4
        pl5 = plot(tt, d4I, title = "d4I", legend=false)
        if extra
            pl6 = plot(tt, extras, title = "extra", legend=false)
            plot(pl1,pl2,pl3,pl4,pl5,pl6, layout = (2,3),size = (1000, 500))
        else
            plot(pl1,pl2,pl3,pl4,pl5, layout = (2,3),size = (1000, 500))
        end
    else
        plot(pl1,pl2,pl3,pl4, layout = (2,2),size = (1000, 500))
    end
end


function vector_str(v)
    return join(string.(v), ",")
end


function SoluzioPeriodikoakGorde(iterkop, filepath)
    z0_list = Any[]
    Zopt_list = Any[]
    n_list = Any[]
    sigma_list = Any[]
    
    z0 = zeros(4)
    counter = 0
    for i in 1:iterkop
        try
            z0 = pi*rand(4)
            u0 = HasierakoEgoera(z0)
            U0 = LCinvFcn(u0)
            h = ThreeBodyRelEnergy(u0)
            while (E_bitarra(U0, rel2abs(u0))>0 || TBP_Bitarra(U0, 1/256, 50., h))
                z0 = pi*rand(4)
                u0 = HasierakoEgoera(z0)
                U0 = LCinvFcn(u0)
                h = ThreeBodyRelEnergy(u0)
            end
            
            tauend = 20.
            dtau = 1/256
            
            #0 FASEA
            z = copy(z0)
            F0 = FLoss0(z,dtau,tauend,h)
            w0 = [z[3], z[4]]
            res0 = optimize(F0,w0)
            w1 = Normalizatu_z(res0.minimizer)


            #1 FASEA
            z1 = copy(F0.z)
            z1[3:4] .= w1
            F1 = FLoss1(dtau, tauend, h)
            res1 = optimize(F1, z1)
            z1_opt = Normalizatu_z(res1.minimizer)

            #2 FASEA
            T, sigma = F1(z1_opt, 2)            
            Z0 = copy(z1_opt)
            push!(Z0, T)
            F2 = FLoss2(sigma, Int64(ceil(T/ F1.dtau)), F1.h)
            res2 = optimize(F2, Z0, g_tol=1e-13, autodiff=:forward, method=BFGS())
            Zopt2 = Normalizatu_z(res2.minimizer)


            W0 = copy(Zopt2)
            push!(W0, h)
            F3 = FLoss3(sigma,Int64(ceil(W0[5]/F1.dtau)))
            res3 = optimize(F3, W0, g_tol=1e-13, autodiff=:forward, method=BFGS())
            Zopt3 = Normalizatu_z(res3.minimizer)
            
            
            push!(z0_list, z0)
            push!(Zopt_list, Zopt3)
            push!(n_list, F3.nsteps)
            push!(sigma_list, sigma)
            
        catch
            println(i, ". iterazioan akats bat egon da.")
            println("z0= ", z0)
            println("__________________________________")
            counter+=1
        end
    end
    
    # Goiburua zehaztu
    header = "z0_1, z0_2, z0_3, z0_4, Zopt_1, Zopt_2, Zopt_3, Zopt_4, Zopt_5, Zopt_6, n,sigma_1, sigma_2, sigma_3, t*"
    
    # Fitxategia existitzen den zehaztu
    file_exists = isfile(filepath)
    
    open(filepath, file_exists ? "a" : "w") do io
        # Fitxategia ez bada existitzen goiburua idatzi
        if !file_exists
            println(io, header)
        end
    
        # Idatzi lerro bakoitza
        for i in 1:length(z0_list)
            row = join([
                vector_str(z0_list[i]),
                vector_str(Zopt_list[i]),
                string(n_list[i]),
                vector_str(sigma_list[i]),
            ], ",")
            println(io, row)
        end
    end
    
    println("CSV fitxategian $(length(z0_list)) lerro berri gehitu dira.")
    println("Guztira $(counter) akats izan dira.")
end


function csvIrakurri(filepath, kodea = false)
    df = CSV.read(filepath, DataFrame)
    # Lortu lista bakoitzaren zutabe izenak
    list1_cols = filter(col -> startswith(string(col), "list1_"), names(df))
    list2_cols = filter(col -> startswith(string(col), "list2_"), names(df))
    list3_cols = filter(col -> startswith(string(col), "list3_"), names(df))
    list4_cols = filter(col -> startswith(string(col), "list4_"), names(df))
    
    
    # Listak berreraiko DataFrame lerroetatik
    list1 = [Vector{Float64}(df[i, list1_cols]) for i in 1:nrow(df)]
    list2 = [Vector{Float64}(df[i, list2_cols]) for i in 1:nrow(df)]
    list3 = [[eval(Meta.parse(df[i,list3_cols[1]])), Int64(df[i,list3_cols[2]]), df[i,list3_cols[3]]] for i in 1:nrow(df)];
    list4 = [Vector{Float64}(df[i, list4_cols]) for i in 1:nrow(df)]

    if !kodea
        return list1, list2, list3, list4
    else
        list5 = eval.(Meta.parse.(df[!, :list5]))
        list6 = eval.(Meta.parse.(df[!, :list6]))
        return list1, list2, list3, list4, list5, list6
    end
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



function filterSoluzioPeriodikoakGorde(sol_filepath, new_filepath)
    z0_list, Zopt_list, n_list, sigma_list = CSVIrakurri(sol_filepath);
    # Goiburua zehaztu
    header = "z0_1, z0_2, z0_3, z0_4, Zopt_1, Zopt_2, Zopt_3, Zopt_4, Zopt_5, Zopt_6, n,sigma_1, sigma_2, sigma_3"
        
    # Fitxategia existitzen den zehaztu
    file_exists = isfile(new_filepath)
    
    open(new_filepath, file_exists ? "a" : "w") do io
        # Fitxategia ez bada existitzen goiburua idatzi
        if !file_exists
            println(io, header)
        end
    
        # Idatzi lerro bakoitza
        for i in 1:length(z0_list)
            T = Zopt_list[i][5]
            n = n_list[i]
            F3 = FLoss3(sigma_list[i], n)
            helb_fun = F3(Zopt_list[i])
            h = Zopt_list[i][6]
            if T>0.1 && helb_fun<1e-5 && (-0.51<h<-0.49)
                row = join([
                    vector_str(z0_list[i]),
                    vector_str(Zopt_list[i]),
                    string(n_list[i]),
                    vector_str(sigma_list[i])
                ], ",")
                println(io, row)
            end
        end
    end
end


function filterEskuzSoluzioPeriodikoakGorde(sol_filepath, new_filepath, i_list)
    z0_list, Zopt_list, n_list, sigma_list = CSVIrakurri(sol_filepath);
    # Goiburua zehaztu
    header = "z0_1, z0_2, z0_3, z0_4, Zopt_1, Zopt_2, Zopt_3, Zopt_4, Zopt_5, Zopt_6, n,sigma_1, sigma_2, sigma_3"
        
    # Fitxategia existitzen den zehaztu
    file_exists = isfile(new_filepath)
    
    open(new_filepath, file_exists ? "a" : "w") do io
        # Fitxategia ez bada existitzen goiburua idatzi
        if !file_exists
            println(io, header)
        end
    
        # Idatzi lerro bakoitza
        for i in i_list
            row = join([
                    vector_str(z0_list[i]),
                    vector_str(Zopt_list[i]),
                    string(n_list[i]),
                    vector_str(sigma_list[i])
                ], ",")
                println(io, row)
        end
    end
end