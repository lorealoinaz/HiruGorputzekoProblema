include("TBP.jl");

struct FLoss0
    z::Vector{Float64}
    dtau::Float64
    taumax::Float64
    h::Float64
end

struct FLoss1
    dtau::Float64
    taumax::Float64
    h::Float64
end

struct FLoss2
    sigma::Tuple
    nsteps::Int64
    h::Float64
end


struct FLoss3
    sigma::Tuple
    nsteps::Int64
end


function Helb_fun_ald!(rws::Vector{Float64}, Uk::Vector{Float64}, uk::Vector{Float64}, σ::NTuple{3, Int}, vv::Vector{Float64}, dtdtaus::Float64)
    s12 = 2*(σ[1]-1)+1
    s23 = 2*(σ[2]-1)+1
    s31 = 2*(σ[3]-1)+1

    rs12 = Uk[s12]^2+Uk[s12+1]^2    
    rs23 = Uk[s23]^2+Uk[s23+1]^2    
    rs31 = Uk[s31]^2+Uk[s31+1]^2
    
    qqdots12 = uk[s12] * vv[s12] + uk[s12+1] * vv[s12+1]
    qqdots23 = uk[s23] * vv[s23] + uk[s23+1] * vv[s23+1]  
    qqdots31 = uk[s31] * vv[s31] + uk[s31+1] * vv[s31+1]
    
    
    ws12 = 2 * qqdots12 * dtdtaus
    ws23 = 2 * qqdots23 * dtdtaus
    ws31 = 2 * qqdots31 * dtdtaus
    
    rws[1] = rs12
    rws[2] = rs23
    rws[3] = rs31
    rws[4] = ws12
    rws[5] = ws23
    rws[6] = ws31
end

function (F0::FLoss0)(w::Vector{Float64}, out = 1, extrakop = 8)  
    F0.z[3] = w[1]
    F0.z[4] = w[2]
    z = F0.z
    dtau = F0.dtau
    taumax = F0.taumax
    h = F0.h
    par = h
    
    minn = Inf
    mink_ = 0.
    mink = Inf
    
    Tmin = -1.

    sigmas = collect(Tuple.(permutations(1:3)))
    
    sigmamin = sigmas[1]
    sigmamink = sigmas[1]

    
    gorde = false

    rw = zeros(6)
    rws = zeros(6)
    vv = zeros(6)

    u0 = HasierakoEgoera(z)
    U0 = LCinvFcn(u0)
    
    n = Int64(ceil(taumax/dtau))
    dUk = similar(U0)
    tk = Base.TwicePrecision(0.)
    uk = copy(u0)
    Uk_hi = copy(U0)
    Uk_lo = zero(U0)

    
    w = similar(U0)
    w_ = similar(U0)
    
    ∇ = [zero(U0) for i in 1:extrakop]    
    ∇_ = zero(U0)   
    ∇_new = zero(U0)
    
    ss = Float64[]
    
    #Abiadura bektore erlatiboak kalkulatu
    vv[1] = uk[9] - uk[7]
    vv[2] = uk[10] - uk[8]
    vv[3] = uk[11] - uk[9]
    vv[4] = uk[12] - uk[10]
    vv[5] = uk[7] - uk[11]
    vv[6] = uk[8] - uk[12]

    rs12 = Uk_hi[1]^2+Uk_hi[2]^2    
    rs23 = Uk_hi[3]^2+Uk_hi[4]^2    
    rs31 = Uk_hi[5]^2+Uk_hi[6]^2

    dtdtaus = sqrt(rs12+rs23+rs31)/(1/rs12 + 1/rs23 + 1/rs31)
    
    Helb_fun_ald!(rw, U0, u0, sigmas[1], vv, dtdtaus)
    r12, r23, r31, w12, w23, w31 = rw

    for j in eachindex(sigmas)
        Helb_fun_ald!(rws, Uk_hi, uk, sigmas[j], vv, dtdtaus)
        rs12, rs23, rs31, ws12, ws23, ws31 = rws 
        helb = (r12-rs12)^2+(r23-rs23)^2+(r31-rs31)^2+(w12-ws12)^2+(w23-ws23)^2+(w31-ws31)^2

        if helb<=mink
            mink = helb
            sigmamink = sigmas[j]
        end            
    end
    
    if out==3
        push!(ss, mink)
    end

    
    for k in 1:n
        dUk .= 0
        for i in 1:min(k-1, extrakop)
            @. dUk +=  ∇[i]
        end
        
        EMIStep!(dUk,Uk_hi,Uk_lo,dtau,par,w,w_)

        for j in eachindex(U0)
            Ukj = Base.TwicePrecision(Uk_hi[j],Uk_lo[j]) + dUk[j]
            Uk_hi[j] = Ukj.hi
            Uk_lo[j] = Ukj.lo
        end
        tk += dtau
        LCFcn!(uk, Uk_hi) 

            
        mink = Inf
        sigmamink = sigmas[1]
        
        #Abiadura bektore erlatiboak kalkulatu
        vv[1] = uk[9] - uk[7]
        vv[2] = uk[10] - uk[8]
        vv[3] = uk[11] - uk[9]
        vv[4] = uk[12] - uk[10]
        vv[5] = uk[7] - uk[11]
        vv[6] = uk[8] - uk[12]
        
        rs12 = Uk_hi[1]^2+Uk_hi[2]^2    
        rs23 = Uk_hi[3]^2+Uk_hi[4]^2    
        rs31 = Uk_hi[5]^2+Uk_hi[6]^2
    
        dtdtaus = sqrt(rs12+rs23+rs31)/(1/rs12 + 1/rs23 + 1/rs31)
        for j in eachindex(sigmas)
            Helb_fun_ald!(rws, Uk_hi, uk, sigmas[j], vv, dtdtaus)
            rs12, rs23, rs31, ws12, ws23, ws31 = rws 
            helb = (r12-rs12)^2+(r23-rs23)^2+(r31-rs31)^2+(w12-ws12)^2+(w23-ws23)^2+(w31-ws31)^2

            if helb<=mink
                mink = helb
                sigmamink = sigmas[j]
            end            
        end

        if out==3
            push!(ss, mink)
        end
        
        if gorde
            if mink<=minn
                minn = mink
                Tmin = tk.hi
                sigmamin = sigmamink
            end
        else
            if mink_> mink 
                gorde = true
            end
            mink_ = mink
        end 
    
        ∇_ .= ∇[1]
        ∇[1] .= dUk
        
        for i in 2:min(k, extrakop)
            for j in eachindex(∇_)
                ∇ji = ∇[i-1][j] - ∇_[j]
                ∇_[j] = ∇[i][j]
                ∇[i][j] = ∇ji
            end
        end
    end
    
    if out == 1
        return minn
    elseif out == 2
        return Tmin, sigmamin
    elseif out == 3
        return minn, ss
    end
end
    

function (F1::FLoss1)(z::Vector{Float64}, out = 1, extrakop = 8)
    dtau = F1.dtau
    taumax = F1.taumax
    h = F1.h
    par = h
    
    minn = Inf
    mink_ = 0.
    mink = Inf
    
    Tmin = -1.

    sigmas = collect(Tuple.(permutations(1:3)))
    
    sigmamin = sigmas[1]
    sigmamink = sigmas[1]

    
    gorde = false

    rw = zeros(6)
    rws = zeros(6)
    vv = zeros(6)

    u0 = HasierakoEgoera(z)
    U0 = LCinvFcn(u0)
    
    n = Int64(ceil(taumax/dtau))
    dUk = similar(U0)
    tk = Base.TwicePrecision(0.)
    uk = copy(u0)
    Uk_hi = copy(U0)
    Uk_lo = zero(U0)

    
    w = similar(U0)
    w_ = similar(U0)
    
    ∇ = [zero(U0) for i in 1:extrakop]    
    ∇_ = zero(U0)   
    ∇_new = zero(U0)
    
    ss = Float64[]
    
    #Abiadura bektore erlatiboak kalkulatu    
    vv[1] = uk[9] - uk[7]
    vv[2] = uk[10] - uk[8]
    vv[3] = uk[11] - uk[9]
    vv[4] = uk[12] - uk[10]
    vv[5] = uk[7] - uk[11]
    vv[6] = uk[8] - uk[12]

    rs12 = Uk_hi[1]^2+Uk_hi[2]^2    
    rs23 = Uk_hi[3]^2+Uk_hi[4]^2    
    rs31 = Uk_hi[5]^2+Uk_hi[6]^2

    dtdtaus = sqrt(rs12+rs23+rs31)/(1/rs12 + 1/rs23 + 1/rs31)
    
    Helb_fun_ald!(rw, U0, u0, sigmas[1], vv, dtdtaus)
    r12, r23, r31, w12, w23, w31 = rw

    for j in eachindex(sigmas)
        Helb_fun_ald!(rws, Uk_hi, uk, sigmas[j], vv, dtdtaus)
        rs12, rs23, rs31, ws12, ws23, ws31 = rws 
        helb = (r12-rs12)^2+(r23-rs23)^2+(r31-rs31)^2+(w12-ws12)^2+(w23-ws23)^2+(w31-ws31)^2

        if helb<=mink
            mink = helb
            sigmamink = sigmas[j]
        end            
    end
    
    if out==3
        push!(ss, mink)
    end

    
    for k in 1:n
        dUk .= 0
        for i in 1:min(k-1, extrakop)
            @. dUk +=  ∇[i]
        end
        
        EMIStep!(dUk,Uk_hi,Uk_lo,dtau,par,w,w_)

        for j in eachindex(U0)
            Ukj = Base.TwicePrecision(Uk_hi[j],Uk_lo[j]) + dUk[j]
            Uk_hi[j] = Ukj.hi
            Uk_lo[j] = Ukj.lo
        end
        tk += dtau
        LCFcn!(uk, Uk_hi) 

            
        mink = Inf
        sigmamink = sigmas[1]
        
        vv[1] = uk[9] - uk[7]
        vv[2] = uk[10] - uk[8]
        vv[3] = uk[11] - uk[9]
        vv[4] = uk[12] - uk[10]
        vv[5] = uk[7] - uk[11]
        vv[6] = uk[8] - uk[12]
        
        rs12 = Uk_hi[1]^2+Uk_hi[2]^2    
        rs23 = Uk_hi[3]^2+Uk_hi[4]^2    
        rs31 = Uk_hi[5]^2+Uk_hi[6]^2
    
        dtdtaus = sqrt(rs12+rs23+rs31)/(1/rs12 + 1/rs23 + 1/rs31)
        for j in eachindex(sigmas)
            Helb_fun_ald!(rws, Uk_hi, uk, sigmas[j], vv, dtdtaus)
            rs12, rs23, rs31, ws12, ws23, ws31 = rws 
            helb = (r12-rs12)^2+(r23-rs23)^2+(r31-rs31)^2+(w12-ws12)^2+(w23-ws23)^2+(w31-ws31)^2

            if helb<=mink
                mink = helb
                sigmamink = sigmas[j]
            end            
        end

        if out==3
            push!(ss, mink)
        end
        
        if gorde
            if mink<=minn
                minn = mink
                Tmin = tk.hi
                sigmamin = sigmamink                
            end
        else
            if mink_> mink 
                gorde = true
            end
            mink_ = mink
        end 
    
        ∇_ .= ∇[1]
        ∇[1] .= dUk
        
        for i in 2:min(k, extrakop)
            for j in eachindex(∇_)
                ∇ji = ∇[i-1][j] - ∇_[j]
                ∇_[j] = ∇[i][j]
                ∇[i][j] = ∇ji
            end
        end
    end
    
    if out == 1
        return minn
    elseif out == 2
        return Tmin, sigmamin
    elseif out == 3
        return minn, ss
    end
end


function (F2::FLoss2)(z)
    sigma = F2.sigma
    nsteps = F2.nsteps
    h = F2.h
    T = z[5]
    
    u0 = HasierakoEgoera(z)
    U0 = LCinvFcn(u0)

    r12 = U0[1]^2+U0[2]^2    
    r23 = U0[3]^2+U0[4]^2    
    r31 = U0[5]^2+U0[6]^2

    qqdot12 = u0[1] * (u0[9]-u0[7]) + u0[2] * (u0[10]-u0[8]) 
    qqdot23 = u0[3] * (u0[11]-u0[9]) + u0[4] * (u0[12]-u0[10])   
    qqdot31 = u0[5] * (u0[7]-u0[11]) + u0[6] * (u0[8]-u0[12]) 

    dtdtau = sqrt(r12+r23+r31)/(1/r12 + 1/r23 + 1/r31)
    
    w12 = 2 * qqdot12*dtdtau
    w23 = 2 * qqdot23*dtdtau
    w31 = 2 * qqdot31*dtdtau

    Utauend = EMI(U0, 0., T, nsteps, h)
    utauend = LCFcn(Utauend)
    
    s12 = 2*(sigma[1]-1)+1
    s23 = 2*(sigma[2]-1)+1
    s31 = 2*(sigma[3]-1)+1
    
    #Abiadura bektore erlatiboak kalkulatu
    vv= [utauend[9]  - utauend[7],
         utauend[10] - utauend[8],
         utauend[11] - utauend[9],
         utauend[12] - utauend[10],
         utauend[7]  - utauend[11],
         utauend[8]  - utauend[12]]

    rs12 = Utauend[s12]^2+Utauend[s12+1]^2    
    rs23 = Utauend[s23]^2+Utauend[s23+1]^2    
    rs31 = Utauend[s31]^2+Utauend[s31+1]^2

    qqdots12 = utauend[s12] * vv[s12] + utauend[s12+1] * vv[s12+1] 
    qqdots23 = utauend[s23] * vv[s23] + utauend[s23+1] * vv[s23+1]   
    qqdots31 = utauend[s31] * vv[s31] + utauend[s31+1] * vv[s31+1]
    
    dtdtaus = sqrt(rs12+rs23+rs31)/(1/rs12 + 1/rs23 + 1/rs31)
    
    ws12 = 2 * qqdots12 * dtdtaus
    ws23 = 2 * qqdots23 * dtdtaus
    ws31 = 2 * qqdots31 * dtdtaus

    
    helb = (r12-rs12)^2+(r23-rs23)^2+(r31-rs31)^2+(w12-ws12)^2+(w23-ws23)^2+(w31-ws31)^2

    return helb
end



function (F3::FLoss3)(z)
    sigma = F3.sigma
    nsteps = F3.nsteps
    T = z[5]
    h = z[6]
    
    u0 = HasierakoEgoera(z)
    U0 = LCinvFcn(u0)

    r12 = U0[1]^2+U0[2]^2    
    r23 = U0[3]^2+U0[4]^2    
    r31 = U0[5]^2+U0[6]^2

    qqdot12 = u0[1] * (u0[9]-u0[7]) + u0[2] * (u0[10]-u0[8]) 
    qqdot23 = u0[3] * (u0[11]-u0[9]) + u0[4] * (u0[12]-u0[10])   
    qqdot31 = u0[5] * (u0[7]-u0[11]) + u0[6] * (u0[8]-u0[12]) 

    dtdtau = sqrt(r12+r23+r31)/(1/r12 + 1/r23 + 1/r31)
    
    w12 = 2 * qqdot12*dtdtau
    w23 = 2 * qqdot23*dtdtau
    w31 = 2 * qqdot31*dtdtau

    Utauend = EMI(U0, 0., T, nsteps, h)
    utauend = LCFcn(Utauend)
    
    s12 = 2*(sigma[1]-1)+1
    s23 = 2*(sigma[2]-1)+1
    s31 = 2*(sigma[3]-1)+1
    
    #Abiadura bektore erlatiboak kalkulatu
    vv= [utauend[9]  - utauend[7],
         utauend[10] - utauend[8],
         utauend[11] - utauend[9],
         utauend[12] - utauend[10],
         utauend[7]  - utauend[11],
         utauend[8]  - utauend[12]]

    rs12 = Utauend[s12]^2+Utauend[s12+1]^2    
    rs23 = Utauend[s23]^2+Utauend[s23+1]^2    
    rs31 = Utauend[s31]^2+Utauend[s31+1]^2

    qqdots12 = utauend[s12] * vv[s12] + utauend[s12+1] * vv[s12+1] 
    qqdots23 = utauend[s23] * vv[s23] + utauend[s23+1] * vv[s23+1]   
    qqdots31 = utauend[s31] * vv[s31] + utauend[s31+1] * vv[s31+1]
    
    dtdtaus = sqrt(rs12+rs23+rs31)/(1/rs12 + 1/rs23 + 1/rs31)
    
    ws12 = 2 * qqdots12 * dtdtaus
    ws23 = 2 * qqdots23 * dtdtaus
    ws31 = 2 * qqdots31 * dtdtaus

    
    helb = (r12-rs12)^2+(r23-rs23)^2+(r31-rs31)^2+(w12-ws12)^2+(w23-ws23)^2+(w31-ws31)^2

    return helb
end
