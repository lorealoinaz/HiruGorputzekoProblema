include("TBP.jl");

function EMIStep!(dU0,U0,U0_lo,dt,par,w,w_) 
    jmax = 1000
    j = 0
    jarraitu = true
    dt_ = dt/2
    difer_ = Inf
    for j in eachindex(w)
        w[j] = U0[j] + (0.5*dU0[j] + U0_lo[j])
    end
    while jarraitu
        j+=1
        ThreeBodyODE_LC!(dU0, w, par)
        for j in eachindex(w)
            w_[j] = w[j]
            w[j] = U0[j] + (dt_*dU0[j] + U0_lo[j])
            w_[j] = w[j]-w_[j]
        end
        difer = norm(w_,Inf)
        jarraitu = (difer<difer_) && (difer!=0) && (j <= jmax)
        difer_ = difer
        #println("tauj=$(t0.hi), j=$j, difer=$difer")
    end   
    @. dU0 = dt * dU0
    return j
end


function EMI_Irudikatzeko(U0, t0, tauend, dt, par, extrakop=8) 
    n = Int64(ceil((tauend-t0)/dt))
    dUk = similar(U0)
    tk = Base.TwicePrecision(t0)
    Uk_hi = copy(U0)
    Uk_lo = zero(U0)
    w = similar(U0)
    w_ = similar(U0)
    
    iter = 0
    iter_lag = 0
    
    ∇ = [zero(U0) for i in 1:extrakop]    
    ∇_ = zero(U0)   

    UU = [U0]
    tt = [t0]   

    for k in 1:n
        dUk .= 0
        for i in 1:min(k-1, extrakop)
            @. dUk +=  ∇[i]
        end
        
        iter_lag = EMIStep!(dUk,Uk_hi,Uk_lo,dt,par,w,w_)
        iter+=iter_lag
        
        for j in eachindex(U0)
            Ukj = Base.TwicePrecision(Uk_hi[j],Uk_lo[j]) + dUk[j]
            Uk_hi[j] = Ukj.hi
            Uk_lo[j] = Ukj.lo
        end

        tk += dt

        push!(UU, copy(Uk_hi))
        push!(tt, tk)
        
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

    return tt, UU
end


function EMI(U0, t0, tauend, n, par, extrakop=8) 
    dt = (tauend-t0)/n
    dUk = similar(U0)
    tk = Base.TwicePrecision(t0)
    Uk_hi = copy(U0)
    Uk_lo = zero(U0)
    w = similar(U0)
    w_ = similar(U0)
        
    ∇ = [zero(U0) for i in 1:extrakop]    
    ∇_ = zero(U0)   
    for k in 1:n
        dUk .= 0
        for i in 1:min(k-1, extrakop)
            @. dUk +=  ∇[i]
        end
        
        EMIStep!(dUk,Uk_hi,Uk_lo,dt,par,w,w_)
        
        tk += dt

        for j in eachindex(U0)
            Ukj = Base.TwicePrecision(Uk_hi[j],Uk_lo[j]) + dUk[j] #berria, zaharra Uk_hi eta Uk_lo
            Uk_hi[j] = Ukj.hi
            Uk_lo[j] = Ukj.lo
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

    return Uk_hi
end



function EMIiter(U0, t0, tauend, dt, par, extrakop=8) 
    n = Int64(ceil((tauend-t0)/dt))
    dUk = similar(U0)
    tk = Base.TwicePrecision(t0)
    Uk_hi = copy(U0)
    Uk_lo = zero(U0)
    w = similar(U0)
    w_ = similar(U0)
    
    iter = 0
    iter_lag = 0
    
    ∇ = [zero(U0) for i in 1:extrakop]    
    ∇_ = zero(U0)   

    for k in 1:n
        dUk .= 0
        for i in 1:min(k-1, extrakop)
            @. dUk +=  ∇[i]
        end
        
        iter_lag = EMIStep!(dUk,Uk_hi,Uk_lo,dt,par,w,w_)
        iter+=iter_lag
        
        for j in eachindex(U0)
            Ukj = Base.TwicePrecision(Uk_hi[j],Uk_lo[j]) + dUk[j]
            Uk_hi[j] = Ukj.hi
            Uk_lo[j] = Ukj.lo
        end

        tk += dt
        
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

    
    return Uk_hi, iter/n
end


function EMIiter0(U0, t0, tauend, dt, par) 
    n = Int64(ceil((tauend-t0)/dt))
    dUk = similar(U0)
    tk = Base.TwicePrecision(t0)
    Uk_hi = copy(U0)
    Uk_lo = zero(U0)
    w = similar(U0)
    w_ = similar(U0)
    
    iter = 0
    iter_lag = 0
    
    for k in 1:n
        iter_lag = EMIStep!(dUk,Uk_hi,Uk_lo,dt,par,w,w_)
        iter+=iter_lag
        
        for j in eachindex(U0)
            Ukj = Base.TwicePrecision(Uk_hi[j],Uk_lo[j]) + dUk[j]
            Uk_hi[j] = Ukj.hi
            Uk_lo[j] = Ukj.lo
        end

        tk += dt
    end

    return Uk_hi, iter/n
end

