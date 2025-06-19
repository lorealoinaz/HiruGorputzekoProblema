function HasierakoEgoera(z)
    x = (1/4) *(1+cos(z[1]))
    
    y = sqrt(1-(x+0.5)^2) * (1+cos(z[2]))/2
    
    q1x_ = x
    q1y_ = y
    q2x_ = 0.5
    q2y_ = 0.
    q3x_ = -0.5
    q3y_ = 0.

    a_ = sqrt((q2x_ - q1x_)^2 + (q2y_ - q1y_)^2) #||q12||
    b_ = sqrt((q1x_ - q3x_)^2 + (q1y_ - q3y_)^2) #||q31||
    c_ = 1
    
    T = (1/4) *(1+cos(z[3]))
    
    λ = 2/(2*T+1)* (1/a_ + 1/b_ + 1/c_) #Eskalatzeko
    βx = -λ*(q1x_+q2x_+q3x_)/3 #Masa zentroa jatorrian kokatzeko
    βy = -λ*(q1y_+q2y_+q3y_)/3
    
    q1x = λ * q1x_ + βx
    q1y = λ * q1y_ + βy
    
    q2x = λ * q2x_ + βx
    q2y = λ * q2y_ + βy
    
    q3x = λ * q3x_ + βx
    q3y = λ * q3y_ + βy


    theta = 2*z[4]
    
    A = cos(theta)*(q1y-q2y) + sin(theta)*(q2x-q1x)
    B = 1/2*((q2x-q1x)*cos(theta)+(q2y-q1y)*sin(theta))
    C = -(2*B*q3x-A*q3y)/(3*(q3x^2+q3y^2))
    D = -(A*q3x + 2*B*q3y)/(3*(q3x^2+q3y^2))

    rho12 = sqrt(4T/(3*(C^2+D^2)+1))
    
    qdot3x = C*rho12
    qdot3y = D*rho12
    qdot12x = rho12*cos(theta)
    qdot12y = rho12*sin(theta)

    qdot1x = -0.5 * (qdot3x + qdot12x)
    qdot1y = -0.5 * (qdot3y + qdot12y)
    qdot2x =  0.5 * (qdot12x - qdot3x)
    qdot2y =  0.5 * (qdot12y - qdot3y)

    return [q2x-q1x, q2y-q1y, q3x-q2x, q3y-q2y, q1x-q3x, q1y-q3y,
            qdot1x, qdot1y, qdot2x, qdot2y, qdot3x, qdot3y]
    
end


function Normalizatu_z(z)
    norm_z = Float64[]
    for i in 1:min(length(z), 4)
        z_mod = mod(z[i], 2 * pi)
        if z_mod <= pi
            push!(norm_z, z_mod)
        else
            push!(norm_z, 2 * pi - z_mod)
        end
    end
    for i in 5:length(z)
        push!(norm_z, z[i])
    end
    return norm_z
end


"""
    AngularMomentum(u)

Hiru gorputzen problemaren 'u' egoerari dagokion momentu angeluarra itzultzen du
"""
function AngularMomentum(u)
    q = [u[1:2], u[3:4], u[5:6]]
    p = [u[7:8], u[9:10], u[11:12]]
    L = 0.
    for i in 1:3
        L += q[i][1]*p[i][2] - q[i][2] * p[i][1]  #Lz =x py−y px
    end
    return L
end



"""
    ThreeBodyRelPotential(u)

Hiru gorputzen problemaren 'u' egoerari dagokion energia potentziala itzultzen du.
#Arguments
    - u::Vector{Float64}: Hiru gorputzeko problemari dagokion egoera, posizioak erlatiboak eta abiadurak absolutuak izanik
"""
function ThreeBodyRelPotential(u)
    q1 = u[1:2]
    q2 = u[3:4]
    q3 = u[5:6]
    U = -1/norm(q1)-1/norm(q2)-1/norm(q3)
    return U
end


"""
    ThreeBodyRelEnergy(u)

Hiru gorputzen problemaren 'u' egoerari dagokion energia totala itzultzen du.
#Arguments
    - u::Vector{Float64}: Hiru gorputzeko problemari dagokion egoera, posizioak erlatiboak eta abiadurak absolutuak izanik
"""
function ThreeBodyRelEnergy(u)
    q1 = u[1:2]
    q2 = u[3:4]
    q3 = u[5:6]
    v = u[7:12]
    U = -1/norm(q1)-1/norm(q2)-1/norm(q3)
    return dot(v,v)/2 + U
end


"""
    ThreeBodyAbsPotential(u)

Hiru gorputzen problemaren 'u' egoerari dagokion energia potentziala itzultzen du.
#Arguments
    - u::Vector{Float64}: Hiru gorputzeko problemari dagokion egoera, posizioak eta abiadurak absolutuak izanik
"""
function ThreeBodyAbsPotential(u)
    q1 = u[1:2]
    q2 = u[3:4]
    q3 = u[5:6]
    U = -1/norm(q2-q1)-1/norm(q3-q2)-1/norm(q1-q3)
    return U
end


"""
    ThreeBodyAbsEnergy(u)

Hiru gorputzen problemaren 'u' egoerari dagokion energia totala itzultzen du.
#Arguments
    - u::Vector{Float64}: Hiru gorputzeko problemari dagokion egoera, posizioak eta abiadurak absolutuak izanik
"""
function ThreeBodyAbsEnergy(u)
    q1 = u[1:2]
    q2 = u[3:4]
    q3 = u[5:6]
    v = u[7:12]
    U = -1/norm(q2-q1)-1/norm(q3-q2)-1/norm(q1-q3)
    return dot(v,v)/2 + U
end



"""
    g(U)

Hiru gorputzen problemaren 'U' egoera erregularizatuari dagokion inertzia momentua itzultzen du.

#Arguments
    - U::Vector{Float64}: Hiru gorputzeko problemari dagokion egoera erregularizatua.
"""
function g(U)
    #LCFcn
    X12 = U[1]
    Y12 = U[2]  
    px12 = U[7]
    py12 = U[8]
    x12, y12, px12, py12 = LC(X12,Y12,px12,py12)
    X23 = U[3]
    Y23 = U[4]  
    px23 = U[9]
    py23 = U[10]
    x23, y23, px23, py23 = LC(X23,Y23,px23,py23)
    X31 = U[5]
    Y31 = U[6]  
    px31 = U[11]
    py31 = U[12]
    x31, y31, px31, py31 = LC(X31,Y31,px31,py31)
    vx1 = px31-px12
    vy1 = py31-py12
    vx2 = px12-px23
    vy2 = py12-py23
    vx3 = px23-px31
    vy3 = py23-py31
    
    
    #rel2abs
    x1 = (x31-x12)/3.
    y1 = (y31-y12)/3.
    
    x2 = (x12-x23)/3.
    y2 = (y12-y23)/3.
    
    x3 = (x23-x31)/3.    
    y3 = (y23-y31)/3.

    #Inertzia momentuaren deribatua
    return 2*(x1*vx1+y1*vy1 + x2*vx2+y2*vy2 + x3*vx3+y3*vy3)
end

"""
    abs2rel(u)

Hiru gorputzen problemaren balio absolutuak dituen 'u' egoerako posizioak balio erlatiboetara pasa eta itzultzen du. 
"""
abs2rel(u) = [u[3]-u[1], u[4]-u[2], u[5]-u[3], u[6]-u[4], u[1]-u[5], u[2]-u[6],
                u[7], u[8], u[9], u[10], u[11], u[12]#=, u[13]=#]
                 #u[7]-u[9], u[8]-u[10], u[9]-u[11], u[10]-u[12], u[11]-u[7], u[12]-u[8], u[13]]


"""
    rel2abs(u)

Hiru gorputzen problemaren posizio erlatiboak dituen 'u' egoera balio absolutuetara pasa eta itzultzen du. 
"""
function rel2abs(u)
    x1 = (u[5]-u[1])/3.
    x2 = (u[1]-u[3])/3.
    x3 = (u[3]-u[5])/3.
    y1 = (u[6]-u[2])/3.
    y2 = (u[2]-u[4])/3.
    y3 = (u[4]-u[6])/3.
    
    vx1 = (u[11]-u[7])/3.
    vx2 = (u[7]-u[9])/3.
    vx3 = (u[9]-u[11])/3.
    vy1 = (u[12]-u[8])/3.
    vy2 = (u[8]-u[10])/3.
    vy3 = (u[10]-u[12])/3.
    
    return [x1, y1, x2, y2, x3, y3, u[7], u[8], u[9], u[10], u[11], u[12]#=, u[13]=#]
end


"""
    LC(X,Y,Px,Py)

Levi Civita-ren aldagai aldaketa jasan duten 'X', 'Y', 'Px' eta 'Py' aldagaiak jatorrizko eran itzultzen ditu.
"""
function LC(X,Y,Px,Py)
    x = X*X-Y*Y
    y = 2*X*Y 
    r = X*X+Y*Y
    px = (X*Px-Y*Py)/2r
    py = (Y*Px+X*Py)/2r
    return x,y,px,py
end    


"""
    LC(X,Y,Px,Py)

Jatorrizko aldagaiak, 'x', 'y', 'px' eta 'py', Levi Civita-ren aldagai aldaketa jasan ostean, erregularizatuta, itzultzen ditu.
"""
function LCinv(x,y,px,py)
    r = sqrt(x^2+y^2)
    if x>0
        fact = 1/sqrt(2*(r+x))
        X = fact*(x+r)
        Y = fact*y
    else
        fact = 1/sqrt(2*(r-x))
        X = fact*y
        Y = fact*(r-x)
    end
    Px = 2*(X*px+Y*py)
    Py = 2*(-Y*px+X*py)
    return X,Y,Px,Py
end


"""
    LCFcn(U) 

Hiru gorputzen problemako egoera erregularizatua, Levi Civita-ren aldagai aldaketa izan duena, erregularizatu gabe itzultzen du.
"""
function LCFcn(U) 
    X12 = U[1]
    Y12 = U[2]  
    Px12 = U[7]
    Py12 = U[8]
    x12, y12, px12, py12 = LC(X12,Y12,Px12,Py12)
    X23 = U[3]
    Y23 = U[4]  
    Px23 = U[9]
    Py23 = U[10]
    x23, y23, px23, py23 = LC(X23,Y23,Px23,Py23)
    X31 = U[5]
    Y31 = U[6]  
    Px31 = U[11]
    Py31 = U[12]
    x31, y31, px31, py31 = LC(X31,Y31,Px31,Py31)
    vx1 = (px31-px12)
    vy1 = (py31-py12)
    vx2 = (px12-px23)
    vy2 = (py12-py23)
    vx3 = (px23-px31)
    vy3 = (py23-py31)

    #return [x12,y12,x23,y23,x31,y31,vx1,vy1,vx2,vy2,vx3,vy3,U[13]]
    return [x12,y12,x23,y23,x31,y31,vx1,vy1,vx2,vy2,vx3,vy3]        
end


"""
    LCFcn(U) 

Hiru gorputzen problemako egoera erregularizatua, Levi Civita-ren aldagai aldaketa izan duena, erregularizatu gabe itzultzen du.
"""
function LCFcn!(u, U) 
    X12 = U[1]
    Y12 = U[2]  
    Px12 = U[7]
    Py12 = U[8]
    x12, y12, px12, py12 = LC(X12,Y12,Px12,Py12)
    X23 = U[3]
    Y23 = U[4]  
    Px23 = U[9]
    Py23 = U[10]
    x23, y23, px23, py23 = LC(X23,Y23,Px23,Py23)
    X31 = U[5]
    Y31 = U[6]  
    Px31 = U[11]
    Py31 = U[12]
    x31, y31, px31, py31 = LC(X31,Y31,Px31,Py31)
    vx1 = (px31-px12)
    vy1 = (py31-py12)
    vx2 = (px12-px23)
    vy2 = (py12-py23)
    vx3 = (px23-px31)
    vy3 = (py23-py31)
    u[1] = x12
    u[2] = y12
    u[3] = x23
    u[4] = y23
    u[5] = x31
    u[6] = y31
    u[7] = vx1
    u[8] = vy1
    u[9] = vx2
    u[10] = vy2
    u[11] = vx3
    u[12] = vy3          
end



"""
    LCinvFcn(U) 

Hiru gorputzen problemako egoera erregularizatu gabea, Levi Civita-ren aldagai aldaketa aplikatu ostean itzultzen du, erregularizatua.
"""
function LCinvFcn(u) 
    x12 = u[1]
    y12 = u[2]  
    vx1 = u[7]
    vy1 = u[8]
    x23 = u[3]
    y23 = u[4]  
    vx2 = u[9]
    vy2 = u[10]
    x31 = u[5]
    y31 = u[6]  
    vx3 = u[11]
    vy3 = u[12]
    
    norm12 = norm([x12,y12])
    norm23 = norm([x23,y23])
    norm31 = norm([x31,y31])

    if(norm12 >= norm23 && norm31 >= norm23) #Masa urrunena m1 denean
        px12 = -vx1/2
        py12 = -vy1/2
        px23 = (vx3-vx2)/2
        py23 = (vy3-vy2)/2
        px31 = vx1/2
        py31 = vy1/2
        
    elseif(norm12 >= norm31 && norm23 >= norm31) #Masa urrunena m2 denean
        px12 = vx2/2
        py12 = vy2/2
        px23 = -vx2/2
        py23 = -vy2/2
        px31 = (vx1-vx3)/2
        py31 = (vy1-vy3)/2
        
    elseif(norm31 >= norm12 && norm23 >= norm12) #Masa urrunena m2 denean        
        px12 = (vx2-vx1)/2
        py12 = (vy2-vy1)/2
        px23 = vx3/2
        py23 = vy3/2
        px31 = -vx3/2
        py31 = -vy3/2 
    else
        println(u)
    end
    
    X12, Y12, Px12, Py12 = LCinv(x12,y12,px12,py12)
    X23, Y23, Px23, Py23 = LCinv(x23,y23,px23,py23)
    X31, Y31, Px31, Py31 = LCinv(x31,y31,px31,py31)
    #return [X12,Y12,X23,Y23,X31,Y31,Px12,Py12,Px23,Py23,Px31,Py31,u[13]]
    return [X12,Y12,X23,Y23,X31,Y31,Px12,Py12,Px23,Py23,Px31,Py31]
end


"""
    ThreeBodyHamiltonianLC(U,h=-0.5) 

Hiru gorputzen problemaren egoera erregularizatuaren funtzio hamiltondarraren emaitza itzultzen du.
"""
function ThreeBodyHamiltonianLC(U,h=-0.5)
    Qx1 = U[1]
    Qy1 = U[2]
    Qx2 = U[3]
    Qy2 = U[4]
    Qx3 = U[5]
    Qy3 = U[6]
    Px1 = U[7]
    Py1 = U[8]
    Px2 = U[9]
    Py2 = U[10]
    Px3 = U[11]
    Py3 = U[12]
    r12 = Qx1^2+Qy1^2
    r23 = Qx2^2+Qy2^2
    r31 = Qx3^2+Qy3^2
    P12_2 = Px1^2+Py1^2
    P23_2 = Px2^2+Py2^2
    P31_2 = Px3^2+Py3^2
    a12 = r23*r31
    a23 = r31*r12
    a31 = r12*r23
    b12 = P12_2*a12
    b23 = P23_2*a23
    b31 = P31_2*a31
    ux12 = Qx1*Px1-Qy1*Py1
    uy12 = Px1*Qy1+Py1*Qx1
    ux23 = Qx2*Px2-Qy2*Py2
    uy23 = Px2*Qy2+Py2*Qx2
    ux31 = Qx3*Px3-Qy3*Py3
    uy31 = Px3*Qy3+Py3*Qx3
    w31 = ux12*ux23+uy12*uy23
    w12 = ux23*ux31+uy23*uy31
    w23 = ux31*ux12+uy31*uy12
    c12 = w12*r12
    c23 = w23*r23
    c31 = w31*r31
    R = r12*r23*r31
    A = (b12 + b23 + b31) - 
        (c12+c23+c31) -
        4*h*R
    Binv = 4*(a12 + a23 + a31)
    B = 1/Binv
    C = A*B - 1
    L = r12 + r23 + r31
    Lsqrt = sqrt(L)
    Ham = Lsqrt*C
    return Ham
end



"""
    ThreeBodyODE_LC!(du,u,par,t) 

Hiru gorputzen problemaren egoera erregularizatuaren ODE kalkulatzen du, denboraren birparametrizazioa aplikatuz. Denbora fiktizioa itzultzen du?
"""
function ThreeBodyODE_LC!(du,u,par)
    h = par[1]
    Qx1 = u[1]
    Qy1 = u[2]
    Qx2 = u[3]
    Qy2 = u[4]
    Qx3 = u[5]
    Qy3 = u[6]
    Px1 = u[7]
    Py1 = u[8]
    Px2 = u[9]
    Py2 = u[10]
    Px3 = u[11]
    Py3 = u[12]
    r12 = Qx1^2+Qy1^2
    r23 = Qx2^2+Qy2^2
    r31 = Qx3^2+Qy3^2
    P12_2 = Px1^2+Py1^2
    P23_2 = Px2^2+Py2^2
    P31_2 = Px3^2+Py3^2
    a12 = r23*r31
    a23 = r31*r12
    a31 = r12*r23
    b12 = P12_2*a12
    b23 = P23_2*a23
    b31 = P31_2*a31
    ux12 = Qx1*Px1-Qy1*Py1
    uy12 = Px1*Qy1+Py1*Qx1
    ux23 = Qx2*Px2-Qy2*Py2
    uy23 = Px2*Qy2+Py2*Qx2
    ux31 = Qx3*Px3-Qy3*Py3
    uy31 = Px3*Qy3+Py3*Qx3
    w31 = ux12*ux23+uy12*uy23
    w12 = ux23*ux31+uy23*uy31
    w23 = ux31*ux12+uy31*uy12
    c12 = w12*r12
    c23 = w23*r23
    c31 = w31*r31
    R = r12*r23*r31
    A = (b12 + b23 + b31) - 
        (c12+c23+c31) -
        4*h*R
    Binv = 4*(a12 + a23 + a31)
    B = 1/Binv
    C = A*B - 1
    L = r12 + r23 + r31
    Lsqrt = sqrt(L)
    Ham = Lsqrt*C
    # Gradientea kalkulatzeko hasieraketak
    hQx1 = 0.
    hQy1 = 0.
    hQx2 = 0.
    hQy2 = 0.
    hQx3 = 0.
    hQy3 = 0.
    hPx1 = 0.
    hPy1 = 0.
    hPx2 = 0.
    hPy2 = 0.
    hPx3 = 0.
    hPy3 = 0.
    hr12 = 0.
    hr23 = 0.
    hr31 = 0.
    hr12inv = 0.
    hr23inv = 0.
    hr31inv = 0.
    hP12_2 = 0.
    hP23_2 = 0.
    hP31_2 = 0.
    ha12 = 0.
    ha23 = 0.
    ha31 = 0.
    hb12 = 0.
    hb23 = 0.
    hb31 = 0.
    hux12 = 0.
    huy12 = 0.
    hux23 = 0.
    huy23 = 0.
    hux31 = 0.
    huy31 = 0.
    hw31 = 0.
    hw12 = 0.
    hw23 = 0.
    hc12 = 0.
    hc23 = 0.
    hc31 = 0.
    hR = 0.
    hA = 0.
    hBinv = 0.
    hB = 0.
    hC = 0.
    hL = 0.
    hLsqrt = 0.
    hHam = 1.
# Gradientearen kalkulua
    #Ham = Lsqrt*C
    hLsqrt = hHam*C
    hC = hHam*Lsqrt
    #Lsqrt = sqrt(L)
    hL = 0.5*hLsqrt/Lsqrt
    #L = r12 + r23 + r31
    hr12 = hL
    hr23 = hL
    hr31 = hL
    #C = A*B - 1
    hA += hC*B
    hB += hC*A
    #B = 1/Binv
    hBinv += -B^2*hB
    #Binv = 4*(a12 + a23 + a31)
    ha12 += 4*hBinv
    ha23 += 4*hBinv
    ha31 += 4*hBinv
    #A = (b12 + b23 + b31) - 
    #    (c12+c23+c31) -
    #    4*h*R
    hb12 += hA
    hb23 += hA
    hb31 += hA
    hc12 += -hA
    hc23 += -hA
    hc31 += -hA
    hR += -4*hA*h
    #R = r12*r23*r31
    hr12 += r23*r31*hR
    hr23 += r12*r31*hR
    hr31 += r12*r23*hR
    #c12 = w12*r12
    #c23 = w23*r23
    #c31 = w31*r31
    hw12 += hc12*r12
    hr12 += hc12*w12
    hw23 += hc23*r23
    hr23 += hc23*w23
    hw31 += hc31*r31
    hr31 += hc31*w31
    #w31 = ux12*ux23+uy12*uy23
    #w12 = ux23*ux31+uy23*uy31
    #w23 = ux31*ux12+uy31*uy12
    hux12 += hw31*ux23
    hux23 += hw31*ux12
    huy12 += hw31*uy23
    huy23 += hw31*uy12
    
    hux23 += hw12*ux31
    hux31 += hw12*ux23
    huy23 += hw12*uy31
    huy31 += hw12*uy23
    
    hux31 += hw23*ux12
    hux12 += hw23*ux31
    huy31 += hw23*uy12
    huy12 += hw23*uy31
    #ux12 = Qx1*Px1-Qy1*Py1
    #uy12 = Px1*Qy1+Py1*Qx1
    #ux23 = Qx2*Px2-Qy2*Py2
    #uy23 = Px2*Qy2+Py2*Qx2
    #ux31 = Qx3*Px3-Qy3*Py3
    #uy31 = Px3*Qy3+Py3*Qx3
    hQx1 += hux12*Px1
    hPx1 += hux12*Qx1
    hQy1 += -hux12*Py1
    hPy1 += -hux12*Qy1
    hQx1 += huy12*Py1
    hPx1 += huy12*Qy1
    hQy1 += huy12*Px1
    hPy1 += huy12*Qx1
    
    hQx2 += hux23*Px2
    hPx2 += hux23*Qx2
    hQy2 += -hux23*Py2
    hPy2 += -hux23*Qy2
    hQx2 += huy23*Py2
    hPx2 += huy23*Qy2
    hQy2 += huy23*Px2
    hPy2 += huy23*Qx2
    
    hQx3 += hux31*Px3
    hPx3 += hux31*Qx3
    hQy3 += -hux31*Py3
    hPy3 += -hux31*Qy3
    hQx3 += huy31*Py3
    hPx3 += huy31*Qy3
    hQy3 += huy31*Px3
    hPy3 += huy31*Qx3
    
    #b12 = P12_2*a12
    #b23 = P23_2*a23
    #b31 = P31_2*a31
    hP12_2 += hb12*a12
    ha12 += hb12*P12_2
    
    hP23_2 += hb23*a23
    ha23 += hb23*P23_2
    
    hP31_2 += hb31*a31
    ha31 += hb31*P31_2
    
    #a12 = r23*r31
    #a23 = r31*r12
    #a31 = r12*r23
    
    hr23 += ha12*r31
    hr31 += ha12*r23
    
    hr31 += ha23*r12
    hr12 += ha23*r31
    
    hr12 += ha31*r23
    hr23 += ha31*r12
    
    #P12_2 = Px1^2+Py1^2
    #P23_2 = Px2^2+Py2^2
    #P31_2 = Px3^2+Py3^2
    hPx1 += 2*Px1*hP12_2
    hPy1 += 2*Py1*hP12_2
    
    hPx2 += 2*Px2*hP23_2
    hPy2 += 2*Py2*hP23_2
    
    hPx3 += 2*Px3*hP31_2
    hPy3 += 2*Py3*hP31_2
    
    #r12 = Qx1^2+Qy1^2
    #r23 = Qx2^2+Qy2^2
    #r31 = Qx3^2+Qy3^2
    hQx1 += 2*Qx1*hr12
    hQy1 += 2*Qy1*hr12
    
    hQx2 += 2*Qx2*hr23
    hQy2 += 2*Qy2*hr23
    
    hQx3 += 2*Qx3*hr31
    hQy3 += 2*Qy3*hr31
    
    du[1] = hPx1
    du[2] = hPy1
    du[3] = hPx2
    du[4] = hPy2
    du[5] = hPx3
    du[6] = hPy3
    du[7] = -hQx1
    du[8] = -hQy1
    du[9] = -hQx2
    du[10] = -hQy2
    du[11] = -hQx3
    du[12] = -hQy3
    #s = 4*R*B*Lsqrt
    #du[13] = s
    #return s
    return nothing
end


function I_dnI_fcn(U, u; extra = false, n = 4)
    x12 = u[1]
    y12 = u[2]
    x23 = u[3]
    y23 = u[4]
    x31 = u[5]
    y31 = u[6]
    vx1 = u[7]
    vy1 = u[8]
    vx2 = u[9]
    vy2 = u[10]
    vx3 = u[11]
    vy3 = u[12]
    vx12 = vx2-vx1
    vy12 = vy2-vy1
    vx23 = vx3-vx2
    vy23 = vy3-vy2
    vx31 = vx1-vx3
    vy31 = vy1-vy3

    X12 = U[1]
    Y12 = U[2]
    X23 = U[3]
    Y23 = U[4]
    X31 = U[5]
    Y31 = U[6]

    r12 = X12^2+Y12^2    
    r23 = X23^2+Y23^2    
    r31 = X31^2+Y31^2

    dtdtau = sqrt(r12+r23+r31)/(1/r12 + 1/r23 + 1/r31)

    qqdot12 = x12 * vx12 + y12 * vy12    
    qqdot23 = x23 * vx23 + y23 * vy23   
    qqdot31 = x31 * vx31 + y31 * vy31 
    
    r3_12 = r12^3
    r3_23 = r23^3
    r3_31 = r31^3

    if extra
        q12qdot3 = x12 * vx3 + y12 *vy3
        q23qdot1 = x23 * vx1 + y23 *vy1
        q31qdot2 = x31 * vx2 + y31 *vy2
        
        valextra = q12qdot3 + q23qdot1 + q31qdot2
        valextra*=dtdtau
    end

    if n==4
        r5_12 = r12^5
        r5_23 = r23^5
        r5_31 = r31^5

        vnorm12 = vx12^2+vy12^2  
        vnorm23 = vx23^2+vy23^2
        vnorm31 = vx31^2+vy31^2
        
        ax12 = -2*x12/r3_12 + x23/r3_23 + x31/r3_31
        ay12 = -2*y12/r3_12 + y23/r3_23 + y31/r3_31
        ax23 = -2*x23/r3_23 + x31/r3_31 + x12/r3_12
        ay23 = -2*y23/r3_23 + y31/r3_31 + y12/r3_12
        ax31 = -2*x31/r3_31 + x12/r3_12 + x23/r3_23
        ay31 = -2*y31/r3_31 + y12/r3_12 + y23/r3_23
    
        qa12 = x12*ax12+y12*ay12   
        qa23 = x23*ax23+y23*ay23
        qa31 = x31*ax31+y31*ay31

        aux1 = 6 * (1/r5_12 * qqdot12^2 + 1/r5_23 * qqdot23^2 + 1/r5_31 * qqdot31^2) 
        aux2 = -2 * (1/r3_12 * vnorm12 + 1/r3_23 * vnorm23 + 1/r3_31 * vnorm31) 
        aux3 = -2 * (1/r3_12 * qa12 + 1/r3_23 * qa23 + 1/r3_31 * qa31)
        
        d4I = aux1+aux2+aux3
        d4I*=dtdtau^4
    end
    
    I = 1/3 * (r12^2+r23^2+r31^2)
    
    dI = 2/3 * (qqdot12+qqdot23+qqdot31)
    dI*=dtdtau

    d2I = 2*(1/r12 + 1/r23 + 1/r31 - 1)
    d2I*=dtdtau^2
    
    d3I = -2 * (1/r3_12 * qqdot12 + 1/r3_23 * qqdot23 + 1/r3_31 * qqdot31)
    d3I*=dtdtau^3
    
    if extra && n==4
        return (I, dI, d2I, d3I, d4I, valextra)
    elseif n==4
        return (I, dI, d2I, d3I, d4I)
    else
        return (I, dI, d2I, d3I)
    end
end

function IELog(U, u)
    r12 = U[1]^2+U[2]^2    
    r23 = U[3]^2+U[4]^2    
    r31 = U[5]^2+U[6]^2

    I = 1/3 * (r12^2+r23^2+r31^2)

    aux1 = min(r23,r31)
    if r12 <= aux1 #kanpokoa 3
        rmin = r12  
        k = 0
    elseif r23<=r31 #kanpokoa 1
        rmin = r23
        k = 1
    else #kanpokoa 2
        rmin = r31
        k = 2
    end
    #rmin
    r_ratio = minimum([r12,r23,r31])/maximum([r12,r23,r31])
    
    i1 = 2*k+1
    i2 = i1 + 1
    r = U[i1]^2 + U[i2]^2
    k_ = mod(k+1,3)
    j1 = 2*k_+1
    j2 = j1 + 1
    #u = LCFcn(U)
    
    E_barne = (0.25*((u[i1+6]-u[j1+6])^2+(u[i2+6]-u[j2+6])^2) - 1/r + 0.5)

    i1 = 2*mod(k+1,3)+1
    i2 = i1 + 1
    b = U[i1]^2 + U[i2]^2
    k_ = mod(k+2,3)
    j1 = 2*k_+1
    j2 = j1 + 1
    c = U[j1]^2 + U[j2]^2
    #u = LCFcn(U)
    k__ = mod(k-1,3)
    l1 = 2*k__+7
    l2 = l1 + 1
    vx31 = u[l1]
    vy31 = u[l2]
    E_kanpo = 3/4 * (vx31^2 + vy31^2) - (1/b + 1/c) #4/(2*b+ a)

    return (I, rmin, r_ratio, E_barne, E_kanpo)
end

function ind_kanpokoa(U)
    r12 = U[1]^2+U[2]^2    
    r23 = U[3]^2+U[4]^2    
    r31 = U[5]^2+U[6]^2
    aux1 = min(r23,r31)
    if r12 <= aux1
        ind = 2   # Kanpokoa 3
    elseif r23<=r31
        ind = 0 # Kanpokoa 1
    else
        ind = 1 # Kanpokoa 2
    end
    return ind
end


"""
    getFarMiddle(U)

Hiru gorputzen problemaren 'U' egoera erregularizatuan urrunen dagoen eta erdian dagoen gorputzen zenbakia itzultzen du.

#Arguments
    - U::Vector{Float64}: Hiru gorputzeko problemari dagokion egoera erregularizatua.
"""
function getFarMiddle(U)
    r12 = U[1]^2+U[2]^2    
    r23 = U[3]^2+U[4]^2    
    r31 = U[5]^2+U[6]^2
    far = 0
    if r12==r23
        middle = 2
    elseif r31==r12
        middle = 1
    elseif r23==r31
        middle = 3
    else
        aux1 = min(r23,r31)
        if r12 <= aux1
            far = 3   
            if r23<r31
                middle = 2
            else 
                middle = 1
            end
        elseif r23<=r31
            far = 1 
                if r12 < r31
                    middle = 2
                else
                    middle = 3
                end
        else
            far = 2 
            if r12<r23
                middle = 1
            else
                middle = 3
            end
        end
    end
    return far, middle
end



function ind_rmin(U)
    r12 = U[1]^2+U[2]^2    
    r23 = U[3]^2+U[4]^2    
    r31 = U[5]^2+U[6]^2
    aux1 = min(r23,r31)
    if r12 <= aux1
        ind = 0   # Kanpokoa 3
    elseif r23<=r31
        ind = 1 # Kanpokoa 1
    else
        ind = 2 # Kanpokoa 2
    end
    return ind
end

function E_bitarra(U, uabs)
    k = ind_kanpokoa(U)
    i1 = 2*k+1
    i2 = i1 + 1
    normq3 = sqrt(uabs[i1]^2 + uabs[i2]^2)
    norm2_qdot3 = uabs[i1+6]^2 + uabs[i2+6]^2
    return 0.75*norm2_qdot3 - 4/3 * normq3
end


function TBP_Bitarra(U0, dtau, taumax, h)
    par = [h]
 
    n = Int64(ceil(taumax/dtau))
    dUk = similar(U0)
    tk = Base.TwicePrecision(0.)
    Uk_hi = copy(U0)
    Uk_lo = zero(U0)

    w = similar(U0)
    w_ = similar(U0)
   
    far_ = ind_kanpokoa(U0)
    far = 0
    
    for k in 1:n
        
        EMIStep!(dUk,Uk_hi,Uk_lo,dtau,par,w,w_)

        for j in eachindex(U0)
            Ukj = Base.TwicePrecision(Uk_hi[j],Uk_lo[j]) + dUk[j]
            Uk_hi[j] = Ukj.hi
            Uk_lo[j] = Ukj.lo
        end
        tk += dtau
        far = ind_kanpokoa(Uk_hi)
        if far_!=far
            return false
        end
        far_ = far
    end
    return true
end


function I_fcn(U)
    X12 = U[1]
    Y12 = U[2]
    X23 = U[3]
    Y23 = U[4]
    X31 = U[5]
    Y31 = U[6]

    r12 = X12^2+Y12^2    
    r23 = X23^2+Y23^2    
    r31 = X31^2+Y31^2
    return 1/3 * (r12^2+r23^2+r31^2)
end




"""
    g_col(u)

Hiru gorputzen problemaren 'u' egoerako gorputzen posizioek sortzen duten paralelogramoaren zeinudun azalera itzultzen du.

#Arguments
    - u::Vector{Float64}: Hiru gorputzeko problemari dagokion egoera, posizioak erlatiboak eta abiadurak absolutuak izanik
"""
function g_col(u)
    x12 = u[1]
    y12 = u[2]  

    x23 = u[3]
    y23 = u[4]  

    x21 = -x12
    y21 = -y12

    cross = x21*y23 - x23*y21
    return cross
end



function getCode(UU, g_col_u)
    code_middle = Int[]
    code_far = Int[] 
    ind_kolinear = Int[]
    signgu_ = 0
    for i in 1:length(g_col_u)
        signgu = sign(g_col_u[i])

        if(g_col_u[i]==0. || (signgu_!=0 && signgu_!= signgu)) 
            far, middle = getFarMiddle(UU[i])
            
            push!(code_middle,middle)
            push!(code_far,far)
            push!(ind_kolinear, i)
        end
        signgu_ = signgu
    end
    return code_middle,code_far, ind_kolinear
end