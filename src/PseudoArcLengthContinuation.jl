


struct solucion{ttype,utype} 
    t::ttype 
    x::utype 
    retcode::Symbol 
end


function PseudoArcLength(F2_fun::Functype, J2_fun::Jtype, x0::utype, Beta0::utype, p::ptype, 
#    t0::ttype, tF::ttype,
    maxiters=100, tolx=1e-10, tolDelta=1e-10, lambda=1, trace=false) where {Functype, Jtype, utype, ptype}

    N=length(x0)
    xj=copy(x0)
    Betaj=copy(Beta0)
    tj=0    
    dtj=1e-4

    Dx=similar(x0)
    F2=similar(x0)
    J2=Array{eltype(x0)}(undef,N,N)
    for i in 1:N
        for j in 1:N
            J2[i,j]=zero(eltype(x0))
        end
    end 

    # Calcular Beta0

    J2_fun(J2,xj,p,Betaj,x0,0.)
    # compute beta
    b=zero(x0)
    b[end]=1
    Betaj=J2\b
    Betaj=Betaj/norm(Betaj)

    # Asegurarse que x0 está en la curva

    iters=0
    dt=zero(eltype(x0))
    F2_fun(F2,xj,p, Betaj,x0,0.)
    if norm(F2)>0 
        Dx=zero(x0)
        Dx[1]=10
        while (norm(Dx)>1e-14 && iters<maxiters)
            iters=iters+1
            J2_fun(J2,xj,p, Betaj,x0,0.)
            lu_fact=lu(J2)
            F2_fun(F2,xj,p, Beta0,x0,0.)
            b=-F2
            Dx=lu_fact\b
            xj=xj+Dx
       end
    end

    if iters==maxiters
        println("error")
        sol = solucion(tj, xj, :Failure) 
        return(sol)
    end 

    println("x0=", x0, ", norm(Dx)=", norm(Dx))

    find=false 


    """
    find=false 
    while !find

        computeBeta!(J, Betaj)
        J2_fun(J2,x,Betaj,xj,dt)

        stop=false

        LU=descomp_LU(J)
        iter=0
        while !stop
        
            iter=iter+1
            b=-F(x0,p,Betaj,xj,0)
            solveLU!(Deltax,LU,b)
            @. xj=xj+Deltax

            if iter==1 
                normDelta0=norm(Deltax) 
            end
                
            stop=norm(Deltax)<tolx || iter>maxiters  

        end

        find=true
    end
    """
    sol = solucion(tj, xj, :Successs) 
    return(sol)

end