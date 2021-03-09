using Plots, LinearAlgebra;

function potential()
    
    Nlat = 31;
    itermax = 500;
    Delta = zeros(ComplexF64, Nlat, Nlat);
    evals = zeros(ComplexF64, 2 * Nlat^2);

    t = 1.0;
    temp = 0.001;
    mu = -1.0;
    Vimp = 5.0;
    Uinter = 1.8;
    tol = 1e-6;
    counter = 0;
    alpha = 0.0;

    middle = convert(Int64,(Nlat - 1) / 2 + 1);

    for iter in 1:itermax
        Ham = zeros(ComplexF64, 2 * Nlat^2, 2 * Nlat^2);
        for ix in 1:Nlat
            for iy in 1:Nlat

                if iter == 1 # initialize random potential
                    Delta[ix,iy] = rand(ComplexF64)
                end
                
                # chemical potential

                # electron
                Ham[(ix - 1) * Nlat + iy,(ix - 1) * Nlat + iy] = -mu

                # hole
                Ham[Nlat^2 + (ix - 1) * Nlat + iy,Nlat^2 + (ix - 1) * Nlat + iy] = mu

                # hopping in y direction
                if iy == Nlat # periodic boundary
                    # electron
                    Ham[(ix - 1) * Nlat + iy,(ix - 1) * Nlat + 1] = -t

                    # hole
                    Ham[Nlat^2 + (ix - 1) * Nlat + iy,Nlat^2 + (ix - 1) * Nlat + 1] = t
                else
                    # electron
                    Ham[(ix - 1) * Nlat + iy,(ix - 1) * Nlat + iy + 1] = -t

                    # hole
                    Ham[Nlat^2 + (ix - 1) * Nlat + iy,Nlat^2 + (ix - 1) * Nlat + iy + 1] = t
                end 
                
                if iy == 1 # periodic boundary
                    # electron
                    Ham[(ix - 1) * Nlat + iy,(ix - 1) * Nlat + Nlat] = -t

                    # hole
                    Ham[Nlat^2 + (ix - 1) * Nlat + iy,Nlat^2 + (ix - 1) * Nlat + Nlat] = t
                else
                    # electron
                    Ham[(ix - 1) * Nlat + iy,(ix - 1) * Nlat + iy - 1] = -t

                    # hole
                    Ham[Nlat^2 + (ix - 1) * Nlat + iy,Nlat^2 + (ix - 1) * Nlat + iy - 1] = t
                end

                if ix == Nlat
                    # electron  
                    Ham[(ix - 1) * Nlat + iy,iy] = -t

                    # hole
                    Ham[Nlat^2 + (ix - 1) * Nlat + iy,Nlat^2 + iy] = t
                else
                    # electron  
                    Ham[(ix - 1) * Nlat + iy,(ix) * Nlat + iy] = -t

                    # hole
                    Ham[Nlat^2 + (ix - 1) * Nlat + iy,Nlat^2 + (ix) * Nlat + iy] = t
                end
                
                if ix == 1
                    # electron  
                    Ham[(ix-1)*Nlat+iy,(Nlat-1)*Nlat+iy] = -t

                    # hole
                    Ham[Nlat^2+(ix-1)*Nlat+iy,Nlat^2+(Nlat-1)*Nlat+iy] = t
                else
                    # electron  
                    Ham[(ix-1)*Nlat+iy,(ix-2)*Nlat+iy] = -t

                    # hole
                    Ham[Nlat^2+(ix-1)*Nlat+iy,Nlat^2+(ix-2)*Nlat+iy] = t
                end

                # pairing
                Ham[(ix - 1) * Nlat + iy,Nlat^2 + (ix - 1) * Nlat + iy] = Delta[ix, iy]

                Ham[Nlat^2 + (ix - 1) * Nlat + iy,(ix - 1) * Nlat + iy] = conj(Delta[ix, iy])

            end
        end
        
       
        # impurity
        Ham[(middle - 1) * Nlat + middle, (middle - 1) * Nlat + middle] = Vimp
        Ham[Nlat^2 + (middle - 1) * Nlat + middle, Nlat^2 + (middle - 1) * Nlat + middle] = -Vimp
    
    

    evals = eigvals(Ham)
    Ham = eigvecs(Ham)

    Deltanew = zeros(ComplexF64, Nlat, Nlat);

        for ix in 1:Nlat
            for iy in 1:Nlat
                for m1 in 1:2*Nlat^2
                    m2 = (ix-1)*Nlat + iy
                    Deltanew[ix,iy] = Deltanew[ix,iy] - Uinter*Ham[m2,m1]*conj(Ham[m2+Nlat^2,m1])*fermi(evals[m1],temp)
                end

                Deltanew[ix,iy] = alpha*Delta[ix,iy] + (1.0 - alpha)*Deltanew[ix,iy]

                if abs(Deltanew[ix,iy] - Delta[ix,iy]) > tol
                    counter += 1
                end
            end
        end

        if counter == 0 && iter > 5 
            print("solution converged")
            break
        else
            
            print([iter, counter, Nlat^2])
            println(Delta[1,1])
            counter = 0

            Delta = Deltanew

        end
    end
    return real(Delta)
end

function fermi(evals,temp)
    if evals > 0.0
        return exp(-evals/temp)/(1.0 + exp(-evals/temp))
    else
        return 1.0/(1.0 + exp(evals/temp))
    end
end

hea