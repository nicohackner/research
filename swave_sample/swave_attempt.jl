Nlat = 31;
itermax = 500;
Delta = ones(ComplexF64, Nlat, Nlat);
Deltanew = ones(ComplexF64, Nlat, Nlat);
Ham = ones(ComplexF64, 2 * Nlat^2, 2 * Nlat^2);

t = 1.0;
temp = 0.001;
mu = -1.0;
Vimp = 5.0;
Uinter = 1.8;

middle = (Nlat - 1) / 2 + 1;

for iter in 1:itermax
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
                Ham[Nlat^2 + (ix - 1) * Nlat + iy,Nlat^2 + (ix - 1) * Nlat + iy + 1] = t
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
            
            # pairing
            Ham[(ix - 1) * Nlat + iy,Nlat^2 + (ix - 1) * Nlat + iy] = Delta(ix, iy)

            Ham[Nlat^2 + (ix - 1) * Nlat + iy,(ix - 1) * Nlat + iy] = conj(Delta(ix, iy))

        end
    end

    # impurity
    Ham[(middle - 1) * Nlat + middle, (middle - 1) * Nlat + middle] = Vimp
    Ham[Nlat^2 + (middle - 1) * Nlat + middle, Nlat^2 + (middle - 1) * Nlat + middle] = -Vimp




