
function gen_scaffold(conc, vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed)
    # Generates scaffold based on experimental data, scaled down by 'scale'
    # Inputs
    #  conc: initial concentration of MSRs in solution (μg/mL)
    #  vol: initial volume of solution (mL)
    #  max_d: maximum number of micro-rods overlapped at any lattice site in the domain
    #  N_cell_2D: number of cells stacked vertically to approximate this 3D domain in 2D
    #  scale: scales down the domain area by the proportion scale ∈ [0,1]
    #  linemethod: method for generating the lines in scaffolds: "original" or "XiaolinWu"
    #  xpos_PDF(M): probability density function (for x) that governs the x positions of the centre of micro-rods - a function of the domain width
    #  ypos_PDF(N): probability density function (for y) that governs the y positions of the centre of micro-rods - a function of the domain height
    #  rot_PDF: probability density function (of θ, in rad) that governs the rotation of micro-rods
    #  blurring: scale for grid size, in which cells represent the average of the pixels within the cells
    #  seed: seed for random number generation (same seed will generate the same scaffold, NaN to use default seed)

    if !isnan(seed) # set seed if not using default seed
        Random.seed!(seed); 
    end

    well = "24well"; # dish size ("24well", "6well", "10cmpetri")

    # parameters
    T_diam = 5; # diameter of T cell (μm) (small estimate)
    A_24well = 1.9; # surface area of a well in a 24 well culture plate (cm^2)   (“Useful Numbers for Cell Culture - AU.” Accessed: Sep. 01, 2024. [Online]. Available: https://www.thermofisher.com/au/en/home/references/gibco-cell-culture-basics/cell-culture-protocols/cell-culture-useful-numbers.html)
    A_6well = 9.6; # surface area of a well in a 6 well culture plate (cm^2)       ^
    A_10cmpetri = 56.7; # surface area of a 10cm petri dish (cm^2)                 ^

    h = T_diam; # distance between nodes in grid (μm)

    MSR_len = 70.4 .+ [-31.7, 31.7]; # range of lengths the MSRs might have (μm)   (D. K. Y. Zhang, A. S. Cheung, and D. J. Mooney, “Activation and expansion of human T cells using artificial antigen-presenting cell scaffolds,” Nat Protoc, vol. 15, no. 3, pp. 773–798, Mar. 2020, doi: 10.1038/s41596-019-0249-0.)
    MSR_diam = 5.1 .+ [-1.2, 1.2]; # range of diameters the MSRs might have (μm)       ^
    MSR_SA = 572*1000000; # surface area of MSR per gram (μm^2/μg)                                       ^
    MSR_m = π*mean(MSR_diam)*mean(MSR_len)/MSR_SA; # mass of individual MSR (μg)

    if well == "24well"
        N = M = Int(floor(sqrt(A_24well)*10000/T_diam*sqrt(scale))); # side length of square simulation domain (based on chosen dish size) assuming each grid cell is the size of one T cell
        h_well = vol/A_24well*10000; # height that the solution reaches in the chosen well/dish (μm)
    elseif well == "6well"
        N = M = Int(floor(sqrt(A_6well)*10000/T_diam*sqrt(scale)));
        h_well = vol/A_6well*10000; # (μm)
    elseif well == "10cmpetri"
        N = M = Int(floor(sqrt(A_10cmpetri)*10000/T_diam*sqrt(scale)));
        h_well = vol/A_10cmpetri*10000; # (μm)
    end

    wid = (M-1)*T_diam; # total domain width and height
    hei = (N-1)*T_diam;

    xpos_PDFM = xpos_PDF(M); # define the probability distribution functions based off the end points in the domain
    ypos_PDFN = ypos_PDF(N);

    # number of rods in the starting domain (input concentration * input volume / MSR mass (gives number of rods in total) scaled by the ratio between N_cell_2D T cell diameters and the height of the well (attempt to make somewhat 2D) * 'scale'
    num_rods = Int(round(conc*vol/MSR_m*(N_cell_2D*T_diam/h_well)*scale));


    rot_arr = zeros(num_rods); # array containing the rotation angles of each micro-rod

    # generate domain
    ρ = zeros(N,M); # initialise scaffold domain
    @showprogress 1 "Generating scaffold..." for rod = 1:num_rods # for each rod, place them randomly in the domain
        #rod_ind = Int.(floor.([rand()*(N-1), rand()*(M-1)].+1)); # random starting indices for rod
        #x0 = rod_ind[2]; y0 = rod_ind[1]; 
        x0 = rand(xpos_PDFM); y0 = rand(ypos_PDFN); # random starting indices for rod (fractional)
        rod_len = max((rand()*(MSR_len[2]-MSR_len[1])+MSR_len[1])/T_diam,1); # random length for current rod in terms of grid lengths (minimum of 1)
        #rod_len_Int = Int(ceil(rod_len));
        rod_wid = max((rand()*(MSR_diam[2]-MSR_diam[1])+MSR_diam[1])/T_diam,1); # random width for current rod in terms of grid lengths (minimum of 1)
        rod_dir = rand(rot_PDF); # current rod direction in grid domain (starting at rod_ind)
        rot_arr[rod] = rod_dir; # record rod direction

        if linemethod == "XiaolinWu" # Xiaolin Wu's line algorithm    (DOESNT ACCOUNT FOR LINE WIDTH)
            x1 = Int(floor(x0+rod_len*cos(rod_dir))); y1 = Int(floor(y0+rod_len*sin(rod_dir))); # end points of line
            
            if abs(y1-y0) < abs(x1-x0) # if the rod is longer in the x direction
                if x1 < x0 # if the rod is being drawn backwards
                    x0, x1 = x1, x0; # swap start and end points
                    y0, y1 = y1, y0;
                end

                dx = x1-x0; dy = y1-y0; # changes in x and y
                m = dy/dx; # line gradient

                for i = 0:dx # for each point along the line
                    x = x0+i; # define x and y coordinates
                    y = y0+i*m; 

                    if x < 1 # 'wrap around' the rods if they exit the domain
                        x += M;
                    elseif x >= M+1
                        x -= M;
                    end
                    if y < 1
                        y += N;
                    elseif y >= N+1
                        y -= N;
                    end

                    x_int = min(max(Int(round(x)),1),M); 
                    y_int = min(max(Int(round(y)),1),N);
                    dist = y-y_int; # define distance between y and its floor 

                    y_up = y_int+1; # second pixel for plotting
                    if y_up == N+1 # if the second pixel is outside the domain, wrap around
                        y_up = 1;
                    end

                    ρ[y_int,x_int] = min(ρ[y_int,x_int]+(1-dist), max_d); # add scaled rod density to this point and the one below it, ensuring cumulative density does not exceed 1
                    ρ[y_up,x_int] = min(ρ[y_up,x_int]+dist, max_d);
                end
            else # the rod is longer in the y direction
                if y1 < y0 # if the rod is being drawn downwards
                    x0, x1 = x1, x0; # swap start and end points
                    y0, y1 = y1, y0;
                end

                dx = x1-x0; dy = y1-y0; # changes in x and y
                m = dx/dy; # line gradient

                for i = 0:dy # for each point along the line
                    x = x0+i*m; # define x and y coordinates
                    y = y0+i; 

                    if x < 1 # 'wrap around' the rods if they exit the domain
                        x += M;
                    elseif x >= M+1
                        x -= M;
                    end
                    if y < 1
                        y += N;
                    elseif y >= N+1 
                        y -= N;
                    end

                    x_int = min(max(Int(round(x)),1),M); 
                    y_int = min(max(Int(round(y)),1),N);
                    dist = x-x_int; # define distance between x and its floor 

                    x_up = x_int+1; # second pixel for plotting
                    if x_up == M+1 # if the second pixel is outside the domain, wrap around
                        x_up = 1;
                    end

                    ρ[y_int,x_int] = min(ρ[y_int,x_int]+(1-dist), max_d); # add scaled rod density to this point and the one below it, ensuring cumulative density does not exceed max_d
                    ρ[y_int,x_up] = min(ρ[y_int,x_up]+dist, max_d);
                end
            end


        elseif linemethod == "original" # original method I implemented
            a = tan(rod_dir); b = -1; c = y0-tan(rod_dir)*x0; # equation of line ax+by+c=0 (tan(rod_dir)*ind[2] - ind[1] + y0-tan(rod_dir)*x0 = 0)

            xind_range = Int(round(x0-rod_len*2)):Int(round(x0+rod_len*2)); # range of indices to check over for rod placement
            yind_range = Int(round(y0-rod_len*2)):Int(round(y0+rod_len*2));

            # for all potential cells the rod may fall on
            for x_ind = xind_range
                for y_ind = yind_range 
                    dist_n = abs(a*x_ind+b*y_ind+c)/sqrt(a^2+b^2); # distance from new indices to the closest point on the line
                    if dist_n <= rod_wid/2 && sqrt((x_ind-x0)^2+(y_ind-y0)^2) <= rod_len/2 # if the point is within half the width and length of the rod
                        x_ind_n = x_ind;
                        if x_ind_n < 1 # 'wrap around' the rods if they exit the domain
                            x_ind_n += M;
                        elseif x_ind_n > M
                            x_ind_n -= M;
                        end
                        y_ind_n = y_ind;
                        if y_ind_n < 1
                            y_ind_n += N;
                        elseif y_ind_n > N
                            y_ind_n -= N;
                        end
                        ρ[y_ind_n,x_ind_n] = min(ρ[y_ind_n,x_ind_n]+1/(rod_len*rod_wid), max_d); # add rod density to this point, ensuring cumulative density does not exceed max_d
                    end
                end
            end
        end
    end

    # convolve the scaffold using a uniform square kernel with side length "blurring"
    if blurring > 1 && isodd(blurring) 
        ρ_new = zeros(N,M); # new scaffold that has been blurring using convolution
        padding = (blurring-1)÷2;
        for i = 1:N
            for j = 1:M
                for ii = i-padding:i+padding
                    if ii < 1; ii += N; elseif ii > N; ii -= N; end # handle edge cases (periodic BC)
                    for jj = j-padding:j+padding
                        if jj < 1; jj += M; elseif jj > M; jj -= M; end
                        ρ_new[i,j] += ρ[ii, jj]/blurring^2;
                    end
                end
            end
        end
        ρ = ρ_new; # replace scaffold with its convolution
    end


    # compute probability density of rotation angles
    bins = length(rot_arr)÷2;
    rot_dens = zeros(bins);
    for rod = eachindex(rot_arr)
        rot_n = rot_arr[rod];
        while rot_n > π
            rot_n -= π;
        end
        while rot_n < 0
            rot_n += π
        end
        rot_dens[Int(round(rot_n/π*(bins-1)))+1] += 1/length(rot_arr);
    end


    Random.seed!(); # reset seed back to default


    # generate other variables that depend on ρ
    result_scale = vol/(wid*hei*N_cell_2D*T_diam)*1e12; # ratio of total number of rods to those in the simulation domain (assumed all results also scale by this amount)
    ρ_max = maximum(ρ)+1*(maximum(ρ)==0); # maximum value of ρ (1 if there is no scaffold)
    V = N*M*result_scale; # full domain volume
    ρ_avg = mean(ρ); # average rod density: rods/mL times ratio mL per lattice site, alternatively N_rods/(N*M*result_scale)
    N_rods = ρ_avg*N*M*result_scale; # APC_ms_conc*APC_ms_vol/m_MSR; # total number of rods in the domain (integral of ρ(x,y))
    
    if approx == "finite" # finite differences to approximate ∇ρ, has some error in sum of probabilities
        if BC == "noflux"
            ρy_d = [zeros(1,M); (ρ[2:N,:]-ρ[1:N-1,:])/h]; # backward finite differences for down, derivative at bottom edge is 0
            ρy_u = [ρy_d[2:N,:]; zeros(1,M)]; # forward finite differences for up
            ρx_l = [zeros(N,1) (ρ[:,2:M]-ρ[:,1:M-1])/h]; # backward finite differences for left
            ρx_r = [ρx_l[:,2:M] zeros(N,1)]; # forward finite differences for right
        elseif BC == "periodic"
            ρy_d = [reshape((ρ[1,:]-ρ[N,:])/h,(1,:)); (ρ[2:N,:]-ρ[1:N-1,:])/h]; # backward finite differences for down, bottom edge wraps around to top edge
            ρy_u = [ρy_d[2:N,:]; reshape(ρy_d[1,:],(1,:))]; # forward finite differences for up
            ρx_l = [(ρ[:,1]-ρ[:,M])/h (ρ[:,2:M]-ρ[:,1:M-1])/h]; # backward finite differences for left
            ρx_r = [ρx_l[:,2:M] ρx_l[:,1]]; # forward finite differences for right
        end
    elseif approx == "central" # central differences to approximate ∇ρ, has no error
        if BC == "noflux"
            ρy_d = ρy_u = [zeros(1,M); (ρ[3:N,:]-ρ[1:N-2,:])/2h; zeros(1,M)]; # central differences for up and down, derivatives at top and bottom edges are 0
            ρx_l = ρx_r = [zeros(N,1) (ρ[:,3:M]-ρ[:,1:M-2])/2h zeros(N,1)]; # central differences for left and right
        elseif BC == "periodic"
            ρy_d = ρy_u = [reshape((ρ[2,:]-ρ[N,:])/2h,(1,:)); (ρ[3:N,:]-ρ[1:N-2,:])/2h; reshape((ρ[1,:]-ρ[N-1,:])/2h,(1,:))]; # central differences for up and down, top and bottom wrap around to each other
            ρx_l = ρx_r = [(ρ[:,2]-ρ[:,M])/2h (ρ[:,3:M]-ρ[:,1:M-2])/2h (ρ[:,1]-ρ[:,M-1])/2h]; # central differences for left and right
        end
    end
    ρ_bar = max(maximum(abs.(ρy_d)), maximum(abs.(ρx_l))); # maximum out of all finite differences
    if ρ_bar == 0
        ρ_bar = 1; # if ρ is a space-invariant field (constant), set max(ρ_xy)=1 to avoid division by 0
    end
    p_d = (1 .- β*ρy_d/ρ_bar)/4; # density-dependent probabilities of moving down/up/left/right (without considering wait probability - this is included in simulation)
    p_u = (1 .+ β*ρy_u/ρ_bar)/4;
    p_l = (1 .- β*ρx_l/ρ_bar)/4;
    p_r = (1 .+ β*ρx_r/ρ_bar)/4;


    return reverse(ρ,dims=1), rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r
end

function move_cell(x, y, wid, hei, ε, BC, p_d, p_u, p_l, p_r)
    # Outputs the resulting position (x,y) of the current cell after movement
    # Inputs:
    #  x: cell's current x position (μm)
    #  y: cell's current y position (μm)
    #  wid: width of domain (x direction) (μm)
    #  hei: height of domain (y direction) (μm)
    #  ε: cell's step length (μm)
    #  BC: boundary condition ("noflux" or "periodic")
    #  p_d: current probability of moving down
    #  p_u: current probability of moving up
    #  p_l: current probability of moving left
    #  p_r: current probability of moving right
    # Outputs
    #  x, y: next x and y coordinates of cell 

    if BC == "noflux" # if BC is no-flux
        if y - ε < 0 # if the cell can't move down (stay within scaffold domain)
            p_d = 0; 
        elseif y + ε > hei # up
            p_u = 0;
        end
        if x - ε < 0 # left
            p_l = 0;
        elseif x + ε > wid # right
            p_r = 0;
        end
    end

    # cell movement
    r = rand(); # sample random number ∈ [0,1]
    if r < p_d # check probabilities of moving in any direction
        y -= ε;
        if y < 0 && BC == "periodic"
            y += hei;
        end
    elseif p_d <= r < p_d+p_u
        y += ε;
        if y > hei && BC == "periodic"
            y -= hei;
        end
    elseif p_d+p_u <= r < p_d+p_u+p_l
        x -= ε;
        if x < 0 && BC == "periodic"
            x += wid;
        end
    elseif p_d+p_u+p_l <= r < p_d+p_u+p_l+p_r
        x += ε;
        if x > wid && BC == "periodic"
            x -= wid;
        end
    end # else the cell waits at the current position

    return x, y
end

function genA_n(N, M, h, D_n, χ_n, R_nn, ρ, BC)
    # Generates the A matrix that holds the coefficients in the ODE system (nodes from FVM) for u_n
    #  original PDE is:  ∂u_n/∂t = ∇•[D_n(ρ)∇u_n - χ_n(ρ)∇ρ*u_n] + R_nn(ρ)u_n + R_na(c)u_a
    # Inputs:
    #  N: width of discretised grid (number of nodes)
    #  M: height of discretised grid (number of nodes)
    #  h: distance between nodes in discretised grid
    #  D_n: scaffold-dependent diffusivity for u_n
    #  χ_n: scaffold-dependent haptotactic sensitivity for u_n (full advection term is -χ_n∇ρ)
    #  R_nn: scaffold-dependent reaction term driven by naive cells for u_n
    #  ρ: matrix containing scaffold density in the domain
    #  BC: boundary condition: "noflux" or "periodic"
    # Output:
    #  A: sparse A matrix containing all coefficients

    ind(i,j) = i + N*(j-1); # equation to determine index in A for the (i,j)th node

    if BC == "noflux"
        nonzero_entries = (N-2)*(M-2)*5+2*(N+M-4)*4+4*3; # number of equations * number of coefficients  for  interior+edges+corners
    elseif BC == "periodic"
        nonzero_entries = (N-2)*(M-2)*5+2*(N+M-4)*5+4*5;
    end
    
    i_vec = zeros(nonzero_entries,1)[:]; # vector containing all i indices of non-zero entries into A
    j_vec = zeros(nonzero_entries,1)[:]; # vector containing all j indices of non-zero entries into A
    A_vec = zeros(nonzero_entries,1)[:]; # vector containing all non-zero entries into A
    count = 1; # counter for elements added to vectors

    # interior nodes
    for i = 2:N-1
        for j = 2:M-1
            i_vec[count:count+4] .= ind(i,j);
            j_vec[count:count+4] = [ind(i,j), ind(i-1,j), ind(i+1,j), ind(i,j-1), ind(i,j+1)];
            A_vec[count:count+4] = 1/h^2*[-(4D_n(ρ[i,j])+D_n(ρ[i-1,j])+D_n(ρ[i+1,j])+D_n(ρ[i,j-1])+D_n(ρ[i,j+1]))/2 + ((χ_n(ρ[i,j])+χ_n(ρ[i-1,j]))*(ρ[i,j]-ρ[i-1,j])-(χ_n(ρ[i,j])+χ_n(ρ[i+1,j]))*(ρ[i+1,j]-ρ[i,j])+(χ_n(ρ[i,j])+χ_n(ρ[i,j-1]))*(ρ[i,j]-ρ[i,j-1])-(χ_n(ρ[i,j])+χ_n(ρ[i,j+1]))*(ρ[i,j+1]-ρ[i,j]))/4 + R_nn(ρ[i,j])*h^2, (D_n(ρ[i,j])+D_n(ρ[i-1,j]))/2+(χ_n(ρ[i,j])+χ_n(ρ[i-1,j]))*(ρ[i,j]-ρ[i-1,j])/4, (D_n(ρ[i,j])+D_n(ρ[i+1,j]))/2-(χ_n(ρ[i,j])+χ_n(ρ[i+1,j]))*(ρ[i+1,j]-ρ[i,j])/4, (D_n(ρ[i,j])+D_n(ρ[i,j-1]))/2+(χ_n(ρ[i,j])+χ_n(ρ[i,j-1]))*(ρ[i,j]-ρ[i,j-1])/4, (D_n(ρ[i,j])+D_n(ρ[i,j+1]))/2-(χ_n(ρ[i,j])+χ_n(ρ[i,j+1]))*(ρ[i,j+1]-ρ[i,j])/4];
            count = count + 5; 
        end
    end

    if BC == "noflux" # if no-flux boundary condition
        # left and right edge
        for j = 2:M-1
            i_vec[count:count+3] .= ind(1,j);
            j_vec[count:count+3] = [ind(1,j), ind(2,j), ind(1,j-1), ind(1,j+1)];
            A_vec[count:count+3] = 2/h^2*[-(2D_n(ρ[1,j])+D_n(ρ[2,j])+D_n(ρ[1,j-1])/2+D_n(ρ[1,j+1])/2)/2 + (-(χ_n(ρ[1,j])+χ_n(ρ[2,j]))*(ρ[2,j]-ρ[1,j])+(χ_n(ρ[1,j])+χ_n(ρ[1,j-1]))*(ρ[1,j]-ρ[1,j-1])/2-(χ_n(ρ[1,j])+χ_n(ρ[1,j+1]))*(ρ[1,j+1]-ρ[1,j])/2)/4 + R_nn(ρ[1,j])*h^2/2, (D_n(ρ[1,j])+D_n(ρ[2,j]))/2-(χ_n(ρ[1,j])+χ_n(ρ[2,j]))*(ρ[2,j]-ρ[1,j])/4, ((D_n(ρ[1,j])+D_n(ρ[1,j-1]))/2+(χ_n(ρ[1,j])+χ_n(ρ[1,j-1]))*(ρ[1,j]-ρ[1,j-1])/4)/2, ((D_n(ρ[1,j])+D_n(ρ[1,j+1]))/2-(χ_n(ρ[1,j])+χ_n(ρ[1,j+1]))*(ρ[1,j+1]-ρ[1,j])/4)/2];
            count = count + 4; 
        
            i_vec[count:count+3] .= ind(N,j);
            j_vec[count:count+3] = [ind(N,j), ind(N-1,j), ind(N,j-1), ind(N,j+1)];
            A_vec[count:count+3] = 2/h^2*[-(2D_n(ρ[N,j])+D_n(ρ[N-1,j])+D_n(ρ[N,j-1])/2+D_n(ρ[N,j+1])/2)/2 + ((χ_n(ρ[N,j])+χ_n(ρ[N-1,j]))*(ρ[N,j]-ρ[N-1,j])+(χ_n(ρ[N,j])+χ_n(ρ[N,j-1]))*(ρ[N,j]-ρ[N,j-1])/2-(χ_n(ρ[N,j])+χ_n(ρ[N,j+1]))*(ρ[N,j+1]-ρ[N,j])/2)/4 + R_nn(ρ[N,j])*h^2/2, (D_n(ρ[N,j])+D_n(ρ[N-1,j]))/2+(χ_n(ρ[N,j])+χ_n(ρ[N-1,j]))*(ρ[N,j]-ρ[N-1,j])/4, ((D_n(ρ[N,j])+D_n(ρ[N,j-1]))/2+(χ_n(ρ[N,j])+χ_n(ρ[N,j-1]))*(ρ[N,j]-ρ[N,j-1])/4)/2, ((D_n(ρ[N,j])+D_n(ρ[N,j+1]))/2-(χ_n(ρ[N,j])+χ_n(ρ[N,j+1]))*(ρ[N,j+1]-ρ[N,j])/4)/2];
            count = count + 4;
        end 

        # bottom and top edges
        for i = 2:N-1
            i_vec[count:count+3] .= ind(i,1);
            j_vec[count:count+3] = [ind(i,1), ind(i-1,1), ind(i+1,1), ind(i,2)];
            A_vec[count:count+3] = 2/h^2*[-(2D_n(ρ[i,1])+D_n(ρ[i-1,1])/2+D_n(ρ[i+1,1])/2+D_n(ρ[i,2]))/2 + ((χ_n(ρ[i,1])+χ_n(ρ[i-1,1]))*(ρ[i,1]-ρ[i-1,1])/2-(χ_n(ρ[i,1])+χ_n(ρ[i+1,1]))*(ρ[i+1,1]-ρ[i,1])/2-(χ_n(ρ[i,1])+χ_n(ρ[i,2]))*(ρ[i,2]-ρ[i,1]))/4 + R_nn(ρ[i,1])*h^2/2, ((D_n(ρ[i,1])+D_n(ρ[i-1,1]))/2+(χ_n(ρ[i,1])+χ_n(ρ[i-1,1]))*(ρ[i,1]-ρ[i-1,1])/4)/2, ((D_n(ρ[i,1])+D_n(ρ[i+1,1]))/2-(χ_n(ρ[i,1])+χ_n(ρ[i+1,1]))*(ρ[i+1,1]-ρ[i,1])/4)/2, (D_n(ρ[i,1])+D_n(ρ[i,2]))/2-(χ_n(ρ[i,1])+χ_n(ρ[i,2]))*(ρ[i,2]-ρ[i,1])/4];
            count = count + 4;
        
            i_vec[count:count+3] .= ind(i,M);
            j_vec[count:count+3] = [ind(i,M), ind(i-1,M), ind(i+1,M), ind(i,M-1)];
            A_vec[count:count+3] = 2/h^2*[-(2D_n(ρ[i,M])+D_n(ρ[i-1,M])/2+D_n(ρ[i+1,M])/2+D_n(ρ[i,M-1]))/2 + ((χ_n(ρ[i,M])+χ_n(ρ[i-1,M]))*(ρ[i,M]-ρ[i-1,M])/2-(χ_n(ρ[i,M])+χ_n(ρ[i+1,M]))*(ρ[i+1,M]-ρ[i,M])/2+(χ_n(ρ[i,M])+χ_n(ρ[i,M-1]))*(ρ[i,M]-ρ[i,M-1]))/4 + R_nn(ρ[i,M])*h^2/2, ((D_n(ρ[i,M])+D_n(ρ[i-1,M]))/2+(χ_n(ρ[i,M])+χ_n(ρ[i-1,M]))*(ρ[i,M]-ρ[i-1,M])/4)/2, ((D_n(ρ[i,M])+D_n(ρ[i+1,M]))/2-(χ_n(ρ[i,M])+χ_n(ρ[i+1,M]))*(ρ[i+1,M]-ρ[i,M])/4)/2, (D_n(ρ[i,M])+D_n(ρ[i,M-1]))/2+(χ_n(ρ[i,M])+χ_n(ρ[i,M-1]))*(ρ[i,M]-ρ[i,M-1])/4];
            count = count + 4;
        end

        # corners
        i_vec[count:count+2] .= ind(1,1);
        j_vec[count:count+2] = [ind(1,1), ind(2,1), ind(1,2)];
        A_vec[count:count+2] = 2/h^2*[-(2D_n(ρ[1,1])+D_n(ρ[2,1])+D_n(ρ[1,2]))/2 + (-(χ_n(ρ[1,1])+χ_n(ρ[2,1]))*(ρ[2,1]-ρ[1,1])-(χ_n(ρ[1,1])+χ_n(ρ[1,2]))*(ρ[1,2]-ρ[1,1]))/4 + R_nn(ρ[1,1])*h^2/2, (D_n(ρ[1,1])+D_n(ρ[2,1]))/2-(χ_n(ρ[1,1])+χ_n(ρ[2,1]))*(ρ[2,1]-ρ[1,1])/4, (D_n(ρ[1,1])+D_n(ρ[1,2]))/2-(χ_n(ρ[1,1])+χ_n(ρ[1,2]))*(ρ[1,2]-ρ[1,1])/4];
        count = count+3; 
        
        i_vec[count:count+2] .= ind(N,1);
        j_vec[count:count+2] = [ind(N,1), ind(N-1,1), ind(N,2)];
        A_vec[count:count+2] = 2/h^2*[-(2D_n(ρ[N,1])+D_n(ρ[N-1,1])+D_n(ρ[N,2]))/2 + ((χ_n(ρ[N,1])+χ_n(ρ[N-1,1]))*(ρ[N,1]-ρ[N-1,1])-(χ_n(ρ[N,1])+χ_n(ρ[N,2]))*(ρ[N,2]-ρ[N,1]))/4 + R_nn(ρ[N,1])*h^2/2, (D_n(ρ[N,1])+D_n(ρ[N-1,1]))/2+(χ_n(ρ[N,1])+χ_n(ρ[N-1,1]))*(ρ[N,1]-ρ[N-1,1])/4, (D_n(ρ[N,1])+D_n(ρ[N,2]))/2-(χ_n(ρ[N,1])+χ_n(ρ[N,2]))*(ρ[N,2]-ρ[N,1])/4];
        count = count+3; 
        
        i_vec[count:count+2] .= ind(1,M);
        j_vec[count:count+2] = [ind(1,M), ind(2,M), ind(1,M-1)];
        A_vec[count:count+2] = 2/h^2*[-(2D_n(ρ[1,M])+D_n(ρ[2,M])+D_n(ρ[1,M-1]))/2 + (-(χ_n(ρ[1,M])+χ_n(ρ[2,M]))*(ρ[2,M]-ρ[1,M])+(χ_n(ρ[1,M])+χ_n(ρ[1,M-1]))*(ρ[1,M]-ρ[1,M-1]))/4 + R_nn(ρ[1,M])*h^2/2, (D_n(ρ[1,M])+D_n(ρ[2,M]))/2-(χ_n(ρ[1,M])+χ_n(ρ[2,M]))*(ρ[2,M]-ρ[1,M])/4, (D_n(ρ[1,M])+D_n(ρ[1,M-1]))/2+(χ_n(ρ[1,M])+χ_n(ρ[1,M-1]))*(ρ[1,M]-ρ[1,M-1])/4];
        count = count+3; 
        
        i_vec[count:count+2] .= ind(N,M);
        j_vec[count:count+2] = [ind(N,M), ind(N-1,M), ind(N,M-1)];
        A_vec[count:count+2] = 2/h^2*[-(2D_n(ρ[N,M])+D_n(ρ[N-1,M])+D_n(ρ[N,M-1]))/2 + ((χ_n(ρ[N,M])+χ_n(ρ[N-1,M]))*(ρ[N,M]-ρ[N-1,M])+(χ_n(ρ[N,M])+χ_n(ρ[N,M-1]))*(ρ[N,M]-ρ[N,M-1]))/4 + R_nn(ρ[N,M])*h^2/2, (D_n(ρ[N,M])+D_n(ρ[N-1,M]))/2+(χ_n(ρ[N,M])+χ_n(ρ[N-1,M]))*(ρ[N,M]-ρ[N-1,M])/4, (D_n(ρ[N,M])+D_n(ρ[N,M-1]))/2+(χ_n(ρ[N,M])+χ_n(ρ[N,M-1]))*(ρ[N,M]-ρ[N,M-1])/4];

    elseif BC == "periodic" # if periodic boundary condition
        # left and right edge
        for j = 2:M-1
            i_vec[count:count+4] .= ind(1,j);
            j_vec[count:count+4] = [ind(1,j), ind(N,j), ind(2,j), ind(1,j-1), ind(1,j+1)];
            A_vec[count:count+4] = 2/h^2*[-(3D_n(ρ[1,j])+D_n(ρ[N,j])+D_n(ρ[2,j])+D_n(ρ[1,j-1])/2+D_n(ρ[1,j+1])/2)/2 + ((χ_n(ρ[1,j])+χ_n(ρ[N,j]))*(ρ[1,j]-ρ[N,j])-(χ_n(ρ[1,j])+χ_n(ρ[2,j]))*(ρ[2,j]-ρ[1,j])+(χ_n(ρ[1,j])+χ_n(ρ[1,j-1]))*(ρ[1,j]-ρ[1,j-1])/2-(χ_n(ρ[1,j])+χ_n(ρ[1,j+1]))*(ρ[1,j+1]-ρ[1,j])/2)/4 + R_nn(ρ[1,j])*h^2/2, (D_n(ρ[1,j])+D_n(ρ[N,j]))/2+(χ_n(ρ[1,j])+χ_n(ρ[N,j]))*(ρ[1,j]-ρ[N,j])/4, (D_n(ρ[1,j])+D_n(ρ[2,j]))/2-(χ_n(ρ[1,j])+χ_n(ρ[2,j]))*(ρ[2,j]-ρ[1,j])/4, ((D_n(ρ[1,j])+D_n(ρ[1,j-1]))/2+(χ_n(ρ[1,j])+χ_n(ρ[1,j-1]))*(ρ[1,j]-ρ[1,j-1])/4)/2, ((D_n(ρ[1,j])+D_n(ρ[1,j+1]))/2-(χ_n(ρ[1,j])+χ_n(ρ[1,j+1]))*(ρ[1,j+1]-ρ[1,j])/4)/2];
            count = count + 5; 
        
            i_vec[count:count+4] .= ind(N,j);
            j_vec[count:count+4] = [ind(N,j), ind(N-1,j), ind(1,j), ind(N,j-1), ind(N,j+1)];
            A_vec[count:count+4] = 2/h^2*[-(3D_n(ρ[N,j])+D_n(ρ[N-1,j])+D_n(ρ[1,j])+D_n(ρ[N,j-1])/2+D_n(ρ[N,j+1])/2)/2 + ((χ_n(ρ[N,j])+χ_n(ρ[N-1,j]))*(ρ[N,j]-ρ[N-1,j])-(χ_n(ρ[N,j])+χ_n(ρ[1,j]))*(ρ[1,j]-ρ[N,j])+(χ_n(ρ[N,j])+χ_n(ρ[N,j-1]))*(ρ[N,j]-ρ[N,j-1])/2-(χ_n(ρ[N,j])+χ_n(ρ[N,j+1]))*(ρ[N,j+1]-ρ[N,j])/2)/4 + R_nn(ρ[N,j])*h^2/2, (D_n(ρ[N,j])+D_n(ρ[N-1,j]))/2+(χ_n(ρ[N,j])+χ_n(ρ[N-1,j]))*(ρ[N,j]-ρ[N-1,j])/4, (D_n(ρ[N,j])+D_n(ρ[1,j]))/2-(χ_n(ρ[N,j])+χ_n(ρ[1,j]))*(ρ[1,j]-ρ[N,j])/4, ((D_n(ρ[N,j])+D_n(ρ[N,j-1]))/2+(χ_n(ρ[N,j])+χ_n(ρ[N,j-1]))*(ρ[N,j]-ρ[N,j-1])/4)/2, ((D_n(ρ[N,j])+D_n(ρ[N,j+1]))/2-(χ_n(ρ[N,j])+χ_n(ρ[N,j+1]))*(ρ[N,j+1]-ρ[N,j])/4)/2];
            count = count + 5;
        end 

        # bottom and top edges
        for i = 2:N-1
            i_vec[count:count+4] .= ind(i,1);
            j_vec[count:count+4] = [ind(i,1), ind(i-1,1), ind(i+1,1), ind(i,M), ind(i,2)];
            A_vec[count:count+4] = 2/h^2*[-(3D_n(ρ[i,1])+D_n(ρ[i-1,1])/2+D_n(ρ[i+1,1])/2+D_n(ρ[i,M])+D_n(ρ[i,2]))/2 + ((χ_n(ρ[i,1])+χ_n(ρ[i-1,1]))*(ρ[i,1]-ρ[i-1,1])/2-(χ_n(ρ[i,1])+χ_n(ρ[i+1,1]))*(ρ[i+1,1]-ρ[i,1])/2+(χ_n(ρ[i,1])+χ_n(ρ[i,M]))*(ρ[i,1]-ρ[i,M])-(χ_n(ρ[i,1])+χ_n(ρ[i,2]))*(ρ[i,2]-ρ[i,1]))/4 + R_nn(ρ[i,1])*h^2/2, ((D_n(ρ[i,1])+D_n(ρ[i-1,1]))/2+(χ_n(ρ[i,1])+χ_n(ρ[i-1,1]))*(ρ[i,1]-ρ[i-1,1])/4)/2, ((D_n(ρ[i,1])+D_n(ρ[i+1,1]))/2-(χ_n(ρ[i,1])+χ_n(ρ[i+1,1]))*(ρ[i+1,1]-ρ[i,1])/4)/2, (D_n(ρ[i,1])+D_n(ρ[i,M]))/2+(χ_n(ρ[i,1])+χ_n(ρ[i,M]))*(ρ[i,1]-ρ[i,M])/4, (D_n(ρ[i,1])+D_n(ρ[i,2]))/2-(χ_n(ρ[i,1])+χ_n(ρ[i,2]))*(ρ[i,2]-ρ[i,1])/4];
            count = count + 5;
        
            i_vec[count:count+4] .= ind(i,M);
            j_vec[count:count+4] = [ind(i,M), ind(i-1,M), ind(i+1,M), ind(i,M-1), ind(i,1)];
            A_vec[count:count+4] = 2/h^2*[-(3D_n(ρ[i,M])+D_n(ρ[i-1,M])/2+D_n(ρ[i+1,M])/2+D_n(ρ[i,M-1])+D_n(ρ[i,1]))/2 + ((χ_n(ρ[i,M])+χ_n(ρ[i-1,M]))*(ρ[i,M]-ρ[i-1,M])/2-(χ_n(ρ[i,M])+χ_n(ρ[i+1,M]))*(ρ[i+1,M]-ρ[i,M])/2+(χ_n(ρ[i,M])+χ_n(ρ[i,M-1]))*(ρ[i,M]-ρ[i,M-1])-(χ_n(ρ[i,M])+χ_n(ρ[i,1]))*(ρ[i,1]-ρ[i,M]))/4 + R_nn(ρ[i,M])*h^2/2, ((D_n(ρ[i,M])+D_n(ρ[i-1,M]))/2+(χ_n(ρ[i,M])+χ_n(ρ[i-1,M]))*(ρ[i,M]-ρ[i-1,M])/4)/2, ((D_n(ρ[i,M])+D_n(ρ[i+1,M]))/2-(χ_n(ρ[i,M])+χ_n(ρ[i+1,M]))*(ρ[i+1,M]-ρ[i,M])/4)/2, (D_n(ρ[i,M])+D_n(ρ[i,M-1]))/2+(χ_n(ρ[i,M])+χ_n(ρ[i,M-1]))*(ρ[i,M]-ρ[i,M-1])/4, (D_n(ρ[i,M])+D_n(ρ[i,1]))/2-(χ_n(ρ[i,M])+χ_n(ρ[i,1]))*(ρ[i,1]-ρ[i,M])/4];
            count = count + 5;
        end

        # corners
        i_vec[count:count+4] .= ind(1,1);
        j_vec[count:count+4] = [ind(1,1), ind(N,1), ind(2,1), ind(1,M), ind(1,2)];
        A_vec[count:count+4] = 2/h^2*[-(4D_n(ρ[1,1])+D_n(ρ[N,1])+D_n(ρ[2,1])+D_n(ρ[1,M])+D_n(ρ[1,2]))/2 + ((χ_n(ρ[1,1])+χ_n(ρ[N,1]))*(ρ[1,1]-ρ[N,1])-(χ_n(ρ[1,1])+χ_n(ρ[2,1]))*(ρ[2,1]-ρ[1,1])+(χ_n(ρ[1,1])+χ_n(ρ[1,M]))*(ρ[1,1]-ρ[1,M])-(χ_n(ρ[1,1])+χ_n(ρ[1,2]))*(ρ[1,2]-ρ[1,1]))/4 + R_nn(ρ[1,1])*h^2/2, (D_n(ρ[1,1])+D_n(ρ[N,1]))/2+(χ_n(ρ[1,1])+χ_n(ρ[N,1]))*(ρ[1,1]-ρ[N,1])/4, (D_n(ρ[1,1])+D_n(ρ[2,1]))/2-(χ_n(ρ[1,1])+χ_n(ρ[2,1]))*(ρ[2,1]-ρ[1,1])/4, (D_n(ρ[1,1])+D_n(ρ[1,M]))/2+(χ_n(ρ[1,1])+χ_n(ρ[1,M]))*(ρ[1,1]-ρ[1,M])/4, (D_n(ρ[1,1])+D_n(ρ[1,2]))/2-(χ_n(ρ[1,1])+χ_n(ρ[1,2]))*(ρ[1,2]-ρ[1,1])/4];
        count = count+5; 
        
        i_vec[count:count+4] .= ind(N,1);
        j_vec[count:count+4] = [ind(N,1), ind(N-1,1), ind(1,1), ind(N,M), ind(N,2)];
        A_vec[count:count+4] = 2/h^2*[-(4D_n(ρ[N,1])+D_n(ρ[N-1,1])+D_n(ρ[1,1])+D_n(ρ[N,M])+D_n(ρ[N,2]))/2 + ((χ_n(ρ[N,1])+χ_n(ρ[N-1,1]))*(ρ[N,1]-ρ[N-1,1])-(χ_n(ρ[N,1])+χ_n(ρ[1,1]))*(ρ[1,1]-ρ[N,1])+(χ_n(ρ[N,1])+χ_n(ρ[N,M]))*(ρ[N,1]-ρ[N,M])-(χ_n(ρ[N,1])+χ_n(ρ[N,2]))*(ρ[N,2]-ρ[N,1]))/4 + R_nn(ρ[N,1])*h^2/2, (D_n(ρ[N,1])+D_n(ρ[N-1,1]))/2+(χ_n(ρ[N,1])+χ_n(ρ[N-1,1]))*(ρ[N,1]-ρ[N-1,1])/4, (D_n(ρ[N,1])+D_n(ρ[1,1]))/2-(χ_n(ρ[N,1])+χ_n(ρ[1,1]))*(ρ[1,1]-ρ[N,1])/4, (D_n(ρ[N,1])+D_n(ρ[N,M]))/2+(χ_n(ρ[N,1])+χ_n(ρ[N,M]))*(ρ[N,1]-ρ[N,M])/4, (D_n(ρ[N,1])+D_n(ρ[N,2]))/2-(χ_n(ρ[N,1])+χ_n(ρ[N,2]))*(ρ[N,2]-ρ[N,1])/4];
        count = count+5; 
        
        i_vec[count:count+4] .= ind(1,M);
        j_vec[count:count+4] = [ind(1,M), ind(N,M), ind(2,M), ind(1,M-1), ind(1,1)];
        A_vec[count:count+4] = 2/h^2*[-(4D_n(ρ[1,M])+D_n(ρ[N,M])+D_n(ρ[2,M])+D_n(ρ[1,M-1])+D_n(ρ[1,1]))/2 + ((χ_n(ρ[1,M])+χ_n(ρ[N,M]))*(ρ[1,M]-ρ[N,M])-(χ_n(ρ[1,M])+χ_n(ρ[2,M]))*(ρ[2,M]-ρ[1,M])+(χ_n(ρ[1,M])+χ_n(ρ[1,M-1]))*(ρ[1,M]-ρ[1,M-1])-(χ_n(ρ[1,M])+χ_n(ρ[1,1]))*(ρ[1,1]-ρ[1,M]))/4 + R_nn(ρ[1,M])*h^2/2, (D_n(ρ[1,M])+D_n(ρ[N,M]))/2+(χ_n(ρ[1,M])+χ_n(ρ[N,M]))*(ρ[1,M]-ρ[N,M])/4, (D_n(ρ[1,M])+D_n(ρ[2,M]))/2-(χ_n(ρ[1,M])+χ_n(ρ[2,M]))*(ρ[2,M]-ρ[1,M])/4, (D_n(ρ[1,M])+D_n(ρ[1,M-1]))/2+(χ_n(ρ[1,M])+χ_n(ρ[1,M-1]))*(ρ[1,M]-ρ[1,M-1])/4, (D_n(ρ[1,M])+D_n(ρ[1,1]))/2-(χ_n(ρ[1,M])+χ_n(ρ[1,1]))*(ρ[1,1]-ρ[1,M])/4];
        count = count+5; 
        
        i_vec[count:count+4] .= ind(N,M);
        j_vec[count:count+4] = [ind(N,M), ind(N-1,M), ind(1,M), ind(N,M-1), ind(N,1)];
        A_vec[count:count+4] = 2/h^2*[-(4D_n(ρ[N,M])+D_n(ρ[N-1,M])+D_n(ρ[1,M])+D_n(ρ[N,M-1])+D_n(ρ[N,1]))/2 + ((χ_n(ρ[N,M])+χ_n(ρ[N-1,M]))*(ρ[N,M]-ρ[N-1,M])-(χ_n(ρ[N,M])+χ_n(ρ[1,M]))*(ρ[1,M]-ρ[N,M])+(χ_n(ρ[N,M])+χ_n(ρ[N,M-1]))*(ρ[N,M]-ρ[N,M-1])-(χ_n(ρ[N,M])+χ_n(ρ[N,1]))*(ρ[N,1]-ρ[N,M]))/4 + R_nn(ρ[N,M])*h^2/2, (D_n(ρ[N,M])+D_n(ρ[N-1,M]))/2+(χ_n(ρ[N,M])+χ_n(ρ[N-1,M]))*(ρ[N,M]-ρ[N-1,M])/4, (D_n(ρ[N,M])+D_n(ρ[1,M]))/2-(χ_n(ρ[N,M])+χ_n(ρ[1,M]))*(ρ[1,M]-ρ[N,M])/4, (D_n(ρ[N,M])+D_n(ρ[N,M-1]))/2+(χ_n(ρ[N,M])+χ_n(ρ[N,M-1]))*(ρ[N,M]-ρ[N,M-1])/4, (D_n(ρ[N,M])+D_n(ρ[N,1]))/2-(χ_n(ρ[N,M])+χ_n(ρ[N,1]))*(ρ[N,1]-ρ[N,M])/4];
    end

    return A = sparse(i_vec,j_vec,A_vec); # generate A matrix containing coefficients in ODE system
end

function genA_a(N, M, h, D_a, χ_aρ, R_aa, ρ, c, BC)
    # Generates the A matrix that holds the coefficients in the ODE system (nodes from FVM) for u_a
    #  original PDE is:  ∂u_a/∂t = ∇•[D_a(ρ,c)∇u_a - χ_aρ(ρ,c)∇ρu_a] + R_aa u_a + R_an(ρ)u_n
    # Inputs:
    #  N: width of discretised grid (number of nodes)
    #  M: height of discretised grid (number of nodes)
    #  h: distance between nodes in discretised grid
    #  D_a: scaffold-dependent diffusivity for u_a
    #  χ_aρ: scaffold- and cytokine-dependent haptotactic sensitivity in the direction of ∇ρ for u_a
    #  R_aa: constant reaction term driven by activated cells for u_a
    #  ρ: matrix containing scaffold density in the domain
    #  c: cytokine density in the domain (condensed into 1D array) at the previous timestep
    #  BC: boundary condition: "noflux" or "periodic"
    # Output:
    #  A: sparse A matrix containing all coefficients

    ind(i,j) = i + N*(j-1); # equation to determine index in A for the (i,j)th node

    if BC == "noflux"
        nonzero_entries = (N-2)*(M-2)*5+2*(N+M-4)*4+4*3; # number of equations * number of coefficients  for  interior+edges+corners
    elseif BC == "periodic"
        nonzero_entries = (N-2)*(M-2)*5+2*(N+M-4)*5+4*5;
    end
    
    i_vec = zeros(nonzero_entries,1)[:]; # vector containing all i indices of non-zero entries into A
    j_vec = zeros(nonzero_entries,1)[:]; # vector containing all j indices of non-zero entries into A
    A_vec = zeros(nonzero_entries,1)[:]; # vector containing all non-zero entries into A
    count = 1; # counter for elements added to vectors

    # compute function at each node
    D_a_ij = zeros(N,M);
    χ_aρ_ij = zeros(N,M);
    for i = 1:N 
        for j = 1:M 
            D_a_ij[i,j] = D_a(ρ[i,j]);
            χ_aρ_ij[i,j] = χ_aρ(ρ[i,j]);
        end
    end

    # interior nodes
    for i = 2:N-1
        for j = 2:M-1
            i_vec[count:count+4] .= ind(i,j);
            j_vec[count:count+4] = [ind(i,j), ind(i-1,j), ind(i+1,j), ind(i,j-1), ind(i,j+1)];
            A_vec[count:count+4] = 1/h^2*[-(4D_a_ij[i,j]+D_a_ij[i-1,j]+D_a_ij[i+1,j]+D_a_ij[i,j-1]+D_a_ij[i,j+1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4 + R_aa*h^2, (D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4, (D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4, (D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4, (D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4];
            count = count + 5; 
        end
    end

    if BC == "noflux" # if no-flux boundary condition
        # left and right edge
        for j = 2:M-1
            i = 1; 
            i_vec[count:count+3] .= ind(1,j);
            j_vec[count:count+3] = [ind(1,j), ind(2,j), ind(1,j-1), ind(1,j+1)];
            A_vec[count:count+3] = 2/h^2*[-(2D_a_ij[i,j]+D_a_ij[i+1,j]+D_a_ij[i,j-1]/2+D_a_ij[i,j+1]/2)/2 + (-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j])/2)/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4, ((D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4)/2, ((D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4)/2];
            count = count + 4; 
            
            i = N;
            i_vec[count:count+3] .= ind(N,j);
            j_vec[count:count+3] = [ind(N,j), ind(N-1,j), ind(N,j-1), ind(N,j+1)];
            A_vec[count:count+3] = 2/h^2*[-(2D_a_ij[i,j]+D_a_ij[i-1,j]+D_a_ij[i,j-1]/2+D_a_ij[i,j+1]/2)/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j])/2)/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4, ((D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4)/2, ((D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4)/2];
            count = count + 4;
        end 

        # bottom and top edges
        for i = 2:N-1
            j = 1;
            i_vec[count:count+3] .= ind(i,1);
            j_vec[count:count+3] = [ind(i,1), ind(i-1,1), ind(i+1,1), ind(i,2)];
            A_vec[count:count+3] = 2/h^2*[-(2D_a_ij[i,j]+D_a_ij[i-1,j]/2+D_a_ij[i+1,j]/2+D_a_ij[i,j+1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4 + R_aa*h^2/2, ((D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4)/2, ((D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4)/2, (D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4];
            count = count + 4;
        
            j = M;
            i_vec[count:count+3] .= ind(i,M);
            j_vec[count:count+3] = [ind(i,M), ind(i-1,M), ind(i+1,M), ind(i,M-1)];
            A_vec[count:count+3] = 2/h^2*[-(2D_a_ij[i,j]+D_a_ij[i-1,j]/2+D_a_ij[i+1,j]/2+D_a_ij[i,j-1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])/2+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4 + R_aa*h^2/2, ((D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4)/2, ((D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4)/2, (D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4];
            count = count + 4;
        end

        # corners
        i = 1; j = 1;
        i_vec[count:count+2] .= ind(1,1);
        j_vec[count:count+2] = [ind(1,1), ind(2,1), ind(1,2)];
        A_vec[count:count+2] = 2/h^2*[-(2D_a_ij[i,j]+D_a_ij[i+1,j]+D_a_ij[i,j+1])/2 + (-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4, (D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4];
        count = count+3; 
        
        i = N; j = 1;
        i_vec[count:count+2] .= ind(N,1);
        j_vec[count:count+2] = [ind(N,1), ind(N-1,1), ind(N,2)];
        A_vec[count:count+2] = 2/h^2*[-(2D_a_ij[i,j]+D_a_ij[i-1,j]+D_a_ij[i,j+1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4, (D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4];
        count = count+3; 
        
        i = 1; j = M;
        i_vec[count:count+2] .= ind(1,M);
        j_vec[count:count+2] = [ind(1,M), ind(2,M), ind(1,M-1)];
        A_vec[count:count+2] = 2/h^2*[-(2D_a_ij[i,j]+D_a_ij[i+1,j]+D_a_ij[i,j-1])/2 + (-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4, (D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4];
        count = count+3; 
        
        i = N; j = M;
        i_vec[count:count+2] .= ind(N,M);
        j_vec[count:count+2] = [ind(N,M), ind(N-1,M), ind(N,M-1)];
        A_vec[count:count+2] = 2/h^2*[-(2D_a_ij[i,j]+D_a_ij[i-1,j]+D_a_ij[i,j-1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4, (D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4];

    elseif BC == "periodic" # if periodic boundary condition
        # left and right edge
        for j = 2:M-1
            i = 1;
            i_vec[count:count+4] .= ind(1,j);
            j_vec[count:count+4] = [ind(1,j), ind(N,j), ind(2,j), ind(1,j-1), ind(1,j+1)];
            A_vec[count:count+4] = 2/h^2*[-(3D_a_ij[i,j]+D_a_ij[N,j]+D_a_ij[i+1,j]+D_a_ij[i,j-1]/2+D_a_ij[i,j+1]/2)/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[N,j])*(ρ[i,j]-ρ[N,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j])/2)/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[N,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[N,j])*(ρ[i,j]-ρ[N,j]))/4, (D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4, ((D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4)/2, ((D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4)/2];
            count = count + 5; 
        
            i = N;
            i_vec[count:count+4] .= ind(N,j);
            j_vec[count:count+4] = [ind(N,j), ind(N-1,j), ind(1,j), ind(N,j-1), ind(N,j+1)];
            A_vec[count:count+4] = 2/h^2*[-(3D_a_ij[i,j]+D_a_ij[i-1,j]+D_a_ij[1,j]+D_a_ij[i,j-1]/2+D_a_ij[i,j+1]/2)/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[1,j])*(ρ[1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j])/2)/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4, (D_a_ij[i,j]+D_a_ij[1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[1,j])*(ρ[1,j]-ρ[i,j]))/4, ((D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4)/2, ((D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4)/2];
            count = count + 5;
        end 

        # bottom and top edges
        for i = 2:N-1
            j = 1;
            i_vec[count:count+4] .= ind(i,1);
            j_vec[count:count+4] = [ind(i,1), ind(i-1,1), ind(i+1,1), ind(i,M), ind(i,2)];
            A_vec[count:count+4] = 2/h^2*[-(3D_a_ij[i,j]+D_a_ij[i-1,j]/2+D_a_ij[i+1,j]/2+D_a_ij[i,M]+D_a_ij[i,j+1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])/2+(χ_aρ_ij[i,j]+χ_aρ_ij[i,M])*(ρ[i,j]-ρ[i,M])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4 + R_aa*h^2/2, ((D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4)/2, ((D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4)/2, (D_a_ij[i,j]+D_a_ij[i,M])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,M])*(ρ[i,j]-ρ[i,M]))/4, (D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4];
            count = count + 5;
        
            j = M;
            i_vec[count:count+4] .= ind(i,M);
            j_vec[count:count+4] = [ind(i,M), ind(i-1,M), ind(i+1,M), ind(i,M-1), ind(i,1)];
            A_vec[count:count+4] = 2/h^2*[-(3D_a_ij[i,j]+D_a_ij[i-1,j]/2+D_a_ij[i+1,j]/2+D_a_ij[i,j-1]+D_a_ij[i,1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])/2-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])/2+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,1])*(ρ[i,1]-ρ[i,j]))/4 + R_aa*h^2/2, ((D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4)/2, ((D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4)/2, (D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4, (D_a_ij[i,j]+D_a_ij[i,1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,1])*(ρ[i,1]-ρ[i,j]))/4];
            count = count + 5;
        end

        # corners
        i = 1; j = 1;
        i_vec[count:count+4] .= ind(1,1);
        j_vec[count:count+4] = [ind(1,1), ind(N,1), ind(2,1), ind(1,M), ind(1,2)];
        A_vec[count:count+4] = 2/h^2*[-(4D_a_ij[i,j]+D_a_ij[N,j]+D_a_ij[i+1,j]+D_a_ij[i,M]+D_a_ij[i,j+1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[N,j])*(ρ[i,j]-ρ[N,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,M])*(ρ[i,j]-ρ[i,M])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[N,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[N,j])*(ρ[i,j]-ρ[N,j]))/4, (D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4, (D_a_ij[i,j]+D_a_ij[i,M])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,M])*(ρ[i,j]-ρ[i,M]))/4, (D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4];
        count = count+5; 
        
        i = N; j = 1;
        i_vec[count:count+4] .= ind(N,1);
        j_vec[count:count+4] = [ind(N,1), ind(N-1,1), ind(1,1), ind(N,M), ind(N,2)];
        A_vec[count:count+4] = 2/h^2*[-(4D_a_ij[i,j]+D_a_ij[i-1,j]+D_a_ij[1,j]+D_a_ij[i,M]+D_a_ij[i,j+1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[1,j])*(ρ[1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,M])*(ρ[i,j]-ρ[i,M])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4, (D_a_ij[i,j]+D_a_ij[1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[1,j])*(ρ[1,j]-ρ[i,j]))/4, (D_a_ij[i,j]+D_a_ij[i,M])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,M])*(ρ[i,j]-ρ[i,M]))/4, (D_a_ij[i,j]+D_a_ij[i,j+1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,j+1])*(ρ[i,j+1]-ρ[i,j]))/4];
        count = count+5; 
        
        i = 1; j = M;
        i_vec[count:count+4] .= ind(1,M);
        j_vec[count:count+4] = [ind(1,M), ind(N,M), ind(2,M), ind(1,M-1), ind(1,1)];
        A_vec[count:count+4] = 2/h^2*[-(4D_a_ij[i,j]+D_a_ij[N,j]+D_a_ij[i+1,j]+D_a_ij[i,j-1]+D_a_ij[i,1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[N,j])*(ρ[i,j]-ρ[N,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,1])*(ρ[i,1]-ρ[i,j]))/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[N,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[N,j])*(ρ[i,j]-ρ[N,j]))/4, (D_a_ij[i,j]+D_a_ij[i+1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i+1,j])*(ρ[i+1,j]-ρ[i,j]))/4, (D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4, (D_a_ij[i,j]+D_a_ij[i,1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,1])*(ρ[i,1]-ρ[i,j]))/4];
        count = count+5; 
        
        i = N; j = M;
        i_vec[count:count+4] .= ind(N,M);
        j_vec[count:count+4] = [ind(N,M), ind(N-1,M), ind(1,M), ind(N,M-1), ind(N,1)];
        A_vec[count:count+4] = 2/h^2*[-(4D_a_ij[i,j]+D_a_ij[i-1,j]+D_a_ij[1,j]+D_a_ij[i,j-1]+D_a_ij[i,1])/2 + ((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j])-(χ_aρ_ij[i,j]+χ_aρ_ij[1,j])*(ρ[1,j]-ρ[i,j])+(χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1])-(χ_aρ_ij[i,j]+χ_aρ_ij[i,1])*(ρ[i,1]-ρ[i,j]))/4 + R_aa*h^2/2, (D_a_ij[i,j]+D_a_ij[i-1,j])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i-1,j])*(ρ[i,j]-ρ[i-1,j]))/4, (D_a_ij[i,j]+D_a_ij[1,j])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[1,j])*(ρ[1,j]-ρ[i,j]))/4, (D_a_ij[i,j]+D_a_ij[i,j-1])/2+((χ_aρ_ij[i,j]+χ_aρ_ij[i,j-1])*(ρ[i,j]-ρ[i,j-1]))/4, (D_a_ij[i,j]+D_a_ij[i,1])/2+(-(χ_aρ_ij[i,j]+χ_aρ_ij[i,1])*(ρ[i,1]-ρ[i,j]))/4];
    end

    return A = sparse(i_vec,j_vec,A_vec); # generate A matrix containing coefficients in ODE system
end

function genA_I_IMEX(N, M, h, D_c, BC)
    # Generates the A matrix that holds the coefficients in the ODE system (nodes from FVM) for c using an IMEX method
    #  original PDE is:  ∂c/∂t = ∇•(D_c∇c) - λ(u_a)c + S(ρ,t)
    # Inputs:
    #  N: width of discretised grid (number of nodes)
    #  M: height of discretised grid (number of nodes)
    #  h: distance between nodes in discretised grid
    #  D_c: constant diffusivity for c
    #  BC: boundary condition: "noflux" or "periodic"
    # Output:
    #  A: sparse A matrix containing all coefficients

    ind(i,j) = i + N*(j-1); # equation to determine index in A for the (i,j)th node
    if BC == "noflux"
        nonzero_entries = (N-2)*(M-2)*5+2*(N+M-4)*4+4*3; # number of equations * number of coefficients  for  interior+edges+corners
    elseif BC == "periodic"
        nonzero_entries = (N-2)*(M-2)*5+2*(N+M-4)*5+4*5;
    end
    i_vec = zeros(nonzero_entries,1)[:]; # vector containing all i indices of non-zero entries into A
    j_vec = zeros(nonzero_entries,1)[:]; # vector containing all j indices of non-zero entries into A
    A_vec = zeros(nonzero_entries,1)[:]; # vector containing all non-zero entries into A
    count = 1; # counter for elements added to vectors

    # interior nodes
    for i = 2:N-1
        for j = 2:M-1
            i_vec[count:count+4] .= ind(i,j);
            j_vec[count:count+4] = [ind(i,j), ind(i-1,j), ind(i+1,j), ind(i,j-1), ind(i,j+1)];
            A_vec[count:count+4] = [-1/h^2*4D_c, 1/h^2*D_c, 1/h^2*D_c, 1/h^2*D_c, 1/h^2*D_c];
            count = count + 5; 
        end
    end

    if BC == "noflux" # if no-flux boundary condition
        # left and right edge
        for j = 2:M-1
            i_vec[count:count+3] .= ind(1,j);
            j_vec[count:count+3] = [ind(1,j), ind(2,j), ind(1,j-1), ind(1,j+1)];
            A_vec[count:count+3] = [-2/h^2*2D_c, 2/h^2*D_c, 1/h^2*D_c, 1/h^2*D_c];
            count = count + 4; 
        
            i_vec[count:count+3] .= ind(N,j);
            j_vec[count:count+3] = [ind(N,j), ind(N-1,j), ind(N,j-1), ind(N,j+1)];
            A_vec[count:count+3] = [-2/h^2*2D_c, 2/h^2*D_c, 1/h^2*D_c, 1/h^2*D_c];
            count = count + 4;
        end 

        # bottom and top edges
        for i = 2:N-1
            i_vec[count:count+3] .= ind(i,1);
            j_vec[count:count+3] = [ind(i,1), ind(i-1,1), ind(i+1,1), ind(i,2)];
            A_vec[count:count+3] = [-2/h^2*2D_c, 1/h^2*D_c, 1/h^2*D_c, 2/h^2*D_c];
            count = count + 4;
        
            i_vec[count:count+3] .= ind(i,M);
            j_vec[count:count+3] = [ind(i,M), ind(i-1,M), ind(i+1,M), ind(i,M-1)];
            A_vec[count:count+3] = [-2/h^2*2D_c, 1/h^2*D_c, 1/h^2*D_c, 2/h^2*D_c];
            count = count + 4;
        end

        # corners
        i_vec[count:count+2] .= ind(1,1);
        j_vec[count:count+2] = [ind(1,1), ind(2,1), ind(1,2)];
        A_vec[count:count+2] = [-4/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
        count = count+3; 
        
        i_vec[count:count+2] .= ind(N,1);
        j_vec[count:count+2] = [ind(N,1), ind(N-1,1), ind(N,2)];
        A_vec[count:count+2] = [-4/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
        count = count+3; 
        
        i_vec[count:count+2] .= ind(1,M);
        j_vec[count:count+2] = [ind(1,M), ind(2,M), ind(1,M-1)];
        A_vec[count:count+2] = [-4/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
        count = count+3; 
        
        i_vec[count:count+2] .= ind(N,M);
        j_vec[count:count+2] = [ind(N,M), ind(N-1,M), ind(N,M-1)];
        A_vec[count:count+2] = [-4/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];


    elseif BC == "periodic" # if periodic boundary condition
        # left and right edge
        for j = 2:M-1
            i_vec[count:count+4] .= ind(1,j);
            j_vec[count:count+4] = [ind(1,j), ind(N,j), ind(2,j), ind(1,j-1), ind(1,j+1)];
            A_vec[count:count+4] = [-2/h^2*3D_c, 2/h^2*D_c, 2/h^2*D_c, 1/h^2*D_c, 1/h^2*D_c];
            count = count + 5; 
        
            i_vec[count:count+4] .= ind(N,j);
            j_vec[count:count+4] = [ind(N,j), ind(N-1,j), ind(1,j), ind(N,j-1), ind(N,j+1)];
            A_vec[count:count+4] = [-2/h^2*3D_c, 2/h^2*D_c, 2/h^2*D_c, 1/h^2*D_c, 1/h^2*D_c];
            count = count + 5;
        end 

        # bottom and top edges
        for i = 2:N-1 
            i_vec[count:count+4] .= ind(i,1);
            j_vec[count:count+4] = [ind(i,1), ind(i-1,1), ind(i+1,1), ind(i,M), ind(i,2)];
            A_vec[count:count+4] = [-2/h^2*3D_c, 1/h^2*D_c, 1/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
            count = count + 5;
        
            i_vec[count:count+4] .= ind(i,M);
            j_vec[count:count+4] = [ind(i,M), ind(i-1,M), ind(i+1,M), ind(i,M-1), ind(i,1)];
            A_vec[count:count+4] = [-2/h^2*3D_c, 1/h^2*D_c, 1/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
            count = count + 5;
        end

        # corners
        i_vec[count:count+4] .= ind(1,1);
        j_vec[count:count+4] = [ind(1,1), ind(N,1), ind(2,1), ind(1,M), ind(1,2)];
        A_vec[count:count+4] = [-4/h^2*2D_c, 2/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
        count = count+5; 
        
        i_vec[count:count+4] .= ind(N,1);
        j_vec[count:count+4] = [ind(N,1), ind(N-1,1), ind(1,1), ind(N,M), ind(N,2)];
        A_vec[count:count+4] = [-4/h^2*2D_c, 2/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
        count = count+5; 
        
        i_vec[count:count+4] .= ind(1,M);
        j_vec[count:count+4] = [ind(1,M), ind(N,M), ind(2,M), ind(1,M-1), ind(1,1)];
        A_vec[count:count+4] = [-4/h^2*2D_c, 2/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
        count = count+5; 
        
        i_vec[count:count+4] .= ind(N,M);
        j_vec[count:count+4] = [ind(N,M), ind(N-1,M), ind(1,M), ind(N,M-1), ind(N,1)];
        A_vec[count:count+4] = [-4/h^2*2D_c, 2/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c, 2/h^2*D_c];
    end

    return A = sparse(i_vec,j_vec,A_vec); # generate A matrix containing coefficients in ODE system
end

function linear_solver(A, b, u0, lsolve, atol, rtol, maxiters, pretype, precond)
    # Compute the solution at the next timestep by solving the linear
    #  system Au=b where u0=u(n-1)
    # Inputs
    #  A: A matrix in linear system Ax=b
    #  b: b vector in linear system
    #  u0: initial guess for the solution u
    #  lsolve: linear system solution method 'backslash', 'FOM' or 'GMRES'
    #  atol: FOM/GMRES absolute error tolerance
    #  rtol: FOM/GMRES relative error tolerance
    #  maxiters: FOM maximum number of iterations
    #  pretype: FOM/GMRES preconditioner type 'none', 'left' or 'right'
    #  precond: FOM/GMRES preconditioner option 'Jacobi' or 'Gauss-Seidel'
    # Outputs
    #  u: next step in solution u(i,t) at node i and time t

    # handle additional FOM or GMRES parameters
    if lsolve == "backslash"
        return x = A\b;
    end
    
    # select preconditioner if not using none
    if pretype != "none"
        if precond == "Jacobi"
            M = diagm(diag(A));
        elseif precond == "Gauss-Seidel"
            M = tril(A);
        end
    end


    # initialise
    N  = size(A,1);
    H = zeros(maxiters+1,maxiters);
    V = zeros(N,maxiters+1);
    rnorm = zeros(maxiters,1);
    y = zeros(maxiters,1); 
    
    # left or right-preconditioned FOM or GMRES:
    r = b - A*u0; 
    if pretype == "right"
        beta = norm(r,2); V[:,1] = r/beta; 
    elseif pretype == "left"
        rtilde = M\r; beta = norm(rtilde,2); V[:,1] = rtilde/beta;
    end
    converged = false;
    m_con = maxiters;
    rnorm0 = norm(r,2);

    for m = 1:maxiters
        # Arnoldi
        if pretype == "left"
            V[:,m+1] = M\(A*V[:,m]);
        elseif pretype == "right"
            V[:,m+1] = A*(M\V[:,m]);
        end
        for j = 1:m
            H[j,m] = V[:,j]'*V[:,m+1];
            V[:,m+1] = V[:,m+1] - H[j,m]*V[:,j];
        end
        H[m+1,m] = norm(V[:,m+1],2);
        
        y = zeros(m,1); 
        # check for breakdown
        if H[m+1,m] < 1e-14
            y = H[1:m,1:m] \ ([beta; zeros(m-1,1)]);
            converged = true;
            break
        else
            V[:,m+1] = V[:,m+1]/H[m+1,m];
        end
        
        # solve small m dimensional linear system for y
        if lsolve == "FOM"
            y = H[1:m,1:m] \ ([beta; zeros(m-1,1)]);
            if pretype == "right"
                rnorm = H[m+1,m]*abs(y[m]);
            elseif pretype == "left"
                rnorm = norm(M*H[m+1,m]*V[:,m+1]*y[m],2);
            end
            rnorm = H[m+1,m]*abs(y[m]);
        else
            y = H[1:m+1,1:m] \ ([beta; zeros(m,1)]);
            if pretype == "right"
                rnorm = norm([beta; zeros(m,1)] - H[1:m+1,1:m]*y,2);
            elseif pretype == "left"
                rnorm = norm(M[m+1,m+1]*([beta; zeros(m,1)] - H[1:m+1,1:m]*y),2);
            end
        end
        
        if rnorm < rtol*rnorm0 + atol
            converged = true;
            m_con = m;
            break
        end
        
    end
    
    # compute approximate solution
    if converged
        if pretype == "left"
            return x = u0 + V[:,1:m_con]*y;
        elseif pretype == "right"
            return x = u0 + M\(V[:,1:m_con]*y);
        end
    else
        return x = [];
    end

end

function IL2_ODE(dy,y,p,t)
    λ, S, N_a, ρ, V = p; # λ is a function of a and I, and S is a function of ρ and t
    dy[1] = -λ(N_a/V,y[1]/V)*V + S(ρ/V,t)*V; # dN_I/dt
end

function IL2_ODE_conc(dy,y,p,t)
    λ, S, a, ρ = p; # λ is a function of a and I, and S is a function of ρ and t
    dy[1] = -λ(a,y[1]) + S(ρ,t); # dN_I/dt
end

function ODE_system(dy,y,p,t)
    R_nn, R_na, R_aa, R_an, λ, S, ρ, V = p; # R_nn and R_an are functions of ρ, R_na and R_aa are functions of I, λ is a function of a and I, and S is a function of ρ and t
    dy[1] = R_nn(ρ/V)*y[1] + R_na(y[3]/V)*y[2]; # dN_n/dt
    dy[2] = R_aa(y[3]/V)*y[2] + R_an(ρ/V)*y[1]; # dN_a/dt
    dy[3] = -λ(y[2]/V,y[3]/V)*V + S(ρ/V,t)*V; # dN_I/dt
end

function simulate_ABM()
    # Simatules ABM ntotal times to get an average

    #use_ODE_for_IL2, cyt_IC, I0_mass, IC_I, cell_IC, halfN0, N0, x_inj, y_inj, inj_ind, cell, ind, start_with_active, ε, τ, N, M, wid, hei, T, Ts_ABM, ntotal, r_d, r_a, r_p, ρ = p;

    start_time = time(); # define start time for averaged ABM solving

    Ts_ABM = T÷τ_ABM; # number of timesteps for IL-2 solving

    n = zeros(N, M, Int(floor(T))+1); # initialise T cell density array for naive T cells, n(x,y,t)
    a = zeros(N, M, Int(floor(T))+1); # initialise T cell density array for activated T cells, a(x,y,t)

    if use_ODE_for_IL2 # if the ODE is being used for IL-2
        # setup ODE solver
        if cyt_IC == "zero" # if starting with some initial cytokine
            I_i = 0
        else
            I_i = I0_mass;
        end
        I = zeros(Ts_ABM+1); I[1] = I_i; # cytokine density for current simulation
        I_ABM = zeros(Ts_ABM+1); # stores the average cytokine mass over time
    else
        I = zeros(N*M, Ts_ABM+1); I[:,1] = IC_I[:]; # cytokine density as a 2-dimensional array (for use in PDE solving)
        I_ABM = zeros(N*M, Ts_ABM+1); # stores the average cytokine solution
    end


    n_tot_sims = zeros(T+1, ntotal); # initialise arrays to hold time-series total populations for all simulations
    a_tot_sims = zeros(T+1, ntotal);
    I_tot_sims = zeros(Ts_ABM+1, ntotal); 
    if start_with_active
        n_tot_sims[1,:] .= halfN0; # both cells start with half the total number of cells
        a_tot_sims[1,:] .= halfN0; 
    else
        n_tot_sims[1,:] .= N0; # all starting cells are naive
    end


    # pre-compute probability of death, activation and waiting for current scaffold 
    P_r_d = r_d*τ;
    P_r_a_ij = r_a.(ρ)*τ;
    w_ij = w.(ρ);

    # average over ntotal simulations
    @showprogress 1 "Averaging simulations... " for sim = 1:ntotal # for each independent simulation     -    can use Threads.@threads 
        n_sim = zeros(N,M,Int(floor(T))+1); # masses of naive and activated T cells in the current simulation
        a_sim = zeros(N,M,Int(floor(T))+1); 
        if cell_IC == "pinpoint"
            if start_with_active
                cells = [cell(x_inj, y_inj, false) for i in 1:halfN0]; # half of total cells are naive and have initial coordinates at the same location
                append!(cells, [cell(x_inj, y_inj, true) for i in 1:halfN0]); # add other half which are activated
                n_sim[inj_ind[1],inj_ind[2],1] = halfN0; # add to the arrays containing the current simulation
                a_sim[inj_ind[1],inj_ind[2],1] = halfN0;
            else
                cells = [cell(x_inj, y_inj, false) for i in 1:N0]; # each cell is naive and has initial coordinates at the same location
                n_sim[inj_ind[1],inj_ind[2],1] = N0;
            end
        elseif cell_IC == "uniform"
            if start_with_active
                cells = [cell(rand()*wid, rand()*hei, false) for i in 1:halfN0]; # half of total cells are naive and have random initial coordinates (sampled from uniform distribution)
                append!(cells, [cell(rand()*wid, rand()*hei, true) for i in 1:halfN0]); # add other half which are activated
            else
                cells = [cell(rand()*wid, rand()*hei, false) for i in 1:N0]; # each cell is naive and has random initial coordinates (sampled from uniform distribution)
            end
            for n = 1:N0 # for each initial T cell
                x = cells[n].x; y = cells[n].y; # x and y positions of current cell
                ind_n = ind(x, y); # current cell's coords in indices (as [i,j] array)
                if cells[n].activated
                    a_sim[ind_n[1],ind_n[2],1] += 1; # contribute to array containing only the current simulation
                else
                    n_sim[ind_n[1],ind_n[2],1] += 1; # contribute to array containing only the current simulation
                end
            end
        end
        if cyt_IC == "zero" # if starting with some initial cytokine
            I_i = 0
        else
            I_i = I0_mass;
        end

        for i = 1:Int(floor(T)) # for each timestep
            #display("Simulation $(sim) out of $(ntotal). Timestep $(i) out of $(T)  =>  $(round((i+(sim-1)*T)/(ntotal*T)*100,digits=2))% complete.") 
            dead_cells = []; # list of indices of cells that die in the current timestep
            if use_ODE_for_IL2
                I_ij = I[i]/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # current IL-2 concentration at each lattice site (spatially constant for ODE)
                P_r_p_n = r_p(I_ij)*τ; # current corresponding probability of proliferation (spatially constant)
            end
            for n in eachindex(cells) # for each current T cell
                cell_alive = true; # initialise each cell to be alive
                # properties at current cell position
                x = cells[n].x; y = cells[n].y; # x and y positions of current cell
                a_n = cells[n].activated; # activation state of current cell

                ind_n = ind(x, y); # current cell's coords in indices (as [i,j] array)
                i_n = ind_n[1]; j_n = ind_n[2];
                
                # cytokine (IL-2) density at current cell's position
                if use_ODE_for_IL2 # if using ODE for IL-2
                    I_n = I[i÷τ_ABM+1*(mod(i,τ_ABM)!=0)]/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # use the current (or most recently calcualted) value and average over space
                else # using PDE for IL-2
                    I_n = I[i_n + N*(j_n-1), i÷τ_ABM+1*(mod(i,τ_ABM)!=0)]; 
                end

                if !use_ODE_for_IL2
                    P_r_p_n = r_p(I_n)*τ; # current cell's probability of proliferation
                end
                
                if s # if multiple actions in one timestep are prohibited
                    r = rand(); # sample random number ∈ [0,1]
                    if a_n # if current cell is activated, it may only die, prolfierate, wait or move
                        if r < P_r_p_n # cell proliferates
                            append!(cells, [cell(x, y, false)]); # create new naive T cell at current cell's location
                            n_sim[i_n, j_n, i+1] += 1; # contribute naive cell at the current location in current timestep
                            n_tot_sims[i+1,sim] += 1; # contribute naive cell to total for this simulation
                        elseif P_r_p_n <= r < P_r_p_n + P_r_d # cell dies
                            append!(dead_cells, n); # record current cell as a dead cell (to be removed after the current loop)
                            cell_alive = false;
                        elseif r >= P_r_p_n + P_r_d + w_ij[i_n,j_n] # cell is not waiting at the current position
                            # movement probabilities
                            p_d_curr = p_d[i_n, j_n]; # current probability of moving down (not considering boundary condition)
                            p_u_curr = p_u[i_n, j_n]; # up
                            p_l_curr = p_l[i_n, j_n]; # left
                            p_r_curr = p_r[i_n, j_n]; # right
                            
                            x, y = move_cell(x, y, wid, hei, ε, BC, p_d_curr, p_u_curr, p_l_curr, p_r_curr); # move the cell based on movement probabilities and other factors
                        end # otherwise cell waits
                    else # else it is naive, and may only die, activate, wait or move
                        P_r_a_n = P_r_a_ij[i_n,j_n]; # current activation probability
                        if r < P_r_a_n # cell activates
                            cells[n].activated = true;
                        elseif P_r_a_n <= r < P_r_a_n + P_r_d # cell dies
                            append!(dead_cells, n); # record current cell as a dead cell (to be removed after the current loop)
                            cell_alive = false;
                        elseif r >= P_r_a_n + P_r_d + w_ij[i_n,j_n] # cell is not waiting at the current position
                            # movement probabilities
                            p_d_curr = p_d[i_n, j_n]; # current probability of moving down (not considering boundary condition)
                            p_u_curr = p_u[i_n, j_n]; # up
                            p_l_curr = p_l[i_n, j_n]; # left
                            p_r_curr = p_r[i_n, j_n]; # right
                            
                            x, y = move_cell(x, y, wid, hei, ε, BC, p_d_curr, p_u_curr, p_l_curr, p_r_curr); # move the cell based on movement probabilities and other factors
                        end # otherwise cell waits
                    end
                else # multiple actions in one timestep are allowed
                    P_r_a_n = P_r_a_ij[i_n,j_n]; # current activation probability
                    r = rand(); # sample random number ∈ [0,1]
                    if r < P_r_a_n && !a_n # cell activates
                        cells[n].activated = true;
                        a_n = true;
                    end
                    r = rand();
                    if r < P_r_p_n && a_n # cell reproduces (re-calculate proliferation probability in case of activation in current step)
                        append!(cells, [cell(x, y, false)]); # create new naive T cell at current cell's location
                        n_sim[i_n, j_n, i+1] += 1; # contribute naive cell at the current location in current timestep
                        n_tot_sims[i+1,sim] += 1; # contribute naive cell to total for this simulation
                    end
                    r = rand();
                    if r < P_r_d # cell dies
                        append!(dead_cells, n); # record current cell as a dead cell (to be removed after the current loop)
                        cell_alive = false;
                    end
                    # move the cell or choose to wait at current position
                    if cell_alive && r < 1-w_ij[i_n,j_n] # cell moves
                        # movement probabilities
                        p_d_curr = p_d[i_n, j_n]; # current probability of moving down (not considering boundary condition)
                        p_u_curr = p_u[i_n, j_n]; # up
                        p_l_curr = p_l[i_n, j_n]; # left
                        p_r_curr = p_r[i_n, j_n]; # right
                        
                        x, y = move_cell(x, y, wid, hei, ε, BC, p_d_curr, p_u_curr, p_l_curr, p_r_curr); # move the cell based on movement probabilities and other factors
                    end
                end

                if cell_alive # if the cell is still alive
                    cells[n].x = x; cells[n].y = y; # record new position of current cell

                    new_ind = ind(x, y); # new cell location index
                    if cells[n].activated # if the current cell is activated
                        a_sim[new_ind[1], new_ind[2], i+1] += 1; # contribute activated cell at the current location in current timestep
                        a_tot_sims[i+1,sim] += 1; # contribute activated cell to total for this simulation
                    else # cell is naive
                        n_sim[new_ind[1], new_ind[2], i+1] += 1; # contribute naive cell at the current location in current timestep
                        n_tot_sims[i+1,sim] += 1; # contribute naive cell to total for this simulation
                    end
                end
            end # end for each cell

            for n in eachindex(dead_cells) # remove dead cells from the population
                deleteat!(cells, dead_cells[n]-(n-1)); # subtract n-1 to account for previously deleted cells
            end


            # simulate cytokine secretion and decay
            if mod(i,τ_ABM)==0
                sum_a = sum(a_sim[:,:,i])*result_scale; # total mass of activated T cells at previous timestep

                if use_ODE_for_IL2 # solve IL-2 ODE
                    yi = [I_i]; # value from previous timestep, [I(i-1)]
                    tspan = [i-1, i]; # previous and current times (separated by one timestep τ) (min)
                    p = (λ, S, sum_a, N_rods, V); # parameters and functions inside the ODE

                    prob = ODEProblem(IL2_ODE, yi, tspan, p); # ODE problem
                    sol = solve(prob, alg, dt=τ, adaptive=false); # solution with fields t and u
                    sol_arr = Array(sol); # solution as an array

                    I_i = sol_arr[2]; # next cytokine mass
                    if supplement_IL2 && in(i+1,supp_IL2_time) # if this is a time at which IL-2 is supplemented
                        I_i = I_i + supp_IL2_mass[findfirst(x->x==i+1,supp_IL2_time)]; # add supppemented mass of IL-2
                    end
                    I[i+1] = I_i/result_scale; # update ODE solution for cytokine

                else # solve IL-2 PDE
                    avg_a = sum_a/V; # average concentration of activated T cells

                    # matrices for cytokine (depends on current activated T cell concentration)
                    #A_I = genA_I(N, M, h, D_I, λ, a[:,:,i][:], BC); # generate A matrix (with dummy matrix for a)

                    for ii = 1:N
                        for jj = 1:M
                            #=
                            conc_scale = 1; # scale a mass to concentration for PDE solving
                            if (ii == 1) || (ii == N)
                                conc_scale = 2; # volume is half normal
                                if (jj == 1) || (jj == M)
                                    conc_scale = 4; # volume is quarter normal
                                end
                            elseif (jj == 1) || (jj == M)
                                conc_scale = 2; # volume is half normal
                            end
                            b_I[ii + N*(jj-1)] = S(ρ[ii,jj],i) - λ(a_sim[ii,jj,i]*conc_scale, I[ii+N*(jj-1),i÷τ_ABM]); # include RHS vector in matrix system dI/dt=AI+b (source term, S(ρ,t), and decay term if using IMEX)
                            =#
                            b_I[ii + N*(jj-1)] = S(ρ[ii,jj],i) - λ(avg_a, I[ii+N*(jj-1),i÷τ_ABM]); # include RHS vector in matrix system dI/dt=AI+b (source term, S(ρ,t), and decay term if using IMEX)
                        end
                    end
                    I[:,i÷τ_ABM+1] = linear_solver(Ã_I,b̃_I*I[:,i÷τ_ABM]+τ_ABM*b_I,I[:,i÷τ_ABM],lsolve,atol,rtol,maxiters_θ,pretype,precond); # solve matrix system Ãc=b̃
                    if supplement_IL2 && in(i+1,supp_IL2_time) # if this is a time at which IL-2 is supplemented
                        I[:,i÷τ_ABM+1] = I[:,i÷τ_ABM+1] .+ supp_IL2_mass[findfirst(x->x==i+1,supp_IL2_time)]/result_scale/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # add supppemented mass of IL-2
                    end
                end
            end

        end # end for each timestep

        
        # contribute simulation to average
        n += n_sim/ntotal;
        a += a_sim/ntotal;
        I_ABM = I_ABM + I/ntotal; # store the cytokine solution in the average

        if use_ODE_for_IL2
            I_mass = I; # ODE solution
        else
            I_reshaped = reshape(I, N, M, Ts_ABM+1); # reshaped PDE solution
            I_mass = sum(I_reshaped[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(I_reshaped[1,2:M-1,:],dims=1)+sum(I_reshaped[N,2:M-1,:],dims=1)+sum(I_reshaped[2:N-1,1,:],dims=1)+sum(I_reshaped[2:N-1,M,:],dims=1))[:]/2+(I_reshaped[1,1,:]+I_reshaped[1,M,:]+I_reshaped[N,1,:]+I_reshaped[N,M,:])/4; # cytokine mass
        end
        I_tot_sims[:,sim] = I_mass;
    end # end for each simulation 

    if !use_ODE_for_IL2
        I_ABM = reshape(I_ABM, N, M, Ts_ABM+1); # reshape spatial cytokine solution into a 3D array: c(x,y,t)
    end

    runtime = time() - start_time; # runtime for averaged ABM solving

    return n, a, I_ABM, n_tot_sims, a_tot_sims, I_tot_sims, runtime
end

function ABM_track_cells()
    # Simatules ABM once (or num_sim times) and track individual cells for trajectory and activation times

    Ts_ABM = T÷τ_ABM; # number of timesteps for IL-2 solving

    act_time = [[],[],[],[],[]]; # activation times for cells up to 5 generations
    act_prop = [0,0,0,0,0]; # average proportion of cells which activate before they die for each generation
    avg_runtime = 0; # average runtime over num_sims simulations

    if !isnan(seed_ABM) # set seed if not using default seed
        Random.seed!(seed_ABM); # apply RNG seed to T cells in ABM 
    end

    n = zeros(N, M, Int(floor(T))+1); # initialise T cell density array for naive T cells, n(x,y,t)
    a = zeros(N, M, Int(floor(T))+1); # initialise T cell density array for activated T cells, a(x,y,t)
    if use_ODE_for_IL2 # if the ODE is being used for IL-2
        I = zeros(Ts_ABM+1); # IL-2 density
    else
        I = zeros(N*M, Ts_ABM+1); # IL-2 density as a 2-dimensional array (for use in PDE solving)
    end
    cells_track = []; # array tracking all cells (dead and alive) [time of birth, time of death (NaN if not), time of activation (NaN if not), (x,y) positions over time, generation]
    
    @showprogress 1 "Simulating ABM..." for sim = 1:num_sims # for each simulation

        if place_cells
            cells = [cell_track(cell_xy[i][1], abs(cell_xy[i][2]-hei), false, i) for i in eachindex(cell_xy)]; # num_plot cells are naive and has chosen initial coordinates
            if length(cell_xy) < N0_track
                append!(cells, [cell_track(rand()*wid, rand()*hei, false, i) for i in length(cell_xy)+1:N0_track]); # the rest are naive with random positions (uniformly distributed)
            end
        else
            cells = [cell_track(rand()*wid, rand()*hei, false, i) for i in 1:N0_track]; # each cell is naive and has random initial coordinates (sampled from uniform distribution)
        end
        cells_track = []; # array tracking all cells (dead and alive) [time of birth, time of death (NaN if not), time of activation (NaN if not), (x,y) positions over time, generation]

        n = zeros(N, M, Int(floor(T))+1); # initialise T cell density array for naive T cells, n(x,y,t)
        a = zeros(N, M, Int(floor(T))+1); # initialise T cell density array for activated T cells, a(x,y,t)
        for cell = 1:N0_track # for each initial T cell
            x = cells[cell].x; y = cells[cell].y; # x and y positions of current cell
            append!(cells_track, [[1, NaN, NaN, zeros(2,Int(floor(T))+1), 1]]); # [time of birth, time of death (NaN if not), time of activation (NaN if not), (x,y) positions over time, generation]
            cells_track[cell][4][:,1] = [x,y];
            ind_n = ind(x, y); # current cell's coords in indices (as [i,j] array)
            if cells[cell].activated
                a[ind_n[1],ind_n[2],1] += 1; # contribute a new initial activated cell
            else
                n[ind_n[1],ind_n[2],1] += 1; # contribute a new initial activated cell
            end
        end

        if use_ODE_for_IL2 # if the ODE is being used for IL-2
            # setup ODE solver
            if cyt_IC == "zero" # if starting with some initial cytokine
                I_i = 0
            else
                I_i = I0_mass;
            end
            I = zeros(Ts_ABM+1); I[1] = I_i; # cytokine density for current simulation
            I_ABM = zeros(Ts_ABM+1); # stores the average cytokine mass over time
        else
            I = zeros(N*M, Ts_ABM+1); I[:,1] = IC_I[:]; # cytokine density as a 2-dimensional array (for use in PDE solving)
            I_ABM = zeros(N*M, Ts_ABM+1); # stores the average cytokine solution
        end

        start_time = time(); # define start time for averaged ABM solving
        # @showprogress 1 "Simulating ABM once..." 
        for i = 1:Int(floor(T)) # for each timestep
            dead_cells = []; # list of indices of cells that die in the current timestep
            for cell_n in eachindex(cells) # for each current T cell
                cell_alive = true; # initialise each cell to be alive
                # properties at current cell position
                x = cells[cell_n].x; y = cells[cell_n].y; # x and y positions of current cell
                a_n = cells[cell_n].activated; # activation state of current cell
                ID = cells[cell_n].ID; # current cell's ID

                ind_n = ind(x, y); # current cell's coords in indices (as [i,j] array)
                ρ_n = ρ[ind_n[1], ind_n[2]]; # scaffold density at current cell's position
                
                # cytokine (IL-2) density at current cell's position
                if use_ODE_for_IL2 # if using ODE for IL-2
                    I_n = I[i÷τ_ABM+1*(mod(i,τ_ABM)!=0)]/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # use the current (or most recently calcualted) value and average over space
                else # using PDE for IL-2
                    I_n = I[ind_n[1] + N*(ind_n[2]-1), i÷τ_ABM+1*(mod(i,τ_ABM)!=0)]; 
                end

                r_a_n = (1-a_n)*r_a(ρ_n)*τ_ABM; # current activation probability
                r_p_n = a_n*r_p(I_n)*τ_ABM; # current proliferation probability
                r_d_n = r_d*τ_ABM; # current death probability


                r = rand(); # sample random number ∈ [0,1]
                if r < r_a_n # cell activates
                    cells[cell_n].activated = true;
                    cells_track[ID][3] = i+1;
                elseif r_a_n <= r < r_a_n + r_p_n # cell proliferates
                    append!(cells, [cell_track(x, y, false, length(cells_track)+1)]); # create new naive T cell at current cell's location
                    n[ind_n[1], ind_n[2], i+1] += 1; # contribute to average naive cell density at the current location in current timestep for one daughter cell only
                    append!(cells_track, [[i+1, NaN, NaN, zeros(2,Int(floor(T))+1), cells_track[ID][5]+1]]); # [time of birth, time of death (NaN if not), time of activation (NaN if not), (x,y) positions over time, generation]
                    cells_track[end][4][:,i+1] = [x,y];
                elseif r_a_n + r_p_n <= r < r_a_n + r_p_n + r_d_n # cell dies
                    append!(dead_cells, cell_n); # record current cell as a dead cell (to be removed after the current loop)
                    cell_alive = false;
                    cells_track[ID][2] = i+1;
                elseif r >= r_a_n + r_p_n + r_d_n + w(ρ_n) # cell is not waiting at the current position
                    # movement probabilities
                    p_d_curr = p_d[ind_n[1], ind_n[2]]; # current probability of moving down (not considering boundary condition)
                    p_u_curr = p_u[ind_n[1], ind_n[2]]; # up
                    p_l_curr = p_l[ind_n[1], ind_n[2]]; # left
                    p_r_curr = p_r[ind_n[1], ind_n[2]]; # right
                    
                    x, y = move_cell(x, y, wid, hei, ε, BC, p_d_curr, p_u_curr, p_l_curr, p_r_curr); # move the cell based on movement probabilities and other factors
                end # otherwise cell waits
                cells_track[ID][4][:,i+1] = [x,y];

                if cell_alive # if the cell is still alive
                    cells[cell_n].x = x; cells[cell_n].y = y; # record new position of current cell

                    new_ind = ind(x, y); # new cell location index
                    if cells[cell_n].activated # if the current cell is activated
                        a[new_ind[1], new_ind[2], i+1] += 1; # contribute to average activated cell density at the current location in current timestep
                    else # cell is naive
                        n[new_ind[1], new_ind[2], i+1] += 1; # contribute to average naive cell density at the current location in current timestep
                    end
                end
            end # end for each cell

            for cell_n in eachindex(dead_cells) # remove dead cells from the population
                deleteat!(cells, dead_cells[cell_n]-(cell_n-1)); # subtract n-1 to account for previously deleted cells
            end


            # simulate cytokine diffusion
            if mod(i,τ_ABM)==0
                sum_a = sum(a[:,:,i])*result_scale; # total mass of activated T cells (mean-field approx.) at previous timestep

                if use_ODE_for_IL2 # solve IL-2 ODE
                    yi = [I_i]; # value from previous timestep, [I(i-1)]
                    tspan = [i-1, i]; # previous and current times (separated by one timestep τ) (min)
                    p = (λ, S, sum_a, N_rods, V); # parameters and functions inside the ODE

                    prob = ODEProblem(IL2_ODE, yi, tspan, p); # ODE problem
                    sol = solve(prob, alg, dt=τ_ABM, adaptive=false); # solution with fields t and u
                    sol_arr = Array(sol); # solution as an array

                    I_i = sol_arr[2]; # next cytokine mass
                    I[i+1] = I_i/result_scale; # update ODE solution for cytokine

                else # solve IL-2 PDE
                    avg_a = sum_a/V; # average concentration of activated T cells

                    # matrices for cytokine (depends on current activated T cell concentration)
                    #A_I = genA_I(N, M, h, D_I, λ, a[:,:,i][:], BC); # generate A matrix (with dummy matrix for a)
                    A_I = genA_I_IMEX(N, M, h, D_I, BC); # generate A matrix, assuming using IMEX method
                    Ã_I = spI-τ_ABM*θ*A_I; b̃_I = spI+τ_ABM*(1-θ)*A_I; # generate matrices used in Theta method

                    for ii = 1:N
                        for jj = 1:M
                            #=
                            conc_scale = 1; # scale a mass to concentration for PDE solving
                            if (ii == 1) || (ii == N)
                                conc_scale = 2; # volume is half normal
                                if (jj == 1) || (jj == M)
                                    conc_scale = 4; # volume is quarter normal
                                end
                            elseif (jj == 1) || (jj == M)
                                conc_scale = 2; # volume is half normal
                            end
                            b_I[ii + N*(jj-1)] = S(ρ[ii,jj],i) - λ(a_sim[ii,jj,i]*conc_scale, I[ii+N*(jj-1),i÷τ_ABM]); # include RHS vector in matrix system dI/dt=AI+b (source term, S(ρ,t), and decay term if using IMEX)
                            =#
                            b_I[ii + N*(jj-1)] = S(ρ[ii,jj],i) - λ(avg_a, I[ii+N*(jj-1),i÷τ_ABM]); # include RHS vector in matrix system dI/dt=AI+b (source term, S(ρ,t), and decay term if using IMEX)
                        end
                    end
                    I[:,i÷τ_ABM+1] = linear_solver(Ã_I,b̃_I*I[:,i÷τ_ABM]+τ_ABM*b_I,I[:,i÷τ_ABM],lsolve,atol,rtol,maxiters_θ,pretype,precond); # solve matrix system Ãc=b̃
                end
            end
        end # end for each timestep
        avg_runtime = time() - start_time; # runtime for averaged ABM solving

        # add each time from birth to activation and contribute to the average proportion of cells that activate
        gen_count = [0,0,0,0,0];
        for cell in cells_track 
            cell_gen = cell[5]; # current cell's generation
            if cell_gen <= 5 # if the cell generation is 5 or below
                gen_count[cell_gen] += 1;
                if !isnan(cell[3]) # if it activated before the end of simulation
                    append!(act_time[cell_gen],(cell[3]-cell[1]-1)*τ) # append time to activation (min)
                    act_prop[cell_gen] += 1; # increment counter for number of activated cells in each generation
                end
            end
        end
        act_prop = act_prop ./ gen_count; # proportion out of the total in each generation

    end # end for each simulation
    avg_runtime = avg_runtime/num_sims;

    Random.seed!(); # reset RNG seed

    if use_ODE_for_IL2
        I_ABM = I;
    else
        I_ABM = reshape(I, N,M,Ts_ABM+1); # reshape cytokine solution into a 3D array: c(x,y,t)
    end

    return n, a, I_ABM, cells_track, act_time, act_prop, avg_runtime
end

function simulate_PDE()
    # Simatules PDE model

    start_time = time(); # define start time for PDE solving

    Ts = T÷τ; # number of timesteps for PDE/ODE solving

    n = zeros(N*M, Ts+1); n[:,1] = IC_u[:]; # naive T cell density as a 2-dimensional array
    a = zeros(N*M, Ts+1); # activated T cell density as a 2-dimensional array
    a_PDE = zeros(N,M,Ts+1); # rearranged 3D array for activated cells: a(x,y,t)
    if start_with_active
        a[:,1] = IC_u[:]; # apply IC to activated cells
        a_PDE[:,:,1] = IC_u; 
    end

    # consider matrix equations for T cell densities du/dt = Au+b, where u=[u(1,1),u(2,1),...,u(N,1),u(1,2),u(2,2),...,u(N,M)]
    A_n = genA_n(N, M, h, D_n, χ_n, R_nn, ρ, BC); # generate A matrix for naive cells
    b_n = zeros(N*M, 1); # constant terms for naive cells
    Ã_n = spI-τ*θ*A_n; b̃_n = spI+τ*(1-θ)*A_n; # generate matrices used in Theta method for naive cells

    b_a = zeros(N*M, 1); # constant terms for activated cells

    if use_ODE_for_IL2 # if the cytokine needs to be updated mid-simulation
        # setup ODE solver
        if cyt_IC == "zero" # if starting with some initial cytokine
            I_i = 0;
        else
            I_i = I0_mass*result_scale;
        end
        I = zeros(Ts+1); I[1] = I_i/result_scale; # IL-2 mass
    else
        I = zeros(N*M, Ts+1); I[:,1] = IC_I[:]; # re-initialise cytokine concentration
    end

    @showprogress 1 "Solving continuum model..." for i = 1:Ts # for each timestep
        # current matrix of IL-2 concentrations 
        if use_ODE_for_IL2
            I_curr = zeros(N*M) .+ I[i]/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4);
        else
            I_curr = I[:,i];
        end

        # matrices for activated cells (depends on current cytokine concentration)
        A_a = genA_a(N, M, h, D_a, χ_aρ, R_aa(0), ρ, I_curr, BC); # generate A matrix for activated cells with current cytokine concentration
        Ã_a = spI-τ*θ*A_a; b̃_a = spI+τ*(1-θ)*A_a; # generate matrices used in Theta method for naive cells

        # compute RHS vectors in matrix systems for T cells dT/dt=AT+b (source term)
        for ii = 1:N
            for jj = 1:M
                b_n[ii + N*(jj-1)] = R_na(I_curr[ii+N*(jj-1)]) * a[ii+N*(jj-1), i]; # source term R_na(I)a
                b_a[ii + N*(jj-1)] = R_an(ρ[ii, jj]) * n[ii+N*(jj-1), i]; # source term R_an(ρ)n
            end
        end
        
        # naive T cells
        n[:,i+1] = linear_solver(Ã_n,b̃_n*n[:,i]+τ*b_n,n[:,i],lsolve,atol,rtol,maxiters_θ,pretype,precond); # solve matrix system Ãn=b̃

        # activated T cells
        a[:,i+1] = linear_solver(Ã_a,b̃_a*a[:,i]+τ*b_a,a[:,i],lsolve,atol,rtol,maxiters_θ,pretype,precond); # solve matrix system Ãa=b̃

        # simulate IL-2
        if use_ODE_for_IL2
            #a_PDE[:,:,i+1] = reshape(a[:,i+1], N, M); # compute mass of activated cells in scaled up domain
            #a_mass = (sum(a_PDE[2:N-1,2:M-1,i],dims=[1,2])[1]+(sum(a_PDE[1,2:M-1,i],dims=1)+sum(a_PDE[N,2:M-1,i],dims=1)+sum(a_PDE[2:N-1,1,i],dims=1)+sum(a_PDE[2:N-1,M,i],dims=1))[1]/2+(a_PDE[1,1,i]+a_PDE[1,M,i]+a_PDE[N,1,i]+a_PDE[N,M,i])/4)*result_scale;
            a_avg = mean(a[:,i]); # average activated cell concentration

            yi = [I_i/V]; # value from previous timestep, [I(i-1)]
            tspan = [(i-1)*τ, i*τ]; # previous and current times (separated by one timestep τ) (min)
            p = (λ, S, a_avg, ρ_avg); # parameters and functions inside the ODE

            prob = ODEProblem(IL2_ODE_conc, yi, tspan, p); # ODE problem
            sol = solve(prob, alg, dt=τ, adaptive=false); # solution with fields t and u
            sol_arr = Array(sol); # solution as an array

            I_i = sol_arr[2]*V; # next cytokine mass
            if supplement_IL2 && in(i+1,supp_IL2_time) # if this is a time at which IL-2 is supplemented
                I_i = I_i + supp_IL2_mass[findfirst(x->x==i+1,supp_IL2_time)]*result_scale; # add supppemented mass of IL-2
            end
            I[i+1] = I_i/result_scale; # update ODE solution for cytokine
        else
            # matrices for cytokine (depends on current activated T cell concentration)
            #A_I = genA_I(N, M, h, D_I, λ, a[:,i], BC); # generate A matrix (with dummy matrix for a)

            # compute RHS vector for cytokine
            #a_curr = reshape(a[:,i+1], N, M); # current spatial array for a
            #avg_a = (sum(a_curr[2:N-1,2:M-1])+(sum(a_curr[1,2:M-1])+sum(a_curr[N,2:M-1])+sum(a_curr[2:N-1,1])+sum(a_curr[2:N-1,M]))/2+(a_curr[1,1]+a_curr[1,M]+a_curr[N,1]+a_curr[N,M])/4)/V;
            for ii = 1:N
                for jj = 1:M
                    b_I[ii + N*(jj-1)] = S(ρ[ii,jj],i*τ) - λ(a[ii+N*(jj-1),i], I[ii+N*(jj-1),i]); # include RHS vector in matrix system dI/dt=AI+b (source term, S(ρ,t), and decay term if using IMEX)
                    #b_I[ii + N*(jj-1)] = S(ρ[ii,jj],i*τ) - λ(avg_a, I[ii+N*(jj-1),i]); # include RHS vector in matrix system dI/dt=AI+b (source term, S(ρ,t), and decay term if using IMEX)
                end
            end
            I[:,i+1] = linear_solver(Ã_I,b̃_I*I[:,i]+τ*b_I,I[:,i],lsolve,atol,rtol,maxiters_θ,pretype,precond); # solve matrix system ÃI=b̃
            if supplement_IL2 && in(i+1,supp_IL2_time) # if this is a time at which IL-2 is supplemented
                I[:,i+1] = I[:,i+1] .+ supp_IL2_mass[findfirst(x->x==i+1,supp_IL2_time)]/result_scale/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # add supppemented mass of IL-2
            end
        end
    end
    a_PDE = reshape(a, N,M,Ts+1); # reshape T cell solutions into 3D arrays: a(x,y,t) and n(x,y,t)
    n_PDE = reshape(n, N,M,Ts+1); 
    if use_ODE_for_IL2
        I_PDE = I;
    else
        I_PDE = reshape(I, N,M,Ts+1); # reshape cytokine solution into a 3D array: I(x,y,t)
    end

    runtime = time() - start_time; # runtime for PDE solving

    return n_PDE, a_PDE, I_PDE, runtime
end

function simulate_ODE()
    # Simatules ODE model

    start_time = time(); # define start time for ODE solving

    V = N*M*result_scale; # full domain volume

    #ρ_avg = APC_ms_conc/m_MSR*mLtoh2; # average rod density: rods/mL times ratio mL per lattice site, alternatively N_rods/(N*M*result_scale)

    #S_ODE(ρ,t) = k*I_loaded*exp(-k*t); # time-dependent source term for cytokine (same as S(N_rods,t) but without K_S)

    function ODE_system(dy,y,p,t)
        R_nn, R_na, R_aa, R_an, λ, S, ρ, V = p; # R_nn and R_an are functions of ρ, R_na and R_aa are functions of I, λ is a function of a and I, and S is a function of ρ and t
        dy[1] = R_nn(ρ/V)*y[1] + R_na(y[3]/V)*y[2]; # dN_n/dt
        dy[2] = R_aa(y[3]/V)*y[2] + R_an(ρ/V)*y[1]; # dN_a/dt
        dy[3] = -λ(y[2]/V,y[3]/V)*V + S(ρ/V,t)*V; # dN_I/dt
    end
    alg = RK4(); # algorithm for solving ODE system
    p = (R_nn, R_na, R_aa, R_an, λ, S, N_rods, V); # parameters and functions inside the ODE system

    if supplement_IL2 # if supplementing IL-2 at specific times, solve ODE iteratively
        I_ODE = zeros(Int(floor(T))+1); # initialise cytokine mass
        if cyt_IC != "zero" # if starting with some initial cytokine
            I_ODE[1] = I0_mass*result_scale; 
        end
        n_ODE = zeros(Int(floor(T))+1); # initialise naive cells
        a_ODE = zeros(Int(floor(T))+1); # initialise activated cells
        if start_with_active # if starting with 50/50 naive/activated cells
            n_ODE[1] = T_loaded/2;
            a_ODE[1] = T_loaded/2;
        else
            n_ODE[1] = T_loaded;
        end

        for i = 1:Int(floor(T)) # for each timestep
            yi = [n_ODE[i], a_ODE[i], I_ODE[i]]; # values from previous timestep, [n((i-1)τ), a((i-1)τ), I((i-1)τ)]
            tspan = [(i-1)*τ, i*τ]; # previous and current times (separated by one timestep τ) (min)
            
            prob = ODEProblem(ODE_system, yi, tspan, p); # ODE problem
            sol = solve(prob, alg, dt=τ, adaptive=false); # solution with fields t and u
            sol_arr = Array(sol); # solution as an array

            I_i = sol_arr[1,2]; # next cytokine mass
            n_ODE[i+1] = sol_arr[1,2];# update ODE solutions
            a_ODE[i+1] = sol_arr[2,2];
            I_ODE[i+1] = sol_arr[3,2];
            if supplement_IL2 && in(i+1,supp_IL2_time) # if this is a time at which IL-2 is supplemented
                I_ODE[i+1] = I_ODE[i+1] + supp_IL2_mass[findfirst(x->x==i+1,supp_IL2_time)]; # add supppemented mass of IL-2
            end
        end
        T_ODE = n_ODE + a_ODE;
    else
        if start_with_active # if starting with 50/50 naive/activated cells
            y0 = [T_loaded/2, T_loaded/2, 0]; # initial condition, [n(0), a(0), I(0)]
        else
            y0 = [T_loaded, 0, 0];
        end
        if cyt_IC != "zero"
            y0[3] = I0_mass*result_scale;
        end
        tspan = [0.0, T]; # timespan (min)

        prob = ODEProblem(ODE_system, y0, tspan, p); # ODE problem
        sol = solve(prob, alg, dt=τ, adaptive=false); # solution with fields t and u
        sol_arr = Array(sol); # solution as an array

        t = sol.t; # time vector
        n_ODE = sol_arr[1,:]; # ODE solution for naive cells
        a_ODE = sol_arr[2,:]; # ODE solution for activated cells
        I_ODE = sol_arr[3,:]; # ODE solution for cytokine
    end

    runtime = time() - start_time; # runtime for ODE solving

    return n_ODE, a_ODE, I_ODE, runtime
end