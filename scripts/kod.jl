using Random
using DelimitedFiles
using LinearAlgebra

# potencjal Andersona; w_i z zakresu [-W/2, W/2]
function gen_pot_Anderson(N, W)
    w = zeros(Float64, N)
    for i in 1:N
        w[i] = W * (rand() - 0.5)
    end
    return w
end

# Hamiltoanin 1D jednnocząstkowy
function gen_H_anderson_1D(L, t, W, PBC::Bool)
    N = L
    H = zeros(Float64, N, N)
    w = gen_pot_Anderson(N, W)

    for i in 1:L
        H[i, i] = w[i] #potencjał Andersona na każdym węźle

        if i < L
            H[i, i + 1] = -t
            H[i + 1, i] = -t
        end
    end

    if PBC == 1 && L > 2
        H[1, L] = -t
        H[L, 1] = -t
    end

    return H
end

# zwraca numerr wezła w przypadku 2D
function site_number_2D(x, y, Lx)
    return x + (y - 1) * Lx
end

# Hamiltoanin 2D jednnocząstkowy
function gen_H_anderson_2D(Lx, Ly, t, W, PBC::Bool)
    N = Lx * Ly
    H = zeros(Float64, N, N)
    w = gen_pot_Anderson(N, W)

    for y in 1:Ly
        for x in 1:Lx
            i = site_number_2D(x, y, Lx)
            # potencjał Andersona na każdym węźle
            H[i, i] = w[i]

        # przeskok w kierunku x
        if x < Lx
            j = site_number_2D(x + 1, y, Lx)
            H[i, j] = -t
            H[j, i] = -t
        elseif PBC == 1 && Lx > 2
            j = site_number_2D(1, y, Lx)
            H[i, j] = -t
            H[j, i] = -t
        end

        # przeskok w kierunku y
        if y < Ly
            j = site_number_2D(x, y + 1, Lx)
            H[i, j] = -t
            H[j, i] = -t
        elseif PBC == 1 && Ly > 2
            j = site_number_2D(x, 1, Lx)
            H[i, j] = -t
            H[j, i] = -t
        end
        end
    end
    return H
end

# zwraca numerr wezła w przypadku 3D
function site_number_3D(x, y, z, Lx, Ly)
    return x + (y - 1) * Lx + (z - 1) * Lx * Ly
end

# Hamiltoanin 3D jednnocząstkowy
function gen_H_anderson_3D(Lx, Ly, Lz, t, W, PBC::Bool)
    N = Lx * Ly * Lz
    H = zeros(Float64, N, N)
    w = gen_pot_Anderson(N, W)

    for z in 1:Lz
    for y in 1:Ly
        for x in 1:Lx
            i = site_number_3D(x, y, z, Lx, Ly)
            H[i, i] = w[i] # potencjal Andersona

            # przeskol w kierunku x
            if x < Lx
                j = site_number_3D(x + 1, y, z, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t
            elseif PBC == 1 && Lx > 2
                j = site_number_3D(1, y, z, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t
            end

            # przeskok w kierunku y
            if y < Ly
                j = site_number_3D(x, y + 1, z, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t

            elseif PBC == 1 && Ly > 2
                j = site_number_3D(x, 1, z, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t
            end

            # przeskok w kierunku z
            if z < Lz
                j = site_number_3D(x, y, z + 1, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t
            elseif PBC == 1 && Lz > 2
                j = site_number_3D(x, y, 1, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t

            end
        end
    end
    end
    return H
end

# Zapis do pliku .data
function save_data(data, filename, W, PBC; seed=nothing, L_x=nothing, L_y=nothing, L_z=nothing)
    folder = raw"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
    filename_full = "$(filename)_W$(W)_PBC$(PBC)"
    if L_x !== nothing
        filename_full *= "_Lx$(L_x)"
    end
    if L_y !== nothing
        filename_full *= "_Ly$(L_y)"
    end
    if L_z !== nothing
        filename_full *= "_Lz$(L_z)"
    end
    if seed !== nothing
        filename_full *= "_seed$(seed)"
    end
    filename_full *= ".data"
    filepath = joinpath(folder, filename_full)

    writedlm(filepath, data)
    println("Zapisano do pliku:")
    println(filepath)
end

####################################################################################################################
function main()
    ### ZAD 1
    # seed = 6767
    # Random.seed!(seed) #losowe ale zawsze takie same póki jest takie same ziarno
    # L_x = 50
    # L_y = 50
    # L_z = 3
    # t = 1.0
    # W = 8.0
    # PBC = true

    # # H1 = gen_H_anderson_1D(L_x, t, W, PBC)
    # H2 = gen_H_anderson_2D(L_x, L_y, t, W, PBC)
    # # H3 = gen_H_anderson_3D(L_x, L_y, L_z, t, W, PBC)
    # # save_data(H2, "Hamiltonian_2D", W, PBC; seed, L_x, L_y)

    # ED = eigen(Symmetric(H2))
    # eigenvalues = ED.values
    # eigenvectors = ED.vectors
    # save_data(eigenvalues, "eigenvalues_2D", W, PBC; seed, L_x, L_y)
    # save_data(eigenvectors, "eigenvectors_2D", W, PBC; seed, L_x, L_y)

    
    ### ZAD 2
    L_x = 1000
    t = 1.0
    W = 0.5
    PBC = true
    R = 100

    seed_start = time_ns()
    eigenvalues_R = zeros(Float64, N, R)
    IPR_R = zeros(Float64, N, R)

    # Pętla po realizacjach nieporządku
    for r in 1:R
        # println("Realizacja $r / $R")

        seed = seed_start + (r - 1)
        Random.seed!(seed)

        H1 = gen_H_anderson_1D(L_x, t, W, PBC)
        ED = eigen(Symmetric(H1))

        eigenvalues = ED.values
        eigenvectors = ED.vectors

        # IPR dla wszystkich wektorów własnych, każda kolumna macierzy to jeden wektor własny
        IPR = vec(sum(abs2.(eigenvectors).^2, dims=1)) # dims -> oznacza sumowanie po wierszach, osobno dla każdej kolumny

        eigenvalues_R[:, r] = eigenvalues 
        IPR_R[:, r] = IPR
        IPR_od_E = hcat(eigenvalues, IPR)

        #save_data(IPR_od_E, "IPR_1D_R$(R)_r$(r)",W, PBC; seed, L_x)
    end

    # Uśrednianie po realizacjach
    eigenvalues_mean = vec(sum(eigenvalues_R, dims=2)) / R
    IPR_mean = vec(sum(IPR_R, dims=2)) / R

    # Zapis średnich: energia, IPR
    IPR_mean_od_E = hcat(eigenvalues_mean, IPR_mean)
    save_data(IPR_mean_od_E, "IPR_mean_1D_R$(R)",W, PBC; seed=seed_start, L_x)

    end

main()