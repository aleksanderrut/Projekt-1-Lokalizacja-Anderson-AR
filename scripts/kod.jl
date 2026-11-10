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

# Liczenie wspołczynników Lapunowa dla całego wektora długości łańcucha
function gamma_rekurencja(w, E, t, L_values)
    k = acos(complex(-E / (2*t)))
    psi0 = exp(im * k)
    psi1 = 1.0 + 0im

    suma_log = 0.0
    gamma = zeros(length(L_values))
    indeks = 1
    for l in 1:maximum(L_values)
        # Rekurencja
        psi2 = (w[l] - E)/t * psi1 - psi0
        norma = sqrt(abs2(psi2) + abs2(psi1))
        suma_log += log(norma)

        psi0 = psi1 / norma
        psi1 = psi2 / norma

        # Zapis gamma dla wybranych długości L
        if l == L_values[indeks]
            gamma[indeks] = suma_log / l
            indeks += 1
        end
    end

    return gamma
end

# Hamiltoanin model Aubry'ego-Andre 1D jednnocząstkowy
function gen_H_aubry_andre_1D(L, t, λ, PBC::Bool)
    N = L
    H = zeros(Float64, N, N)
    β = (sqrt(5) - 1) / 2 
    for i in 1:L
        H[i, i] = λ * cos(2*π*β*i) # potencjał Aubry'ego-Andre

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
        # seed = 1234
        # Random.seed!(seed) #losowe ale zawsze takie same póki jest takie same ziarno
        # L_x = 2500
        # L_y = 50
        # L_z = 3
        # t = 1.0
        # W = 0.5
        # PBC = false

        # H1 = gen_H_anderson_1D(L_x, t, W, PBC)
        # # H2 = gen_H_anderson_2D(L_x, L_y, t, W, PBC)
        # # H3 = gen_H_anderson_3D(L_x, L_y, L_z, t, W, PBC)
        # # save_data(H2, "Hamiltonian_2D", W, PBC; seed, L_x, L_y)

        # ED = eigen(Symmetric(H1))
        # eigenvalues = ED.values
        # eigenvectors = ED.vectors
        # save_data(eigenvalues, "eigenvalues_1D", W, PBC; seed, L_x)
        # save_data(eigenvectors, "eigenvectors_1D", W, PBC; seed, L_x)


#     ### ZAD 2
#     L_x = 1000
#     t = 1.0
#     W = 8.0
#     PBC = true
#     R = 500

#     rd = RandomDevice()
#     eigenvalues_R = zeros(Float64, L_x, R)
#     IPR_R = zeros(Float64, L_x, R)

#     # Pętla po realizacjach nieporządku
#     for r in 1:R
#         # println("Realizacja $r / $R")

#         seed = rand(rd, UInt64)
#         Random.seed!(seed)

#         H1 = gen_H_anderson_1D(L_x, t, W, PBC)
#         ED = eigen(Symmetric(H1))

#         eigenvalues = ED.values
#         eigenvectors = ED.vectors

#         # IPR dla wszystkich wektorów własnych, każda kolumna macierzy to jeden wektor własny
#         IPR = vec(sum(abs2.(eigenvectors).^2, dims=1)) # dims -> oznacza sumowanie po wierszach, osobno dla każdej kolumny

#         eigenvalues_R[:, r] = eigenvalues 
#         IPR_R[:, r] = IPR
#     end

#     # Uśrednianie po realizacjach
#     eigenvalues_mean = vec(sum(eigenvalues_R, dims=2)) / R # dims=2 oznacza sumowanie po każdej kolumie, osbno dla każdego wiersza 
#     IPR_mean = vec(sum(IPR_R, dims=2)) / R

#     # Zapis średnich: energia, IPR
#     #IPR_mean_od_E = hcat(eigenvalues_mean, IPR_mean)
#     #save_data(IPR_mean_od_E, "IPR_mean_1D_R$(R)",W, PBC; L_x)

#     # Binowanie
#     delta_E = 0.1
#     E_min = -6
#     E_max = 6
#     N_binow = round(Int, (E_max - E_min) / delta_E)

#     IPR_suma = zeros(N_binow)
#     liczba_stanow = zeros(Int, N_binow)

#     # Przypisywanie IPR do odpowiednich koszyków ze względu na średnią wartość enrgi przypisaną do IPR'ów
#     for i in eachindex(eigenvalues_mean)
#         E = eigenvalues_mean[i]
#         IPR = IPR_mean[i]
#         numer_binu = floor(Int, (E - E_min) / delta_E) + 1
#         IPR_suma[numer_binu] += IPR
#         liczba_stanow[numer_binu] += 1
#     end

#     energie_binow = zeros(N_binow)
#     IPR_bin = zeros(N_binow)
#     for i in 1:N_binow
#         energie_binow[i] = E_min + (i - 1/2) * delta_E # śreodek koszyka
#         if liczba_stanow[i] > 0
#             IPR_bin[i] = IPR_suma[i] / liczba_stanow[i]
#         end
#     end

#     IPR_data = hcat( energie_binow, IPR_bin, liczba_stanow)
#     save_data(IPR_data, "IPR_1D_R$(R)", W, PBC; L_x)

    
#     ### ZAD 3
#     L_x = 10
#     L_y = 10
#     L_z = 10
#     t = 1.0
#     PBC = true
#     R = 70

#     W_values = 1.0:0.5:20.0
#     N_W = length(W_values)

#     N = L_x * L_y * L_z
#     delta_E = 0.1
#     E_min = -16.0
#     E_max = 16.0
#     N_binow = round(Int, (E_max - E_min) / delta_E)

#     energie_binow = zeros(N_binow)
#     for i in 1:N_binow
#         energie_binow[i] = E_min + (i - 0.5) * delta_E
#     end

#     # Wiersze -> W
#     # Kolumny -> biny energetyczne
#     IPR_matrix = fill(NaN, N_W, N_binow)
#     liczba_stanow_matrix = zeros(Int, N_W, N_binow)

#     rd = RandomDevice()

#     for (w_idx, W) in enumerate(W_values)
#         eigenvalues_R = zeros(Float64, N, R)
#         IPR_R = zeros(Float64, N, R)
#         for r in 1:R
#             seed = rand(rd, UInt64)
#             Random.seed!(seed)
#             H3 = gen_H_anderson_3D(L_x, L_y, L_z, t, W, PBC)
#             ED = eigen(Symmetric(H3))
#             eigenvalues = ED.values
#             eigenvectors = ED.vectors

#             # IPR dla wszystkich wektorów własnych
#             IPR = vec(sum(abs2.(eigenvectors).^2, dims=1))
#             eigenvalues_R[:, r] = eigenvalues
#             IPR_R[:, r] = IPR
#         end

#         eigenvalues_mean = vec(sum(eigenvalues_R, dims=2)) / R
#         IPR_mean = vec(sum(IPR_R, dims=2)) / R

#         IPR_suma = zeros(N_binow)
#         liczba_stanow = zeros(Int, N_binow)

#         for i in eachindex(eigenvalues_mean)
#             E = eigenvalues_mean[i]
#             IPR = IPR_mean[i]
#             numer_binu = floor(Int, (E - E_min) / delta_E) + 1
#             IPR_suma[numer_binu] += IPR
#             liczba_stanow[numer_binu] += 1
#         end

#         for i in 1:N_binow
#             if liczba_stanow[i] > 0
#                 IPR_matrix[w_idx, i] = IPR_suma[i] / liczba_stanow[i]
#             end
#             liczba_stanow_matrix[w_idx, i] = liczba_stanow[i]
#         end
#     end

#     save_data(IPR_matrix, "IPR_matrix_3D_R$(R)", 0, PBC; L_x, L_y, L_z)




    # ### zad 5 od L
    # t = 1.0
    # #E_values = [-3.0,-2.0,-1.0, 0.0, 1.0, 2.0, 3.0]
    # E_values = vcat(-3.0:0.05:3.0)
    # W = 2.0
    # PBC = false
    # R = 100

    # #L_values = vcat(1, 10:10:1000, 1100:100:10000, 11000:1000:1000000, 1100000:100000:10000000)
    # L_values = [100000]
    # L_max = maximum(L_values)
    # gamma_R = zeros(length(L_values), R)

    # for E in E_values
    #     gamma_R = zeros(length(L_values), R)

    #     # Pętla po realizacjach nieporządku
    #     for r in 1:R
    #         rd = RandomDevice()
    #         seed = rand(rd, UInt64)
    #         Random.seed!(seed)
    #         w = gen_pot_Anderson(L_max, W) 
    #         println("w[1] = $(w[1])")#sprawdzałem czy serio loswe
    #         gamma_R[:, r] = gamma_rekurencja(w, E, t, L_values)
    #     end

    #     gamma_mean = vec(sum(gamma_R, dims=2)) / R # Uśrednianie po realizacjach
    #     gamma_data = hcat(L_values, gamma_mean) # Pierwsza kolumna: L, druga kolumna: gamma
    #     save_data(gamma_data, "gamma_1D_E$(E)_R$(R)", W, PBC; L_x=L_max) # Zapis wyników
    # end 

    # ### zad 5 od E
    # t = 1.0
    # E_values = -2.0:0.02:2.0
    # W = 5.0
    # PBC = false
    # R = 100

    # L_values = [100000]
    # L_max = maximum(L_values)
    # gamma_mean_values = zeros(length(E_values))
    # for (i, E) in enumerate(E_values)
    #     gamma_R = zeros(R)
    #     # Pętla po realizacjach nieporządku
    #     for r in 1:R
    #         rd = RandomDevice()
    #         seed = rand(rd, UInt64)
    #         Random.seed!(seed)
    #         w = gen_pot_Anderson(L_max, W)
    #         #println("w[1] = $(w[1])")#sprawdzałem czy serio loswe
    #         gamma_R[r] = gamma_rekurencja(w, E, t, L_values)[1]
    #     end
    #     gamma_mean_values[i] = sum(gamma_R) / R # uśrednianie po realizacjach
    #     println("E=$E, gamma=$(gamma_mean_values[i])")
    # end
    # gamma_data = hcat(collect(E_values), gamma_mean_values) # Pierwsza kolumna: E, druga kolumna: średnie gamma
    # save_data(gamma_data, "gamma_vs_E_R$(R)", W, PBC; L_x=L_max) 
    

    # ### zad 5 od W 
    # t = 1.0 
    # E_values = [0.0, 1.0, 1.5]
    # W_values = 0.1:0.02:10.0
    # PBC = false
    # R = 100

    # L_values = [100000]
    # L_max = maximum(L_values)
    # for E in E_values
    #     gamma_mean_values = zeros(length(W_values))
    #     for (i, W) in enumerate(W_values)
    #         gamma_R = zeros(R)
    #         # Pętla po realizacjach nieporządku
    #         for r in 1:R
    #             rd = RandomDevice()
    #             seed = rand(rd, UInt64)
    #             Random.seed!(seed)
    #             w = gen_pot_Anderson(L_max, W)
    #             gamma_R[r] = gamma_rekurencja(w, E, t, L_values)[1]
    #         end
    #         gamma_mean_values[i] = sum(gamma_R) / R # Uśrednianie po realizacjach
    #         println("W=$W, gamma=$(gamma_mean_values[i])")
    #     end
    #     gamma_data = hcat(collect(W_values), gamma_mean_values) # Pierwsza kolumna: W, druga kolumna: średnie gamma
    #     save_data(gamma_data, "gamma_vs_W_E$(E)_R$(R)", "0.1-10.0", PBC; L_x=L_max)
    #     end

    ### ZAD 6
    L_values = [377, 610, 987, 1597, 2584]
    t = 1.0
    λ_values = 0.0:0.1:4.0
    PBC = true
    for L_x in L_values
        println("Aktualnie liczone L=$L_x")
        IPR_mean = zeros(length(λ_values))
        sigma_rel = zeros(length(λ_values))
        # Środkowa część widma: 0.35L-0.65L 
        indeks_min = ceil(Int, 0.35 * L_x)
        indeks_max = floor(Int, 0.65 * L_x)

        for (i, λ) in enumerate(λ_values)
            H_aa = gen_H_aubry_andre_1D(L_x, t, λ, PBC)
            ED = eigen(Symmetric(H_aa))
            #eigenvalues = ED.values
            eigenvectors = ED.vectors

            # IPR dla wszystkich wektorów własnych, każda kolumna macierzy to jeden wektor własny
            IPR = vec(sum(abs2.(eigenvectors).^2, dims=1)) # dims -> oznacza sumowanie po wierszach, osobno dla każdej kolumny
            #Liczenie średniej IPR i sigma_rel dla środkowej części widma
            IPR_srodek = IPR[indeks_min:indeks_max]
            IPR_mean[i] = sum(IPR_srodek) / length(IPR_srodek)
            IPR2_mean = sum(IPR_srodek.^2) / length(IPR_srodek)
            sigma_rel[i] = (IPR2_mean - IPR_mean[i]^2) / IPR_mean[i]^2
        end
        IPR_aa_data = hcat(λ_values, IPR_mean, sigma_rel)
        save_data(IPR_aa_data, "IPR_sigma_rel_AA", "lambda_0.2-4.0", PBC; L_x=L_x)
    end

end
main()