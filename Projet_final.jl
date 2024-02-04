## Importation des packages nécessaires. 
using Plots
using LaTeXStrings
using Statistics
using Random
using LsqFit
using Trapz
using EasyFit

gr() # Choix du backend pour le package Plots.

## Paramètres initiaux.

taille_Réseau = 50 # Réseau carré.
nb_particules = 400 # Nombre de particules initial.

itermax = 100000 # Itération temporelle maximale.

beta = 0.5 # Coefficients à spécifier pour la friction dynamique.
alpha = 0.5

## Initialisation des paramètres aléatoires.

seed = 1649267719 # Germe fixé pour les calculs du rapport. On peut prendre seed = #Integer(round(time())) pour que le germe corresponde au nombre de seconde écoulé depuis le 1er janvier 1970.

print("$(seed)\n") # On imprime le seed pour pouvoir reproduire le résultat ensuite.

rng = MersenneTwister(seed) # Initialisation de la distribution Mersenne Twister

# Fonctions utilisées dans la routine principale.

function prob_ag(positionx, positiony, masse, beta, alpha, i, j) # Fonction qui calcule la probabilité d'aggrégation.

    F = Vector{Float64}() # Initialisation des vecteurs.
    r = Vector{Float64}()

    for k in 1:size(positionx)[1] # On itère sur toutes les particules.

        if (k != i) && (k != j) # On veut toutes les particules sauf celles qui entre en collision, soit i et j. 

            diff = round(sqrt((positionx[k] - positionx[i])^2 + (positiony[k] - positiony[i])^2), digits=1) # On prend la différence entre les vecteurs de position.
            fk = exp(-diff / (masse[k]^beta)) # On calcule l'effet de la particule sur la collision.

            push!(F, fk) # On rajoute les valeurs pertinentes aux vecteurs F et r.
            push!(r, round(sqrt(positionx[k]^2 + positiony[k]^2), digits=1))

        end

    end

    perm = sortperm(r) # On organise r du plus petit au plus grand et on organise F avec le même ordre.
    r = r[perm]
    F = F[perm]

    P = trapz(r, F)^alpha # On intègre F sur r et on met à la puissance alpha pour obtenir la probabilité selon la friction dynamique.

    return P # On retourne la probabilité.

end

function correl(positionx, positiony, masse, taille_Réseau) # Fonction qui calcule la corrélation spatiale (avec la masse) entre les particules.

    r = vec(0:0.1:round(sqrt(2 * ((taille_Réseau + 1)^2)), digits=1)) # On définit un vecteur r qui peut prendre toutes les valeurs possibles de r dans notre réseau.
    ri = round.(sqrt.(positionx .^ 2 .+ positiony .^ 2), digits=1) # ri correspond à la position des particules.
    perm = sortperm(ri) # On veut réorganiser ri du plus petit au plus grand et le tableau de masse avec la même permutation.

    ri = ri[perm]
    masse = masse[perm]

    Aire = (taille_Réseau + 1)^2 * π # L'aire totale du réseau en polaire.
    rho_tot = zeros(size(ri)[1], size(r)[1]) # On initialise un tableau de fonctions de densité.

    for i in 1:size(ri)[1] # On itère sur chaque particule.

        rho_0 = masse[i] # rho_0 est la densité à la position de la particule.

        rho_0_plus_r = [] # On veut calculer rho à ri + r pour tous les r.

        for j in 1:size(r)[1] # On itère sur toutes les positions r possibles.

            b = findall(x -> x == ri[i] + r[j], ri) # On trouve les particules à la position ri + r.

            if b != []

                rho_iter = sum(masse[b]) / size(b)[1] # On fait la moyenne des densités à cette position.

            else

                rho_iter = 0 # S'il n'y a pas de particules on met 0 pour la densité.

            end

            push!(rho_0_plus_r, rho_iter)

        end

        rho_tot[i, :] .= rho_0 .* rho_0_plus_r # Une ligne de rho_tot correspond à la densité par rapport à une particule pour tous les r.

    end

    C = Vector{Float64}()

    for k in 1:size(r)[1] # On veut calculer la moyenne d'ensemble sur toutes les particules.

        C_iter = trapz(ri, rho_tot[:, k]) / Aire # On intègre sur ri et on divise par l'aire.

        push!(C, C_iter) # On obtient C(r) pour tout r à la fin.

    end

    return C, r # On retourne C et r.

end

function initialisation(nb_particules, taille_Réseau, rng) # Fonction pour initialiser le réseau.

    direction = ones(nb_particules) # On initialise tous les tableaux de position, direction, vitesse et masse.
    masse = ones(nb_particules)

    positionx = ones(nb_particules)
    positiony = ones(nb_particules)
    vitesse = ones(nb_particules)

    for i in 1:nb_particules # On assigne des positions aléatoires en x et y sur le réseau.

        positionx[i] = rand(rng, 1:1:taille_Réseau)[1]
        positiony[i] = rand(rng, 1:1:taille_Réseau)[1]

    end

    return masse, direction, vitesse, positionx, positiony # On retourne tous les tableaux initialisés.

end

function changement_position(positionx, positiony, direction, vitesse, taille_Réseau, rng) # Fonctin qui calcule le changement de position à chaque itération temporelle.

    for i in 1:size(positionx)[1]

        direction[i] = rand(rng, 0:0.1:2π)[1] # On définit des directions aléatoires pour chaque particule.

    end

    for i in 1:size(positionx)[1]

        positionx[i] = round(positionx[i] + vitesse[i] * cos(direction[i]), digits=1) # On calcule la nouvelle position pour chaque particule.
        positiony[i] = round(positiony[i] + vitesse[i] * sin(direction[i]), digits=1)

        if positionx[i] < 0 # Conditions pour les particules qui sortent du réseau, elles collent aux frontières jusqu'à la prochaine itération (frontières à taille_Réseau+1).

            positionx[i] = 0

        end

        if positionx[i] > taille_Réseau + 1

            positionx[i] = taille_Réseau + 1

        end

        if positiony[i] < 0

            positiony[i] = 0

        end

        if positiony[i] > taille_Réseau + 1

            positiony[i] = taille_Réseau + 1

        end

    end

    return positionx, positiony, direction # On retourne les nouveaux tableaux.

end

function collision(positionx, positiony, direction, vitesse, masse, t, beta, alpha, rng) # Fonction qui identifie les collisions et fait le changement de vitesse. Fonctionne pour des collisions de deux particules seulement.

    collision = [] # Tableau qui prendra les indices des particules en collision.

    for j in 1:size(positionx)[1] # On itère sur toutes les positions 2 fois pour avoir chaque combinaison i et j.

        for i in 1:size(positionx)[1]

            if (i != j) # On ne veut pas que ce soit la même particule.

                if (positionx[i] == positionx[j]) && (positiony[i] == positiony[j]) # Si les particules sont à la même position, on regarde s'il y a une collision.

                    ra = rand(rng)[1] # Nombre aléatoire entre 0 et 1.
                    pa = prob_ag(positionx, positiony, masse, beta, alpha, i, j) # La probabilité de collision est définie selon la fonction prob_ag.


                    if ra < pa # Si le nombre aléatoire ra est plus petit que la probabilité, il y a collision.

                        push!(collision, sort([i, j])) # On rajoute la paire d'indices de la collision au tableau de collision.

                        vitesse[i] += -sqrt((2 * pa) / masse[i]) # La perte d'énergie est proportionnelle à prob_ag, on calcule la perte de vitesse avec E=0.5mv^2.
                        vitesse[j] += -sqrt((2 * pa) / masse[j])

                        pv = 1 - ((masse[i] * vitesse[i]) + (masse[j] * vitesse[j])) / (masse[i] + masse[j]) # Probabilité que la nouvelle particule reste stationnaire après la collision.
                        rv = rand(rng)[1] # Nombre aléatoire entre 0 et 1.

                        if rv < pv # Si rv est plus petit que la probabilité, les vitesses sont mises à 0 pour les particules en collision.

                            vitesse[i] = 0
                            vitesse[j] = 0

                        end


                    end

                end

            end


        end

    end


    collision = unique(sort, collision) # On veut supprimer les doubles [i,j] et [j,i] et ordonner les paires d'indices.
    print("$(collision), $(t)\n") # On imprime les particules en collision avec l'indice temporel.

    supprimer = [] # Tableau pour les indices à supprimer.

    if collision != []

        for i in 1:size(collision)[1] # On itère sur toutes les collisions.

            masse[collision[i][1]] += masse[collision[i][2]] # On ajoute la masse de la particule j à celle de la particule i.
            push!(supprimer, collision[i][2])

        end


        deleteat!(positionx, sort(supprimer)) # On supprime l'indice j de chaque collision dans tous les tableaux, i est la nouvelle particule.
        deleteat!(positiony, sort(supprimer))
        deleteat!(direction, sort(supprimer))
        deleteat!(masse, sort(supprimer))
        deleteat!(vitesse, sort(supprimer))

    end

    return positionx, positiony, direction, vitesse, masse # On retourne les tableaux modifiés.

end

## Fonction pour la routine principale.

function routine_principale(nb_particules, taille_Réseau, itermax, rng, beta, alpha) # Fonction qui fait le calcul principal et utilise les autres fonctions.

    masse, direction, vitesse, positionx, positiony = initialisation(nb_particules, taille_Réseau, rng) # On initialise les tableaux.

    scatter(positionx, positiony, markershape=:circle, markercolor=:black, markersize=masse, dpi=300, legend=false) # On crée un graphique qui montre les particules sur le réseau initialement.
    display(scatter!(xlab=L"$x$", ylab=L"$y$"))
    savefig("/home/laurent/Documents/PHY 3075/Projet final/position_initiale.png")

    for t in 1:itermax # On itère sur le temps.

        positionx, positiony, direction = changement_position(positionx, positiony, direction, vitesse, taille_Réseau, rng) # On calcule le changement de position.

        positionx, positiony, direction, vitesse, masse = collision(positionx, positiony, direction, vitesse, masse, t, beta, alpha, rng) # On calcule les collisions.

        if t == itermax # Si on est à la dernière itération temporelle, on veut tracer la corrélation et trouver la dimension fractale.

            c, rayon = correl(positionx, positiony, masse, taille_Réseau) # On obtient les tableaux pour r et C(r).


            @. model(x, p) = p[1] * x^(p[2]) # Fonction pour fitter la corrélation, C=Ar^(k).
            p0 = [1.0, -1.0] # Guess initial pour les paramètres.

            c_smooth = movavg(c, 10).x # On smooth la corrélation avec un moving average.

            fit = curve_fit(model, rayon[2:100], c_smooth[2:100], p0) # On fit la fonction sur des points qui ne divergent pas trop de la loi de puissance, indices trouvés en regardant le graphique de C(r) en fonction de r.
            param = fit.param # Paramètres fittés.
            rfit = 0:0.1:rayon[size(rayon)[1]]
            yfit = model(rfit, param) # Fonction à tracer pour le fit.
            erreur = stderror(fit) # Erreur sur les paramètres.

            scatter(rayon, c_smooth, markershape=:circle, markercolor=:black, markersize=2, dpi=300, legend=false) # On trace C(r) (smooth) en fonction de r. 
            plot!(rfit, yfit, c=:black) # On trace la fonction fitté.
            display(plot!(xlab=L"$r$", ylab=L"$C(r)$"))
            savefig("/home/laurent/Documents/PHY 3075/Projet final/fit.png")
            print("$(-(-param[2]-2)) +/- $(erreur[2]), $(itermax)\n") # On imprime la dimension fractale, soit -(-k-2) où 2 est la dimension physique de 2.

        end

    end

    scatter(positionx, positiony, markershape=:circle, markercolor=:black, markersize=masse, dpi=300, legend=false) # On crée un graphique qui montre les particules sur le réseau à la fin avec les masses qui correspondent à la taille des points.
    display(scatter!(xlab=L"$x$", ylab=L"$y$"))
    savefig("/home/laurent/Documents/PHY 3075/Projet final/position_finale.png")

    return

end

## Fonction pour créer un gif des particules.

function gif_temporel(nb_particules, taille_Réseau, itermax, rng, beta, alpha) # Fonction qui fait un gif des particules dans le temps.

    anim = @animate for t in 0:itermax # On itère sur le temps pour avoir chaque frame du gif.

        if t == 0

            global masse, direction, vitesse, positionx, positiony = initialisation(nb_particules, taille_Réseau, rng) # On initialise les tableaux.

            scatter(positionx, positiony, markershape=:circle, markercolor=:black, markersize=masse, dpi=300, legend=false, xlab=L"$x$", ylab=L"$y$") # On crée un graphique qui montre les particules sur le réseau initialement.

        else

            positionx, positiony, direction = changement_position(positionx, positiony, direction, vitesse, taille_Réseau, rng) # On calcule le changement de position.

            positionx, positiony, direction, vitesse, masse = collision(positionx, positiony, direction, vitesse, masse, t, beta, alpha, rng) # On calcule les collisions.

            scatter(positionx, positiony, markershape=:circle, markercolor=:black, markersize=masse, dpi=300, legend=false, xlab=L"$x$", ylab=L"$y$") # On crée un graphique qui montre les particules sur le réseau à chaque itération temporelle.      

        end


    end

    gif(anim, "/home/laurent/Documents/PHY 3075/Projet final/animation.gif", fps=15) # On sauve le gif.

    return

end

## Lancer la routine principale.


routine_principale(nb_particules, taille_Réseau, itermax, rng, beta, alpha) # On peut lancer la routine pour divers paramètres alpha et beta et différents temps maximals pour voir l'effet sur la dimension fractale et la distribution des particules.


## Créer un gif 

#=
gif_temporel(nb_particules, taille_Réseau, itermax, rng, beta, alpha) # Appel de la fonction qui créé un gif.
=#

## Lancer la routine principale avec plusieurs valeurs de itermax pour voir comment la dimension évolue dans le temps.

#=
temps = [1, 10, 50, 100, 300, 500, 1000, 3000, 5000, 8000, 10000, 30000, 50000, 100000] # Tableau d'itérations maximales.


for i in temps # On itère sur les différents itermax

    rng = MersenneTwister(seed) # On redéfinit le rng pour avoir la même progressions temporelle à chaque fois.

    routine_principale(nb_particules, taille_Réseau, i, rng, beta, alpha) # On peut lancer la routine pour différents temps maximals pour voir l'effet sur la dimension fractale. # Commenter les traçages de graphique et l'imprimage des collisions dans la routine principale et la fonction collision pour avoir juste les  dimensions qui sont imprimées et ne pas avoir de graphiques.  

end
=#

## Tracer un graphique de la dimension en fonction du temps. 

#=
temps = [1, 10, 50, 100, 300, 500, 1000, 3000, 5000, 8000, 10000, 30000, 50000, 100000] # Tableau d'itérations maximales.

Dim = [1.9922400119558688, 1.9888186940460557, 1.9766458123573671, 1.974689193890168, 1.983197150198585, 1.9855233476875844, 1.9595568034991022, 1.7322321008566621, 1.7029610988386812, 1.686577304143917, 1.6911329093927612, 1.7109215808912175, 1.7094184244249708, 1.708908401366789] # Dim et err_Dim ont été déterminés par les dimensions imprimés par la section précédente.
err_Dim = [0.00391184752266812, 0.004291377346208982, 0.00408262111222954, 0.003299238880132864, 0.005279461404815343, 0.005288591699887982, 0.008650668197550353, 0.016785251646763505, 0.01961820211885456, 0.022735713505804587, 0.025143579834917145, 0.0255924005171586, 0.025988939470755954, 0.025972744471327476]

scatter(temps, Dim, xscale=:log10, yerror=err_Dim, markershape=:circle, markersize=3, legend=false, dpi=300, markercolor=:black, xlab=L"$t$", ylab=L"$D$")
savefig("/home/laurent/Documents/PHY 3075/Projet final/Dimension_.png")
=#

##
