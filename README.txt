README - Projet 2 : Modélisation de Matrices de Comptage (PFM)
Objectif du projet
L'objectif de ce projet est de modéliser des classes Python permettant de manipuler des matrices de comptage (PFM, Position Frequency Matrix), qui représentent des motifs dans des séquences d'ADN. Ce projet se concentre sur la construction d'une PFM à partir d'un ensemble de séquences d'ADN, en comptabilisant les occurrences de chaque nucléotide à chaque position dans les séquences.

#######################Fonctionnalités principales##########################
Classe Multifasta :

Représente un fichier multifasta contenant des séquences d'ADN.
Fournit une méthode sequences() qui retourne la liste des séquences sous forme de chaînes de caractères.

Classe PFM :

Permet de créer et manipuler une matrice de comptage à partir d'une liste de sites (séquences d'ADN).
Calcul de la longueur de la matrice.
Vérification des colonnes conservées (colonnes où tous les nucléotides sont identiques).
Calcul du consensus (nucléotide le plus fréquent pour chaque position).
Création d'un consensus mou qui indique les possibilités à chaque position en utilisant une notation entre crochets ([ ]).
Possibilité d'ajouter des sites à la PFM et de concaténer plusieurs PFM.
Recherche de motif (search) :

Recherche des positions où le consensus mou de la PFM apparaît dans une séquence d'ADN donnée.
Fichiers inclus
Le dépôt contiendra les fichiers suivants :

main_PFM.py : Contient la fonction search et la partie principale du programme qui prend un fichier multifasta de sites et une séquence ADN en entrée, et qui affiche la PFM, les sites conservés, le consensus mou, ainsi que les positions où la PFM apparaît dans la séquence.
multifasta.py : Contient la classe Multifasta pour gérer les fichiers multifasta.
PFM.py : Contient la classe PFM pour manipuler les matrices de comptage et effectuer les calculs associés.
exemple.fasta : Un fichier d'exemple contenant un ensemble de sites test au format multifasta.
README.txt : Ce fichier, expliquant le projet, les fonctionnalités et comment exécuter le programme.

#####################Étapes pour exécuter le programme ###################

python3 main_PFM.py <file_input.fasta> <sequence.fasta>


#####################Détails d'implémentation#############################

Classe Multifasta
Méthode sequences() : Cette méthode retourne la liste des séquences présentes dans le fichier multifasta.
Classe PFM
Méthode __init__(self, sites) : Initialise une PFM à partir d'une liste de sites.
Méthode is_conserved(self, col_num) : Vérifie si une colonne est conservée.
Méthode consensus(self) : Retourne la séquence de consensus basée sur la PFM.
Méthode weak_consensus(self) : Retourne le consensus mou avec la notation [ ] pour les ambiguïtés.
Méthode append(self, site) : Ajoute un site à la PFM.
Méthode __add__(self, other) : Permet de concaténer deux PFM.
Fonction search
La fonction search recherche les positions dans une séquence où le consensus mou de la PFM apparaît.

#######################Remarques############################################

Le programme suppose que toutes les séquences sont de la même longueur.
Les colonnes conservées dans la PFM sont celles où le même nucléotide apparaît sur toutes les lignes.


##################### resultat ###################

Liste des séquences des sites sous forme de chaînes de caractères :
['ATTGCGGTC', 'ATTGCGGTT', 'TTTGCGGTC', 'ATTGTGGTT', 'ACTGCGGTT', 'TTTGTGGTC', 'CTTGTGGTC', 'ACTGTGGTT', 'CCTGCGGTT', 'ATTGCGGTA', 'TATGCGGTT', 'GCTGCGGTT', 'AATGTGGTA', 'TTTGTGGAC', 'AAAGTGGTT', 'TTCTTGGTT', 'TCAGTGGGT']
---
Matrice de comptage (PFM) :
A 8 3 2 0 0 0 0 1 2
T 6 9 14 1 9 0 0 15 10
C 2 5 1 0 8 0 0 0 5
G 1 0 0 16 0 17 17 1 0
---
Longueur des séquences (nombre de colonnes) :
9
---
La colonne 5 est non conservée.
La colonne 6 est conservée.
La colonne 7 est conservée.
La colonne 8 est non conservée.
---
Colonnes conservées :
[6, 7]
---
Nucléotide le plus fréquent dans la colonne 1 :
A
---
Consensus :
ATTGTGGTT
---
Consensus mou :
[ATCG][ATC][ATC][TG][TC]GG[ATG][ATC]
---
Matrice après ajout d'un nouveau site :
A 9 3 2 1 0 0 1 1 3
T 6 10 15 1 9 0 0 16 10
C 2 5 1 0 8 0 0 0 5
G 1 0 0 16 1 18 17 1 0
---
Colonnes conservées après ajout :
[6]
---
Concaténation de la PFM avec elle-même :
A 9 3 2 1 0 0 1 1 3 9 3 2 1 0 0 1 1 3
T 6 10 15 1 9 0 0 16 10 6 10 15 1 9 0 0 16 10
C 2 5 1 0 8 0 0 0 5 2 5 1 0 8 0 0 0 5
G 1 0 0 16 1 18 17 1 0 1 0 0 16 1 18 17 1 0
---
Séquence ADN lue : GCCAGCGGAAGGAACTTGCTCAT
Positions d'occurrence du consensus mou dans la séquence :
[2]
###################fonction qui ne fonctionne pas ######################

La fonction search est censée rechercher des correspondances entre un consensus mou généré à partir d'une PFM (Position Frequency Matrix) et une séquence donnée. Cependant, la fonction ne semble pas renvoyer les résultats attendus pour certaines positions dans la séquence cible. Plus précisément, elle renvoie seulement [2] au lieu de [1, 12], ce qui suggère un problème dans la comparaison entre la séquence et le consensus mou.

Raisons possibles du dysfonctionnement :
Consensus mou mal formé : Si le consensus mou n'est pas correctement généré, la fonction pourrait ne pas effectuer la bonne comparaison. Par exemple, une position de la séquence pourrait ne pas correspondre à ce qui est attendu en raison de la mauvaise gestion des ambiguïtés (bases possibles à chaque position).

Indexation incorrecte : Il se pourrait qu'il y ait un problème dans la façon dont vous vérifiez les positions. Par exemple, l'indexation dans la séquence pourrait être décalée, ce qui rendrait la correspondance erronée.

Problème de comparaison : La fonction peut ne pas gérer correctement la comparaison entre la séquence et le consensus mou, en particulier pour les bases ambiguës (comme [A,C,T]), ce qui fait que certaines correspondances sont ignorées.