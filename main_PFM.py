import argparse
from Multifasta_parser import Multifasta
from PFM_parser import PFM

def search(pfm, sequence):
    """
    Recherche les positions où le consensus mou de la PFM peut apparaître dans une séquence d'ADN.
    
    Arguments :
    pfm -- un objet de la classe PFM, représentant la matrice de fréquence des motifs.
    sequence -- la séquence d'ADN dans laquelle on cherche des motifs.

    Retourne :
    positions -- une liste contenant les indices des positions où le consensus mou est trouvé dans la séquence.
    """
    # Création de la PFM à partir des séquences passées
    m = PFM(pfm.sequences)
    
    # Obtenir le consensus mou
    weak_consensus = m.weak_consensus()
    positions = []
    motif_length = len(pfm)  # Longueur du motif à rechercher dans la séquence d'ADN.

    # Recherche dans la séquence d'ADN pour toutes les positions possibles
    for i in range(len(sequence) - motif_length + 1):
        window = sequence[i:i + motif_length]  # Fenêtre de la séquence où l'on cherche le motif
        match = True
        j = 0  # Index pour le consensus mou
        k = 0  # Index pour la fenêtre de la séquence

        # Vérification de la correspondance avec le consensus mou
        while j < len(weak_consensus) and k < motif_length:
            if weak_consensus[j] == '[':
                # Gestion des ambiguïtés dans le consensus mou (bases multiples possibles)
                end_bracket = weak_consensus.index(']', j)
                possible_bases = weak_consensus[j + 1:end_bracket]
                if sequence[i + k] not in possible_bases:
                    match = False
                    break
                j = end_bracket + 1
            else:
                # Correspondance de base unique
                if sequence[i + k] != weak_consensus[j]:
                    match = False
                    break
                j += 1
            k += 1

        if match:
            positions.append(i + 1)  # Ajouter la position trouvée (indexé à partir de 1)

    return positions

if __name__ == "__main__":
    # Analyse des arguments de la ligne de commande
    parser = argparse.ArgumentParser(description="Multifasta and PFM processing")
    parser.add_argument("input_file", help="Path to the multifasta input file")
    parser.add_argument("sequence_file", help="Path to the fasta file containing the DNA sequence")
    args = parser.parse_args()

    # Création d'un objet Multifasta pour lire le fichier multifasta
    f = Multifasta(args.input_file)
    sequences = f.sequences()  # Obtenir les séquences du fichier multifasta

    # Afficher les séquences lues
    print("Liste des séquences des sites sous forme de chaînes de caractères :")
    print(sequences)
    print("---")

    # Création d'un objet PFM à partir des séquences
    m = PFM(sequences)

    # Lecture de la séquence ADN à partir du fichier fasta
    with open(args.sequence_file, "r") as fasta_file:
        fasta_lines = fasta_file.readlines()
    dna_sequence = "".join([line.strip() for line in fasta_lines if not line.startswith(">")])

    # Afficher la matrice de comptage (PFM)
    print("Matrice de comptage (PFM) :")
    print(m)
    print("---")

    # Afficher la longueur de la PFM
    print("Longueur des séquences (nombre de colonnes) :")
    print(len(m))
    print("---")

    # Vérifier la conservation pour des colonnes spécifiques
    for column in [5, 6, 7, 8]:
        is_conserved = m.is_conserved(column)
        status = "conservée" if is_conserved else "non conservée"
        print(f"La colonne {column} est {status}.")
    print("---")

    # Afficher les colonnes conservées
    print("Colonnes conservées :")
    print(m.conserved_columns())
    print("---")

    # Afficher le nucléotide le plus fréquent dans la première colonne
    print("Nucléotide le plus fréquent dans la colonne 1 :")
    print(m.most_frequent_base(1))
    print("---")

    # Afficher la séquence consensus
    print("Consensus :")
    print(m.consensus())
    print("---")

    # Afficher le consensus mou
    print("Consensus mou :")
    print(m.weak_consensus())
    print("---")

    # Ajouter un nouveau site à la PFM et afficher la mise à jour
    m.append("ATTAGGATA")
    print("Matrice après ajout d'un nouveau site :")
    print(m)
    print("---")

    # Afficher les colonnes conservées après mise à jour
    print("Colonnes conservées après ajout :")
    print(m.conserved_columns())
    print("---")

    # Concaténer la PFM avec elle-même et afficher le résultat
    print("Concaténation de la PFM avec elle-même :")
    print(m + m)
    print("---")

    # Afficher la séquence ADN lue
    print("Séquence ADN lue :", dna_sequence)

    # Recherche des positions d'occurrence du consensus mou dans la séquence ADN
    positions = search(m, dna_sequence)
    print("Positions d'occurrence du consensus mou dans la séquence :")
    print(positions if positions else "Aucune occurrence trouvée.")
