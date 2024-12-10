class PFM:
    def __init__(self, sequences):
        """
        Initialise la PFM à partir d'une liste de séquences.
        :param sequences: Liste des séquences d'ADN à analyser.
        """
        self.sequences = sequences  # Liste des séquences
        self.matrix = self._build_matrix()  # Matrice de comptage construite à partir des séquences

    def _build_matrix(self):
        """
        Construit la matrice de comptage à partir des séquences.
        La matrice compte la fréquence des bases A, T, C et G à chaque position (colonne) des séquences.
        :return: Dictionnaire représentant la matrice de fréquence des bases.
        """
        matrix = {base: [0] * len(self.sequences[0]) for base in "ATCG"}  # Initialisation de la matrice de comptage
        for seq in self.sequences:  # Pour chaque séquence
            for i, base in enumerate(seq):  # Pour chaque base dans la séquence
                if base in matrix:  # Si la base est valide (A, T, C, G)
                    matrix[base][i] += 1  # Incrémenter le comptage pour la base à la position i
        return matrix  # Retourner la matrice construite

    def __str__(self):
        """
        Représentation textuelle de la PFM.
        Affiche chaque base (A, T, C, G) et son comptage à chaque position de la séquence.
        :return: Représentation en chaîne de caractères de la matrice.
        """
        return "\n".join(
            f"{base} {' '.join(map(str, counts))}" for base, counts in self.matrix.items()
        )

    def __len__(self):
        """
        Retourne la longueur des séquences (nombre de colonnes dans la matrice).
        :return: Longueur des séquences.
        """
        return len(self.sequences[0])

    def __add__(self, other):
        """
        Combine deux matrices PFM en ajoutant les mêmes colonnes pour chaque PFM.
        La matrice résultante aura toutes les séquences des deux PFM.
        :param other: Une autre instance de la classe PFM à combiner.
        :return: Une nouvelle instance de PFM avec les séquences combinées.
        """
        if len(self) != len(other):
            raise ValueError("Les deux PFM doivent avoir le même nombre de colonnes.")  # Vérification des dimensions

        # Dupliquer les séquences et ajouter les colonnes
        combined_sequences = self.sequences + other.sequences
        combined_matrix = {base: counts + self.matrix[base] for base, counts in self.matrix.items()}

        # Retourner une nouvelle instance de PFM avec la matrice combinée
        new_pfm = PFM(combined_sequences)
        new_pfm.matrix = combined_matrix  # Mise à jour de la matrice combinée
        return new_pfm

    def append(self, new_site):
        """
        Ajoute un nouveau site à la PFM (vérifie la longueur de la séquence).
        :param new_site: La nouvelle séquence (site) à ajouter.
        """
        if len(new_site) != len(self):
            raise ValueError("Le nouveau site doit avoir la même longueur que les séquences existantes.")  # Vérification de la longueur
        self.sequences.append(new_site)  # Ajouter la nouvelle séquence
        self.matrix = self._build_matrix()  # Recalculer la matrice après l'ajout

    def is_conserved(self, column):
        """
        Vérifie si une colonne est conservée, c'est-à-dire si toutes les séquences ont la même base à cette position.
        :param column: Numéro de la colonne à vérifier (1-indexé).
        :return: True si la colonne est conservée, False sinon.
        """
        counts = [self.matrix[base][column - 1] for base in "ATCG"]  # Comptage des bases dans la colonne donnée
        return max(counts) == sum(counts)  # Si une seule base est présente, la colonne est conservée

    def conserved_columns(self):
        """
        Retourne la liste des colonnes conservées, c'est-à-dire les positions où toutes les séquences ont la même base.
        :return: Liste des indices des colonnes conservées.
        """
        return [col + 1 for col in range(len(self)) if self.is_conserved(col + 1)]  # Colonnes conservées (1-indexées)

    def most_frequent_base(self, column):
        """
        Retourne le nucléotide le plus fréquent dans une colonne donnée.
        :param column: Numéro de la colonne (1-indexé).
        :return: Le nucléotide le plus fréquent à la position donnée.
        """
        column_index = column - 1  # Convertir en index 0-indexé
        max_count = 0
        most_frequent = None
        # Parcourir les bases (A, T, C, G) et identifier la base la plus fréquente
        for base in "ATCG":
            if self.matrix[base][column_index] > max_count:
                max_count = self.matrix[base][column_index]
                most_frequent = base
        return most_frequent  # Retourner la base la plus fréquente

    def consensus(self):
        """
        Retourne la séquence consensus en prenant la base la plus fréquente à chaque position.
        :return: La séquence consensus sous forme de chaîne de caractères.
        """
        return "".join(self.most_frequent_base(col + 1) for col in range(len(self)))  # Séquence consensus

    def weak_consensus(self):
        """
        Retourne le consensus mou, où les positions ambiguës sont représentées entre crochets.
        :return: Le consensus mou sous forme de chaîne de caractères.
        """
        consensus = []
        for col in range(len(self)):
            bases = [base for base in "ATCG" if self.matrix[base][col] > 0]  # Bases présentes dans la colonne
            if len(bases) == 1:
                consensus.append(bases[0])  # Si une seule base est présente, l'ajouter
            else:
                consensus.append(f"[{''.join(bases)}]")  # Si plusieurs bases sont présentes, les mettre entre crochets
        return "".join(consensus)  # Retourner le consensus mou
