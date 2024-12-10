import argparse

class Multifasta:
    def __init__(self, filepath):
        """
        Initialise la classe Multifasta avec le chemin d'accès au fichier.
        Charge les séquences du fichier au format multifasta.
        """
        self.filepath = filepath
        self._sequences = []
        self._load_sequences()
    
    def _load_sequences(self):
        """
        Charge les séquences à partir du fichier multifasta.
        Ignore les noms des sites (lignes commençant par '>').
        """
        try:
            # Vérifier si le fichier existe avant de tenter de l'ouvrir
            with open(self.filepath, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line.startswith(">"):  # Ignore les noms des sites
                        # Vérifier si la ligne est vide (séquence invalide)
                        if line:
                            self._sequences.append(line.upper())  # Convertit en majuscules
                        else:
                            raise ValueError("Erreur : La séquence dans le fichier est vide ou mal formatée.")
        except FileNotFoundError:
            print(f"Erreur : Le fichier {self.filepath} n'existe pas.")
            raise  # Relancer l'exception pour indiquer l'échec de l'opération
        except ValueError as e:
            print(e)
            raise  # Relancer l'exception après l'avoir capturée
        except Exception as e:
            print(f"Erreur inattendue : {e}")
            raise  # Relancer l'exception générique pour tout autre problème

    def sequences(self):
        """
        Retourne la liste des séquences sous forme de chaînes de caractères.
        """
        if not self._sequences:
            raise ValueError("Erreur : Aucune séquence n'a été chargée. Vérifiez le format du fichier multifasta.")
        return self._sequences