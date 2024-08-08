# Essai
Ce dépôt contient tous les documents et fichiers de code R utilisés pour la réalisation de mon essai.
# Résumé
Cet essai vise à comprendre et à quantifier le risque d'incendie, en se concentrant particulièrement sur l'assurance construction pour les chantiers en bois massif. Pour ce faire, on utilise diverses bases de données de sinistres incendie, notamment celles de la ville de Toronto et la *National Fire Information Database*. À l'aide de ces bases de données, on crée plusieurs modèles probabilistes en utilisant différentes familles paramétriques. On examine l'impact des types de construction sur les pertes totales en sinistre incendie en calculant la part allouée à chaque type en se basant sur la théorie du partage des risques (*risk sharing*).

# Description des dossiers
## Data
Ce dossier contient la base de données des périls d'incendie de la ville de Toronto. En raison de sa taille et des restrictions liées aux droits, la National Fire Information Database ne peut pas être incluse dans ce dépôt. Les autres bases de données fournies sont des bases de données bien connues qui ont permis de confirmer les résultats.

## Articles
Ce dossier regroupe l'ensemble des articles consultés au cours de ma maîtrise pour la rédaction de mon essai.

## Code
Ce dossier contient tout le code développé pour mon essai.

#### POT_Exemple
Ce dossier contient le fichier RMarkdown illustrant la méthode *Peak-over-Threshold* (POT). Les bases de données nécessaires se trouvent dans le dossier Data.

#### Modelisation
Ce dossier contient tous les codes relatifs à la modélisation du péril incendie pour différentes bases de données.
- `Fonction.R` : contiens toutes les fonctions nécessaires pour la modélisation.
- `Documentation.R` : contiens la documentation des fonctions qui se trouve dans le fichier Fonction.R.
- `Toronto.R` : contient le code pour la modélisation du péril incendie pour la base de données de la ville de Toronto, qui représente une portion importante de l'essai.
- `NFID.R` : contient le code pour la modélisation du péril incendie pour la NFID, qui représente une portion importante de l'essai.
- `MixedErlang.R` : contient le code qui teste la loi de mélange d'Erlang pour la portion sous le seuil.

Les autres fichiers de ce dossier sont similaires aux fichiers `Toronto.R` et `NFID.R`, mais utilisent d'autres bases de données.

#### RiskSharing
Ce dossier contient le code nécessaire pour la section sur le partage de risque de l'essai.
- `RiskSharing_LN.R` : contient le code pour le partage de risque en utilisant la loi LN-PaG, tel que utilisé dans l'essai.
- `RiskSharing_GB2.R` : contiens le code pour le partage de risque en utilisant la loi GB2-PaG.

