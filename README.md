OMEGA — Modélisation biophysique de la plasticité synaptique

Projet de recherche : Axiomes Base Labs
Version : Omega V1
Date : 16 novembre 2026
Langage principal : C++17

1. Présentation

OMEGA est un projet de modélisation computationnelle consacré à la dynamique biophysique de la plasticité synaptique.

Son objectif est de construire une architecture modulaire permettant de représenter et d'étudier les interactions moléculaires impliquées dans la régulation de la plasticité synaptique, notamment les mécanismes associés à la calmoduline (CaM), à la CaMKII, à la phosphatase PP1 et à l'inhibiteur-1 (I-1).

Le projet vise à développer un outil de modélisation reproductible, documenté et progressivement confrontable aux données expérimentales publiées.

2. Architecture d'Omega V1

La version actuelle est organisée autour de modules interconnectés.

Module| Fonction
M1 — Calmoduline| Représentation simplifiée des états de liaison du calcium à CaM
M2 — CaMKII| Dynamique réactionnelle réduite de la kinase
M3 — PP1/I-1| Régulation de PP1 par phosphorylation et inhibition via I-1
M3.2 — Déphosphorylation| Interaction de PP1 avec les états phosphorylés de CaMKII
M2–M3–M3.2| Intégration des mécanismes kinase/phosphatase
M1–M2–M3–M3.2| Réseau intégré avec partage des pools de CaM et de PP1

L'architecture permet d'étudier séparément les mécanismes, puis leurs interactions dans un système commun.

3. Principes de modélisation

Omega V1 repose sur les principes suivants :

- Décomposition des mécanismes en modules indépendants.
- Formalisation des réactions et de leurs dérivées.
- Conservation des quantités moléculaires définies par le modèle.
- Contrôle de la positivité et de la finitude des états.
- Analyse de la convergence numérique.
- Intégration progressive des sous-systèmes.
- Séparation entre validation mathématique, validation numérique et validation biologique.

Les paramètres actuellement utilisés dans les tests de développement sont synthétiques. Ils servent à vérifier le fonctionnement des équations et des implémentations.

Ils ne doivent pas être interprétés comme des paramètres biologiques calibrés.

4. Validation numérique

La suite de développement comprend des tests portant notamment sur :

- Les invariants de conservation.
- La dynamique des états moléculaires.
- La positivité des concentrations.
- La stabilité numérique.
- La convergence temporelle.
- Le comportement du réseau intégré.

Les tests réalisés sur les modules et le réseau intégré ont permis de vérifier plusieurs propriétés numériques, notamment la conservation des pools représentés et la convergence de l'intégration d'Euler dans les configurations testées.

Ces résultats établissent un niveau de cohérence numérique interne. Ils ne constituent pas, à eux seuls, une validation expérimentale ou biologique du modèle.

5. Limites actuelles

Plusieurs dimensions restent à développer :

1. Constitution d'un registre de paramètres biologiques avec références bibliographiques.
2. Vérification de la compatibilité entre les constantes publiées et les équations du modèle.
3. Calibration des paramètres à partir de données expérimentales.
4. Analyses de sensibilité et d'identifiabilité.
5. Comparaison avec des modèles publiés et des observations expérimentales.
6. Documentation des hypothèses et des domaines de validité.

Le calcium et l'activité PKA sont actuellement représentés par des entrées imposées dans le réseau intégré. Leur dynamique complète n'est pas encore modélisée.

6. Travaux antérieurs

Le dépôt comprend également des travaux biophysiques antérieurs portant notamment sur :

- La dynamique synaptique couplée à CaMKII et PP1.
- La régulation homéostatique.
- L'adaptation à des entrées pré- et postsynaptiques.
- L'intégration avec la bibliothèque FANN.
- Les mécanismes de bruit, de stabilité et de rupture.

Ces travaux constituent une base historique distincte d'Omega V1. Leur relation avec l'architecture actuelle doit être examinée avant toute intégration.

Les résultats annoncés dans les documents historiques doivent être reproduits à partir de leurs programmes, paramètres et sorties avant d'être considérés comme établis.

7. Orientation scientifique

L'objectif à terme est de produire un outil computationnel utilisable par des chercheurs pour :

- Explorer des hypothèses mécanistiques.
- Tester des scénarios de plasticité synaptique.
- Étudier les interactions entre voies enzymatiques.
- Identifier les paramètres influençant les dynamiques simulées.
- Concevoir des expériences permettant de départager des hypothèses concurrentes.

OMEGA est conçu comme un outil d'exploration et de formulation d'hypothèses, et non comme un substitut aux validations expérimentales.

8. Environnement technique

- C++17
- Compilation avec Clang ou GCC
- Tests exécutables indépendamment
- Architecture modulaire
- Gestion des versions avec Git

9. Statut du projet

Omega V1 : architecture intégrée et validation numérique interne en cours de consolidation.

La prochaine étape prioritaire est la constitution d'un registre de paramètres biologiques documentés et la confrontation du modèle aux données publiées.

---

Axiomes Base Labs
Architecture des systèmes vivants et noétique appliquée.

Auteur : Xioranne Lessard
