# OMEGA — Modélisation biophysique de la plasticité synaptique

**Projet de recherche :** Axiomes Base Labs  
**Version :** Omega V1  
**Date :** 16 novembre 2026  
**Langage principal :** C++17  
**Auteure :** Xioranne Lessard

## 1. Présentation

OMEGA est un projet de modélisation computationnelle consacré à la dynamique biophysique de la plasticité synaptique.

Son objectif est de construire une architecture modulaire permettant de représenter et d’étudier les interactions moléculaires impliquées dans la régulation de la plasticité synaptique, notamment les mécanismes associés à :

- la calmoduline (CaM) ;
- la CaMKII ;
- la phosphatase PP1 ;
- l’inhibiteur-1 (I-1).

Le projet vise à développer un outil de modélisation reproductible, documenté et progressivement confrontable aux données expérimentales publiées.

## 2. Architecture d’Omega V1

La version actuelle est organisée autour de modules interconnectés.

| Module | Fonction |
|---|---|
| **M1 — Calmoduline** | Représentation simplifiée des états de liaison du calcium à la calmoduline. |
| **M2 — CaMKII** | Dynamique réactionnelle réduite de la kinase. |
| **M3 — PP1/I-1** | Phosphorylation de I-1 et inhibition réversible de PP1 par I1p. |
| **M3.2 — Déphosphorylation** | Déphosphorylation de l’état KpC de CaMKII par PP1. |
| **M2–M3–M3.2** | Intégration des mécanismes kinase/phosphatase. |
| **M1–M2–M3–M3.2** | Réseau complet avec partage des pools de CaM et de PP1. |

L’architecture permet d’étudier séparément les mécanismes, puis leurs interactions dans un système commun.

## 3. Principes de modélisation

Omega V1 repose sur les principes suivants :

- Décomposition des mécanismes en modules indépendants.
- Formalisation explicite des réactions et de leurs dérivées temporelles.
- Conservation des quantités moléculaires définies par le modèle.
- Contrôle de la positivité et de la finitude des états.
- Analyse de la convergence numérique.
- Intégration progressive des sous-systèmes.
- Séparation entre validation mathématique, validation numérique et validation biologique.

Les paramètres actuellement utilisés dans les tests de développement sont synthétiques. Ils servent à vérifier le fonctionnement des équations et des implémentations.

Ils ne doivent pas être interprétés comme des paramètres biologiques calibrés.

## 4. Validation numérique

Les vérifications numériques portent notamment sur :

- Les invariants de conservation.
- La dynamique des états moléculaires.
- La positivité des concentrations.
- La stabilité numérique.
- La convergence temporelle.
- Le comportement du réseau intégré.

Les résultats de conservation et de convergence doivent être documentés pour chaque configuration testée, avec :

- les paramètres utilisés ;
- les conditions initiales ;
- les pas de temps ;
- les sorties numériques correspondantes.

Ces vérifications visent à établir la cohérence numérique interne du modèle. Elles ne constituent pas, à elles seules, une validation expérimentale ou biologique.

## 5. Limites actuelles

Plusieurs dimensions restent à développer :

1. Constitution d’un registre de paramètres biologiques avec références bibliographiques.
2. Vérification de la compatibilité entre les constantes publiées et les équations du modèle.
3. Calibration des paramètres à partir de données expérimentales.
4. Analyses de sensibilité et d’identifiabilité.
5. Comparaison avec des modèles publiés et des observations expérimentales.
6. Documentation des hypothèses et des domaines de validité.

Le calcium et l’activité de la PKA sont actuellement représentés par des entrées imposées dans le réseau complet. Leur dynamique propre n’est pas modélisée.

## 6. Travaux antérieurs

Le dépôt comprend également des travaux biophysiques antérieurs portant notamment sur :

- La dynamique synaptique couplée à CaMKII et PP1.
- La régulation homéostatique.
- L’adaptation à des entrées pré- et postsynaptiques.
- L’intégration avec la bibliothèque FANN.
- Les mécanismes de bruit, de stabilité et de rupture.

Ces travaux constituent une base historique distincte d’Omega V1. Leur relation avec l’architecture actuelle doit être examinée avant toute intégration.

Les résultats annoncés dans les documents historiques doivent être reproduits à partir de leurs programmes, paramètres et sorties avant d’être considérés comme établis.

## 7. Orientation scientifique

L’objectif à terme est de produire un outil computationnel utilisable par des chercheurs pour :

- Explorer des hypothèses mécanistiques.
- Tester des scénarios de plasticité synaptique.
- Étudier les interactions entre voies enzymatiques.
- Identifier les paramètres influençant les dynamiques simulées.
- Concevoir des expériences permettant de départager des hypothèses concurrentes.

OMEGA est conçu comme un outil d’exploration et de formulation d’hypothèses, et non comme un substitut aux validations expérimentales.

## 8. Environnement technique

- C++17.
- Compilateurs visés : Clang et GCC.
- Tests exécutables indépendamment.
- Architecture modulaire.
- Gestion des versions avec Git.

## 9. Statut du projet

**Omega V1 : architecture intégrée et validation numérique interne en cours de consolidation.**

La prochaine étape prioritaire est la constitution d’un registre de paramètres biologiques documentés et la confrontation du modèle aux données publiées.

---

**Axiomes Base Labs**  
Architecture des systèmes vivants et noétique appliquée.

**Auteure :** Xioranne Lessard
