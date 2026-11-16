# OMEGA — Modélisation biophysique de la plasticité synaptique

Projet de recherche : Axiomes Base Labs  
Version : Omega V1  
Date : 16 novembre 2026  
Langage principal : C++17  
Auteur : Xioranne Lessard

## 1. Présentation

OMEGA est un projet de modélisation computationnelle consacré à la dynamique biophysique de la plasticité synaptique.

Son objectif est de construire une architecture modulaire permettant de représenter et d'étudier les interactions moléculaires impliquées dans sa régulation, notamment les mécanismes associés à la calmoduline (CaM), à la kinase CaMKII, à la phosphatase PP1 et à l'inhibiteur-1 (I-1).

Le projet vise à développer un outil reproductible, documenté et progressivement confrontable aux données expérimentales publiées.

Omega V1 représente un réseau réactionnel réduit. Il ne constitue pas, à lui seul, un modèle complet de la plasticité synaptique.

## 2. Architecture d'Omega V1

La version actuelle est organisée autour de modules interconnectés.

| Module | Fonction |
| --- | --- |
| M1 — Calmoduline | Représentation simplifiée des états de liaison du calcium à CaM. |
| M2 — CaMKII | Dynamique réactionnelle réduite comprenant la liaison de CaM, la dimérisation et la phosphorylation. |
| M3 — PP1/I-1 | Phosphorylation de I-1 sous une activité PKA imposée et inhibition réversible de PP1 par I1p. |
| M3.2 — Déphosphorylation | Association de PP1 à KpC et conversion catalytique du complexe en PP1 libre et KC. |
| M1–M2 | Couplage de la calmoduline et de CaMKII avec partage du pool de M4 libre. |
| M2–M3.2 | Couplage de CaMKII avec sa déphosphorylation par PP1. |
| M2–M3–M3.2 | Intégration des mécanismes kinase/phosphatase et de l'inhibition de PP1. |
| M1–M2–M3–M3.2 | Réseau complet avec partage des pools de CaM libre et de PP1 libre. |

### Fichiers principaux

| Fichiers | Rôle |
| --- | --- |
| `calmodulin.hpp`, `calmodulin.cpp` | Module M1. |
| `camkii.hpp`, `camkii.cpp` | Module M2. |
| `pp1_i1.hpp`, `pp1_i1.cpp` | Module M3. |
| `camkii_dephosphorylation.hpp`, `camkii_dephosphorylation.cpp` | Module M3.2. |
| `coupled.hpp`, `coupled.cpp` | Couplage M1–M2. |
| `coupled_dephosphorylation.hpp`, `coupled_dephosphorylation.cpp` | Couplage M2–M3.2. |
| `integrated_network.hpp`, `integrated_network.cpp` | Réseau M2–M3–M3.2. |
| `omega_full.hpp`, `omega_full.cpp` | Réseau complet M1–M2–M3–M3.2. |

Cette architecture permet d'étudier séparément les mécanismes, puis leurs interactions dans un système commun.

## 3. Principes de modélisation

Omega V1 repose sur les principes suivants :

- Décomposition des mécanismes en modules.
- Formalisation explicite des réactions et de leurs dérivées.
- Prise en compte de la stœchiométrie dans les bilans moléculaires.
- Partage explicite des pools entre les modules couplés.
- Vérification de la finitude et de la non-négativité des entrées des fonctions réactionnelles.
- Analyse de la conservation, de la positivité et de la convergence numérique.
- Intégration progressive des sous-systèmes.
- Séparation entre validation mathématique, validation numérique et validation biologique.

Les paramètres utilisés dans les tests de développement sont synthétiques. Ils servent à vérifier le fonctionnement des équations et des implémentations.

Ils ne doivent pas être interprétés comme des paramètres biologiques calibrés.

Les unités des concentrations, du temps, des constantes cinétiques et de l'activité PKA doivent être documentées de manière cohérente.

## 4. Bilans moléculaires

Les quantités totales suivantes sont définies pour le réseau complet :

| Quantité | Expression |
| --- | --- |
| Calmoduline totale | `M0 + M2N + M2C + M4 + KC + 2 D + KpC + PP1KpC` |
| CaMKII totale | `K0 + KC + 2 D + KpC + Kp + PP1KpC` |
| PP1 totale | `PP1 + PP1I1p + PP1KpC` |
| I-1 total | `I1 + I1p + PP1I1p` |

Le dimère `D` compte pour deux unités de kinase et deux unités de calmoduline liée.

Le complexe `PP1KpC` contient une unité de PP1, une unité de kinase et une unité de calmoduline liée.

Le calcium étant une entrée imposée, le réseau ne représente pas un système fermé pour le calcium total.

## 5. Validation numérique

La validation numérique doit porter notamment sur :

- Les invariants de conservation.
- La dynamique des états moléculaires.
- La positivité des concentrations.
- La finitude des flux et des dérivées calculés.
- Le rejet des entrées négatives ou non finies.
- La stabilité numérique.
- La convergence temporelle.
- Le comportement du réseau intégré.

Les fonctions de calcul des dérivées définissent le système réactionnel. Elles ne réalisent pas, à elles seules, l'intégration temporelle.

La fonction `advance_calmodulin` utilise Euler explicite avec un sous-découpage déterminé par les taux de sortie des états. Ce mécanisme vise à préserver la positivité du module de calmoduline isolé ; il ne fournit pas d'estimation de l'erreur locale ni de garantie générale pour le réseau complet.

Les résultats des tests doivent être accompagnés des paramètres, conditions initiales, commandes, versions du code et sorties permettant leur reproduction.

Les affirmations de conservation, de stabilité et de convergence doivent être rattachées aux configurations effectivement testées. Ce README ne remplace pas un rapport d'exécution des tests.

Une cohérence numérique interne ne constitue pas, à elle seule, une validation expérimentale ou biologique du modèle.

### Convention de mise à jour de M4

Dans `CoupledDerivative`, `M4_total` contient la dérivée totale de M4 libre, incluant les contributions de M1 et de M2.

Lors de l'intégration, ce champ doit être utilisé pour mettre à jour `CoupledState::cam.M4`, sans ajouter une seconde fois la contribution de `cam.M4`.

Dans `OmegaFullDerivative`, la contribution de CaMKII est déjà ajoutée directement à `cam.M4`.

## 6. Limites actuelles

Plusieurs dimensions restent à développer :

1. Constitution d'un registre de paramètres biologiques avec références bibliographiques.
2. Vérification des unités et des conventions cinétiques.
3. Vérification de la compatibilité entre les constantes publiées et les équations du modèle.
4. Calibration des paramètres à partir de données expérimentales.
5. Analyses de sensibilité et d'identifiabilité.
6. Comparaison avec des modèles publiés et des observations expérimentales.
7. Documentation des hypothèses et des domaines de validité.

Le calcium et l'activité PKA sont actuellement représentés par des entrées imposées dans le réseau complet. Leur dynamique propre n'est pas modélisée.

Le module de calmoduline ne représente pas explicitement les états à un ou trois ions calcium liés.

Le module PP1/I-1 ne représente pas la déphosphorylation de I1p.

Le module de déphosphorylation de CaMKII agit sur KpC et ne représente pas la déphosphorylation directe de Kp.

Les résultats doivent être interprétés dans les limites des réactions effectivement implémentées.

## 7. Travaux antérieurs

Le dépôt comprend également des travaux biophysiques antérieurs portant notamment sur :

- La dynamique synaptique couplée à CaMKII et PP1.
- La régulation homéostatique.
- L'adaptation à des entrées pré- et postsynaptiques.
- L'intégration avec la bibliothèque FANN.
- Les mécanismes de bruit, de stabilité et de rupture.

Ces travaux constituent une base historique distincte d'Omega V1. Leur relation avec l'architecture actuelle doit être examinée avant toute intégration.

Les résultats annoncés dans les documents historiques doivent être reproduits à partir de leurs programmes, paramètres et sorties avant d'être considérés comme établis.

## 8. Orientation scientifique

L'objectif à terme est de produire un outil computationnel utilisable par des chercheurs pour :

- Explorer des hypothèses mécanistiques.
- Tester des scénarios de plasticité synaptique.
- Étudier les interactions entre voies enzymatiques.
- Identifier les paramètres influençant les dynamiques simulées.
- Concevoir des expériences permettant de départager des hypothèses concurrentes.

OMEGA est conçu comme un outil d'exploration et de formulation d'hypothèses, et non comme un substitut aux validations expérimentales.

## 9. Environnement technique

- Langage : C++17.
- Compilateurs visés : GCC et Clang.
- Architecture modulaire.
- Gestion des versions : Git.

Les commandes de compilation et d'exécution des tests doivent être documentées selon l'arborescence réelle du dépôt.

Les chemins des directives `#include` doivent correspondre à cette arborescence.

## 10. Statut du projet

Omega V1 : architecture intégrée et validation numérique interne en cours de consolidation.

Les priorités sont :

1. Vérifier la compilation à partir de l'arborescence réelle du dépôt.
2. Documenter les tests et rendre leurs résultats reproductibles.
3. Vérifier les bilans de conservation, la positivité et la convergence du réseau complet.
4. Documenter les unités et les hypothèses réactionnelles.
5. Constituer un registre de paramètres biologiques sourcés.
6. Confronter progressivement le modèle aux données publiées.

---

Axiomes Base Labs  
Architecture des systèmes vivants et noétique appliquée.

Auteur : Xioranne Lessard
