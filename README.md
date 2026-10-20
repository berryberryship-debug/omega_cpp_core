# Omega C++ Core : Synapse Biophysique à Verrou CaMKII/I-1

Ce projet implémente un modèle de plasticité synaptique bio-inspiré basé sur la compétition enzymatique entre la CaMKII (Kinase) et la PP1 (Phosphatase), régulée par le verrou moléculaire de l'Inhibiteur-1 (I-1).

## État de validation — réconciliation locale et branche publique

La branche `main` n'est pas encore réconciliée avec la version locale exécutée sous Termux. Les modifications de l'interface de l'adaptateur et le protocole de test avec graines distinctes restent à intégrer et à compiler dans une copie propre du dépôt.

Un balayage exploratoire de la force homéostatique c, de 5 à 30, est prévu sur 100 essais par valeur. Les résultats, les conditions d'exécution et les intervalles de confiance seront publiés après exécution effective et vérification.

Aucun seuil critique ni aucune robustesse générale ne sont considérés comme établis avant cette validation.

## 📊 Résultats Précédents (à titre indicatif)
- **Bifurcation et Hystérésis** : Transition de phase nette (LTD vers LTP) à la concentration critique de calcium.
- **Intégration FANN Native** : Le modèle est couplé avec succès à la bibliothèque C `libfann` (test de compilation réussi).

## 🛠️ Compilation et Exécution
Nécessite `libfann-dev`.
```bash
g++ -O3 -std=c++17 test_fann_integration.cpp biophysical_synapse.cpp biophysical_fann_adapter.cpp -o test_fann -lfann -lm
./test_fann
