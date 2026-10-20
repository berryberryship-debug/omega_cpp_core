# Omega C++ Core : Synapse Biophysique à Verrou CaMKII/I-1

Ce projet implémente un modèle de plasticité synaptique bio-inspiré basé sur la compétition enzymatique entre la CaMKII (Kinase) et la PP1 (Phosphatase), régulée par le verrou moléculaire de l'Inhibiteur-1 (I-1).

## 📊 Résultats Clés
1. **Bifurcation et Hystérésis** : Transition de phase nette (LTD vers LTP) à la concentration critique de calcium.
2. **Intégration FANN Native** : Le modèle est couplé avec succès à la bibliothèque C `libfann`. Le poids synaptique biophysique pilote directement les connexions du réseau de neurones.
3. **Seuil de Rupture Mesuré** : Le verrou résiste à une attaque homéostatique linéaire jusqu'à c=10 s^-1. Il ne cède qu'à une attaque massive sur la Kinase (f* ≈ 20 s^-1).

## 🛠️ Compilation et Exécution
Nécessite `libfann-dev`.
```bash
g++ -O3 -std=c++17 test_fann_integration.cpp biophysical_synapse.cpp biophysical_fann_adapter.cpp -o test_fann -lfann -lm
./test_fann
