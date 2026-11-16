# Omega C++ Core : Synapse Biophysique à Verrou CaMKII/I-1

Ce projet implémente un modèle de plasticité synaptique bio-inspiré basé sur la compétition enzymatique entre la CaMKII (Kinase) et la PP1 (Phosphatase), régulée par le verrou moléculaire de l'Inhibiteur-1 (I-1).

## 🔬 Découverte : Le mécanisme de l'oubli catastrophique

Grâce aux tests de contrôle, nous avons identifié le mécanisme exact de l'oubli dans ce modèle :
- **Sans apprentissage de B** : La mémoire de A survit à 100% (même sous forte attaque homéostatique).
- **Avec apprentissage de B** : La mémoire de A s'effondre (0% de survie).
- **Intervention causale** : Si l'on désactive l'homéostasie de A *uniquement* pendant l'apprentissage de B, la mémoire de A survit à 100%.

**Conclusion** : L'oubli catastrophique n'est pas causé par le temps qui passe ou une perturbation externe, mais par l'activité post-synaptique générée lors de l'apprentissage d'une nouvelle tâche, qui déclenche la dépression homéostatique de l'ancienne.

## 📊 Résultats Clés
- **Bifurcation et Hystérésis** : Transition de phase nette (LTD vers LTP) à la concentration critique de calcium.
- **Intégration FANN Native** : Le modèle est couplé avec succès à la bibliothèque C `libfann` (test de compilation réussi avec `-lfann`).

## 🛠️ Compilation et Exécution
Nécessite `libfann-dev`.

```bash
# Compilation du test principal
g++ -O3 -std=c++17 test_fann_integration.cpp biophysical_synapse.cpp biophysical_fann_adapter.cpp -o test_fann -lfann -lm
./test_fann

# Compilation du test de causalité
g++ -O3 -std=c++17 test_causalite.cpp biophysical_synapse.cpp biophysical_fann_adapter.cpp -o test_causalite
./test_causalite
