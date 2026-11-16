<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <title>Simulateur Omega V1</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body { font-family: Arial; background: #1e1e1e; color: white; text-align: center; }
        .container { width: 90%; max-width: 800px; margin: auto; }
        canvas { background: #2a2a2a; margin-top: 20px; border-radius: 10px; }
        button { padding: 10px 20px; font-size: 16px; cursor: pointer; background: #00bcd4; border: none; border-radius: 5px; margin-top: 20px;}
        .stats { margin-top: 15px; font-family: monospace; color: #00ff00; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Simulateur Biophysique Omega V1</h1>
        <p>Lancement d'un protocole de stimulation (Train de Calcium)</p>
        
        <button onclick="lancerSimulation()">Démarrer la simulation</button>

        <div class="stats" id="stats"></div>
        <canvas id="graphique" width="800" height="400"></canvas>
    </div>

    <script type="module">
        import { simulate, DEFAULT_PROTOCOL } from './omegaModel.js';

        let monGraphique = null;

        window.lancerSimulation = function() {
            // On simule 100 secondes
            const resultat = simulate({ 
                protocol: DEFAULT_PROTOCOL, 
                duration: 100 
            });

            // On prépare les données pour les courbes
            const temps = resultat.rows.map(r => r.t);
            const calcium = resultat.rows.map(r => r.ca);
            const camkiiActive = resultat.rows.map(r => r.Kp); // Kp = CaMKII phosphorylée
            const pp1Libre = resultat.rows.map(r => r.PP1);

            // Affichage des preuves de conservation (Invariants)
            document.getElementById('stats').innerHTML = 
                `Preuve de conservation de masse :<br>` +
                `CaM: ${resultat.invariants.cam.toFixed(3)} | ` +
                `CaMKII: ${resultat.invariants.camkii.toFixed(3)} | ` +
                `PP1: ${resultat.invariants.pp1.toFixed(3)} | ` +
                `I1: ${resultat.invariants.i1.toFixed(3)}`;

            // Création du graphique
            if (monGraphique) monGraphique.destroy();
            const ctx = document.getElementById('graphique').getContext('2d');
            
            monGraphique = new Chart(ctx, {
                type: 'line',
                data: {
                    labels: temps,
                    datasets: [
                        { label: 'Calcium (µM)', data: calcium, borderColor: 'cyan', borderWidth: 1, pointRadius: 0 },
                        { label: 'CaMKII Active (Kp)', data: camkiiActive, borderColor: 'yellow', borderWidth: 2, pointRadius: 0 },
                        { label: 'PP1 Libre', data: pp1Libre, borderColor: 'red', borderWidth: 2, pointRadius: 0 }
                    ]
                },
                options: {
                    responsive: false,
                    plugins: { legend: { labels: { color: 'white' } } },
                    scales: {
                        x: { ticks: { color: 'white' }, title: { display: true, text: 'Temps (s)', color: 'white' } },
                        y: { ticks: { color: 'white' }, title: { display: true, text: 'Concentration (µM)', color: 'white' } }
                    }
                }
            });
        };
    </script>
</body>
</html>
