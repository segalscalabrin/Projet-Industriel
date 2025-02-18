# Projet Industriel
Mise en œuvre de la résolution de Navier-Stockes sur maillage quadtree via la bibliothèque Neos.

## Prérequis
Avant d'exécuter le projet, vous devez télécharger la bibliothèque Neos. Le plus simple est de l'utiliser avec Docker :

### Installation de Docker
```bash
sudo apt-get update
sudo apt-get install -y docker.io
```
Démarrez et activez Docker :
```bash
sudo systemctl start docker
sudo systemctl enable docker
```
Vérifiez l'installation :
```bash
docker --version
```
### Téléchargement et installation de Neos
Clonez le dépôt GitLab de Neos :
```bash
git clone https://gitlab.inria.fr/memphis/neos.git
cd neos
```
Construisez l'image Docker de Neos :
```bash
sudo docker build -t neos -f docker/Dockerfile .
```
### Lancement de l'environnement Docker
Placez-vous dans le dossier du projet (Projet-Volume-Fini), puis exécutez la commande suivante :
```bash
sudo docker run -it -v $(pwd):/builds/workspace neos
```
**Note** : Cette commande devra être exécutée à chaque fois que vous souhaitez travailler dans le dossier `Projet-Industriel`.
## Compilation et exécution

Lancez la compilation et l'exécution avec la commande suivante :
```bash
./main.sh run 
```

### Nettoyage du répertoire
Pour nettoyer le répertoire, utilisez la commande suivante :
```bash
./main.sh clean
```

