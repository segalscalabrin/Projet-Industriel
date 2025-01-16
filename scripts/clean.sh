#!/bin/bash
# Script pour nettoyer les fichiers générés

# Supprimer le répertoire builds
rm -rf builds

# Supprimer le contenu du répertoire solutions
rm -rf solutions

# Supprimer les logs
rm -rf logs

# Supprimer l'exécutable
rm -f levelset

echo "Nettoyage terminé."
