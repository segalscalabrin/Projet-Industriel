#!/bin/bash
# Script pour exécuter l'exécutable

# Définir les chemins relatifs
ROOT_DIR=$(pwd)
LOG_DIR="${ROOT_DIR}/logs"

# Vérifier si l'exécutable existe
if [[ ! -f ./levelset ]]; then
    echo "L'exécutable 'levelset' n'existe pas. Veuillez construire le projet d'abord."
    exit 1
fi

# Exécuter l'exécutable et enregistrer les logs
./navierstokes > "${LOG_DIR}/run_log.txt" 2>&1

# Déplacer les fichiers log générés par la bibliothèque dans le répertoire logs
mv *.log "${LOG_DIR}/" 2>/dev/null

echo "Exécution terminée. Les logs sont dans ${LOG_DIR}/."
