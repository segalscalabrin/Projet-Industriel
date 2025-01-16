#!/bin/bash
# Script principal pour gérer le projet

# Vérifier si le mot "clean" est passé en argument
if [[ "$1" == "clean" ]]; then
    ./scripts/clean.sh
    exit 0
fi

# Créer les répertoires nécessaires
./scripts/create_dirs.sh

# Construire le projet
./scripts/build.sh

# Exécuter le projet si aucun autre argument n'est donné
if [[ "$1" == "run" ]]; then
    ./scripts/run.sh
fi
