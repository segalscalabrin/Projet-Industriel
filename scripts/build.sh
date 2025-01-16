#!/bin/bash
# Script pour construire le projet

# Définir les chemins relatifs
ROOT_DIR=$(pwd)
BUILD_DIR="${ROOT_DIR}/builds"
SRC_DIR="${ROOT_DIR}/src"
LOG_DIR="${ROOT_DIR}/logs"

# Crée les répertoires nécessaires
mkdir -p "${BUILD_DIR}"
mkdir -p "${LOG_DIR}"

# Crée dynamiquement le fichier CMakeLists.txt dans le répertoire builds
cat <<EOL > "${BUILD_DIR}/CMakeLists.txt"
cmake_minimum_required(VERSION 3.1)

project(navierstokes)

# Définir le standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Définir le chemin vers le répertoire source
set(SRC_DIR \${CMAKE_SOURCE_DIR}/../src)

# Récupérer tous les fichiers .cpp du répertoire source
file(GLOB_RECURSE SRCS
    \${SRC_DIR}/*.cpp
    \${SRC_DIR}/**/*.cpp
)

# Vérifier si des fichiers .cpp ont été trouvés
if (NOT SRCS)
    message(FATAL_ERROR "Aucun fichier .cpp trouvé dans le répertoire \${SRC_DIR}")
endif()

# Inclure le répertoire contenant les headers
include_directories(
    \${SRC_DIR}
    \${SRC_DIR}/**
)

# Trouver la bibliothèque NEOS
find_package(NEOS REQUIRED)

# Ajouter l'exécutable avec tous les fichiers trouvés
add_executable(navierstokes \${SRCS})

# Lier la bibliothèque NEOS
target_link_libraries(navierstokes neos::neos)

# Définir le répertoire de sortie de l'exécutable
set_target_properties(navierstokes PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY \${CMAKE_SOURCE_DIR}/..
)

EOL

# Aller dans le répertoire builds
cd "${BUILD_DIR}" || exit

# Exécuter cmake et enregistrer les logs
cmake . > "${LOG_DIR}/cmake_log.txt" 2>&1
if [[ $? -ne 0 ]]; then
    echo "CMake a échoué. Consultez les logs dans ${LOG_DIR}/cmake_log.txt."
    exit 1
fi

# Exécuter make et enregistrer les logs
make > "${LOG_DIR}/make_log.txt" 2>&1
if [[ $? -ne 0 ]]; then
    echo "Make a échoué. Consultez les logs dans ${LOG_DIR}/make_log.txt."
    exit 1
fi

# Revenir au répertoire principal
cd "${ROOT_DIR}" || exit

echo "Build terminé. Les logs sont dans ${LOG_DIR}."
