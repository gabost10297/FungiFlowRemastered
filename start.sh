#!/usr/bin/env bash
# Start the FungiFlow app at http://127.0.0.1:3838
# Usage: bash start.sh             builds only the images that are missing
#        bash start.sh --rebuild   rebuilds all three images
# Stop it with: docker compose -f docker/compose.yml down

set -euo pipefail
cd "$(dirname "$0")"

build() {
    if [[ "${1}" == --rebuild ]] || ! docker image inspect "$2" > /dev/null 2>&1; then
        echo "Building $2"
        docker build -t "$2" "$3"
    fi
}
build "${1:-}" fungiflow-remastered docker/pipeline
build "${1:-}" fungiflow-ngspeciesid docker/ngspeciesid
build "${1:-}" fungiflow-viz docker/viz

export HOST_PROJECT="$PWD"
export RUN_AS="$(id -u):$(id -g)"
# The socket's group as the containers see it; the app joins it to start pipeline containers.
export DOCKER_GID="$(docker run --rm -v /var/run/docker.sock:/var/run/docker.sock \
    fungiflow-viz stat -c %g /var/run/docker.sock)"

mkdir -p runs
docker compose -f docker/compose.yml up -d --force-recreate
echo "FungiFlow is running at http://127.0.0.1:3838"
