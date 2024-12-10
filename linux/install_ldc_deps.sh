#!/bin/bash
set -eEo pipefail
shopt -s failglob
trap 'echo "${BASH_SOURCE[0]}{${FUNCNAME[0]}}:${LINENO}: Error: command \`${BASH_COMMAND}\` failed with exit code $?"' ERR
cd "$(dirname "$0")"

echo "Installing dependencies..."
sudo apt update
sudo apt install -y \
	python3-pil.imagetk

echo "Application dependencies (except for OpenCV) has been successfully installed."
