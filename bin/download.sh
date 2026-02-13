#!/bin/bash

# Download tractography results from remote server
# Usage: ./download.sh <filename>
# Example: ./download.sh tracks_2025-09-22_01-58-34.mat

if [ $# -eq 0 ]; then
    echo "Usage: $0 <filename>"
    echo "Example: $0 tracks_2025-09-22_01-58-34.mat"
    exit 1
fi

FILENAME="$1"
REMOTE_PATH="toosalty0000@olea.uic.yonsei.ac.kr:~/hinec/hinec_runs/${FILENAME}"
LOCAL_PATH="./hinec_runs/${FILENAME}"

# Create local directory if it doesn't exist
mkdir -p ./hinec_runs

echo "Downloading ${FILENAME}..."
scp -r "${REMOTE_PATH}" "${LOCAL_PATH}"

if [ $? -eq 0 ]; then
    echo " Successfully downloaded: ${LOCAL_PATH}"
else
    echo " Download failed"
    exit 1
fi