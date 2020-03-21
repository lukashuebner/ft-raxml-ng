#!/usr/bin/env bash

COMMAND=$1
FILES=("work_by_rank.csv" "*.err" "*.out")

if [ "$COMMAND" == "" ]; then
    echo "Usage: $0 <run,clear,collect>"
    exit 1
fi

are_you_sure () {
    read -p "Are you sure you want to $1?" -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 2
    fi
}

if [ "$COMMAND" == "run" ]; then
    for batch_file in *.sbatch; do
        short_code=$(echo "$batch_file" | grep -o "[AD][0-9][a]*")
        sbatch -J "$short_code" "$batch_file"
    done
elif [ "$COMMAND" == "clear" ]; then
    are_you_sure "clear all files"
    for pattern in "${FILES[@]}"; do
        find . -type f -name "$pattern" -delete
    done
elif [ "$COMMAND" == "collect" ]; then
    for pattern in "${FILES[@]}"; do
        find . -type f -name "$pattern" -print0 |
            while IFS= read -r -d '' f; do
                dir=$(dirname "$f")
                newName="$dir.workByRank.csv"
                cp "$f" "$TARGET_DIR/$newName"
             done
    done
fi

