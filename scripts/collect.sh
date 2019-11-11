#!/bin/bash

find . -type f -name profile.csv -print0 |
    while IFS=read -r -d '' f; do
        dir=$(diname "f");
        newName="$dir.csv";
        mv "f" "collection/$newName";
    done

find . -type f -name out.txt -print0 |
    while IFS=read -r -d '' f; do
        dir=$(diname "f");
        newName="$dir.out";
        mv "f" "collection/$newName";
    done

find . -type f -name err.txt -print0 |
    while IFS=read -r -d '' f; do
        dir=$(diname "f");
        newName="$dir.err";
        mv "f" "collection/$newName";
    done
