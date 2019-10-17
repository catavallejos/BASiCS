#!/usr/bin/env bash

dv=$(awk '/Version/ {print $2}' DESCRIPTION)
nv=$(grep version NEWS | tail -1 | awk '{print $2}')

if [ "$dv" == "$nv" ]; then
    if [ -z "$(git status --porcelain)" ]; then 
        git tag "v$dv"
    else
        echo "Some changes have not been committed!"
        exit 1
    fi
else
    echo "DESCRIPTION and NEWS versions do not match!"
    exit 1
fi
