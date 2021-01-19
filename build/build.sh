#!/bin/bash

if [[ -z "${CG_VERSION}" ]]; then
    echo "CG_VERSION undefined"
    exit 1
else
    echo "Building version: ${CG_VERSION}"
fi

gcloud builds submit \
    --config build/cloudbuild_gisaid.yaml \
    --substitutions=TAG_NAME="${CG_VERSION}" . && \
gcloud builds submit \
    --config build/cloudbuild_gisaid_private.yaml \
    --substitutions=TAG_NAME="${CG_VERSION}" . && \
gcloud builds submit \
    --config build/cloudbuild_genbank.yaml \
    --substitutions=TAG_NAME="${CG_VERSION}" .
