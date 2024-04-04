#!/bin/bash

if [[ -z "${CG_VERSION}" ]]; then
    echo "CG_VERSION undefined"
    exit 1
else
    echo "Deploying version: ${CG_VERSION}"
fi

if [[ -z "${PROJECT_ID}" ]]; then
    echo "PROJECT_ID undefined"
    exit 1
else
    echo "Project: ${PROJECT_ID}"
fi

gcloud run deploy cg         --image "gcr.io/${PROJECT_ID}/cg:${CG_VERSION}" && \
gcloud run deploy cg-private --image "gcr.io/${PROJECT_ID}/cg-private:${CG_VERSION}" && \
gcloud run deploy cg-genbank --image "gcr.io/${PROJECT_ID}/cg-genbank:${CG_VERSION}" && \
gcloud run deploy cg-alpha   --image "gcr.io/${PROJECT_ID}/cg-alpha:${CG_VERSION}"

# RSV

gcloud run deploy rsv-genbank --image "gcr.io/${PROJECT_ID}/rsv:${CG_VERSION}"

# FLU

gcloud run deploy flu-gisaid --image "gcr.io/${PROJECT_ID}/flu-gisaid:${CG_VERSION}"
gcloud run deploy flu-genbank --image "gcr.io/${PROJECT_ID}/flu-genbank:${CG_VERSION}"
