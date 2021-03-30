#!/bin/bash

if [[ -z "${CG_VERSION}" ]]; then
    echo "CG_VERSION undefined"
    exit 1
else
    echo "Building version: ${CG_VERSION}"
fi

gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="cg",_CONFIGFILE="config/config_gisaid.yaml",_TAG_NAME="${CG_VERSION}", . && \
gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="cg-genbank",_CONFIGFILE="config/config_genbank.yaml",_TAG_NAME="${CG_VERSION}" . && \
gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="cg-private",_CONFIGFILE="config/config_gisaid_private.yaml",_TAG_NAME="${CG_VERSION}" .
gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="cg-alpha",_CONFIGFILE="config/config_alpha.yaml",_TAG_NAME="${CG_VERSION}", .
