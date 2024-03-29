#!/bin/bash
# Build application in Google Cloud's docker context
# !! RUN FROM REPO ROOT !!
# TODO: run these builds in parallel?

if [[ -z "${CG_VERSION}" ]]; then
    echo "CG_VERSION undefined"
    exit 1
else
    echo "Building version: ${CG_VERSION}"
fi

# SARS2

gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="cg",_CONFIGFILE="config/config_sars2_gisaid.yaml",_TAG_NAME="${CG_VERSION}", . && \
gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="cg-genbank",_CONFIGFILE="config/config_sars2_genbank.yaml",_TAG_NAME="${CG_VERSION}" . && \
gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="cg-private",_CONFIGFILE="config/config_sars2_gisaid_private.yaml",_TAG_NAME="${CG_VERSION}" . && \
gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="cg-alpha",_CONFIGFILE="config/config_sars2_alpha.yaml",_TAG_NAME="${CG_VERSION}", .

# RSV

gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="rsv",_CONFIGFILE="config/config_rsv_genbank.yaml",_TAG_NAME="${CG_VERSION}" .

# FLU

gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="flu-gisaid",_CONFIGFILE="config/config_flu_gisaid.yaml",_TAG_NAME="${CG_VERSION}" .
gcloud builds submit --config build/cloudbuild.yaml --substitutions=_TARGET="flu-genbank",_CONFIGFILE="config/config_flu_genbank.yaml",_TAG_NAME="${CG_VERSION}" .