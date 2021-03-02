FROM node:12-alpine

ARG CONFIGFILE

WORKDIR /app

ADD $CONFIGFILE ./config.yaml
ENV CONFIGFILE /app/config.yaml

ADD package.json ./package.json
ADD package-lock.json ./package-lock.json
RUN npm install

# RUN npm run build-only
