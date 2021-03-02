FROM node:12-alpine

ARG CONFIGFILE

WORKDIR /app

ADD $CONFIGFILE ./config.yaml
ENV CONFIGFILE /app/config.yaml

# Add root files
ADD package.json ./package.json
ADD package-lock.json ./package-lock.json

ADD .browserslistrc ./.browserslistrc
ADD .eslintignore ./.eslintignore
ADD .npmrc ./.npmrc
ADD .watchmanconfig ./.watchmanconfig
ADD babel.config.js ./babel.config.js
ADD tsconfig.json ./tsconfig.json
ADD jsconfig.json ./jsconfig.json
ADD webpack.config.dev.js ./webpack.config.dev.js
ADD webpack.config.prod.js ./webpack.config.prod.js

ADD tools ./tools

RUN npm install

# RUN npm run build-only
