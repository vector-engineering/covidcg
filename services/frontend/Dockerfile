# --------------------
# FOR DEVELOPMENT ONLY
# --------------------

FROM node:12-alpine

ARG CONFIGFILE

WORKDIR /app

# We need to manually add all the files on the root,
# since we're only mounting the src/ folder as a bind
# volume.
#
# Why not just mount the entire project directory?
#
# We only want to do the npm install at build-time,
# and volume mounting is not done at build-time but at
# run-time. So in order to run `npm install` prior to
# run-time, we have to move all the required files for
# the install over here first

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

RUN npm ci
