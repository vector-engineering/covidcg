# build environment
FROM node:12-alpine as react-build
WORKDIR /app
COPY . ./
ENV CONFIGFILE config/config_gisaid.yaml
RUN npm install
RUN npm run build-only

# server environment
FROM nginx:alpine
COPY server/htpasswd /etc/nginx/htpasswd
COPY server/cg_gcr_auth.conf /etc/nginx/conf.d/configfile.template
ENV PORT 8080
ENV HOST 0.0.0.0
RUN sh -c "envsubst '\$PORT'  < /etc/nginx/conf.d/configfile.template > /etc/nginx/conf.d/default.conf"
COPY --from=react-build /app/dist /usr/share/nginx/prod
EXPOSE 8080
CMD ["nginx", "-g", "daemon off;"]