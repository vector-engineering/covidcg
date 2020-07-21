COVID-19 CG SERVER
------------------

Since the app is all front-end, this server is just an NGINX file server with HTTP authentication

The following files were written and run on Debian 10 (buster) on a Google Compute Engine VM (1 vCPU, 600MB RAM).

Setup
-----

All the following steps require root privileges

1. Install requirements

  > sudo apt update
  > sudo apt install nginx apache2-utils
  # Make sure you also have gzip, gunzip

2. Move COVID-19 CG files

  > cp -r dist/ /var/www/covid

3. Set up NGINX

  # Remove default configuration
  > rm -f /etc/nginx/sites-enabled/default

  # Edit covid.conf, replace the server_name with your own
  # And edit the path to the correct path in /var/www
  > vim covid.conf

  # Move and enable the covid.conf file
  > cp covid.conf /etc/nginx/sites-available/covid
  > ln -s /etc/nginx/sites-available/covid /etc/nginx/sites-enabled/covid

4. Create user passwords for authentication

  # Create user passwords
  > mkdir -p /etc/apache2

  # Replace "user1" with your desired username
  # Remove the "-c" flag for subsequent user entries
  > sudo htpasswd -c /etc/apache2/.htpasswd user1

  # Confirm addition of users
  > cat /etc/apache2/.htpasswd

5. Finalize

  # Confirm NGINX configuration is valid
  > nginx -t

  # Restart NGINX service
  > systemctl restart nginx

6. Set up SSL via. LetsEncrypt (optional)

  Follow instructions: https://certbot.eff.org/lets-encrypt/debianbuster-nginx

Attributions
------------

- https://docs.nginx.com/nginx/admin-guide/security-controls/configuring-http-basic-authentication/
- https://medium.com/@jgefroh/a-guide-to-using-nginx-for-static-websites-d96a9d034940
