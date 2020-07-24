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


Accession ID Decryption
-----------------------

Data served to clients has Accession IDs hashed to prevent dataset reconstruction, so to get the actual Accession IDs of the displayed data, the users must have to hit a server which takes in hashed IDs and return real IDs

For covidcg.org this is deployed as a Google Cloud Function, and the code for this function is in the decrypt_function/ folder. See the test_script.sh file for instructions on running it locally on your machine and a test request to see if it works. This function requires a map of hash ID -> real ID, and we have this map stored as a CSV in a Google Cloud Storage bucket. If you want, you can define the map directly in the python file or load it in as a file along with the function python script. Keep in mind that this map must be updated when the data is updated, or else hashed IDs will return as null values.


Attributions
------------

- https://docs.nginx.com/nginx/admin-guide/security-controls/configuring-http-basic-authentication/
- https://medium.com/@jgefroh/a-guide-to-using-nginx-for-static-websites-d96a9d034940
